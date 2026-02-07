#!/usr/bin/env python3

from itertools import accumulate
import os
import sys
import json
import argparse
import numpy as np
from datetime import datetime
import platform

# Import external potential functions
from external_potential import create_potential_function

# Add current directory to Python path for importing pymcpolymer module
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

# Add python directory to Python path for importing our data processing modules
python_dir = os.path.join(script_dir, 'python')
sys.path.append(python_dir)

# Import pymcpolymer module after adding paths
release_path = os.path.join(script_dir,"./Release")
sys.path.append(release_path)

import pymcpolymer
print("pymcpolymer_path:", pymcpolymer.__file__)
# Import our data processing modules
from system_analyzer import (
    PolymerAnalyzer,
    DistributionAnalyzer,
    GlobalPropertyAnalyzer,
    cal_polymer_rg,
    cal_density_profile,
    cal_mono_density_profile,
    cal_W,
    cal_G_profile,
    cal_Wz_profile
)
from data_recorder import DataRecorder
#from trace_recorder import TraceRecorder


def read_config(config_file):
    """Read configuration file"""
    with open(config_file, 'r') as f:
        return json.load(f)


def run_simulation(config):

    project_root = os.path.dirname(os.path.abspath(__file__))
    # Extract configuration parameters
    input_params = config['input_params']
    simulation_params = config['simulation_params']
    output_params = config['output_params']
    

    
    # Convert output directory to absolute path based on project root
    # Ensure all output goes to the project folder
    output_dir =  output_params['output_dir']

    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"Output directory: {output_dir}")
    
    # Initialize simulation object
    if input_params['polymer_type'] == 'ring':
        mc_sys = pymcpolymer.MuVT_MC_RingPolymer(
            input_params['configuration'],
            input_params['mu_b'],
            input_params['M'],
            input_params['init_N'],
            input_params['rho_b'],
            input_params['box_xy'],
            input_params['H'],
            input_params['rcut'],
            input_params['max_N']
        )
    else:
        mc_sys = pymcpolymer.MuVT_MC_LinearPolymer(
            input_params['configuration'],
            input_params['mu_b'],
            input_params['M'],
            input_params['init_N'],
            input_params['rho_b'],
            input_params['box_xy'],
            input_params['H'],
            input_params['rcut'],
            input_params['max_N']
        )
    
    # Initialize system
    print("Initializing system...")
    mc_sys.init_second()
    
    # Set simulation parameters
    mc_sys.set_sim_parameters(
        simulation_params['EPS_TRANS'],
        simulation_params['ROT_RATIO'],
        simulation_params['K_MAX']
    )
    
    # Set external potential if specified
    external_potential = input_params['external_potential']
    print(f"External potential type: {external_potential}")
    
    if external_potential == 'custom':
        # Check if Vext_params exists
        if 'Vext_params' in config:
            Vext_params = config['Vext_params']
           
            # Extract parameters using try-except
            try:
                An = Vext_params['An']
                phi_n = Vext_params['phi_n']
                Vlin_par = Vext_params['Vlin_par']
                x_tar = Vext_params['x_tar']
                C = Vext_params['C']
                box_size_z = input_params['H']  # Use H as box size in z-direction
            except KeyError as e:
                raise ValueError(f"Missing required external potential parameter: {e}")
            
            # Create potential function
            potential_func = create_potential_function(An, phi_n, box_size_z, Vlin_par, x_tar, C)
           
            # Set external potential
            mc_sys.set_external_potential(potential_func, 'custom_potential')
            print("Custom external potential set successfully")
        else:
            print("Warning: external_potential is 'custom' but Vext_params is not provided")
            print("Using no external potential instead")
    elif external_potential == 'None':
        mc_sys.set_external_potential("hs_wall")
        print("No external potential will be used")
    else:
        print(f"Unknown external potential type: {external_potential}")
        print("Using no external potential instead")
    
    # Print all parameters
    mc_sys.print_all_parameters()
    
    
    # Initialize data recorder
    recorder = DataRecorder()
    
    # Initialize analyzers
    dz = simulation_params['dz']  # z-axis interval for density distribution
    n_bins = int(input_params['H'] / dz)  # Number of bins for density distribution

    """
    polymer_analyzer = PolymerAnalyzer(
        M=input_params['M'],
        N=input_params['max_N'],
        trace_or_not=np.False_,
        max_samples=0
    )
    """

    rho_profile = DistributionAnalyzer(name = "rho_two_profile", dz=dz , bins=n_bins)
    
    G_profile = DistributionAnalyzer(name = "GLGR" , dz = dz , bins = n_bins)
    
    W_insert = GlobalPropertyAnalyzer(name = "mu")

    Wz_insert = DistributionAnalyzer(name = "Wz_profile", dz=dz , bins=n_bins)
        # Initialize W_insert accumulation variables


    #trace_file = f"{output_params['output_dir']}/{output_params['output_prefix']}_trajectory.xyz"
    #trace_recorder = TraceRecorder(trace_file)
    


    # Initialize current simulation parameters
    current_eps_trans = simulation_params['EPS_TRANS']
    current_k_max = simulation_params['K_MAX']
    current_rot_ratio = simulation_params['ROT_RATIO']
    

    
    # Main simulation loop
    print(f"Starting simulation, {simulation_params['sample_block']} blocks, {simulation_params['sample_time']} steps per block")
    
    static_mono = 2
    insert_time = 800
    V = mc_sys.get_box_xy() * mc_sys.get_box_xy() * mc_sys.get_H()
    for block in range(simulation_params['sample_block']):
        print(f"\n--- Block {block} started ---")
        
        # Reset simulation records
        mc_sys.reset_mc_record()
        
        # Intra-block simulation loop
        for step in range(simulation_params['sample_time']):
            if step % 10000 == 0:
                print(f"  Block {block}, step {step}/{simulation_params['sample_time']}")
            
            mc_sys.mc_one_step_NVT()
            mc_sys.mc_one_step_MuVT()
            
            # Sample data periodically
            if step % simulation_params['sample_interval'] == 0:
                
                # Record trajectory (every 100 sample points)
                #if step % (simulation_params['sample_interval'] * 100) == 0:
                #    trace_recorder.write_frame(mc_sys.r_total, mc_sys.get_N_now(), block * simulation_params['sample_time'] + step)
                
                rho_profile.accumulate(mc_sys,lambda s,b,d :cal_mono_density_profile(s,b,d,static_mono))
                #rho_profile.accumulate(mc_sys,cal_density_profile)
                
                Wz_insert.accumulate(
                    mc_sys,
                    lambda s,d,b :cal_Wz_profile(s,d,b,static_mono,simulation_params['K_MAX'], insert_time)
                )
                
                """
                G_profile.accumulate(
                    mc_sys,
                    lambda s,d,b :cal_G_profile(s,d,b,static_mono,simulation_params['K_MAX'], insert_time)
                    )
                """
                W_insert.accumulate(mc_sys,lambda s:cal_W(s,simulation_params['K_MAX']) )
                W_insert.accumulate(mc_sys,lambda s:cal_W(s,simulation_params['K_MAX']) )
                W_insert.accumulate(mc_sys,lambda s:cal_W(s,simulation_params['K_MAX']) )
                W_insert.accumulate(mc_sys,lambda s:cal_W(s,simulation_params['K_MAX']) )
                #distribution_analyzer.accumulate(mc_sys, calculate_density_profile)

                
        
        # End block
        mc_sys.end_block(block)
        
        
        
        

        # Get average results for current block
        # block_avg_rg_values = polymer_analyzer.get_average_values()
        block_avg_density_profile = rho_profile.average
        mu_block = -np.log(W_insert.average)
        Wz_insert_average = -np.log(Wz_insert.average,where=(Wz_insert.average>0))
        n_now = mc_sys.get_N_now()

        #block_overall_avg_rg = np.mean(block_avg_rg_values[:n_now])
        
        # Save current block results, only keep actual N_now rows
        # recorder.save(f"{output_params['output_dir']}/block_{block}_rg_values.dat", block_avg_rg_values[:n_now])
        recorder.save(f"{output_params['output_dir']}/block_{block}_{rho_profile.name}.dat", block_avg_density_profile)
        recorder.save(f"{output_params['output_dir']}/block_{block}_{Wz_insert.name}.dat", Wz_insert_average)

        
        # Save current block overall average radius of gyration, chemical potential, and simulation parameters
        trans_acceptance = mc_sys.get_trans_acceptance()
        rot_acceptance = mc_sys.get_rot_acceptance()
        insert_acceptance = mc_sys.get_insert_acceptance()
        delete_acceptance = mc_sys.get_delete_acceptance()
        with open(f"{output_params['output_dir']}/block_{block}_ensemble_data.dat", 'w') as f:
            # f.write(f"average_rg {block_overall_avg_rg:.6f}\n")
            f.write(f"end_rho_avg {n_now/V:.6f}\n")
            f.write(f"average_mu {mu_block+np.log(n_now/V):.6f}\n")
            f.write(f"eps_trans {current_eps_trans:.6f}\n")
            f.write(f"K_MAX {current_k_max}\n")
            f.write(f" Block {block} trans acceptance: {trans_acceptance:.4f}, rot acceptance: {rot_acceptance:.4f}\n")
            f.write(f" Block {block} insert acceptance: {insert_acceptance:.4f}, delete acceptance: {delete_acceptance:.4f}\n")
        
        

        
       

        print(f"  Parameters before adjustment: EPS_TRANS={current_eps_trans}, K_MAX={current_k_max}")
        if trans_acceptance > 0.5:
            current_eps_trans *= 2
            print(f"  Trans acceptance ratio above 50%, increasing EPS_TRANS")
        elif trans_acceptance < 0.2:
            current_eps_trans *= 0.2
            print(f"  Trans acceptance ratio below 20%, decreasing EPS_TRANS")
        
        # Set EPS_TRANS bounds
        current_eps_trans = max(0.001, min(current_eps_trans, 1.0))
        
        # Output adjusted parameters
        print(f"  Parameters after adjustment: EPS_TRANS={current_eps_trans:.6f}")
        
        # Update simulation parameters
        mc_sys.set_sim_parameters(current_eps_trans, current_rot_ratio, current_k_max)
        
        # Reset analyzers to accumulate data from scratch for next block
        G_profile.reset()
        rho_profile.reset()
        W_insert.reset()
        Wz_insert.reset()        

        
        # Save block results
        print(f"--- Block {block} ended ---")
    
    # Close trajectory recorder
    # trace_recorder.close()
    
    
    # Save simulation configuration to file
    with open(f"{output_params['output_dir']}/config.dat", 'w') as f:
        f.write(f"# Simulation Configuration\n")
        f.write(f"polymer_type: {input_params['polymer_type']}\n")
        f.write(f"M: {input_params['M']}\n")
        f.write(f"init_N: {input_params['init_N']}\n")
        f.write(f"rho_b: {input_params['rho_b']:.6f}\n")
        f.write(f"box_xy: {input_params['box_xy']:.6f}\n")
        f.write(f"H: {input_params['H']:.6f}\n")
        f.write(f"rcut: {input_params['rcut']:.6f}\n")
        f.write(f"max_N: {input_params['max_N']}\n")
    print(f"Simulation configuration saved to {output_params['output_dir']}/config.dat")
    
    print("\nSimulation completed!")
    print(f"Output files saved in: {output_params['output_dir']}")
    


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(description="Monte Carlo polymer simulation")
    parser.add_argument('--config', '-c', type=str, default='config.json',
                        help='Path to configuration file')
    parser.add_argument('--monomer-index', '-m', type=int, default=-1,
                        help='Index of the monomer to analyze density for (-1 for all monomers)')
    
    args = parser.parse_args()
    
    # Read configuration file
    
    
    try:
        config = read_config(args.config)
    except FileNotFoundError:
        # 重新抛出带有自定义信息的异常
        raise FileNotFoundError(f"Configuration file {args.config} not found")
    except Exception as e:
        # (可选) 捕获文件存在但格式不对等其他错误
        raise RuntimeError(f"Error reading configuration file: {e}")
    
    # Run simulation
    run_simulation(config)


if __name__ == "__main__":
    main()
