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
from external_potential import setup_external_potential

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
# Import our data processing modules (ring polymer specific)
from system_analyzer import (
    DistributionAnalyzer,
    GlobalPropertyAnalyzer,
    cal_mono_density_profile,
    cal_density_profile,
    # cal_W_ring,
    # cal_W_ring_z
)
from data_recorder import DataRecorder


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

    # Initialize simulation object for ring polymers
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

    # Initialize system
    print("Initializing system...")
    mc_sys.init_second()

    # Set simulation parameters
    mc_sys.set_sim_parameters(
        simulation_params['EPS_TRANS'],
        simulation_params['ROT_RATIO'],
        simulation_params['K_MAX']
    )

    # 默认行为：外势beta=1.0
    mc_sys.set_external_beta(1.0)
    print("Using default beta_ext=1.0")

    # Set external potential using the unified function
    setup_external_potential(mc_sys, input_params, config, output_dir)

    # Print all parameters
    mc_sys.print_all_parameters()


    # Initialize data recorder
    recorder = DataRecorder()

    # Initialize analyzers
    dz = simulation_params['dz']  # z-axis interval for density distribution
    n_bins = int(input_params['H'] / dz)  # Number of bins for density distribution


    rho_profile = DistributionAnalyzer(name = "rho_profile", dz=dz , bins=n_bins)


    # W_insert = GlobalPropertyAnalyzer(name = "mu")
    #
    #
    # num_wz_points = 3  # number of uniformly distributed points
    # H = mc_sys.get_H()
    # z_points = np.array([1.5,3,4.5])#np.linspace(0.5*dz, H-0.5*dz, num_wz_points)
    # Wz_fix_z_analyzers = [GlobalPropertyAnalyzer(name=f"wz_fix_z_{i}") for i in range(num_wz_points)]
    

    trace_file = f"{output_params['output_dir']}/{output_params['output_prefix']}_trajectory.xyz"




    # Initialize current simulation parameters
    current_eps_trans = simulation_params['EPS_TRANS']
    current_k_max = simulation_params['K_MAX']
    current_rot_ratio = simulation_params['ROT_RATIO']



    # Main simulation loop
    print(f"Starting simulation, {simulation_params['sample_block']} blocks, {simulation_params['sample_time']} steps per block")

    # static_mono = 2
    # insert_time = 50
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

                #rho_profile.accumulate(mc_sys,lambda s,b,d :cal_mono_density_profile(s,b,d,static_mono))
                rho_profile.accumulate(mc_sys,cal_density_profile)
                
                # for time in range(0,insert_time):
                #     W_insert.accumulate(mc_sys,lambda s:cal_W_ring(s,simulation_params['K_MAX']) )
                #     for analyzer, z_pos in zip(Wz_fix_z_analyzers, z_points):
                #         analyzer.accumulate(mc_sys,lambda s,z=z_pos:cal_W_ring_z(s,z,simulation_params['K_MAX']) )



        # End block
        mc_sys.end_block(block)

        # Get average results for current block
        block_avg_density_profile = rho_profile.average
        # mu_block = -np.log(W_insert.average)
        
        n_now = mc_sys.get_N_now()

        # Calculate w(z) for each fixed point
        # wz_fix_z_values = [-np.log(analyzer.average) if analyzer.average > 0 else 0.0 for analyzer in Wz_fix_z_analyzers]
        

        # Calculate density at fixed z positions and corresponding mu_z_i = wz_value + log(rho(z))
        
        # mu_z_values = []
        # for i, (z_pos, wz_value) in enumerate(zip(z_points, wz_fix_z_values)):
        #     # Convert z position to bin index
        #     bin_index = int(z_pos / dz)
        #     rho_z = block_avg_density_profile[bin_index]/input_params['M']
        #     mu_z = wz_value + np.log(rho_z) if rho_z > 0 else np.inf
        #     mu_z_values.append(mu_z)
        

        # Save current block results, only keep actual N_now rows
        recorder.save(f"{output_params['output_dir']}/block_{block}_{rho_profile.name}.dat", block_avg_density_profile)


        # Save current block overall average radius of gyration, chemical potential, and simulation parameters
        trans_acceptance = mc_sys.get_trans_acceptance()
        rot_acceptance = mc_sys.get_rot_acceptance()
        insert_acceptance = mc_sys.get_insert_acceptance()
        delete_acceptance = mc_sys.get_delete_acceptance()

        with open(f"{output_params['output_dir']}/block_{block}_ensemble_data.dat", 'w') as f:
            f.write(f"end_rho_avg {n_now/V:.6f}\n")
            # f.write(f"average_mu {mu_block+np.log(n_now/V):.6f}\n")
            # Output mu_z_i = wz_value + log(rho(z)) for each fixed point
            # for i, (z_pos, mu_z) in enumerate(zip(z_points, mu_z_values)):
            #     f.write(f"mu_z_{i} z={z_pos:.3f} value={mu_z:.6f}\n")
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
        rho_profile.reset()
        # W_insert.reset()
        # Reset w(z) analyzers for each point
        # for analyzer in Wz_fix_z_analyzers:
        #     analyzer.reset()


        # Save block results
        print(f"--- Block {block} ended ---")



    # Save simulation configuration to file
    with open(f"{output_params['output_dir']}/config.dat", 'w') as f:
        f.write(f"# Simulation Configuration\n")
        f.write(f"polymer_type: ring\n")
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
    Main function for ring polymer simulation
    """
    parser = argparse.ArgumentParser(description="Monte Carlo ring polymer simulation")
    parser.add_argument('--config', '-c', type=str, 
                        help='Path to configuration file')

    args = parser.parse_args()

    # Read configuration file
    config = read_config(args.config)

    # Run simulation
    run_simulation(config)


if __name__ == "__main__":
    main()