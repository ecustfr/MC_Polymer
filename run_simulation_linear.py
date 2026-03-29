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
# Import our data processing modules (linear polymer specific)
from system_analyzer import (
    DistributionAnalyzer,
    GlobalPropertyAnalyzer,
    cal_mono_density_profile,
    cal_insert_wall
)
from data_recorder import DataRecorder


def read_config(config_file):
    """Read configuration file"""
    with open(config_file, 'r') as f:
        return json.load(f)
    
def cal_all_monomer_density(r,dz:float,n_bins:int,M:int,bin_volume:float)->np.ndarray:
    
    profile = np.zeros((n_bins, M))
    
    z_coords = r[:, 2]
    
    for mono in range(M):
        counts, _ = np.histogram(z_coords[mono:-1:M], bins=n_bins, range=(0, n_bins * dz))
        profile[:,mono] = counts.astype(float)/bin_volume 
    
    return profile


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

    # Initialize simulation object for linear polymers
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

    # 默认行为：外势beta=1.0


    # Set external potential using the unified function
    setup_external_potential(mc_sys, input_params, config, output_dir)

    # Print all parameters
    mc_sys.print_all_parameters()


    # Initialize data recorder
    recorder = DataRecorder()

    # Initialize analyzers
    dz = simulation_params['dz']  # z-axis interval for density distribution
    n_bins = int(input_params['H'] / dz)  # Number of bins for density distribution

    H = mc_sys.get_H()
    M = mc_sys.get_M()
    # rho_profile = DistributionAnalyzer(name = "rho_two_profile", dz=dz , bins=n_bins)
    
    z_list_1 = np.arange(-0.5,0.5,dz)
    z_list_2 = np.arange(H-0.5,H+0.5,dz)
    z_list = np.concatenate( (z_list_1,z_list_2),axis=0)
    wall_bins = z_list.size

    mu_ex_wall = DistributionAnalyzer(name = 'mu_ex_wall',dz=dz, bins = wall_bins)


    rho_all_profile = np.zeros([n_bins,M])

    # Initialize current simulation parameters
    current_eps_trans = simulation_params['EPS_TRANS']
    current_k_max = simulation_params['K_MAX']
    current_rot_ratio = simulation_params['ROT_RATIO']



    # Main simulation loop
    print(f"Starting simulation, {simulation_params['sample_block']} blocks, {simulation_params['sample_time']} steps per block")

    # static_mono = 2
    # insert_time = 50
    V = mc_sys.get_box_xy() * mc_sys.get_box_xy() * mc_sys.get_H()
    bin_volume = mc_sys.get_box_xy()*mc_sys.get_box_xy()*dz

    for block in range(simulation_params['sample_block']):
        print(f"\n--- Block {block} started ---")
        n_sample = 0
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
                r = mc_sys.r_total
                rho_all_profile += cal_all_monomer_density(r,dz,n_bins,M,bin_volume)
                n_sample = n_sample+1 
                mu_ex_wall.accumulate(mc_sys,lambda sim,dz,wall_bins: cal_insert_wall(sim,z_list))



        # End block
## ---------------------------------------------------------------------------------------------------------------------------------------------------
        mc_sys.end_block(block)



        n_now = mc_sys.get_N_now()
        rho_all_profile = rho_all_profile/n_sample
        n_sample = 0
        
## ---------------------------------------------------------------------------------------------------------------------------------------------------
# Save to only-data file
        # recorder.save(f"{output_params['output_dir']}/block_{block}_rg_values.dat", block_avg_rg_values[:n_now])
        # recorder.save(f"{output_params['output_dir']}/block_{block}_{rho_profile.name}.dat", block_avg_density_profile)
        #recorder.save(f"{output_params['output_dir']}/block_{block}_{Wz_insert.name}.dat", Wz_insert_average)
        recorder.save(f"{output_params['output_dir']}/block_{block}_rho_profile.dat",rho_all_profile)
        rho_all_profile.fill(0)

        recorder.save(f"{output_params['output_dir']}/block_{block}_mu_ex_wall.dat",-np.log(mu_ex_wall.average))
        mu_ex_wall.reset()
## ---------------------------------------------------------------------------------------------------------------------------------------------------
# save to ensemble average file
        trans_acceptance = mc_sys.get_trans_acceptance()
        rot_acceptance = mc_sys.get_rot_acceptance()
        insert_acceptance = mc_sys.get_insert_acceptance()
        delete_acceptance = mc_sys.get_delete_acceptance()


        with open(f"{output_params['output_dir']}/block_{block}_ensemble_data.dat", 'w') as f:
            # f.write(f"average_rg {block_overall_avg_rg:.6f}\n")
            f.write(f"end_rho_avg {n_now/V:.6f}\n")
            f.write(f"eps_trans {current_eps_trans:.6f}\n")
            f.write(f"K_MAX {current_k_max}\n")
            # f.write(f"current_beta_ext {current_beta_ext:.6f}\n")
            f.write(f" Block {block} trans acceptance: {trans_acceptance:.4f}, rot acceptance: {rot_acceptance:.4f}\n")
            f.write(f" Block {block} insert acceptance: {insert_acceptance:.4f}, delete acceptance: {delete_acceptance:.4f}\n")
            # f.write(f" Block {block} reptation acceptance: {reptation_acceptance:.4f}\n")



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# reset simulation parameters and 

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
        #G_profile.reset()
        #rho_profile.reset()
        #W_insert.reset()
        # Reset w(z) analyzers for each point
        #for analyzer in Wz_fix_z_analyzers:
        #    analyzer.reset()
        #Wz_insert.reset()


        # Save block results
        print(f"--- Block {block} ended ---")

    # Close trajectory recorder
    # trace_recorder.close()

    # Close energy trace file
    #energy_trace_fp.close()

##---------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Save ending of simulation 
    with open(f"{output_params['output_dir']}/config.dat", 'w') as f:
        f.write(f"# Simulation Configuration\n")
        f.write(f"polymer_type: linear\n")
        f.write(f"M: {input_params['M']}\n")
        f.write(f"init_N: {input_params['init_N']}\n")
        f.write(f"rho_b: {input_params['rho_b']:.6f}\n")
        f.write(f"box_xy: {input_params['box_xy']:.6f}\n")
        f.write(f"H: {input_params['H']:.6f}\n")
        f.write(f"rcut: {input_params['rcut']:.6f}\n")
        f.write(f"max_N: {input_params['max_N']}\n")
    print(f"Simulation configuration saved to {output_params['output_dir']}/config.dat")

    # 保存最后一帧所有聚合物坐标
    try:
        # 获取坐标数组和单体总数
        positions = mc_sys.r_total  # 形状为 (MN_now, 3) 的numpy数组
        MN_now = mc_sys.get_MN_now()

        # 构建输出文件名
        coords_file = os.path.join(output_params['output_dir'], f"{output_params['output_prefix']}_final_coordinates.dat")

        # 保存坐标，每行x y z（无头信息）
        with open(coords_file, 'w') as f:
            for i in range(MN_now):
                f.write(f"{positions[i, 0]:.8f} {positions[i, 1]:.8f} {positions[i, 2]:.8f}\n")

        print(f"最终聚合物坐标已保存到: {coords_file}")
        print(f"单体总数: {MN_now}")
    except Exception as e:
        print(f"警告: 保存最终坐标时出错: {e}")

    print("\nSimulation completed!")
    print(f"Output files saved in: {output_params['output_dir']}")



def main():
    """
    Main function for linear polymer simulation
    """
    parser = argparse.ArgumentParser(description="Monte Carlo linear polymer simulation")
    parser.add_argument('--config', '-c', type=str, default='config_0000_M8_Trivial_H8.0_mu1.23.json',
                        help='Path to configuration file')

    args = parser.parse_args()

    # Read configuration file
    config = read_config(args.config)

    # Run simulation
    run_simulation(config)


if __name__ == "__main__":
    main()