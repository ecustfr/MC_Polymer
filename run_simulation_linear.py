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
    cal_W,
    cal_wz_fix_z
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


    rho_profile = DistributionAnalyzer(name = "rho_two_profile", dz=dz , bins=n_bins)


    W_insert = GlobalPropertyAnalyzer(name = "mu")
    # Analyzer for fixed z w(z) calculation at multiple points
    num_wz_points = 3  # number of uniformly distributed points
    H = mc_sys.get_H()
    # Generate uniformly distributed points (avoiding boundaries)
    z_points = np.array([1.5,3,4.5])#np.linspace(0.5*dz, H-0.5*dz, num_wz_points)
    # Create analyzers for each point
    Wz_fix_z_analyzers = [GlobalPropertyAnalyzer(name=f"wz_fix_z_{i}") for i in range(num_wz_points)]

    # Energy trajectory file
    # energy_trace_file = f"{output_params['output_dir']}/{output_params['output_prefix']}_energy_trace.dat"

    # Wz_insert = DistributionAnalyzer(name = "Wz_profile", dz=dz , bins=n_bins)
    # Initialize W_insert accumulation variables



    trace_file = f"{output_params['output_dir']}/{output_params['output_prefix']}_trajectory.xyz"

    # Open energy trace file
    #energy_trace_file = f"{output_params['output_dir']}/{output_params['output_prefix']}_energy_trace.dat"
    #energy_trace_fp = open(energy_trace_file, 'w')
    #energy_trace_fp.write("# Step Block External_Energy\n")  # Header



    # Initialize current simulation parameters
    current_eps_trans = simulation_params['EPS_TRANS']
    current_k_max = simulation_params['K_MAX']
    current_rot_ratio = simulation_params['ROT_RATIO']



    # Main simulation loop
    print(f"Starting simulation, {simulation_params['sample_block']} blocks, {simulation_params['sample_time']} steps per block")

    static_mono = 2
    insert_time = 50
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


            # Optional polymer-specific moves (uncomment if needed)
            # mc_sys.insert_move(simulation_params['K_MAX'])
            # mc_sys.delete_move(simulation_params['K_MAX'], np.random.randint(0, mc_sys.get_N_now())   )


            # Sample data periodically
            if step % simulation_params['sample_interval'] == 0:

                # Record trajectory (every 100 sample points)
                #if step % (simulation_params['sample_interval'] * 100) == 0:
                  #trace_recorder.write_frame(mc_sys.r_total, mc_sys.get_N_now(), block * simulation_params['sample_time'] + step)

                rho_profile.accumulate(mc_sys,lambda s,b,d :cal_mono_density_profile(s,b,d,static_mono))
                #rho_profile.accumulate(mc_sys,cal_density_profile)

                # Record external potential energy to trace file
                # current_energy = mc_sys.calculate_external_energy()
                # total_step = block * simulation_params['sample_time'] + step
                # energy_trace_fp.write(f"{total_step} {block} {current_energy:.6f}\n")

                """
                Wz_insert.accumulate(
                    mc_sys,
                    lambda s,d,b :cal_Wz_profile(s,d,b,static_mono,simulation_params['K_MAX'], insert_time)
                )

                G_profile.accumulate(
                    mc_sys,
                    lambda s,d,b :cal_G_profile(s,d,b,static_mono,simulation_params['K_MAX'], insert_time)
                    )
                """


                for time in range(0,insert_time):
                    W_insert.accumulate(mc_sys,lambda s:cal_W(s,simulation_params['K_MAX']) )
                    for analyzer, z_pos in zip(Wz_fix_z_analyzers, z_points):
                        analyzer.accumulate(mc_sys,lambda s,z=z_pos:cal_wz_fix_z(s,z,static_mono,simulation_params['K_MAX']) )



        # End block
        mc_sys.end_block(block)

        # Get average results for current block
        # block_avg_rg_values = polymer_analyzer.get_average_values()
        block_avg_density_profile = rho_profile.average
        mu_block = -np.log(W_insert.average)
        # Wz_insert_average = -np.log(Wz_insert.average,where=(Wz_insert.average>0))
        # Calculate w(z) for each fixed point
        wz_fix_z_values = [-np.log(analyzer.average) if analyzer.average > 0 else 0.0 for analyzer in Wz_fix_z_analyzers]
        n_now = mc_sys.get_N_now()

        # Calculate density at fixed z positions and corresponding mu_z_i = wz_value + log(rho(z))
        mu_z_values = []
        for i, (z_pos, wz_value) in enumerate(zip(z_points, wz_fix_z_values)):
            # Convert z position to bin index
            bin_index = int(z_pos / dz)
            rho_z = block_avg_density_profile[bin_index]
            mu_z = wz_value + np.log(rho_z) if rho_z > 0 else wz_value
            mu_z_values.append(mu_z)

        #block_overall_avg_rg = np.mean(block_avg_rg_values[:n_now])

        # Save current block results, only keep actual N_now rows
        # recorder.save(f"{output_params['output_dir']}/block_{block}_rg_values.dat", block_avg_rg_values[:n_now])
        recorder.save(f"{output_params['output_dir']}/block_{block}_{rho_profile.name}.dat", block_avg_density_profile)
        #recorder.save(f"{output_params['output_dir']}/block_{block}_{Wz_insert.name}.dat", Wz_insert_average)


        # Save current block overall average radius of gyration, chemical potential, and simulation parameters
        trans_acceptance = mc_sys.get_trans_acceptance()
        rot_acceptance = mc_sys.get_rot_acceptance()
        insert_acceptance = mc_sys.get_insert_acceptance()
        delete_acceptance = mc_sys.get_delete_acceptance()
        # reptation_acceptance = mc_sys.get_reptation_acceptance()
        # Get current external beta value
        # current_beta_ext = mc_sys.get_current_external_beta()

        with open(f"{output_params['output_dir']}/block_{block}_ensemble_data.dat", 'w') as f:
            # f.write(f"average_rg {block_overall_avg_rg:.6f}\n")
            f.write(f"end_rho_avg {n_now/V:.6f}\n")
            f.write(f"average_mu {mu_block+np.log(n_now/V):.6f}\n")
            # Output mu_z_i = wz_value + log(rho(z)) for each fixed point
            for i, (z_pos, mu_z) in enumerate(zip(z_points, mu_z_values)):
                f.write(f"mu_z_{i} z={z_pos:.3f} value={mu_z:.6f}\n")
            f.write(f"eps_trans {current_eps_trans:.6f}\n")
            f.write(f"K_MAX {current_k_max}\n")
            # f.write(f"current_beta_ext {current_beta_ext:.6f}\n")
            f.write(f" Block {block} trans acceptance: {trans_acceptance:.4f}, rot acceptance: {rot_acceptance:.4f}\n")
            f.write(f" Block {block} insert acceptance: {insert_acceptance:.4f}, delete acceptance: {delete_acceptance:.4f}\n")
            # f.write(f" Block {block} reptation acceptance: {reptation_acceptance:.4f}\n")






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
        rho_profile.reset()
        W_insert.reset()
        # Reset w(z) analyzers for each point
        for analyzer in Wz_fix_z_analyzers:
            analyzer.reset()
        #Wz_insert.reset()


        # Save block results
        print(f"--- Block {block} ended ---")

    # Close trajectory recorder
    # trace_recorder.close()

    # Close energy trace file
    #energy_trace_fp.close()


    # Save simulation configuration to file
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
    parser.add_argument('--monomer-index', '-m', type=int, default=-1,
                        help='Index of the monomer to analyze density for (-1 for all monomers)')

    args = parser.parse_args()

    # Read configuration file
    config = read_config(args.config)

    # Run simulation
    run_simulation(config)


if __name__ == "__main__":
    main()