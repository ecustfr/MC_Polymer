#!/usr/bin/env python3
"""
Script to read input JSON file and plot custom external potential distribution
"""

import os
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
import traceback

# Import external potential functions
from external_potential import (
    calculate_custom_Vext,
    calculate_step_Vext,
    plot_potential,
    plot_step_potential,
    create_potential_function_from_params
)

def plot_custom_potential(json_file):
    """
    Plot external potential from JSON file (supports both custom and step potentials)
    """
    try:
        # 检查文件是否存在
        if not os.path.exists(json_file):
            print(f"Error: File not found: {json_file}")
            return

        # 读取 JSON 文件
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # 获取基础参数
        input_params = data.get('input_params', {})
        external_potential = input_params.get('external_potential', '')
        box_size_z = input_params.get('H', 10.0)

        if external_potential not in ['custom', 'step']:
            print(f"External potential is not custom or step: {external_potential}")
            return

        # 获取 Vext 参数
        vext_params = data.get('Vext_params', {})

        # 创建势能函数
        potential_func = create_potential_function_from_params(vext_params, box_size_z)

        # 根据势能类型生成文件名
        if external_potential == 'custom':
            An = vext_params.get('An', [0.0] * 4)
            phi_n = vext_params.get('phi_n', [0.0] * 4)
            Vlin_par = vext_params.get('Vlin_par', [[0.0] * 4, [0.0] * 4])
            x_tar = vext_params.get('x_tar', [[0.0] * 4, [0.0] * 4])
            C = vext_params.get('C', 1.0)
            print(f"Loaded custom parameters: An={An}, C={C}, H={box_size_z}")

            # 使用plot_potential函数绘制势能分布图
            output_file = os.path.join(os.path.dirname(json_file), 'custom_potential_plot.png')
            plot_potential(An, phi_n, box_size_z, Vlin_par, x_tar, C=C, output_file=output_file)

        elif external_potential == 'step':
            boundaries = vext_params.get('boundaries', [0.0, box_size_z])
            potentials = vext_params.get('potentials', [0.0])
            C = vext_params.get('C', 1.0)
            print(f"Loaded step parameters: boundaries={boundaries}, potentials={potentials}, C={C}, H={box_size_z}")

            # 使用plot_step_potential函数绘制势能分布图
            output_file = os.path.join(os.path.dirname(json_file), 'step_potential_plot.png')
            plot_step_potential(boundaries, potentials, box_size_z, C=C, output_file=output_file)

        # 额外输出沿z方向的势能值列表，硬编码dz=0.05
        dz = 0.05
        z_values = np.arange(0, box_size_z + dz/2, dz)  # 包含端点
        potential_values = []
        for z in z_values:
            V = potential_func(z)
            potential_values.append(V)

        # 输出到文本文件
        values_file = os.path.join(os.path.dirname(json_file), f'{external_potential}_potential_values.txt')
        with open(values_file, 'w', encoding='utf-8') as f:
            f.write(f"# {external_potential.capitalize()} External Potential Values\n")
            f.write(f"# z\tpotential\n")
            for z, V in zip(z_values, potential_values):
                f.write(f"{z:.6f}\t{V:.6f}\n")

        print(f"Potential values saved to: {values_file}")
        print(f"Number of points: {len(z_values)}, dz={dz}")
        print(f"z range: [{z_values[0]:.2f}, {z_values[-1]:.2f}]")

        # 可选：打印前几个值到控制台
        print("\nFirst 5 values:")
        for i in range(min(5, len(z_values))):
            print(f"  z={z_values[i]:.2f}, V={potential_values[i]:.6f}")

    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot custom external potential from JSON file')
    parser.add_argument('json_file', type=str, help='Path to input JSON file')
    args = parser.parse_args()

    plot_custom_potential(args.json_file)
    #json_file = "config_0001_M6_Linear_H6.0_mu0.83.json"
    # plot_custom_potential(json_file)