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
from external_potential import calculate_Vext, plot_potential

def plot_custom_potential(json_file):
    """
    Plot custom external potential from JSON file
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
        
        if external_potential != 'custom':
            print(f"External potential is not custom: {external_potential}")
            return
        
        # 获取 Vext 参数
        vext_params = data.get('Vext_params', {})
        An = vext_params.get('An', [0.0] * 4)
        phi_n = vext_params.get('phi_n', [0.0] * 4)
        Vlin_par = vext_params.get('Vlin_par', [[0.0] * 4, [0.0] * 4])
        x_tar = vext_params.get('x_tar', [[0.0] * 4, [0.0] * 4])
        C = vext_params.get('C', 1.0)
        
        print(f"Loaded parameters: An={An}, C={C}, H={box_size_z}")
        
        # 使用新的plot_potential函数绘制势能分布图
        output_file = os.path.join(os.path.dirname(json_file), 'custom_potential_plot.png')
        plot_potential(An, phi_n, box_size_z, Vlin_par, x_tar, C=C, output_file=output_file)

    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot custom external potential from JSON file')
    parser.add_argument('json_file', type=str, help='Path to input JSON file')
    args = parser.parse_args()
    
    plot_custom_potential(args.json_file)