#!/usr/bin/env python3
"""
后处理程序：从input/Ring_configs文件夹中的JSON文件出发，
找到对应的output文件夹，读取除了第一个block外的密度分布数据，
计算平均密度分布，并与H, dz, mu_b, M参数一起使用numpy保存。
"""

import json
import numpy as np
import os
import glob
import re
from pathlib import Path
from external_potential import create_potential_function_from_params

def natural_sort_key(s):
    """
    自然排序键函数，用于按数字顺序排序block文件。
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]

def find_block_files(output_dir):
    """
    在output_dir中找到所有block_*_rho_profile.dat文件，
    并按block编号自然排序。
    """
    pattern = os.path.join(output_dir, 'block_*_rho_profile.dat')
    files = glob.glob(pattern)
    # 按block编号排序
    files.sort(key=natural_sort_key)
    return files

def load_density_profile(file_path):
    """
    加载密度分布文件，每行一个浮点数。
    """
    return np.loadtxt(file_path)

def process_json(json_path):
    """
    处理单个JSON配置文件。
    """
    print(f"处理文件: {json_path}")

    # 读取JSON文件
    with open(json_path, 'r') as f:
        config = json.load(f)

    # 提取参数
    H = config['input_params']['H']
    mu_b = config['input_params']['mu_b']
    M = config['input_params']['M']
    dz = config['simulation_params']['dz']
    rho_b = config['input_params']['rho_b']
    sigma = config['input_params'].get('rcut', 1.0)  # 使用rcut作为sigma，默认为1.0

    # 获取Vext参数并计算Vext(z)值
    vext_params = config.get('Vext_params', {})
    vext_values = None
    if vext_params:
        try:
            # 创建势能函数
            potential_func = create_potential_function_from_params(vext_params, H)

            # 计算z坐标：与密度分布对齐，z = [0.5*dz, 1.5*dz, ..., H-0.5*dz]
            n_bins = int(round(H / dz))
            z_values = np.array([(i + 0.5) * dz for i in range(n_bins)])

            # 计算每个z点的Vext值
            vext_values = np.array([potential_func(z) for z in z_values])
            print(f"计算了Vext(z)值，形状: {vext_values.shape}")
        except Exception as e:
            print(f"警告: 无法计算Vext值: {e}")
            vext_values = None
    else:
        print("警告: JSON文件中没有Vext_params，无法计算Vext值")

    # 获取output目录
    output_dir = config['output_params']['output_dir']
    # 确保output_dir是绝对路径（相对于项目根目录）
    if not os.path.isabs(output_dir):
        # 假设JSON文件路径相对于项目根目录，计算绝对路径
        project_root = Path(json_path).parent.parent.parent
        output_dir = os.path.join(project_root, output_dir)

    if not os.path.exists(output_dir):
        print(f"警告: output目录不存在: {output_dir}")
        return None

    # 查找所有block文件
    block_files = find_block_files(output_dir)
    if len(block_files) < 2:
        print(f"警告: 至少需要2个block文件，但只找到 {len(block_files)} 个")
        return None

    print(f"找到 {len(block_files)} 个block文件")

    # 跳过第一个block（block_0），加载剩余block的数据
    densities = []
    for block_file in block_files[1:]:  # 跳过第一个
        density = load_density_profile(block_file)
        densities.append(density)

    # 转换为numpy数组
    densities_array = np.array(densities)  # 形状: (n_blocks, n_bins)

    # 计算平均密度分布
    avg_density = np.mean(densities_array, axis=0)

    # 准备保存的数据
    result = {
        'avg_density': avg_density,
        'H': H,
        'dz': dz,
        'mu_b': mu_b,
        'M': M,
        'rho_b': rho_b,
        'sigma': sigma,
        'n_blocks_used': len(densities_array),
        'n_bins': len(avg_density),
        'config_file': os.path.basename(json_path),
        'output_dir': output_dir
    }

    # 如果计算了Vext值，添加到结果中
    if vext_values is not None:
        result['vext_values'] = vext_values
        result['z_values'] = np.array([(i + 0.5) * dz for i in range(len(avg_density))])

    return result

def save_result(result, output_path):
    """
    将结果保存为numpy .npz文件。
    """
    # 准备保存的数据字典
    save_dict = {
        'avg_density': result['avg_density'],
        'H': result['H'],
        'dz': result['dz'],
        'mu_b': result['mu_b'],
        'M': result['M'],
        'rho_b': result['rho_b'],
        'sigma': result['sigma'],
        'n_blocks_used': result['n_blocks_used'],
        'n_bins': result['n_bins'],
        'config_file': result['config_file'],
        'output_dir': result['output_dir']
    }

    # 如果有Vext值，也保存
    if 'vext_values' in result:
        save_dict['vext_values'] = result['vext_values']
        save_dict['z_values'] = result['z_values']
        print(f"保存Vext值，形状: {result['vext_values'].shape}")

    # 使用savez保存多个数组
    np.savez(output_path, **save_dict)
    print(f"结果已保存到: {output_path}")

def main():
    """
    主函数：遍历input/Ring_configs中的所有JSON文件。
    """
    # 设置路径
    input_dir = os.path.join('input', 'Ring_configs')
    if not os.path.exists(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        return

    # 查找所有JSON文件
    json_pattern = os.path.join(input_dir, '*.json')
    json_files = glob.glob(json_pattern)

    if not json_files:
        print(f"未找到JSON文件: {json_pattern}")
        return

    print(f"找到 {len(json_files)} 个JSON文件")

    # 处理每个JSON文件
    for json_file in json_files:
        try:
            result = process_json(json_file)
            if result is None:
                continue

            # 从result中获取output_dir
            output_dir = result['output_dir']
            # 确保output_dir存在
            os.makedirs(output_dir, exist_ok=True)

            # 生成输出文件名
            config_name = os.path.splitext(os.path.basename(json_file))[0]
            output_file = os.path.join(output_dir, f'avg_density_{config_name}.npz')

            # 保存结果
            save_result(result, output_file)

        except Exception as e:
            print(f"处理文件 {json_file} 时出错: {e}")
            import traceback
            traceback.print_exc()

if __name__ == '__main__':
    main()