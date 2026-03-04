#!/usr/bin/env python3
"""
后处理程序：从input/Ring_configs文件夹中的JSON文件出发，
找到对应的output文件夹，读取除了第一个block外的密度分布数据，
计算平均密度分布，并与H, dz, mu_b, M参数一起使用numpy保存。
支持并行处理以加速大量文件的后处理。
"""

import json
import numpy as np
import os
import glob
import sys
import argparse
import multiprocessing
from external_potential import create_potential_function_from_params
from . import utils


def pre_process_json_npz(json_path):
    """
    处理单个JSON配置文件。
    """
    # print(f"处理文件: {json_path}")

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
    mu_ex_0 = config['input_params']['mu_ex_0']
    # 计算mu_ex_b（使用CSV文件中的三次样条插值）
    mu_ex_b = None
    try:
        interpolator = utils.load_mu_ex_interpolator()
        if interpolator is not None:
            mu_ex_b = utils.calculate_mu_ex_b(mu_b, interpolator)
            if mu_ex_b is not None:
                pass
                # print(f"计算mu_ex_b: mu_b={mu_b}, mu_ex_b={mu_ex_b}")
            else:
                print(f"警告: 无法计算mu_ex_b (mu_b={mu_b})")
        else:
            print("警告: 无法创建插值器，跳过mu_ex_b计算")
    except Exception as e:
        print(f"计算mu_ex_b时出错: {e}")

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
            # print(f"计算了Vext(z)值，形状: {vext_values.shape}")
        except Exception as e:
            print(f"警告: 无法计算Vext值: {e}")
            vext_values = None
    else:
        print("警告: JSON文件中没有Vext_params，无法计算Vext值")

    # 获取output目录（使用utils中的函数解析为绝对路径）
    output_dir = utils.resolve_output_dir(config, json_path)

    if not os.path.exists(output_dir):
        print(f"警告: output目录不存在: {output_dir}")
        return None

    # 查找所有block文件
    block_files = utils.find_block_files(output_dir)
    if len(block_files) < 2:
        print(f"警告: 至少需要2个block文件，但只找到 {len(block_files)} 个")
        return None

    # print(f"找到 {len(block_files)} 个block文件")

    # 跳过第一个block（block_0），加载剩余block的数据
    densities = []
    for block_file in block_files[1:]:  # 跳过第一个
        density = utils.load_density_profile(block_file)
        densities.append(density)

    # 转换为numpy数组
    densities_array = np.array(densities)  # 形状: (n_blocks, n_bins)

    # 计算平均密度分布
    avg_density = np.mean(densities_array, axis=0)

    # 准备保存的数据
    result = {
        'rho_profile': avg_density,
        'H': H,
        'dz': dz,
        'mu_b': mu_b,
        'mu_ex_b': mu_ex_b,
        'mu_ex_0':mu_ex_0,
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

"""
def save_result_npz(result, output_path):
    
    save_dict = {
        'avg_density': result['avg_density'],
        'H': result['H'],
        'dz': result['dz'],
        'mu_b': result['mu_b'],
        'mu_ex_b': result['mu_ex_b'] if ('mu_ex_b' in result and result['mu_ex_b'] is not None) else np.nan,
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

    # 使用utils中的函数保存数据
    utils.save_npz_data(output_path, save_dict)
    print(f"结果已保存到: {output_path}")
"""

def pre_single_json_to_npz(json_file):
    """
    处理单个JSON文件，计算平均密度并保存结果。
    返回 (json_file, success, error_message)
    """
    try:
        save_paras = pre_process_json_npz(json_file) # 什么时候返回 None
        if save_paras is None:
            return (json_file, False, "process_json返回None，可能缺少数据")

        # 生成输出文件名
        config_name = utils.extract_config_name(json_file)
        output_file = os.path.join(save_paras['output_dir'], f'avg_density_{config_name}.npz')

        # 保存结果
        utils.save_npz_data(output_file, save_paras)
        # save_result_npz(save_paras, output_file)
        return (json_file, True, "")

    except Exception as e:
        error_msg = f"{e}"
        import traceback
        traceback_str = traceback.format_exc()
        print(f"处理文件 {json_file} 时出错: {error_msg}")
        print(traceback_str)
        return (json_file, False, error_msg)


def end_recal_npz():
    pass

def main():
    pass

if __name__ == '__main__':
    main()