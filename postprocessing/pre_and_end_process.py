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

def _get_mu_ex_b_ring(mu_b):
    # now only for M = 8
    try:
        interpolator = utils.load_mu_ex_interpolator()
        if interpolator is not None:
            mu_ex_b = utils.calculate_mu_ex_b(mu_b, interpolator)
        else:
            print("警告: 无法创建插值器，跳过mu_ex_b计算")
    except Exception as e:
        print(f"计算mu_ex_b时出错: {e}")

def _get_rho_total_ring(output_dir):

    block_files = utils.find_block_files(output_dir)
    if len(block_files) < 2:
        print(f"警告: 至少需要2个block文件，但只找到 {len(block_files)} 个")
        return None

    densities = []
    for block_file in block_files[1:]:  # ignore the first block   
        density = utils.load_total_density_profile(block_file)
        densities.append(density)

    # 转换为numpy数组
    densities_array = np.array(densities)  # 形状: (n_blocks, n_bins)

    # 计算平均密度分布
    # avg_density = np.mean(densities_array, axis=0)
    avg_density = np.mean(densities_array, axis=0)

    return densities_array


def _get_rho_total_linear(output_dir):
    rho_array = utils._get_block_average(output_dir=output_dir,pattern='block_*_rho_profile.dat')
    return rho_array



def pre_process_json_npz(json_path):
    """
    处理单个JSON配置文件。
    """
    # print(f"处理文件: {json_path}")

    # 读取JSON文件
    with open(json_path, 'r') as f:
        config = json.load(f)

    # 提取参数
    polymer_type = config['input_params']['polymer_type'] # Ring or Linear

    H = config['input_params']['H']
    mu_b = config['input_params']['mu_b']
    M = config['input_params']['M']
    dz = config['simulation_params']['dz']
    rho_b = config['input_params']['rho_b']
    sigma = config['input_params'].get('rcut', 1.0)  # 使用rcut作为sigma，默认为1.0
    polymer_type = config['input_params']['polymer_type']



    output_dir = utils.resolve_output_dir(config, json_path)
    if not os.path.exists(output_dir):
        print(f"警告: output目录不存在: {output_dir}")
        return None


    result={  # the same result for linear polymer and ring polymer
        'H' : H,
        'dz': dz,
        'mu_b' : mu_b,
        'M' : M,
        'sigma': sigma,
        'config_file' : os.path.basename(json_path),
        'output_dir' : output_dir,
        'polymer_type' :polymer_type
    }


    # 如果计算了Vext值，添加到结果中
    # 获取Vext参数并计算Vext(z)值    
    vext_params = config.get('Vext_params', {})
    vext_values = None
    if vext_params:
        try:
            potential_func = create_potential_function_from_params(vext_params, H)

            # 计算z坐标：与密度分布对齐，z = [0.5*dz, 1.5*dz, ..., H-0.5*dz]
            n_bins = int(round(H / dz))
            z_values = np.array([(i + 0.5) * dz for i in range(n_bins)])

            # 计算每个z点的Vext值
            vext_values = np.array([potential_func(z) for z in z_values])
        
            #result['vext_values'] = vext_values
            #result['z_values'] = z_values
            result.update({"vext_values":vext_values,"z_values":z_values})
            
        except Exception as e:
            print(f"警告: 无法计算Vext值: {e}")
            vext_values = None
    else:
        print("警告: JSON文件中没有Vext_params，无法计算Vext值")



   

    if polymer_type == "Ring":
        mu_ex_0 = config['input_params']['mu_ex_0']
        
        # 计算mu_ex_b（使用CSV文件中的三次样条插值）
        mu_ex_b = _get_mu_ex_b_ring(mu_b = mu_b)

        result.update({"mu_ex_b":mu_ex_b,"mu_ex_0":mu_ex_0})

        avg_density = _get_rho_total_ring(output_dir)  # 形状: (n_blocks, n_bins)

        
    elif polymer_type == "Linear":

        avg_density = _get_rho_total_linear(output_dir) # densities_array [rho_profile(s=1);rho_profile(s=2);,...,rho_profile(s=M)]
        result.update({"rho_profile":avg_density})

    
    else: 
        pass


    return result



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
