#!/usr/bin/env python3
"""
处理output文件夹中内容并再计算的整合功能。

利用postprocess_ring_density.py和recal_equ.py两个库，自动化完成以下流程：
1. 扫描input/Ring_configs目录下的所有JSON配置文件。
2. 对每个JSON文件，运行后处理生成平均密度分布NPZ文件（如果尚未生成）。
3. 对每个NPZ文件，运行recal_equ.py中的重新计算。
4. 保存重新计算结果。

支持并行处理以加速大量文件的处理。

使用方法：
    python process_outputs.py [--add-extra] [--overwrite] [--max-workers] [--no-ml]

选项：
    --add-extra: 在JSON文件中添加额外字段
    --overwrite: 即使NPZ文件已存在也重新生成后处理文件
    --max-workers: 并行工作进程数，默认为CPU核心数，设为1则顺序处理
    --no-ml: 不应用于机器学习，只生成NPZ文件，不进行recal和数据集分配
"""

import os
import sys
import json
import argparse
import numpy as np
import datetime
from pathlib import Path
import concurrent.futures
import multiprocessing
import random

# 添加父目录到sys.path以便导入postprocessing模块
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# 导入现有库的功能
from postprocessing import utils
from postprocessing.pre_and_end_process import pre_single_json_to_npz
from postprocessing.recal_equ import recal_simulation


def add_extra_fields(json_path, extra_fields=None):
    """
    向 JSON 文件中合并新的字段。支持嵌套字典的深度合并，不会覆盖原有其他结构。
    """
    # === 定义一个内部的深度合并工具函数 ===
    def deep_update(original_dict, new_dict):
        for key, value in new_dict.items():
            # 如果传进来的值是一个字典
            if isinstance(value, dict):
                # 如果原配置中这个键存在，且也是个字典，就钻进去继续合并
                if key in original_dict and isinstance(original_dict[key], dict):
                    deep_update(original_dict[key], value)
                else:
                    # 如果原配置没有这个键，或者原来不是字典，直接赋值
                    original_dict[key] = value
            else:
                # 如果是普通的值（数字、字符串、列表等），直接更新或新增
                original_dict[key] = value

    try:
        # 1. 读取原有的 JSON 文件
        with open(json_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
            

        if extra_fields:
            # 2. 执行深度合并
            deep_update(config, extra_fields)

            # 3. 保存更新后的 JSON
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(config, f, ensure_ascii=False, indent=2)
            
            # print(f"已成功合并新字段到文件: {json_path}")
        else:
            print(f"无需更新: extra_fields 为空 ({json_path})")

    except Exception as e:
        raise RuntimeError(f"更新JSON文件失败 {json_path}: {e}")

def ensure_npz_file(json_path, overwrite=False):
    """
    确保JSON文件对应的NPZ文件存在。

    参数:
        json_path: JSON配置文件路径
        overwrite: 如果为True，即使NPZ文件已存在也重新生成

    返回:
        npz_file: NPZ文件路径，生成失败时抛出异常
    """
    # 获取预期的NPZ文件路径
    npz_file = utils.get_npz_path_from_json(json_path)

    # 检查NPZ文件是否存在且overwrite为False
    if os.path.exists(npz_file) and not overwrite:
        print(f"NPZ文件已存在: {npz_file}")
        return npz_file

    # print(f"生成NPZ文件: {npz_file}")

    # 运行后处理
    pre_single_json_to_npz(json_path)
    

    # 确保输出目录存在
    #output_dir = save_result['output_dir']
    # os.makedirs(output_dir, exist_ok=True)
    #print("预处理结束")
    # 保存结果
    #utils.save_npz_data(output_dir, save_result)
    #print(f"成功生成NPZ文件: {npz_file}")

    return npz_file


def process_single_json_file(json_file, overwrite=False, add_json_or_not=False, add_json=None, ml=True):
    """
    处理单个JSON文件（不包含绘图）。
    用于并行处理。

    参数:
        json_file: JSON文件路径
        overwrite: 是否覆盖已存在的NPZ文件
        add_extra: 是否在JSON中添加额外字段
        extra_fields: 要添加的额外字段字典，格式为{"section.key": value}
        ml: 是否应用于机器学习，默认True。如果为False，只生成NPZ文件，不进行recal

    返回:
        如果ml为True: 返回recal_simulation的结果（字典或None）
        如果ml为False: 返回包含npz_file路径的字典{'npz_file': npz_file}
    """
    # 步骤1: 可选 - 在JSON中添加额外字段
    if add_json_or_not:
        add_extra_fields(json_file, add_json)

    # 步骤2: 必要 - 确保NPZ文件存在
    npz_file = ensure_npz_file(json_file, overwrite=overwrite)

    # 步骤3: 根据ml标志决定是否运行重新计算
    if ml:
        recal_npz_file = recal_simulation(npz_file)
        return recal_npz_file
    else:
        # 只生成NPZ文件，不进行recal
        return {'npz_file': npz_file}    
    



def find_output_json_from_input(input_json_path):
    """
    从input目录的JSON文件找到对应的output目录中的JSON文件。

    参数:
        input_json_path: input目录中的JSON文件路径

    返回:
        str: output目录中的JSON文件路径，如果找不到则返回None
    """
    try:
        with open(input_json_path, 'r') as f:
            config = json.load(f)

        # 使用utils中的resolve_output_dir函数获取output目录
        output_dir = utils.resolve_output_dir(config, input_json_path)

        # 获取JSON文件名
        json_filename = os.path.basename(input_json_path)

        # 构建output目录中的JSON文件路径
        output_json_path = os.path.join(output_dir, json_filename)

        if os.path.exists(output_json_path):
            return output_json_path
        else:
            print(f"警告: output JSON文件不存在: {output_json_path}")
            return None
    except Exception as e:
        print(f"查找output JSON文件失败 {input_json_path}: {e}")
        return None



def main():
    parser = argparse.ArgumentParser(description='处理output文件夹中内容并再计算的整合功能')
    parser.add_argument('--input-dir', default=os.path.join('input', 'bulk_Linear_M6'), #INPUT2
                       help='input目录中的JSON配置文件目录')
    parser.add_argument('--overwrite', action='store_true',default=True,
                       help='即使NPZ文件已存在也重新生成后处理文件')
    parser.add_argument('--max-workers', type=int, default=6,
                       help='并行工作进程数，默认为CPU核心数，设为1则顺序处理')
    parser.add_argument('--error-log', type=str, default=None,
                       help='失败文件日志输出路径（JSON格式），如果不指定则不保存')
    parser.add_argument('--extra-field', action='append', default=[],
                       help='要添加的额外字段，格式为"key=value"或"section.key=value"，可多次使用')
    parser.add_argument('--no-ml', action='store_false', dest='ml', default=False,
                       help='不应用于机器学习，只生成NPZ文件，不进行recal和数据集分配')

    args = parser.parse_args()

    # 解析extra-field参数为字典
    
    add_json = None # {"input_params":{"mu_ex_0":0}}   {"input_params":{"sigma":1,"mu_ex_0":4.848994294}}

    add_json_or_not = bool(add_json)  # 如果extra_fields非空，则启用add_extra

    # 检查输入目录是否存在

    # 查找input目录中的所有JSON文件
    import glob
    input_json_pattern = os.path.join(args.input_dir, '*.json')
    input_json_files = glob.glob(input_json_pattern)

    # 查找对应的output目录中的JSON文件
    output_json_files = []
    for input_json in input_json_files:
        output_json = find_output_json_from_input(input_json)
        if output_json:
            output_json_files.append(output_json)

    json_files = output_json_files  # 要处理的是output目录中的JSON文件

    print(f"找到 {len(json_files)} 个JSON文件")

    # 确定工作进程数

    max_workers = args.max_workers

    if max_workers > 1:
        print(f"使用 {max_workers} 个工作进程进行并行处理")
    else:
        print("使用顺序处理模式")

    success_files = []
    failed_files = []

    need_to_sum_files = [] # 后续处理的文件
    #need_to_sum_indexs = [] # 后续处理的主键值 
    # 如果只有一个文件或指定单进程，则顺序处理（便于调试）

    print("使用并行处理模式")

        # 使用ProcessPoolExecutor进行并行处理
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
        future_to_json = {}
        for json_file in json_files:
            future = executor.submit(
                process_single_json_file,
                json_file, args.overwrite, add_json_or_not, add_json, args.ml
            )
            future_to_json[future] = json_file

            # 收集结果并显示进度
        completed = 0
        total = len(json_files)

        for future in concurrent.futures.as_completed(future_to_json):
            json_file = future_to_json[future]
            completed += 1

            try:
                sim_result = future.result() # future.result() 返回值对应的
                if sim_result is not None:
                    success_files.append(json_file)
                    if args.ml:
                        # 只有机器学习模式才收集到need_to_sum_files
                        #need_to_sum_indexs.append(sim_result['sim_index'])
                        #need_to_sum_files.append(sim_result['path'])
                        need_to_sum_files.append(sim_result)
                        # print(f"[{completed}/{total}] 成功处理: {json_file}")
            except Exception as e:
                error_info = {
                    'json_file': json_file,
                    'error_message': str(e),
                    'exception': type(e).__name__
                }
                failed_files.append(error_info)
                print(f"[{completed}/{total}] 处理失败: {json_file} - {e}")

    print("\n" + "="*60)
    print("处理完成!")
    total_files = len(success_files) + len(failed_files)
    print(f"总计文件: {total_files}")
    print(f"成功处理: {len(success_files)} 个配置文件")
    print(f"处理失败: {len(failed_files)} 个配置文件")

    if args.ml:
        # 只有机器学习模式才进行数据集分配和保存
        percent = [0.8,0.2]
        train_num = int(percent[0]*len(success_files))
        test_num = int(percent[1]*len(success_files))
        # valid_num = len(success_files) - train_num - test_num
        
        simData = {
            'train':{},
            'test':{}
            }
        shuffle_or_not = True
        if shuffle_or_not:
            random.shuffle(need_to_sum_files)
        # result = {'sim_index':config_index, 'path':npy_path}  
    #-----------------------------------------------------------------------------------------------------------
        for (i,result) in enumerate(need_to_sum_files):
            key = "sim_"+result['sim_index']
            value = np.load(result['path']) # ,mmap_mode='r'
    
            if i<train_num:
                simData['train'][key] = value
            else:
                simData['test'][key] = value
    
            
    
        np.save('simData2.npy',simData,allow_pickle = True)
    

#-----------------------------------------------------------------------------------------------------------
    """
    if success_files:
        print("\n成功处理的文件:")
        for json_file in success_files:
            print(f"  - {json_file}")
    """

    if failed_files:
        print(f"\n失败文件列表 ({len(failed_files)} 个):")
        for failed in failed_files:
            json_file = failed.get('json_file', 'unknown')
            # 简单输出文件路径
            print(f"  - {json_file}")

    # 保存错误日志到文件（如果指定了--error-log）
    if args.error_log:
        error_log_path = args.error_log
        try:
            # 确保目录存在
            dir_name = os.path.dirname(error_log_path)
            if dir_name:
                os.makedirs(dir_name, exist_ok=True)

            # 保存失败文件路径列表（每行一个）
            with open(error_log_path, 'w', encoding='utf-8') as f:
                for failed in failed_files:
                    json_file = failed.get('json_file', '')
                    if json_file:
                        f.write(json_file + '\n')
            print(f"\n失败文件列表已保存到: {error_log_path}")

            # 同时保存JSON格式的详细信息（可选）
            json_path = error_log_path.replace('.txt', '.json') if error_log_path.endswith('.txt') else error_log_path + '.json'
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(failed_files, f, ensure_ascii=False, indent=2)
            print(f"详细错误日志已保存到: {json_path}")
        except Exception as e:
            print(f"保存错误日志失败: {e}")

    print("="*60)


if __name__ == '__main__':
    main()