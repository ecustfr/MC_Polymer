#!/usr/bin/env python3
"""
处理output文件夹中内容并再计算的整合功能。

利用postprocess_ring_density.py和recal_equ.py两个库，自动化完成以下流程：
1. 扫描input/Ring_configs目录下的所有JSON配置文件。
2. 对每个JSON文件，运行后处理生成平均密度分布NPZ文件（如果尚未生成）。
3. 对每个NPZ文件，运行recal_equ.py中的重新计算。
4. 保存重新计算结果，并生成图表。

使用方法：
    python process_outputs.py [--plot] [--overwrite]

选项：
    --plot: 生成并保存图表（PNG格式）
    --overwrite: 即使NPZ文件已存在也重新生成后处理文件
"""

import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 导入现有库的功能
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from postprocess_ring_density import process_json, save_result
except ImportError as e:
    print(f"错误: 无法导入postprocess_ring_density: {e}")
    sys.exit(1)

try:
    from recal_equ import recal_simulation
except ImportError as e:
    print(f"错误: 无法导入recal_equ: {e}")
    sys.exit(1)


def get_npz_file_from_json(json_path):
    """
    根据JSON文件路径获取对应的npz文件路径
    假设npz文件位于JSON中指定的output_dir目录中

    参数:
        json_path: JSON配置文件路径

    返回:
        npz_file: NPZ文件路径
    """
    # 读取JSON文件获取output_dir
    with open(json_path, 'r') as f:
        config = json.load(f)

    # 获取output_dir
    output_dir = config['output_params']['output_dir']

    # 确保output_dir是绝对路径（相对于项目根目录）
    if not os.path.isabs(output_dir):
        # 假设JSON文件路径相对于项目根目录，计算绝对路径
        # 项目根目录是JSON文件父目录的父目录的父目录
        project_root = Path(json_path).parent.parent.parent
        output_dir = os.path.join(project_root, output_dir)

    # 从JSON文件名提取配置名
    config_name = os.path.splitext(os.path.basename(json_path))[0]

    # 构建两种可能的npz文件路径（新格式和旧格式）
    npz_file_new = os.path.join(output_dir, f"avg_density_{config_name}.npz")
    npz_file_old = os.path.join(output_dir, f"{config_name}_avg_density.npz")

    # 检查文件是否存在（优先新格式）
    if os.path.exists(npz_file_new):
        return npz_file_new
    elif os.path.exists(npz_file_old):
        return npz_file_old
    else:
        # 如果都不存在，返回新格式路径（用于创建新文件）
        return npz_file_new


def ensure_npz_file(json_path, overwrite=False):
    """
    确保JSON文件对应的NPZ文件存在。

    参数:
        json_path: JSON配置文件路径
        overwrite: 如果为True，即使NPZ文件已存在也重新生成

    返回:
        npz_file: NPZ文件路径，如果生成失败则返回None
    """
    # 获取预期的NPZ文件路径
    npz_file = get_npz_file_from_json(json_path)

    # 检查NPZ文件是否存在且overwrite为False
    if os.path.exists(npz_file) and not overwrite:
        print(f"NPZ文件已存在: {npz_file}")
        return npz_file

    print(f"生成NPZ文件: {npz_file}")

    # 运行后处理
    result = process_json(json_path)
    if result is None:
        print(f"警告: 处理JSON文件失败: {json_path}")
        return None

    # 确保输出目录存在
    output_dir = result['output_dir']
    os.makedirs(output_dir, exist_ok=True)

    # 保存结果
    try:
        save_result(result, npz_file)
        print(f"成功生成NPZ文件: {npz_file}")
        return npz_file
    except Exception as e:
        print(f"保存NPZ文件时出错: {e}")
        return None


def run_recalculation(npz_file, sigma=1.0, **kwargs):
    """
    运行重新计算。

    参数:
        npz_file: NPZ文件路径
        sigma: 粒子直径
        **kwargs: 传递给recal_simulation的其他参数

    返回:
        result_dict: 包含计算结果和元数据的字典
    """
    print(f"运行重新计算: {npz_file}")

    try:
        G_new, exp_lambda_cal, rho_profile = recal_simulation(
            npz_file=npz_file,
            sigma=sigma,
            **kwargs
        )

        # 从NPZ文件名提取配置名（支持新旧格式）
        base_name = os.path.basename(npz_file)

        # 尝试匹配新格式: avg_density_config_XXXX.npz
        if base_name.startswith('avg_density_') and base_name.endswith('.npz'):
            config_name = base_name[len('avg_density_'):-len('.npz')]
        # 尝试匹配旧格式: config_XXXX_avg_density.npz
        elif base_name.endswith('_avg_density.npz'):
            config_name = base_name[:-len('_avg_density.npz')]
        else:
            config_name = os.path.splitext(base_name)[0]

        return {
            'config_name': config_name,
            'G_new': G_new,
            'exp_lambda_cal': exp_lambda_cal,
            'rho_profile': rho_profile,
            'npz_file': npz_file
        }
    except Exception as e:
        print(f"重新计算失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def save_recal_results(result, output_dir=None):
    """
    保存重新计算结果。

    参数:
        result: 包含计算结果的字典
        output_dir: 输出目录（如果为None，则使用NPZ文件所在目录）
    """
    config_name = result['config_name']

    if output_dir is None:
        output_dir = os.path.dirname(result['npz_file'])

    os.makedirs(output_dir, exist_ok=True)

    # 保存为NPZ文件（只保存实部）
    npz_output = os.path.join(output_dir, f'recal_results_{config_name}.npz')
    np.savez(npz_output,
             G_new_real=np.real(result['G_new']),
             exp_lambda_cal_real=np.real(result['exp_lambda_cal']),
             rho_profile=result['rho_profile'],
             config_name=result['config_name'])

    print(f"重新计算结果保存到: {npz_output}")

    # 同时保存为文本文件（CSV格式）以便于查看
    txt_output = os.path.join(output_dir, f'recal_results_{config_name}.csv')
    n_points = len(result['G_new'])
    grid_indices = np.arange(n_points)

    data = np.column_stack([
        grid_indices,
        np.real(result['G_new']),
        np.real(result['exp_lambda_cal']),
        result['rho_profile']
    ])

    header = "grid_index,G_new_real,exp_lambda_cal_real,rho_profile"
    np.savetxt(txt_output, data, delimiter=',', header=header, comments='')

    print(f"文本格式结果保存到: {txt_output}")

    return npz_output, txt_output


def plot_recal_results(result, save_dir=None):
    """
    绘制重新计算结果并保存为PNG文件。

    参数:
        result: 包含计算结果的字典
        save_dir: 保存图表的目录（如果为None，则使用NPZ文件所在目录）
    """
    config_name = result['config_name']
    G_new = result['G_new']
    exp_lambda_cal = result['exp_lambda_cal']
    rho_profile = result['rho_profile']

    if save_dir is None:
        save_dir = os.path.dirname(result['npz_file'])

    os.makedirs(save_dir, exist_ok=True)

    fig = plt.figure(figsize=(15, 10))

    # 子图1: G_new (实部)
    ax1 = plt.subplot(2, 2, 1)
    grid_indices = np.arange(len(G_new))
    G_new_real = np.real(G_new)
    ax1.plot(grid_indices, G_new_real, '--', label='G_new (real)', color='blue')
    ax1.set_xlabel('Grid index')
    ax1.set_ylabel('G_new (real)')
    ax1.set_title(f"{config_name}\nG_new (real part)")
    ax1.legend()
    ax1.grid(True)

    # 子图2: exp_lambda_cal (实部)
    ax2 = plt.subplot(2, 2, 2)
    exp_lambda_real = np.real(exp_lambda_cal)
    ax2.plot(grid_indices, exp_lambda_real, '-', label='exp_lambda_cal (real)', color='green')
    ax2.set_title(f"{config_name}\neffective exp_lambda_cal (real part)")
    ax2.set_xlabel('Grid index')
    ax2.set_ylabel('exp_lambda_cal (real)')
    ax2.legend()
    ax2.grid(True)

    # 子图3: 密度分布 (rho_profile)
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(grid_indices, rho_profile, '-', label='Density profile', color='red')
    ax3.set_xlabel('Grid index')
    ax3.set_ylabel('Density')
    ax3.set_title(f"{config_name}\nDensity profile")
    ax3.legend()
    ax3.grid(True)

    # 子图4: 处理后的结果 - exp_lambda_cal 的对数图
    ax4 = plt.subplot(2, 2, 4)
    # 只处理正数部分取对数，避免对数域错误
    exp_lambda_pos = np.where(exp_lambda_real > 0, exp_lambda_real, np.nan)
    ax4.plot(grid_indices, np.log(exp_lambda_pos), '-', label='log(exp_lambda_cal)', color='purple')
    ax4.set_title(f"{config_name}\nlog(exp_lambda_cal)")
    ax4.set_xlabel('Grid index')
    ax4.set_ylabel('log(exp_lambda_cal)')
    ax4.legend()
    ax4.grid(True)

    plt.tight_layout()

    # 保存图表
    plot_file = os.path.join(save_dir, f'recal_plot_{config_name}.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"图表保存到: {plot_file}")


def process_all_json_files(input_dir, plot=True, overwrite=False, sigma=1.0, **kwargs):
    """
    处理指定目录中的所有JSON文件。

    参数:
        input_dir: 包含JSON文件的目录
        plot: 是否生成图表
        overwrite: 是否覆盖已存在的NPZ文件
        sigma: 粒子直径
        **kwargs: 传递给recal_simulation的其他参数
    """
    # 查找所有JSON文件
    json_pattern = os.path.join(input_dir, '*.json')
    import glob
    json_files = glob.glob(json_pattern)

    if not json_files:
        print(f"未找到JSON文件: {json_pattern}")
        return []

    print(f"找到 {len(json_files)} 个JSON文件")

    results = []
    for json_file in json_files:
        print("\n" + "="*60)
        print(f"处理文件: {json_file}")
        print("="*60)

        # 步骤1: 确保NPZ文件存在
        npz_file = ensure_npz_file(json_file, overwrite=overwrite)
        if npz_file is None:
            print(f"跳过 {json_file}，无法生成NPZ文件")
            continue

        # 步骤2: 运行重新计算
        result = run_recalculation(npz_file, sigma=sigma, **kwargs)
        if result is None:
            print(f"跳过 {json_file}，重新计算失败")
            continue

        # 步骤3: 保存结果
        npz_out, txt_out = save_recal_results(result)
        result['output_npz'] = npz_out
        result['output_txt'] = txt_out

        # 步骤4: 生成图表（如果需要）
        if plot:
            plot_recal_results(result)

        results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(description='处理output文件夹中内容并再计算的整合功能')
    parser.add_argument('--input-dir', default='input/Ring_configs',
                       help='JSON配置文件目录 (默认: input/Ring_configs)')
    parser.add_argument('--plot', action='store_true',default=True,
                       help='生成并保存图表（PNG格式）')
    parser.add_argument('--overwrite', action='store_true',
                       help='即使NPZ文件已存在也重新生成后处理文件')
    parser.add_argument('--sigma', type=float, default=1.0,
                       help='粒子直径 (默认: 1.0)')
    parser.add_argument('--max-iters', type=int, default=1000,
                       help='最大迭代次数 (默认: 1000)')
    parser.add_argument('--mix-rate', type=float, default=0.02,
                       help='混合比例 (默认: 0.02)')
    parser.add_argument('--init-val', type=float, default=0.1,
                       help='G_new初始值 (默认: 0.1)')
    parser.add_argument('--tol', type=float, default=1e-6,
                       help='收敛容差 (默认: 1e-6)')

    args = parser.parse_args()

    # 检查输入目录是否存在
    if not os.path.exists(args.input_dir):
        print(f"错误: 输入目录不存在: {args.input_dir}")
        sys.exit(1)

    # 运行处理流程
    results = process_all_json_files(
        input_dir=args.input_dir,
        plot=args.plot,
        overwrite=args.overwrite,
        sigma=args.sigma,
        max_iters=args.max_iters,
        mix_rate=args.mix_rate,
        init_val=args.init_val,
        tol=args.tol
    )

    print("\n" + "="*60)
    print("处理完成!")
    print(f"成功处理 {len(results)} 个配置文件")

    if results:
        print("\n处理结果:")
        for result in results:
            print(f"  - {result['config_name']}:")
            print(f"     结果文件: {result['output_npz']}")
            print(f"     文本文件: {result['output_txt']}")
            if args.plot:
                plot_file = os.path.join(os.path.dirname(result['npz_file']),
                                       f"recal_plot_{result['config_name']}.png")
                print(f"     图表文件: {plot_file}")

    print("="*60)


if __name__ == '__main__':
    main()