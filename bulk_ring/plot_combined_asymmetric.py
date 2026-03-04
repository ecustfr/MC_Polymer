#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
合并图形（使用块统计信息）：
1. 左图：平均密度 vs 化学势（带对称误差条）
2. 右图：偏差 (mu_b - mu_id) = mu_ex 与理论过剩化学势比较（带对称误差条）

使用块统计信息（block_mean 和 block_se）从 bulk_density_results_mu_*.txt 文件。
假设 M = 8（每个聚合物的单体数）。
"""
import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端，避免显示图形
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 尝试设置控制台编码以避免GBK问题
try:
    import locale
    locale.setlocale(locale.LC_ALL, '')
except:
    pass

# 尝试设置合适的编码用于标准输出
if sys.platform.startswith('win'):
    try:
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    except:
        pass

def parse_result_file(filepath):
    """
    解析单个结果文件，提取相关统计信息。
    返回字典，包含 mu_b, block_mean_density, block_se_density,
           block_mean_mu_ex, block_se_mu_ex, block_mean_mu_id, block_se_mu_id 等。
    """
    data = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') and ':' in line:
                # 移除开头的 '# '
                key_value = line[2:] if line.startswith('# ') else line[1:]
                key, value = key_value.split(':', 1)
                key = key.strip()
                value = value.strip()

                # 尝试转换为浮点数
                try:
                    if key in ['mu_b', 'ave_density', 'std_density', 'block_mean_density',
                              'block_se_density', 'block_std_density', 'mu_ex', 'mu_id',
                              'block_mean_mu_ex', 'block_se_mu_ex', 'block_std_mu_ex',
                              'block_mean_mu_id', 'block_se_mu_id', 'block_std_mu_id']:
                        data[key] = float(value)
                    elif key == 'block_mu_ex_values' or key == 'block_mu_id_values':
                        # 解析多个值
                        values = [float(x) for x in value.split()]
                        data[key] = values
                except ValueError:
                    # 如果转换失败，保持字符串格式
                    data[key] = value
    return data

def read_block_statistics(bulk_dir):
    """
    读取 bulk_ring 目录中的所有 bulk_density_results_mu_*.txt 文件，
    提取块统计信息。

    返回排序后的数组：
        mu_b, block_mean_density, block_std_density, block_se_density,
        block_mean_mu_ex, block_std_mu_ex, block_se_mu_ex,
        block_mean_mu_id, block_std_mu_id, block_se_mu_id
    """
    pattern = os.path.join(bulk_dir, 'bulk_density_results_mu_*.txt')
    files = glob.glob(pattern)

    if not files:
        raise ValueError(f"未找到匹配的文件: {pattern}")

    print(f"找到 {len(files)} 个结果文件")

    data_list = []
    for filepath in files:
        try:
            data = parse_result_file(filepath)
            if 'mu_b' in data:
                data_list.append(data)
                print(f"  解析: {os.path.basename(filepath)} (μ_b={data['mu_b']})")
        except Exception as e:
            print(f"  警告: 解析文件 {filepath} 时出错: {e}")

    if not data_list:
        raise ValueError("没有成功解析任何数据文件")

    # 按 mu_b 排序
    data_list.sort(key=lambda x: x['mu_b'])

    # 提取数组
    mu_b = np.array([d['mu_b'] for d in data_list])
    block_mean_density = np.array([d.get('block_mean_density', d.get('ave_density', 0)) for d in data_list])
    block_std_density = np.array([d.get('block_std_density', d.get('std_density', 0)) for d in data_list])
    block_se_density = np.array([d.get('block_se_density', 0) for d in data_list])
    block_mean_mu_ex = np.array([d.get('block_mean_mu_ex', d.get('mu_ex', 0)) for d in data_list])
    block_std_mu_ex = np.array([d.get('block_std_mu_ex', 0) for d in data_list])
    block_se_mu_ex = np.array([d.get('block_se_mu_ex', 0) for d in data_list])
    block_mean_mu_id = np.array([d.get('block_mean_mu_id', d.get('mu_id', 0)) for d in data_list])
    block_std_mu_id = np.array([d.get('block_std_mu_id', 0) for d in data_list])
    block_se_mu_id = np.array([d.get('block_se_mu_id', 0) for d in data_list])

    return mu_b, block_mean_density, block_std_density, block_se_density, \
           block_mean_mu_ex, block_std_mu_ex, block_se_mu_ex, \
           block_mean_mu_id, block_std_mu_id, block_se_mu_id

def main():
    try:
        # 确定 bulk_ring 目录路径
        polymer_type_ad = 'M8'

        script_dir = os.path.dirname(os.path.abspath(__file__))
        bulk_dir = os.path.join(script_dir, polymer_type_ad)

        # 读取块统计信息
        print("正在读取块统计信息...")
        mu_b, block_mean_density, block_std_density, block_se_density, \
        block_mean_mu_ex, block_std_mu_ex, block_se_mu_ex, \
        block_mean_mu_id, block_std_mu_id, block_se_mu_id = read_block_statistics(bulk_dir)

        # 设置参数
        M = 8.0  # 每个聚合物的单体数（硬编码，与C++代码一致）

        # 检查数据有效性
        if len(mu_b) == 0:
            print("错误: 没有读取到有效数据")
            return

        print(f"成功读取 {len(mu_b)} 个数据点")
        print(f"化学势范围: {mu_b[0]:.2f} 到 {mu_b[-1]:.2f}")

        # ====================================================================
        # 计算 derived quantities
        # ====================================================================

        # 计算 mu_b - mu_id（用于过剩化学势比较）
        calculated_excess = mu_b - block_mean_mu_id

        # 计算总化学势 mu_total = mu_id + mu_ex
        mu_total = block_mean_mu_id + block_mean_mu_ex

        # ====================================================================
        # 計算 derived quantities
        # ====================================================================

        # 計算 mu_b - mu_id（用於過剩化學勢比較）
        calculated_excess = mu_b - block_mean_mu_id

        # 計算總化學勢 mu_total = mu_id + mu_ex
        mu_total = block_mean_mu_id + block_mean_mu_ex

        # ====================================================================
        # 保存圖片中使用的數據到 .npz 和 .csv
        # ====================================================================
        print("正在保存繪圖數據到 .npz 和 .csv 檔案...")
        
        # 1. 保存為 .npz
        npz_output = os.path.join(bulk_dir,  polymer_type_ad+'_mu_ex.npz')
        np.savez(npz_output, 
                 mu_b=mu_b, 
                 block_mean_density=block_mean_density, 
                 block_std_density=block_std_density,
                 calculated_excess=calculated_excess,
                 block_std_mu_id=block_std_mu_id,
                 block_mean_mu_ex=block_mean_mu_ex,
                 block_std_mu_ex=block_std_mu_ex)
        print(f"  .npz 檔案已保存至: {npz_output}")

        # 2. 保存為 .csv
        csv_output = os.path.join(bulk_dir, polymer_type_ad+'_mu_ex.csv')
        # 將一維陣列組合成二維陣列的列 (columns)
        csv_data = np.column_stack((
            mu_b, 
            block_mean_density, 
            block_std_density,
            calculated_excess,
            block_std_mu_id,
            block_mean_mu_ex,
            block_std_mu_ex
        ))
        # 設定 CSV 標頭
        csv_header = "mu_b,block_mean_density,block_std_density,calculated_excess,block_std_mu_id,block_mean_mu_ex,block_std_mu_ex"
        np.savetxt(csv_output, csv_data, delimiter=',', header=csv_header, comments='', fmt='%.6f')
        print(f"  .csv 檔案已保存至: {csv_output}")

        # ====================================================================
        # 創建合併圖形（1行2列）
        # ====================================================================
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # ====================================================================
        # 创建合并图形（1行2列）
        # ====================================================================
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # ====================================================================
        # 左图：平均密度 vs 化学势（带对称误差条）
        # ====================================================================
        ax1.errorbar(mu_b, block_mean_density, yerr=block_std_density, fmt='o-',
                    capsize=4, capthick=1, elinewidth=1,
                    linewidth=2, markersize=6, color='tab:blue')

        ax1.set_xlabel('Chemical Potential (μ_b)', fontsize=12)
        ax1.set_ylabel('Average Density', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.set_title('Average Density vs Chemical Potential (Block Statistics, ±σ)', fontsize=14)

        # ====================================================================
        # 右图：过剩化学势对比（带对称误差条）
        # ====================================================================
        # 绘制 calculated μ_b - μ_id（使用块统计的 μ_id）
        # 误差棒使用 block_std_mu_id（因为 μ_b 是精确值）
        ax2.errorbar(mu_b, calculated_excess, yerr=block_std_mu_id, fmt='o-',
                     capsize=4, capthick=1, elinewidth=1,
                     linewidth=2, markersize=6, color='tab:blue',
                     label='Calculated μ_b - μ_id (with block std error bars)')

        # 绘制直接计算的过剩化学势 μ_ex（来自块统计）
        ax2.errorbar(mu_b, block_mean_mu_ex, yerr=block_std_mu_ex, fmt='s--',
                     capsize=4, capthick=1, elinewidth=1,
                     linewidth=2, markersize=6, color='tab:red',
                     label='Directly Calculated μ_ex (with block std error bars)')

        ax2.set_xlabel('Chemical Potential (μ_b)', fontsize=12)
        ax2.set_ylabel('Excess Chemical Potential', fontsize=12)
        # ax2.set_ylim(4.0,7.0)
        ax2.grid(True, alpha=0.3)
        ax2.set_title('Excess Chemical Potential Comparison (Block Statistics, ±σ)', fontsize=14)
        ax2.legend()

        plt.tight_layout()

        # 保存合并图形
        output_plot = os.path.join(bulk_dir, 'combined_plots_block_stats.png')
        plt.savefig(output_plot, dpi=300, bbox_inches='tight')
        print(f"合并图形（块统计）已保存到: {output_plot}")

        # ====================================================================
        # 输出统计信息
        # ====================================================================
        print("\n=== 块统计信息 ===")
        print(f"化学势范围: {mu_b[0]:.2f} 到 {mu_b[-1]:.2f}")
        print(f"平均密度范围: {block_mean_density[0]:.4f} 到 {block_mean_density[-1]:.4f}")
        print(f"过剩化学势范围: {block_mean_mu_ex[0]:.4f} 到 {block_mean_mu_ex[-1]:.4f}")
        print(f"理想气体化学势范围: {block_mean_mu_id[0]:.4f} 到 {block_mean_mu_id[-1]:.4f}")

        # 计算块标准差统计
        print(f"\n块标准差统计 (block_std):")
        print(f"  密度块标准差相对值 (block_std/block_mean):")
        for i in range(len(mu_b)):
            if block_mean_density[i] > 0:
                rel_err = block_std_density[i] / block_mean_density[i] * 100
                print(f"    μ_b={mu_b[i]:.2f}: {rel_err:.2f}%")

        print(f"\n  μ_id 块标准差绝对大小:")
        for i in range(len(mu_b)):
            print(f"    μ_b={mu_b[i]:.2f}: ±{block_std_mu_id[i]:.4f}")

        print(f"\n  μ_ex 块标准差绝对大小:")
        for i in range(len(mu_b)):
            print(f"    μ_b={mu_b[i]:.2f}: ±{block_std_mu_ex[i]:.4f}")

        # 计算化学势计算偏差
        deviation = mu_b - mu_total
        print(f"\n化学势计算偏差统计 (μ_b - (μ_id + μ_ex)):")
        print(f"  平均偏差: {np.mean(deviation):.6f}")
        print(f"  标准差: {np.std(deviation):.6f}")
        print(f"  最大绝对偏差: {np.max(np.abs(deviation)):.6f}")

        # 计算两种过剩化学势的差异
        diff_excess = calculated_excess - block_mean_mu_ex
        print(f"\n过剩化学势差异统计 (calculated μ_b - μ_id vs direct μ_ex):")
        print(f"  平均差异: {np.mean(diff_excess):.6f}")
        print(f"  标准差: {np.std(diff_excess):.6f}")
        print(f"  最大绝对差异: {np.max(np.abs(diff_excess)):.6f}")

        # 检查差异是否在误差范围内
        within_error = []
        for i in range(len(mu_b)):
            error_sum = np.sqrt(block_std_mu_id[i]**2 + block_std_mu_ex[i]**2)
            if abs(diff_excess[i]) <= 2 * error_sum:  # 2 sigma
                within_error.append(True)
            else:
                within_error.append(False)

        if within_error:
            within_count = sum(within_error)
            total_count = len(within_error)
            print(f"\n差异在2倍标准差范围内的数据点: {within_count}/{total_count} ({within_count/total_count*100:.1f}%)")

    except Exception as e:
        print(f"运行脚本时发生错误: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())