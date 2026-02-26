import numpy as np
import matplotlib.pyplot as plt
import json
import os
import argparse

def weight2_slit_fourier(rho, H, dz, sigma):
    """
    计算傅里叶权重的辅助函数
    假设rho定义在z = [0.5*dz, 1.5*dz, ..., H-0.5*dz]上，共N=H/dz个点
    """
    N = len(rho)
    R = sigma / 2.0

    F = np.fft.fft(rho)
    w2 = np.pi * sigma * np.ones(N)

    # 计算半径对应的网格点数
    idx = int(np.round(R / dz))

    # 对于周期边界条件，设置w2的值
    if idx > 0 and idx + 1 < N - idx:
        w2[idx + 1 : N - idx] = 0
        w2[idx] = w2[idx] / 2.0
        w2[N - idx] = w2[N - idx] / 2.0

    A = np.fft.fft(w2)
    A[0] = 4 * np.pi * R**2 / dz

    n2 = np.fft.ifft(A * F) * dz
    return n2

def load_data_from_npz(npz_file, sigma=1.0):
    """
    从后处理的npz文件加载数据

    参数:
    npz_file: 后处理的npz文件路径
    sigma: 粒子直径（如果npz文件中没有sigma，则使用此值）

    返回:
    data_dict: 包含计算所需参数的字典
    """
    npz_data = np.load(npz_file)

    # 从npz获取基本参数
    rho_profile = npz_data['avg_density']
    H = float(npz_data['H'])
    dz = float(npz_data['dz'])

    # 从npz获取其他参数
    mu_b = float(npz_data.get('mu_b', 0.0))
    M = float(npz_data.get('M', 8.0))
    rho_b = float(npz_data.get('rho_b', 0.1))
    # 如果npz文件中有sigma，则使用文件中的值，否则使用传入的sigma参数
    sigma_file = npz_data.get('sigma')
    if sigma_file is not None:
        sigma = float(sigma_file)

    # 计算C参数
    C = 1.0 / (4.0 * np.pi * sigma**2)

    return {
        'rhob_seg': rho_b,
        'H': H,
        'dz': dz,
        'sigma': sigma,
        'C': C,
        'rho_profile': rho_profile,
        'mu_b': mu_b,
        'M': M
    }


def recal_simulation(npz_file, sigma=1.0, max_iters=400, mix_rate=0.02,
                   init_val=0.1, tol=1e-6):
    """
    运行主计算流程

    参数:
    npz_file: 输入的npz文件路径
    sigma: 粒子直径（默认为1.0，如果npz文件中有sigma则被覆盖）
    max_iters: 最大迭代次数
    mix_rate: 每次迭代新计算结果的混合比例
    init_val: G_new 的初始值
    tol: 提前终止的收敛容差 (Tolerance)
    """
    # 1. 读取数据
    data_dict = load_data_from_npz(npz_file, sigma)
    rhob_seg = data_dict['rhob_seg']
    H = data_dict['H']
    dz = data_dict['dz']
    sigma = data_dict['sigma']
    C = data_dict['C']
    rho_profile = data_dict['rho_profile']

    print(f"从NPZ文件加载数据:")
    print(f"  H={H}, dz={dz}, sigma={sigma}")
    print(f"  C={C}, rhob_seg={rhob_seg}")
    print(f"  rho_profile shape: {rho_profile.shape}")

    # 2. 动态获取数据维度
    data_length = len(rho_profile)
    G_new = np.ones(data_length) * init_val

    # 3. 带有提前终止的迭代计算
    error = float('inf')  # 初始化误差
    for i in range(max_iters):
        term_rho = (rho_profile / rhob_seg) * (G_new ** -1)
        G_temp = weight2_slit_fourier(term_rho, H, dz, sigma * 2) * C

        # 计算本次迭代的新值
        G_next = G_new * (1.0 - mix_rate) + G_temp * mix_rate

        # 计算更新前后的误差 (采用最大绝对值误差)
        error = np.max(np.abs(G_next - G_new))

        # 更新 G_new
        G_new = G_next

        # 判断是否满足提前终止条件
        if error < tol:
            print(f"计算收敛！在第 {i + 1} 次迭代时提前终止 (最大误差: {error:.2e} < {tol})")
            break

    else:
        # 如果 for 循环跑满了 max_iters 都没有被 break 打断，就会执行到这里
        print(f"警告: 已达到最大迭代次数 {max_iters}，可能尚未完全收敛 (最终误差: {error:.2e})")
        
    # 4. 计算最终的 exp_lambda_cal
    exp_effect_lambda_cal = (rho_profile / rhob_seg) / (G_new ** 2)

    return G_new, exp_effect_lambda_cal, rho_profile




if __name__ == "__main__":
    """
    主函数：使用后处理的NPZ文件
    """
    # 设置matplotlib后端为非交互式，避免显示窗口
    import matplotlib
    matplotlib.use('Agg')

    parser = argparse.ArgumentParser(description='重新计算模拟数据')
    parser.add_argument('--npz-file', type=str,
                        default='output/OUT_config_0000_M8_Trivial_H8.0_mu0.50/avg_density_config_0000_M8_Trivial_H8.0_mu0.50.npz',
                        help='后处理的NPZ文件路径（默认：默认示例文件）')
    parser.add_argument('--sigma', type=float, default=1.0,
                        help='粒子直径（如果NPZ文件中有sigma则被覆盖）')
    parser.add_argument('--max-iters', type=int, default=1000,
                        help='最大迭代次数')
    parser.add_argument('--mix-rate', type=float, default=0.02,
                        help='混合比例')
    parser.add_argument('--init-val', type=float, default=0.1,
                        help='G_new初始值')
    parser.add_argument('--tol', type=float, default=1e-6,
                        help='收敛容差')
    parser.add_argument('--output-plot', type=str, default=None,
                        help='绘图输出文件路径（默认：与NPZ文件同目录，文件名添加recal_plot_前缀）')
    parser.add_argument('--dpi', type=int, default=150,
                        help='输出图像分辨率（默认：150）')
    parser.add_argument('--no-plot', action='store_true',
                        help='不生成绘图')
    args = parser.parse_args()

    npz_file = args.npz_file
    if not os.path.exists(npz_file):
        print(f"NPZ文件不存在: {npz_file}")
        print("请先运行 postprocess_ring_density.py 生成后处理文件")
        exit(1)

    print(f"使用NPZ文件: {npz_file}")

    # 运行模拟
    G_new, exp_lambda_cal, rho_profile = recal_simulation(
        npz_file=npz_file,
        sigma=args.sigma,
        max_iters=args.max_iters,
        mix_rate=args.mix_rate,
        init_val=args.init_val,
        tol=args.tol
    )

    # ==========================================
    # 绘图部分 - 四个子图
    # ==========================================
    if not args.no_plot:
        # 确定输出文件路径
        if args.output_plot is None:
            # 自动生成：基于npz文件路径，生成recal_plot_{config_name}.png
            # 首先提取配置名（支持新旧格式）
            npz_basename = os.path.basename(npz_file)
            config_name = None

            # 尝试匹配新格式: avg_density_config_XXXX.npz
            if npz_basename.startswith('avg_density_') and npz_basename.endswith('.npz'):
                config_name = npz_basename[len('avg_density_'):-len('.npz')]
            # 尝试匹配旧格式: config_XXXX_avg_density.npz
            elif npz_basename.endswith('_avg_density.npz'):
                config_name = npz_basename[:-len('_avg_density.npz')]
            else:
                # 否则使用不带扩展名的文件名
                config_name = os.path.splitext(npz_basename)[0]

            # 生成输出文件名
            output_dir = os.path.dirname(npz_file)
            output_plot = os.path.join(output_dir, f'recal_plot_{config_name}.png')
        else:
            output_plot = args.output_plot

        plt.figure(figsize=(15, 10))

        grid_indices = np.arange(len(G_new))
        G_new_real = np.real(G_new)
        exp_lambda_real = np.real(exp_lambda_cal)

        # 子图1: G_new (实部)
        plt.subplot(2, 2, 1)
        plt.plot(grid_indices, G_new_real, '--', label='G_new (real)', color='blue')
        plt.xlabel('Grid index')
        plt.ylabel('G_new (real)')
        plt.title("G_new (real part)")
        plt.legend()
        plt.grid(True)

        # 子图2: exp_lambda_cal (实部)
        plt.subplot(2, 2, 2)
        plt.plot(grid_indices, exp_lambda_real, '-o', label='exp_lambda_cal (real)', color='green')
        plt.title("effective exp_lambda_cal (real part)")
        plt.xlabel('Grid index')
        plt.ylabel('exp_lambda_cal (real)')
        plt.legend()
        plt.grid(True)

        # 子图3: 密度分布 (rho_profile)
        plt.subplot(2, 2, 3)
        plt.plot(grid_indices, rho_profile, '-s', label='Density profile', color='red')
        plt.xlabel('Grid index')
        plt.ylabel('Density')
        plt.title("Density profile")
        plt.legend()
        plt.grid(True)

        # 子图4: 处理后的结果 - exp_lambda_cal 的对数图
        plt.subplot(2, 2, 4)
        # 只处理正数部分取对数，避免对数域错误
        exp_lambda_pos = np.where(exp_lambda_real > 0, exp_lambda_real, np.nan)
        plt.plot(grid_indices, np.log(exp_lambda_pos), '-^', label='log(exp_lambda_cal)', color='purple')
        plt.title("log(exp_lambda_cal)")
        plt.xlabel('Grid index')
        plt.ylabel('log(exp_lambda_cal)')
        plt.legend()
        plt.grid(True)

        plt.tight_layout()

        # 保存图像
        plt.savefig(output_plot, dpi=args.dpi)
        plt.close()
        print(f"绘图已保存至: {output_plot}")