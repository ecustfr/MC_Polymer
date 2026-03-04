import numpy as np
import matplotlib.pyplot as plt
import json
import os
import argparse
from . import utils

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



def recal_simulation(npz_file):
    """
    运行主计算流程

    参数:
    npz_file: 输入的npz文件路径
    sigma: 粒子直径（默认值，如果npz文件和json文件中都没有sigma，则使用此值）
    json_path: JSON配置文件路径（可选，用于从JSON读取sigma）

    硬编码参数:
    max_iters=400: 最大迭代次数
    mix_rate=0.02: 每次迭代新计算结果的混合比例
    init_val=0.1: G_new 的初始值
    tol=1e-6: 提前终止的收敛容差 (Tolerance)
    """
    # 1. 读取数据
    # data_dict = load_data_from_npz(npz_file)
    data_dict = utils.load_npz_data(npz_file)
    np.seterr(divide='ignore', invalid='ignore')
    # 硬编码参数
    max_iters = 1000      # 最大迭代次数
    mix_rate = 0.02       # 每次迭代新计算结果的混合比例
    init_val = 0.1        # G_new 的初始值
    tol = 1e-6            # 提前终止的收敛容差
    # rhob_seg = data_dict['rhob_seg']
    H = data_dict['H']
    dz = data_dict['dz']
    sigma = data_dict['sigma']
    rho_profile = data_dict['rho_profile']
    mu_ex_b = data_dict['mu_ex_b']
    vext = data_dict['vext_values']
    mu_ex_0 = data_dict['mu_ex_0']
    mu_b = data_dict['mu_b']-mu_ex_0
    M = data_dict["M"]
    

    rhob_seg = np.exp(mu_b-mu_ex_b)*M
    
    # print(f'rhob_seg:{rhob_seg},mu_ex_b:{mu_ex_b}')
    C = 1.0 / (4.0 * np.pi * sigma**2)


    # 2. 动态获取数据维度
    data_length = len(rho_profile)
   
    z = np.array([(i + 0.5) * dz for i in range(data_length)])
    G_new = np.ones(data_length) * init_val

    # 3. 带有提前终止的迭代计算
    cal_finish = False


    # for save
    config_name = utils.extract_config_name(npz_file)
    
    config_index = config_name[7:11]


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
            #print(f"计算收敛！在第 {i + 1} 次迭代时提前终止 (最大误差: {error:.2e} < {tol})")
            cal_finish = True
            break

    else:
        # 如果 for 循环跑满了 max_iters 都没有被 break 打断，就会执行到这里
        
        print(f"警告: 已达到最大迭代次数 {max_iters}，可能尚未完全收敛 (最终误差: {error:.2e})")
        

    # 
    if not cal_finish:
        return None
    
    if vext is not None:
        # custom_potential是数组，与rho_profile尺寸相同
        bz_vext = np.exp(-vext)
    else:
        bz_vext = 1.0
  
    mu_ex_cal_1 = np.log((rhob_seg/rho_profile  ) / (G_new ** 2) *bz_vext*np.exp(mu_ex_b/M))
    mu_ex_cal_2 = np.log((rhob_seg/rho_profile  ) * (G_new ** 2) *bz_vext*np.exp(mu_ex_b/M))
# *******************************文件保存部分*********************************************

    npz_dir = os.path.dirname(npz_file)

    output_dir = npz_dir
    recal_npz = os.path.join(output_dir,f'recal_result_{config_name}.npz')
    
    dt = np.dtype([('z','f8'),('G','f8'),('rho_profile','f8'),('mu_ex','f8'),('bz_vext','f8') ])
    G_new = np.real(G_new)
    rho_profile = np.real(rho_profile)
    mu_ex_cal_1 = np.real(mu_ex_cal_1)
    mu_ex_cal_2 = np.real(mu_ex_cal_2)
    bz_vext = np.real(bz_vext)

    
    ml_data = np.zeros(data_length,dtype=dt)


    ml_data['rho_profile'] = rho_profile
    ml_data['bz_vext'] = bz_vext
    ml_data['mu_ex'] = mu_ex_cal_2
    ml_data['G'] = G_new
    ml_data['z'] = z 

    npy_path = os.path.join(output_dir,'simData'+config_index+'.npy')
    # save_data = {'G':G_new,'rho_profile':rho_profile,'mu_ex_1':mu_ex_cal_1,'mu_ex_2':mu_ex_cal_2,'bz_vext':bz_vext}

    np.save(npy_path , ml_data)    
    # utils.save_npz_data(recal_npz,save_data)

#   csv文件 
    csv_output = os.path.join(output_dir, f'recal_results_{config_name}.csv')
    n_points = len(G_new)
    grid_indices = np.arange(n_points)
    csv_data = np.column_stack([
            grid_indices,
            G_new,
            mu_ex_cal_1,
            mu_ex_cal_2,
            rho_profile,
            bz_vext
            ])
    
    header = "grid_index,G,mu_ex_1,mu_ex_2,rho_profile,bz_Vext"
    np.savetxt(csv_output, csv_data, delimiter=',', header=header, comments='')

    return {'sim_index':config_index, 'path':npy_path}  



     


if __name__ == "__main__":
    pass