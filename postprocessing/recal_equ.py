import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm

import json
import os
import argparse
from pathlib import Path
from . import utils




try :
    from . import utils
except ImportError:
    import utils
try:
    import scipy.signal
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


    # 1. 设置全局字体为支持中文的字体，'SimHei' 是黑体，'Microsoft YaHei' 是微软雅黑
plt.rcParams['font.sans-serif'] = ['SimHei'] 

# 2. 解决更改中文字体后，坐标轴负号 '-' 变成方块的问题（这也是个必踩的坑）
plt.rcParams['axes.unicode_minus'] = False
PADDING = True
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


def _get_alpha(iteration: int, initial_alpha: float = 0.00001) -> float:

    alpha = initial_alpha
    if iteration > 10:
        alpha = 0.001
    if iteration > 20:
        alpha = 0.01
    if iteration > 50:
        alpha = 0.02
    if iteration > 100:
        alpha = 0.05
    if iteration > 300:
        alpha = 0.02
    if iteration > 5000:
        alpha = 0.01
    return alpha

def _check_convergence(rho_profile: np.ndarray, rho_temp_profile: np.ndarray,
                       max_delta: float = 1e-3) -> tuple[bool, float, float]:

    delta = np.max(np.abs(rho_temp_profile - rho_profile))
    relative_error = delta / np.max(rho_profile)
    converged = delta < max_delta or relative_error < max_delta
    return converged, delta, relative_error


def reflect_aug(rho_profile_mat):
    # 线性高分子 对称增强
    rho_profile_mat_flipped = np.fliplr(rho_profile_mat)
    rho_profile_augmented = (rho_profile_mat + rho_profile_mat_flipped)/2
    return rho_profile_augmented


def _get_padding_mu_ex(output_dir):
    mu_ex = utils._get_block_average(output_dir=output_dir,pattern="block_*_mu_ex_wall.dat") 
    return mu_ex

def filter_signal(data, cutoff, fs, order=4):
    """
    对信号进行低通滤波（Butterworth滤波器，零相位）

    参数:
    data: 输入信号数组
    cutoff: 截止频率（与fs相同单位）
    fs: 采样频率
    order: 滤波器阶数（默认4）

    返回:
    filtered_data: 滤波后的信号
    """
    data = np.squeeze(data)
    if not HAS_SCIPY:
        # 回退到简单移动平均（窗口大小由截止频率推算）
        window_size = int(fs / cutoff)
        if window_size % 2 == 0:
            window_size += 1
        if window_size < 3:
            window_size = 3
        window = np.ones(window_size) / window_size
        filtered = np.convolve(data, window, mode='same')
        return filtered

    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    if normal_cutoff >= 1.0:
        normal_cutoff = 0.99
    elif normal_cutoff <= 0.0:
        normal_cutoff = 0.01

    b, a = scipy.signal.butter(order, normal_cutoff, btype='low', analog=False)
    filtered = scipy.signal.filtfilt(b, a, data)
    return filtered

def recal_simulation(npz_file):
    data_dict = utils.load_npz_data(npz_file)
    PADDING = False
    np.seterr(divide='ignore', invalid='ignore')
    polymer_type = data_dict['polymer_type']
    # print(polymer_type)
    if polymer_type=="Ring":
        result = recal_simulation_ring(npz_file)
    elif polymer_type == "Linear":
        if PADDING:
            result = recal_simulation_linear_padding( npz_file , mu_ex_0 = 0)
        else:
            result = recal_simulation_linear(npz_file)
    else:
        result = None
        print(polymer_type+"polymer is not included in postprocessing.")
    return result

def recal_simulation_linear(npz_file):
    import numpy as np
    data_dict = utils.load_npz_data(npz_file)
    np.seterr(divide='ignore', invalid='ignore')
    max_iters = 20000      # 最大迭代次数
    mix_rate = 0.02       # 每次迭代新计算结果的混合比例
    tol = 1e-5            # 提前终止的收敛容差

    H = data_dict['H']
    dz = data_dict['dz'].astype(float)
    sigma = data_dict['sigma']
    
    vext = data_dict['vext_values']

    vext = np.expand_dims(vext,axis=1)

    # mu_ex_0 = data_dict['mu_ex_0']
    mu_ex_0 =  3
    mu_b = data_dict['mu_b']-mu_ex_0
    M = data_dict["M"].astype(int)
    rho_profile_mat = data_dict['rho_profile'] #*dz 由于之前代码的失误


    C = 1.0 / (4.0 * np.pi * sigma**2)

    

    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]

    
    data_length = rho_profile_mat.shape[0]
    rho_profile_mat = reflect_aug(rho_profile_mat)


    z = (np.arange(data_length) + 0.5) * dz
    G = np.ones_like(rho_profile_mat) 
    mu_ex_loc_mat = np.zeros(shape=(data_length,M))
    lambda_mat = mu_ex_loc_mat  + vext

    
    W = np.exp(-lambda_mat)
    W_new = np.copy(W)
    exp_mub = np.exp(mu_b)

    z_in = np.where( (z>sigma/2-dz/2) & ( z < H-sigma/2 ) )[0]

    # print(f"H-sigma/2:{H-sigma/2},z_end:{z[z_in[-1]]}")
    if np.abs(z[z_in[-1]]-(H-sigma/2))<1e-5:
        print("OK")
        rho_profile_mat[z_in[-1],:] *= 2  
    
    mu_ex_loc_avg = np.zeros(data_length)

    mu_ex_loc_avg_temp = np.zeros_like(mu_ex_loc_avg)     
    mu_ex_loc_avg_temp[z_in] = 1
    i = 0
    while True:
        i = i + 1
        alpha = _get_alpha(i)
        G[:,0] = 1
        for mono in range(1,M):
            G[:,mono] =np.real(weight2_slit_fourier( G[:,mono-1]*W[:,mono-1],H,dz,2*sigma)*C)

        W_new = rho_profile_mat/(G*G[:,::-1])/exp_mub
                
        W = W_new*alpha + W*(1-alpha)

        mu_ex_loc_mat[z_in,:] = -np.log(W[z_in,:])-vext[z_in,:]
        
        if i>100:
            mu_ex_loc_avg = np.mean(mu_ex_loc_mat, axis=1)
            
            converged,delta,_ = _check_convergence( mu_ex_loc_avg , mu_ex_loc_avg_temp ,tol)

            mu_ex_loc_avg_temp = mu_ex_loc_avg

            # print(f"delta:{delta}")

            if converged:
                print(f"Converged after {i} iterations (delta = {delta:.6e})")

                break

            if i>max_iters:
                print(f"iter time reach max_iters, Not converged")
                break

    
    
    mu_ex_loc_mat[z_in,:] = -np.log(W[z_in,:])-vext[z_in,:]
    renormal_C = np.zeros(M)
    for mono in range(M):
        renormal_C[mono] = np.mean(mu_ex_loc_mat[z_in,mono].squeeze() - mu_ex_loc_avg[z_in] )
        mu_ex_loc_mat[z_in,mono] -= renormal_C[mono] 
    
    #print(renormal_C)
    #print(np.sum(renormal_C))
    
    G = np.real(G)
    bz_vext = np.exp(-vext).squeeze()
    mu_b = np.ones_like(z)*mu_b

    
    mu_ex_filtered = np.zeros_like(mu_ex_loc_mat)
    
    fs = 1.0 / dz
    cutoff = 2

    # valid = ~np.isnan(mu_ex)
    
    for mono in range(M):
        mu_ex_filtered[z_in,mono] = filter_signal(mu_ex_loc_mat[z_in,mono], cutoff, fs, order=4).reshape(-1)
    


    dt = np.dtype([('z','f8'),
            ('G','f8',(M,)),
            ('rho_profile_mat','f8',(M,)),
            ('mu_ex','f8',(M,)),
            ('mu_ex_f','f8',(M,)),
            ('bz_vext','f8'),
            ('mu_b','f8') ])
    
    ml_data = np.zeros(data_length,dtype=dt)
    ml_data['rho_profile_mat'] = rho_profile_mat
    ml_data['bz_vext'] = bz_vext
    ml_data['mu_ex']   = mu_ex_loc_mat
    ml_data['mu_ex_f'] = mu_ex_filtered
    ml_data['G'] = G
    ml_data['z'] = z 
    ml_data['mu_b'] = mu_b
    
    

    
    npz_dir = os.path.dirname(npz_file)
    output_dir = npz_dir
    
    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]
    
    recal_npz = os.path.join(output_dir,f'recal_result_{config_name}.npy')
    np.save(recal_npz , ml_data)    

    plot_polymer_profiles(
        z = z, 
        rho_profile=rho_profile_mat, 
        mu_ex = mu_ex_loc_mat, 
        mu_ex_f= mu_ex_filtered, 
        output_dir = output_dir, 
        filename = f"{config_name}_linear_profiles.png",
        title_suffix = f"(Linear, M={M})"
        )
    
    return {'sim_index':config_index, 'path':recal_npz}  


def recal_simulation_ring(npz_file):
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
  
    #mu_ex_cal_1 = np.log((rhob_seg/rho_profile  ) / (G_new ** 2) *bz_vext*np.exp(mu_ex_b/M))
    mu_ex = np.log((rhob_seg/rho_profile  ) * (G_new ** 2) *bz_vext*np.exp(mu_ex_b/M))
# *******************************文件保存部分*********************************************

    npz_dir = os.path.dirname(npz_file)

    output_dir = npz_dir
    recal_npz = os.path.join(output_dir,f'recal_result_{config_name}.npz')
    
    dt = np.dtype([('z','f8'),('G','f8'),('rho_profile','f8'),('mu_ex','f8'),('bz_eff_vext','f8'),('rhob','f8'),('mu_ex_f','f8') ])
    G_new = np.real(G_new)
    rho_profile = np.real(rho_profile)
    # mu_ex_cal_1 = np.real(mu_ex_cal_1)
    mu_ex = np.real(mu_ex)
    bz_vext = np.real(bz_vext)

    fs = 1.0 / dz
    mu_ex_filtered = np.zeros_like(mu_ex)
    valid = ~np.isnan(mu_ex)
    cutoff = 2
    mu_ex_filtered_temp = filter_signal(mu_ex[valid], cutoff, fs, order=4)
    mu_ex_filtered[valid] = mu_ex_filtered_temp
    npy_path = os.path.join(output_dir,'simData'+config_index+'.npy')

    ml_data = np.zeros(data_length,dtype=dt)
    ml_data['rho_profile'] = rho_profile
    ml_data['bz_eff_vext'] = bz_vext*np.exp(mu_ex_b/M)
    ml_data['mu_ex'] = mu_ex
    ml_data['mu_ex_f'] = mu_ex_filtered
    ml_data['G'] = G_new
    ml_data['z'] = z 
    ml_data['rhob'] =  np.ones(data_length)*rhob_seg

    np.save(npy_path , ml_data)    

    plot_polymer_profiles(
        z=z,
        rho_profile=rho_profile,
        mu_ex = mu_ex,
        mu_ex_f = mu_ex_filtered,
        output_dir=output_dir,
        filename=f"{config_name}_ring_profiles.png",
        title_suffix="(Ring)"
    )


    return {'sim_index':config_index, 'path':npy_path}  


def padding_linear(npz_file,z, mu_ex_mat ,rho_profile_mat):
    
    dir_path = Path(npz_file).parent
    
    _,M = mu_ex_mat.shape

    temp_z = np.ones(M)

    mu_ex_wall = _get_padding_mu_ex(output_dir=dir_path)
    
    mu_ex_wall_mat = mu_ex_wall[:,None]*temp_z

    mu_ex_mat = np.concatenate((mu_ex_mat,mu_ex_wall_mat),axis=0)


    z_list_file = dir_path/"z_list_mu_ex_wall.dat"
    
    z_list_wall = np.loadtxt(z_list_file)


    z = np.concatenate((z,z_list_wall),axis=0)

    rho_profile_sup = np.zeros_like(z_list_wall)[:,None] * temp_z

    rho_profile_mat = np.concatenate((rho_profile_mat,rho_profile_sup),axis=0)

    sort_indices = np.argsort(z)
    z,rho_profile_mat,mu_ex_mat = z[sort_indices],rho_profile_mat[sort_indices,:],mu_ex_mat[sort_indices,:]
    return z,rho_profile_mat,mu_ex_mat



def recal_simulation_linear_padding( npz_file , mu_ex_0 ):
    import numpy as np
    data_dict = utils.load_npz_data(npz_file)
    np.seterr(divide='ignore', invalid='ignore')
    max_iters = 20000      # 最大迭代次数
    mix_rate = 0.02       # 每次迭代新计算结果的混合比例
    tol = 1e-5            # 提前终止的收敛容差

    H = data_dict['H']
    dz = data_dict['dz'].astype(float)
    sigma = data_dict['sigma']
    
    vext = data_dict['vext_values']
    vext = np.expand_dims(vext,axis=1)

    # mu_ex_0 = data_dict['mu_ex_0']
    mu_b = data_dict['mu_b']-mu_ex_0
    M = data_dict["M"].astype(int)
    rho_profile_mat = data_dict['rho_profile'] #*dz 由于之前代码的失误

    C = 1.0 / (4.0 * np.pi * sigma**2)

    

    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]

    
    data_length = rho_profile_mat.shape[0]
    rho_profile_mat = reflect_aug(rho_profile_mat)


    z = (np.arange(data_length) + 0.5) * dz

    G = np.ones_like(rho_profile_mat) 
    mu_ex_loc_mat = np.zeros(shape=(data_length,M))
    lambda_mat = mu_ex_loc_mat  + vext

    
    W = np.exp(-lambda_mat)
    W_new = np.copy(W)
    exp_mub = np.exp(mu_b)

    z_in = np.where( (z>sigma/2-dz/2) & ( z < H-sigma/2+dz/2 ) )[0]

    mu_ex_loc_avg = np.zeros(data_length)

    mu_ex_loc_avg_temp = np.zeros_like(mu_ex_loc_avg)     
    mu_ex_loc_avg_temp[z_in] = 1
    i = 0
    while True:
        i = i + 1
        alpha = _get_alpha(i)
        G[:,0] = 1
        for mono in range(1,M):
            G[:,mono] =np.real(weight2_slit_fourier( G[:,mono-1]*W[:,mono-1],H,dz,2*sigma)*C)

        W_new = rho_profile_mat/(G*G[:,::-1])/exp_mub
                
        W = W_new*alpha + W*(1-alpha)

        mu_ex_loc_mat[z_in,:] = -np.log(W[z_in,:])-vext[z_in,:]
        
        if i>100:
            mu_ex_loc_avg = np.mean(mu_ex_loc_mat, axis=1)
            
            converged,delta,_ = _check_convergence( mu_ex_loc_avg , mu_ex_loc_avg_temp ,tol)

            mu_ex_loc_avg_temp = mu_ex_loc_avg

            if converged:
                print(f"Converged after {i} iterations (delta = {delta:.6e})")
                break

            if i>max_iters:
                print(f"iter time reach max_iters, Not converged")
                break

    
    
    mu_ex_loc_mat[z_in,:] = -np.log(W[z_in,:])-vext[z_in,:]
    renormal_C = np.zeros(M)
    for mono in range(M):
        renormal_C[mono] = np.mean(mu_ex_loc_mat[z_in,mono].squeeze() - mu_ex_loc_avg[z_in] )
        mu_ex_loc_mat[z_in,mono] -= renormal_C[mono] 
    
    #print(renormal_C)
    #print(np.sum(renormal_C))
    


    # valid = ~np.isnan(mu_ex)
    z,rho_profile_mat,mu_ex_loc_mat = padding_linear(npz_file,z=z[z_in],mu_ex_mat = mu_ex_loc_mat[z_in,:],rho_profile_mat=rho_profile_mat[z_in,:])
    data_length_padding = z.shape
    
    bz_vext = np.zeros_like(z)

    z_no_pad_indices = np.where((z<H+1e-5) &(z>1e-3))[0]
    
    bz_vext[z_no_pad_indices] = np.exp(-vext).squeeze()
    
    mu_ex_filtered = np.zeros_like(mu_ex_loc_mat)
    
    
    mu_b = np.ones_like(z)*mu_b
    

    fs = 1.0 / dz
    cutoff = 2

    for mono in range(M):
        mu_ex_filtered[:,mono] = filter_signal(mu_ex_loc_mat[:,mono], cutoff, fs, order=4).reshape(-1)
    



    dt = np.dtype([('z','f8'),
            ('rho_profile_mat','f8',(M,)),
            ('mu_ex','f8',(M,)),
            ('mu_ex_f','f8',(M,)),
            ('bz_vext','f8'),
            ('mu_b','f8') ])
    
    ml_data = np.zeros(data_length_padding,dtype=dt)
    ml_data['rho_profile_mat'] = rho_profile_mat
    ml_data['bz_vext'] = bz_vext
    ml_data['mu_ex']   = mu_ex_loc_mat
    ml_data['mu_ex_f'] = mu_ex_filtered
    ml_data['z'] = z 
    ml_data['mu_b'] = mu_b
    
    
    npz_dir = os.path.dirname(npz_file)
    output_dir = npz_dir
    
    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]
    
    recal_npz = os.path.join(output_dir,f'pad_recal_result_{config_name}.npy')
    np.save(recal_npz , ml_data)    
    


    plot_polymer_profiles(
    z = z, 
    rho_profile=rho_profile_mat, 
    mu_ex = mu_ex_loc_mat, 
    mu_ex_f= mu_ex_filtered, 
    output_dir = output_dir, 
    filename = f"{config_name}_linear_profiles.png",
    title_suffix = f"(Linear, M={M})"
    )

    return {'sim_index':config_index, 'path':recal_npz}  





def plot_polymer_profiles(z, rho_profile, mu_ex, mu_ex_f=None, output_dir=".", filename="profiles.png", title_suffix=""):
    """
    通用高分子密度与过剩化学势绘图函数 (兼容单链结 Ring 和多链结 Linear)

    参数:
    - z: 一维数组，Z轴坐标 (N,)
    - rho_profile: 密度分布。可以是一维数组 (N,) 或 二维矩阵 (N, M)
    - mu_ex: 原始过剩化学势。可以是一维数组 (N,) 或 二维矩阵 (N, M)
    - mu_ex_f: 滤波后的过剩化学势 (可选)。如果不为None，将绘制3个子图进行对比。
    - output_dir: 图像保存目录
    - filename: 保存的文件名
    - title_suffix: 标题后缀，用于区分 Ring 或 Linear 等信息
    """
    
    # 1. 数据维度统一化：如果是一维数组 (单链结)，强制转换为 (N, 1) 的二维矩阵，方便后续统一用 for 循环处理
    if rho_profile.ndim == 1:
        rho_profile = rho_profile[:, np.newaxis]
        mu_ex = mu_ex[:, np.newaxis]
        if mu_ex_f is not None:
            mu_ex_f = mu_ex_f[:, np.newaxis]
            
    M = rho_profile.shape[1] # 获取链结数量
    has_filter = mu_ex_f is not None
    
    # 2. 初始化画布：根据是否有滤波数据决定画 2 个还是 3 个子图
    fig_height = 14 if has_filter else 10
    n_subplots = 3 if has_filter else 2
    fig, axes = plt.subplots(n_subplots, 1, figsize=(12, fig_height), sharex=True)
    
    ax_rho = axes[0]
    ax_mu = axes[1]
    ax_cmp = axes[2] if has_filter else None

    # 获取颜色映射 (如果只有一个链结用黑色，多个链结用彩虹色)
    colors = cm.viridis(np.linspace(0, 1, M)) if M > 1 else ['#1f77b4']

    # 3. 开始绘图
    for mono in range(M):
        label = f'链结 {mono+1}' if M > 1 else '高分子'
        color = colors[mono]
        
        # 子图 1: 密度分布
        ax_rho.plot(z, rho_profile[:, mono], '-o', markersize=3, linewidth=1.5, color=color, label=label)
        
        # 子图 2: 原始过剩化学势
        alpha_mu = 0.5 if has_filter else 1.0 # 如果下面有对比图，这张图就半透明弱化显示
        ax_mu.plot(z, mu_ex[:, mono], '-o', markersize=3, linewidth=1.5, alpha=alpha_mu, color=color, label=label)
        
        # 子图 3: 滤波前后对比 (仅当提供 mu_ex_f 时)
        if has_filter:
            # 原始数据用虚线
            legend_label_orig = f'原始 {label}' if (M <= 5 and mono == 0) or M == 1 else ""
            ax_cmp.plot(z, mu_ex[:, mono], '--', color=color, linewidth=1.0, alpha=0.5, label=legend_label_orig)
            
            # 滤波数据用实线
            legend_label_filt = f'滤波后 {label}' if (M <= 5 and mono == 0) or M == 1 else ""
            ax_cmp.plot(z, mu_ex_f[:, mono], '-', color=color, linewidth=2.0, alpha=0.9, label=legend_label_filt)

    # 4. 样式设置与标注 (消除重复代码)
    ax_rho.set_ylabel(r'$\rho(z)$', fontsize=12)
    ax_rho.set_title(f'密度分布 {title_suffix}')
    
    ax_mu.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
    ax_mu.set_title(f'原始过剩化学势分布 {title_suffix}')
    
    if has_filter:
        ax_cmp.set_xlabel(r'$z$', fontsize=12)
        ax_cmp.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax_cmp.set_title(f'过剩化学势滤波前后对比 {title_suffix}')
        
        # 添加图例说明线型区别
        custom_lines = [Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, alpha=0.7),
                        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.0, alpha=0.9)]
        ax_cmp.legend(custom_lines, ['原始数据 (虚线)', '滤波后数据 (实线)'], loc='best', fontsize=10)
    else:
        ax_mu.set_xlabel(r'$z$', fontsize=12)

    # 统一给所有子图加网格和图例
    for ax in axes:
        ax.grid(True, alpha=0.3)
        if M > 1: # 只有多链结时才显示单图图例，防止太拥挤
            ax.legend(loc='best', fontsize=10, ncol=min(M, 8))

    # 5. 保存图像
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, filename)
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)

    return plot_path

if __name__ == "__main__":
    npz_file = "D://Desktop//NewMuVTPolymer//output//OUT_config_0002_M6_Linear_H20.0_mu-1.41//avg_density_config_0002_M6_Linear_H20.0_mu-1.41.npz"
    recal_simulation_linear(npz_file)