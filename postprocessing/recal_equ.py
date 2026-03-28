import numpy as np
import matplotlib.pyplot as plt
import json
import os
import argparse

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

    np.seterr(divide='ignore', invalid='ignore')
    polymer_type = data_dict['polymer_type']
    # print(polymer_type)
    if polymer_type=="Ring":
        result = recal_simulation_ring(npz_file)
    elif polymer_type == "Linear":
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

        # converged,delta,_ = _check_convergence(-np.log(W[z_in,:]) ,-np.log(W_new[z_in,:]),tol)
        
        # converged,delta,_ = _check_convergence(W ,W_new,tol)
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
    
    i += 1

    
    npz_dir = os.path.dirname(npz_file)
    output_dir = npz_dir
    
    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]
    
    recal_npz = os.path.join(output_dir,f'recal_result_{config_name}.npy')
    np.save(recal_npz , ml_data)    
    plot_linear_monomer_profiles(recal_npz, use_filtered=True, output_dir=output_dir)
    return {'sim_index':config_index, 'path':recal_npz}  

def recal_simulation_linear_anderson_try(npz_file):

    data_dict = utils.load_npz_data(npz_file)
    np.seterr(divide='ignore', invalid='ignore')
    max_iters = 20000      # 引入了Anderson，不需要20000步这么多，通常几百步就能搞定
    tol = 1e-4            # 提前终止的收敛容差

    H = data_dict['H']
    dz = data_dict['dz'].astype(float)
    sigma = data_dict['sigma']
    
    vext = data_dict['vext_values']
    vext = np.expand_dims(vext,axis=1)

    mu_ex_0 =  0
    mu_b = data_dict['mu_b']-mu_ex_0
    M = data_dict["M"].astype(int)
    rho_profile_mat = data_dict['rho_profile']

    C = 1.0 / (4.0 * np.pi * sigma**2)

    config_name = utils.extract_config_name(npz_file)
    config_index = config_name[7:11]
    
    data_length = rho_profile_mat.shape[0]
    rho_profile_mat = reflect_aug(rho_profile_mat)

    z = (np.arange(data_length) + 0.5) * dz

    G = np.ones_like(rho_profile_mat) 
    mu_ex_loc_mat = np.zeros(shape=(data_length,M))
    lambda_mat = mu_ex_loc_mat + vext  # 直接在 lambda 层面上操作
    
    W = np.exp(-lambda_mat)
    exp_mub = np.exp(mu_b)

    z_in = np.where( (z>sigma/2-dz/2) & ( z < H-sigma/2+dz/2 ) )[0]

# ================= Anderson Mixing 初始化 =================
    alpha = 0.02
    m_hist = 5
    
    # 【关键修正 1】：只提取物理有效区域的网格点，彻底排除边界垃圾数据干扰！
    valid_mask = np.zeros((data_length, M), dtype=bool)
    valid_mask[z_in, :] = True
    vec_len = np.sum(valid_mask) # 现在 vec_len 变小了，全是真实有效的数据
    
    X_hist = np.zeros((vec_len, m_hist))
    F_hist = np.zeros((vec_len, m_hist))
    hist_count = 0
    # ==========================================================

    i = 0
    while True:
        i = i + 1
        
        # 1. 正向传播计算 G
        G[:,0] = 1
        for mono in range(1,M):
            G[:,mono] = np.real(weight2_slit_fourier(G[:,mono-1]*W[:,mono-1],H,dz,2*sigma)*C)

        # 2. 计算新的目标权重与防溢出保护
        denom = G * G[:,::-1] * exp_mub
        denom[denom < 1e-50] = 1e-50
        W_new = rho_profile_mat / denom
        W_new[W_new < 1e-50] = 1e-50
        
        # 3. 计算残差
        lambda_new = -np.log(W_new)
        f_val = lambda_new - lambda_mat
        
        # 【关键修正 2】：只在有效区域内评估收敛误差
        delta = np.max(np.abs(f_val[valid_mask]))
        
        if i % 10 == 0:
            print(f"Iter: {i:4d}, delta: {delta:.6e}")

        if delta < tol:
            print(f"Converged after {i} iterations (delta = {delta:.6e})")
            break

        if i > max_iters:
            print(f"iter time reach max_iters, Not converged")
            break

        # ================= Anderson Mixing 核心计算 =================
        hist_count = min(hist_count + 1, m_hist)
        
        # 【关键修正 3】：只把有效区域的数据拍平成向量交给 Anderson
        x_vec = lambda_mat[valid_mask]
        f_vec = f_val[valid_mask]
        
        # 历史记录往后推一格
        if hist_count > 1:
            X_hist[:, 1:hist_count] = X_hist[:, 0:hist_count-1]
            F_hist[:, 1:hist_count] = F_hist[:, 0:hist_count-1]
            
        X_hist[:, 0] = x_vec
        F_hist[:, 0] = f_vec
        
        if hist_count == 1:
            lambda_next_vec = X_hist[:, 0] + alpha * F_hist[:, 0]
        else:
            k_cols = hist_count - 1
            Delta_X = X_hist[:, 0:1] - X_hist[:, 1:hist_count]
            Delta_F = F_hist[:, 0:1] - F_hist[:, 1:hist_count]
            
            # 【关键修正 4】：抛弃有偏的阻尼正则化，改用奇异值截断的伪逆 (pinv)
            # rcond=1e-5 会自动忽略那些导致发散和报错的极小平坦方向，只保留有物理意义的梯度
            gamma = np.linalg.pinv(Delta_F, rcond=1e-5) @ F_hist[:, 0]
            
            lambda_next_vec = X_hist[:, 0] + alpha * F_hist[:, 0] - (Delta_X + alpha * Delta_F) @ gamma
            
        # 【关键修正 5】：把预测出的新势场精确填回有效区域，死区保持原样！
        lambda_mat[valid_mask] = lambda_next_vec
        W = np.exp(-lambda_mat)
        # ==========================================================

    # 收敛后提取有效过剩化学势
    mu_ex_loc_mat = lambda_mat - vext
    
    # 按照原代码逻辑，将无效区域的值强行归零（可选，保证后续滤波干净）
    mask = np.ones(data_length, dtype=bool)
    mask[z_in] = False
    mu_ex_loc_mat[mask, :] = 0

    G = np.real(G)
    bz_vext = np.exp(-vext).squeeze()
    mu_b = np.ones_like(z)*mu_b

    mu_ex_filtered = np.zeros_like(mu_ex_loc_mat)
    
    fs = 1.0 / dz
    cutoff = 4

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
    
    recal_npz = os.path.join(output_dir,f'recal_result_{config_name}.npy')
    np.save(recal_npz , ml_data)    
    plot_linear_monomer_profiles(recal_npz,data=ml_data, use_filtered=True, output_dir=output_dir)
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



    # save_data = {'G':G_new,'rho_profile':rho_profile,'mu_ex_1':mu_ex_cal_1,'mu_ex_2':mu_ex_cal_2,'bz_vext':bz_vext}

    np.save(npy_path , ml_data)    
    # utils.save_npz_data(recal_npz,save_data)

    # plot_one_data(output_dir,np.arange(len(G_new)),rho_profile,mu_ex_filtered)
    # plot_filter_comparison(output_dir, np.arange(len(G_new)), rho_profile, mu_ex, dz=0.02, cutoff=10)
#   csv文件 
    """
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
    """

    return {'sim_index':config_index, 'path':npy_path}  


def plot_one_data(output_dir, grid_indice ,rho_profile, mu_ex,dz = 0.02 ):
    """
    绘制密度和化学势分布图并保存

    参数:
    output_dir: 输出目录路径
    grid_indice: 网格索引数组
    rho_profile: 密度分布数组
    mu_ex: 过剩化学势分布数组
    dz: 网格间距，默认0.02
    """
    import matplotlib.pyplot as plt

    # 计算x坐标：z = dz * (grid_indice + 0.5)
    z = dz * (grid_indice + 0.5)

    # 创建图形，上下两个子图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # 子图1：密度分布
    ax1.plot(z, rho_profile, 'b-', linewidth=2)
    ax1.set_ylabel(r'$\rho(z)$', fontsize=12)
    ax1.set_title('Density Profile')
    ax1.grid(True, alpha=0.3)

    # 子图2：过剩化学势分布
    ax2.plot(z, mu_ex, 'r-', linewidth=2)
    ax2.set_xlabel(r'$z$', fontsize=12)
    ax2.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
    ax2.set_title('Excess Chemical Potential Profile')
    ax2.grid(True, alpha=0.3)

    # 调整布局
    plt.tight_layout()

    # 保存图像
    import os
    plot_path = os.path.join(output_dir, 'density_muex_profile.png')
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)

    #print(f"Plot saved to: {plot_path}")

     


def plot_filter_comparison(output_dir, grid_indice, rho_profile, mu_ex, dz=0.02, cutoff=None):
    """
    绘制密度分布和滤波前后的化学势分布对比图

    参数:
    output_dir: 输出目录路径
    grid_indice: 网格索引数组
    rho_profile: 密度分布数组
    mu_ex: 过剩化学势分布数组（原始）
    dz: 网格间距，默认0.02
    cutoff: 截止频率（单位：1/长度单位）。如果为None，则不进行滤波。
    """
    import matplotlib.pyplot as plt
    import os

    # 计算x坐标：z = dz * (grid_indice + 0.5)
    z = dz * (grid_indice + 0.5)

    # 创建图形，三个子图
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # 子图1：密度分布
    ax1.plot(z, rho_profile, 'b-', linewidth=2)
    ax1.set_ylabel(r'$\rho(z)$', fontsize=12)
    ax1.set_title('Density Profile')
    ax1.grid(True, alpha=0.3)

    # 子图2：原始化学势分布
    ax2.plot(z, mu_ex, 'r-', linewidth=2, label='Original')
    ax2.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
    ax2.set_title('Excess Chemical Potential Profile (Original)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 子图3：滤波前后对比（如果cutoff不为None）
    if cutoff is not None:
        # 计算采样频率：fs = 1/dz（每单位长度的点数）
        fs = 1.0 / dz
        mu_ex_filtered = np.zeros_like(mu_ex)
        valid = ~np.isnan(mu_ex)
        mu_ex_filtered_temp = filter_signal(mu_ex[valid], cutoff, fs, order=4)
        mu_ex_filtered[valid] = mu_ex_filtered_temp
        ax3.plot(z, mu_ex, 'r-', linewidth=1, alpha=0.5, label='Original')
        ax3.plot(z, mu_ex_filtered, 'b-', linewidth=2, label='Filtered')
        ax3.set_xlabel(r'$z$', fontsize=12)
        ax3.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax3.set_title(f'Excess Chemical Potential Profile Comparison (Cutoff={cutoff:.3f})')
        ax3.grid(True, alpha=0.3)
        ax3.legend()

        # 保存滤波后的数据到文本文件
        data_save_path = os.path.join(output_dir, 'filtered_muex_data.txt')
        np.savetxt(data_save_path, np.column_stack([z, mu_ex, mu_ex_filtered]),
                  header='z, mu_ex_original, mu_ex_filtered', comments='')
        print(f"Filtered data saved to: {data_save_path}")
    else:
        # 如果没有滤波，只显示原始数据
        ax3.plot(z, mu_ex, 'r-', linewidth=2)
        ax3.set_xlabel(r'$z$', fontsize=12)
        ax3.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax3.set_title('Excess Chemical Potential Profile (Original)')
        ax3.grid(True, alpha=0.3)

    # 调整布局
    plt.tight_layout()

    # 保存图像
    if cutoff is not None:
        plot_path = os.path.join(output_dir, f'density_muex_filtered_cutoff_{cutoff:.3f}.png')
    else:
        plot_path = os.path.join(output_dir, 'density_muex_comparison.png')
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"Plot saved to: {plot_path}")
    return plot_path


def plot_linear_monomer_profiles(npz_file=None, data=None, use_filtered=True, output_dir=None, plot_filename=None):
    """
    绘制线性高分子各个链结的密度分布和过剩化学势分布

    参数:
    npz_file: 重新计算结果的npz文件路径（由recal_simulation_linear生成）
    data: 包含绘图数据的字典，如果提供则忽略npz_file。必须包含以下键:
          'z': 形状 (data_length,)
          'rho_profile_mat': 形状 (data_length, M)
          'mu_ex': 形状 (data_length, M)
          'mu_ex_f': 形状 (data_length, M)
    use_filtered: 是否使用滤波后的过剩化学势数据（mu_ex_f），默认为True
                  当use_filtered=False时，仅显示原始数据；当use_filtered=True时，同时显示滤波和非滤波数据（对比模式）
    output_dir: 输出目录路径。如果为None且提供npz_file，则使用npz_file所在目录；
                如果为None且提供data，则使用当前目录
    plot_filename: 输出图像文件名。如果为None，则根据npz_file或默认名称生成

    返回:
    plot_path: 保存的图像文件路径
    """
    import matplotlib.pyplot as plt
    import os
    import numpy as np

    # 确定数据来源
    if data is not None:
        # 使用传入的数据
        z = data['z']
        rho_profile_mat = data['rho_profile_mat']
        mu_ex = data['mu_ex']
        mu_ex_f = data['mu_ex_f']
        data_source = "data"
    elif npz_file is not None:
        # 从文件加载数据
        data_from_file = np.load(npz_file)
        z = data_from_file['z']
        rho_profile_mat = data_from_file['rho_profile_mat']
        mu_ex = data_from_file['mu_ex']
        mu_ex_f = data_from_file['mu_ex_f']
        data_source = "file"
    else:
        raise ValueError("必须提供 npz_file 或 data 参数")

    # 获取链结数量
    M = rho_profile_mat.shape[1]

    # 设置输出目录
    if output_dir is None:
        if data_source == "file":
            output_dir = os.path.dirname(npz_file)
        else:
            output_dir = "."

    # 根据use_filtered参数决定绘图模式
    if not use_filtered:
        # 原始模式：仅显示原始数据，保持原有行为
        mu_ex_data = mu_ex
        data_label = "原始"

        # 创建图形：上下两个子图，分别显示密度分布和过剩化学势分布
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

        # 子图1：各个链结的密度分布
        for mono in range(M):
            ax1.plot(z, rho_profile_mat[:, mono], linewidth=1.5, label=f'链结 {mono+1}')

        ax1.set_ylabel(r'$\rho(z)$', fontsize=12)
        ax1.set_title(f'线性高分子密度分布（M={M}个链结）')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='best', fontsize=10, ncol=min(M, 5))

        # 子图2：各个链结的过剩化学势分布（原始）
        for mono in range(M):
            ax2.plot(z, mu_ex_data[:, mono], linewidth=1.5, label=f'链结 {mono+1}')

        ax2.set_xlabel(r'$z$', fontsize=12)
        ax2.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax2.set_title(f'线性高分子过剩化学势分布（{data_label}，M={M}个链结）')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best', fontsize=10, ncol=min(M, 5))

        # 确定输出文件名
        if plot_filename is None:
            if data_source == "file":
                basename = os.path.splitext(os.path.basename(npz_file))[0]
                plot_filename = f'{basename}_monomer_profiles_original.png'
            else:
                plot_filename = 'linear_monomer_profiles_original.png'

    else:
        # 对比模式：同时显示滤波和非滤波数据，学习环状聚合物的绘图模式
        # 创建图形：三个子图，分别显示密度分布、原始化学势分布、滤波前后对比
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)

        # 子图1：各个链结的密度分布
        for mono in range(M):
            ax1.plot(z, rho_profile_mat[:, mono], linewidth=1.5, label=f'链结 {mono+1}')

        ax1.set_ylabel(r'$\rho(z)$', fontsize=12)
        ax1.set_title(f'线性高分子密度分布（M={M}个链结）')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='best', fontsize=10, ncol=min(M, 5))

        # 子图2：各个链结的原始过剩化学势分布
        for mono in range(M):
            ax2.plot(z, mu_ex[:, mono], linewidth=1.5, alpha=0.7, label=f'链结 {mono+1}')

        ax2.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax2.set_title(f'线性高分子原始过剩化学势分布（M={M}个链结）')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best', fontsize=10, ncol=min(M, 5))

        # 子图3：各个链结的滤波前后对比（学习环状聚合物的plot_filter_comparison模式）
        # 为所有链结绘制对比，使用颜色映射区分不同链结，线型区分原始/滤波数据
        import matplotlib.cm as cm
        import numpy as np

        # 获取颜色映射
        colors = cm.viridis(np.linspace(0, 1, M)) if M > 1 else [(0, 0, 0, 1)]

        for mono in range(M):
            color = colors[mono]
            # 原始数据（虚线，半透明）
            ax3.plot(z, mu_ex[:, mono], '--', color=color, linewidth=0.8, alpha=0.5,
                     label=f'原始 链结 {mono+1}' if M <= 10 and mono == 0 else "")
            # 滤波后数据（实线）
            ax3.plot(z, mu_ex_f[:, mono], '-', color=color, linewidth=1.5, alpha=0.8,
                     label=f'滤波后 链结 {mono+1}' if M <= 10 and mono == 0 else "")

        # 添加图例说明线型区别
        from matplotlib.lines import Line2D
        custom_lines = [Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, alpha=0.7),
                        Line2D([0], [0], color='gray', linestyle='-', linewidth=1.5, alpha=0.9)]
        ax3.legend(custom_lines, ['原始数据 (虚线)', '滤波后数据 (实线)'], loc='best', fontsize=10)

        ax3.set_xlabel(r'$z$', fontsize=12)
        ax3.set_ylabel(r'$\mu_{\mathrm{ex}}(z)$', fontsize=12)
        ax3.set_title(f'线性高分子过剩化学势滤波前后对比（M={M}个链结，颜色区分不同链结）')
        ax3.grid(True, alpha=0.3)

        # 确定输出文件名
        if plot_filename is None:
            if data_source == "file":
                basename = os.path.splitext(os.path.basename(npz_file))[0]
                plot_filename = f'{basename}_monomer_profiles_filtered.png'
            else:
                plot_filename = 'linear_monomer_profiles_filtered.png'

    # 调整布局
    plt.tight_layout()

    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"线性高分子各链结分布图已保存至: {plot_path}")
    return plot_path


if __name__ == "__main__":
    npz_file = "D://Desktop//NewMuVTPolymer//output//OUT_config_0002_M6_Linear_H20.0_mu-1.41//avg_density_config_0002_M6_Linear_H20.0_mu-1.41.npz"
    recal_simulation_linear(npz_file)