import numpy as np

def change_volumn(r_total, M, N, box_size, diameter2, expan_or_shrink, mc_cond_eps_volumn):
    """
    对聚合物系统进行体积缩放并检查是否重叠
    
    参数：
    - r_total: 所有单体位置，形状为 (MN, dim)
    - M: 每个聚合物的单体数
    - N: 聚合物数量
    - box_size: 盒子大小，形状为 (dim,)
    - diameter2: 直径平方，用于判断重叠
    - expan_or_shrink: 缩放方向，正数为缩小，负数为扩大
    - mc_cond_eps_volumn: 体积缩放因子
    
    返回：
    - 0: 无重叠
    - 1: 有重叠
    """
    dim = len(box_size)
    MN = M * N
    
    # 1. 计算聚合物质心和盒子中心
    r_center = np.zeros((N, dim))
    for i in range(N):
        for j in range(dim):
            r_center[i, j] = np.mean(r_total[i*M:(i+1)*M, j])
    
    box_center = box_size / 2
    
    # 2. 计算缩放因子和新盒子大小
    temp_eps_volumn = mc_cond_eps_volumn * expan_or_shrink
    new_box_size = box_size / (1 + temp_eps_volumn)
    
    # 3. 更新单体位置
    r_total_volumn = np.zeros_like(r_total)
    
    for i in range(N):
        # 计算移动量
        delta_move = np.zeros(dim)
        for j in range(dim):
            delta_move[j] = (1/(1 + temp_eps_volumn) * (r_center[i, j] - box_center[j]) + 
                           box_center[j] - r_center[i, j])
        
        # 更新位置
        for par_index in range(i*M, (i+1)*M):
            r_total_volumn[par_index] = r_total[par_index] + delta_move
    
    # 4. 检查是否重叠
    for par_index in range(MN - M):
        next_polymer_fir_monomer = ((par_index // M) + 1) * M
        
        for other_par_index in range(next_polymer_fir_monomer, MN):
            temp_dis = 0.0
            
            # 计算周期性距离平方
            for j in range(dim):
                dx = abs(r_total_volumn[par_index, j] - r_total_volumn[other_par_index, j])
                dx = min(dx, new_box_size[j] - dx)
                temp_dis += dx ** 2
                
                # 如果距离已经大于 diameter2，跳过
                if temp_dis > diameter2:
                    break
            
            # 检查是否重叠
            if temp_dis < diameter2:
                return 1  # 有重叠
    
    return 0  # 无重叠

# 示例使用
