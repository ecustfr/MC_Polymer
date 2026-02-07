#include "MuVT_MC_LinearPolymer_Ext.h"
#include "Utils.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>

MuVT_MC_LinearPolymer_Ext::MuVT_MC_LinearPolymer_Ext(std::string configuration,
    const double mu_b,
    const int M,
    int init_N,
    double rho_b,
    double box_xy,
    double H,
    double rcut,
    int max_N,
    const std::vector<double>& An,
    const std::vector<double>& phi_n,
    const std::vector<std::vector<double>>& Vlin_par,
    const std::vector<std::vector<double>>& x_tar,
    double C)
    : MuVT_MC_LinearPolymer(configuration, mu_b, M, init_N, rho_b, box_xy, H, rcut, max_N),
      An(An), phi_n(phi_n), Vlin_par(Vlin_par), x_tar(x_tar), C(C)
{
    // 构造函数体
}

MuVT_MC_LinearPolymer_Ext::~MuVT_MC_LinearPolymer_Ext()
{
    // 析构函数体
}

// 计算单个位置的势场能量（移除边界检查）
double MuVT_MC_LinearPolymer_Ext::calculate_potential(const double* pos) const
{
    double z = pos[2];
    double box_z = box_size[2];
    
    // 正弦分量 Vext1
    double Vext1 = 0.0;
    for (size_t n = 0; n < An.size(); ++n)
    {
        double arg = 2 * M_PI * (n + 1) * z / box_z + phi_n[n];
        Vext1 += An[n] * sin(arg);
    }
    
    // 线性分段分量 Vext2
    double Vext2 = 0.0;
    for (size_t n = 0; n < 4; ++n)
    {
        if (z >= x_tar[0][n] && z <= x_tar[1][n])
        {
            double slope = (Vlin_par[1][n] - Vlin_par[0][n]) / (x_tar[1][n] - x_tar[0][n]);
            Vext2 += Vlin_par[0][n] + slope * (z - x_tar[0][n]);
        }
    }
    
    return C * (Vext1 + Vext2);
}

// 计算玻尔兹曼因子 exp(-beta*U)
double MuVT_MC_LinearPolymer_Ext::calculate_boltzmann(const double* pos) const
{
    double U = calculate_potential(pos);
    return exp(-beta * U);
}

// 重写平移移动
void MuVT_MC_LinearPolymer_Ext::trans_move(int polymer_index)
{
    num_trans++;
    double delta_move[3];
    
    for (int j = 0; j < 3; j++)
    {
        delta_move[j] = real_dis(gen) * box_size[j] * eps_trans;
    }
    
    // 保存原位置
    std::vector<std::array<double, 3>> r_old(M);
    int start = polymer_index * M;
    for (int i = 0; i < M; i++)
    {
        r_old[i] = { r_total[start + i][0], r_total[start + i][1], r_total[start + i][2] };
    }
    
    // 检查重叠
    bool overlap = false;
    for (int i = 0; i < M; i++)
    {
        double r_try[3] = { r_old[i][0] + delta_move[0], r_old[i][1] + delta_move[1], r_old[i][2] + delta_move[2] };
        if (overlap_all_other_polymer(r_try))
        {
            overlap = true;
            break;
        }
    }
    
    if (!overlap)
    {
        // 计算势场能量差
        double delta_U = 0.0;
        for (int i = 0; i < M; i++)
        {
            double r_new[3] = { r_old[i][0] + delta_move[0], r_old[i][1] + delta_move[1], r_old[i][2] + delta_move[2] };
            delta_U += calculate_potential(r_new) - calculate_potential(r_old[i].data());
        }
        
        // 计算玻尔兹曼因子
        double boltz_factor = exp(-beta * delta_U);
        
        if (boltz_factor > uni_dis(gen))
        {
            acc_trans++;
            // 更新位置
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    r_total[start + i][j] += delta_move[j];
                }
                // update_particle_in_cell_list(start + i);
            }
        }
    }
}

// 重写中间单体旋转（生成2k_max-1个位置）
void MuVT_MC_LinearPolymer_Ext::rot_mid_move(int monomer_index, int polymer_index)
{
    num_rot++;
    
    // 保存原位置（当前位置）
    double r_current[3] = { r_total[monomer_index][0], r_total[monomer_index][1], r_total[monomer_index][2] };
    
    // 获取相邻单体位置
    double r_prev[3] = { r_total[monomer_index - 1][0], r_total[monomer_index - 1][1], r_total[monomer_index - 1][2] };
    double r_next[3] = { r_total[monomer_index + 1][0], r_total[monomer_index + 1][1], r_total[monomer_index + 1][2] };
    
    // 生成2k_max-1个候选位置
    std::vector<std::array<double, 3>> r_try_list(2 * k_max - 1);
    std::vector<double> phi_try(2 * k_max - 1, 1.0);
    
    // 生成前k_max个新位置
    for (int k = 0; k < k_max; k++)
    {
        // 生成随机旋转
        r_try_list[k] = {r_prev[0], r_prev[1], r_prev[2]};
        add_random_unit_vector(r_try_list[k].data(), gen);
        
        // 检查重叠
        if (overlap_other_polymer(r_try_list[k].data(), polymer_index))
        {
            phi_try[k] = 0.0;
            continue;
        }
        
        // 计算势场玻尔兹曼因子
        phi_try[k] = calculate_boltzmann(r_try_list[k].data());
    }
    
    // 生成后k_max-1个位置（老位置方向）
    for (int k = 0; k < k_max - 1; k++)
    {// 生成随机旋转
        r_try_list[k_max + k] = {r_next[0], r_next[1], r_next[2]};
        add_random_unit_vector(r_try_list[k_max + k].data(), gen);
        
        // 检查重叠
        if (overlap_other_polymer(r_try_list[k_max + k].data(), polymer_index))
        {
            phi_try[k_max + k] = 0.0;
            continue;
        }
        
        // 计算势场玻尔兹曼因子
        phi_try[k_max + k] = calculate_boltzmann(r_try_list[k_max + k].data());
    }
    
    // 添加当前位置
    r_try_list[2 * k_max - 2] = { r_current[0], r_current[1], r_current[2] };
    if (overlap_other_polymer(r_current, polymer_index))
    {
        phi_try[2 * k_max - 2] = 0.0;
    }
    else
    {
        phi_try[2 * k_max - 2] = calculate_boltzmann(r_current);
    }
    
    // 轮盘赌选择
    int select_index = RouteWheel(phi_try.data(), 2 * k_max - 1, gen);
    
    // 更新位置
    r_total[monomer_index][0] = r_try_list[select_index][0];
    r_total[monomer_index][1] = r_try_list[select_index][1];
    r_total[monomer_index][2] = r_try_list[select_index][2];
    
    acc_rot++;    // update_particle_in_cell_list(monomer_index);
}

// 重写末端单体旋转（生成2k_max-1个位置）
void MuVT_MC_LinearPolymer_Ext::rot_end_move(int monomer_index, int polymer_index, int dirct)
{
    num_rot++;
    
    // 保存原位置（当前位置）
    double r_current[3] = { r_total[monomer_index][0], r_total[monomer_index][1], r_total[monomer_index][2] };
    
    // 获取相邻单体位置
    int neighbor_index = dirct == 0 ? monomer_index + 1 : monomer_index - 1;
    double r_neighbor[3] = { r_total[neighbor_index][0], r_total[neighbor_index][1], r_total[neighbor_index][2] };
    
    // 生成2k_max-1个候选位置
    std::vector<std::array<double, 3>> r_try_list(2 * k_max - 1);
    std::vector<double> phi_try(2 * k_max - 1, 1.0);
    
    // 生成前k_max个新位置
    for (int k = 0; k < k_max; k++)
    {
        // 生成随机旋转
        r_try_list[k] = {r_neighbor[0], r_neighbor[1], r_neighbor[2]};
        add_random_unit_vector(r_try_list[k].data(), gen);
        
        // 检查重叠
        if (overlap_other_polymer(r_try_list[k].data(), polymer_index))
        {
            phi_try[k] = 0.0;
            continue;
        }
        
        // 计算势场玻尔兹曼因子
        phi_try[k] = calculate_boltzmann(r_try_list[k].data());
    }
    
    // 生成后k_max-1个位置（老位置方向）
    for (int k = 0; k < k_max - 1; k++)
    {// 生成随机旋转
        r_try_list[k_max + k] = {r_current[0], r_current[1], r_current[2]};
        add_random_unit_vector(r_try_list[k_max + k].data(), gen);
        
        // 检查重叠
        if (overlap_other_polymer(r_try_list[k_max + k].data(), polymer_index))
        {
            phi_try[k_max + k] = 0.0;
            continue;
        }
        
        // 计算势场玻尔兹曼因子
        phi_try[k_max + k] = calculate_boltzmann(r_try_list[k_max + k].data());
    }
    
    // 添加当前位置
    r_try_list[2 * k_max - 2] = { r_current[0], r_current[1], r_current[2] };
    if (overlap_other_polymer(r_current, polymer_index))
    {
        phi_try[2 * k_max - 2] = 0.0;
    }
    else
    {
        phi_try[2 * k_max - 2] = calculate_boltzmann(r_current);
    }
    
    // 轮盘赌选择
    int select_index = RouteWheel(phi_try.data(), 2 * k_max - 1, gen);
    
    // 更新位置
    r_total[monomer_index][0] = r_try_list[select_index][0];
    r_total[monomer_index][1] = r_try_list[select_index][1];
    r_total[monomer_index][2] = r_try_list[select_index][2];
    
    acc_rot++;    // update_particle_in_cell_list(monomer_index);
}

// 重写插入移动（仅计算CBMC权重W）
void MuVT_MC_LinearPolymer_Ext::insert_move(int k_max)
{
    num_insert++;
    
    // 生成起始位置
    r_insert_temp[0][0] = box_size[0] * uni_dis(gen);
    r_insert_temp[0][1] = box_size[1] * uni_dis(gen);
    r_insert_temp[0][2] = box_size[2] * uni_dis(gen);
    
    if (overlap_all_other_polymer(r_insert_temp[0]))
    {
        return;
    }
    
    double W = 1.0;
    double r_now[3] = { r_insert_temp[0][0], r_insert_temp[0][1], r_insert_temp[0][2] };
    
    std::vector<double> Phi_HS_try(k_max, 1.0);  // 硬球贡献
    std::vector<double> Phi_ext_try(k_max, 1.0); // 势场玻尔兹曼因子
    
    std::vector<std::array<double, 3>> r_temp(k_max, std::array<double, 3>{ r_now[0], r_now[1], r_now[2] });
    
    for (int monomer_index = 1; monomer_index < M; monomer_index++)
    {
        std::fill(Phi_HS_try.begin(), Phi_HS_try.end(), 1.0);
        std::fill(Phi_ext_try.begin(), Phi_ext_try.end(), 1.0);
        
        // 生成k_max个随机方向
        for (int k = 0; k < k_max; k++)
        {
            add_random_unit_vector(r_temp[k].data(), gen);
        }
        
        // 检查重叠并计算玻尔兹曼因子
        for (int k = 0; k < k_max; k++)
        {
            if (overlap_all_other_polymer(r_temp[k].data()))
            {
                Phi_HS_try[k] = 0.0;
                Phi_ext_try[k] = 0.0;
                continue;
            }
            
            for (int old_monomer = 0; old_monomer < monomer_index - 1; old_monomer++)
            {
                if (overlap_other_monomer_one(r_insert_temp[old_monomer], r_temp[k].data()))
                {
                    Phi_HS_try[k] = 0.0;
                    Phi_ext_try[k] = 0.0;
                    break;
                }
            }
            
            // 计算势场玻尔兹曼因子
            if (Phi_HS_try[k] > 0.0)
            {
                Phi_ext_try[k] = calculate_boltzmann(r_temp[k].data());
            }
        }
        
        if (std::all_of(Phi_HS_try.begin(), Phi_HS_try.end(), [](double v) { return v == 0.0; }))
        {
            return; // 所有候选位置都重叠，插入失败
        }
        
        // 合并硬球和势场贡献
        std::vector<double> Phi_total_try(k_max);
        for (int k = 0; k < k_max; k++)
        {
            Phi_total_try[k] = Phi_HS_try[k] * Phi_ext_try[k];
        }
        
        // 更新权重W
        W *= sum_1d(Phi_total_try.data(), k_max) / k_max;
        
        // 轮盘赌选择
        int select_index = RouteWheel(Phi_total_try.data(), k_max, gen);
        r_insert_temp[monomer_index][0] = r_temp[select_index][0];
        r_insert_temp[monomer_index][1] = r_temp[select_index][1];
        r_insert_temp[monomer_index][2] = r_temp[select_index][2];
        
        // 更新候选位置
        std::fill(r_temp.begin(), r_temp.end(), std::array<double, 3>{ r_temp[select_index][0], r_temp[select_index][1], r_temp[select_index][2] });
    }
    
    // 接受率计算 - 仅使用W，不计算总能量
    double acc_ratio = exp_mu_b * W / (N_now + 1) * V;
    if (acc_ratio > uni_dis(gen))
    {
        acc_insert++;
        // 插入聚合物
        for (int monomer_index = 0; monomer_index < M; monomer_index++)
        {
            for (int j = 0; j < 3; j++)
            {
                r_total[MN_now + monomer_index][j] = r_insert_temp[monomer_index][j];
            }
        }
        MN_now += M;
        N_now++;
        
        if (MN_now > cl_threshold)
        {
            build_cell_list();
        }
    }
}

// 重写删除移动（CBMC实现，包含当前位置权重）
void MuVT_MC_LinearPolymer_Ext::delete_move(int k_max, int delete_polymer_index)
{
    num_delete++;
    
    int delete_first_monomer_index = delete_polymer_index * M;
    int last_polymer_start = MN_now - M;
    
    // 交换要删除的聚合物到末尾
    for (int monomer_index = 0; monomer_index < M; monomer_index++)
    {
        std::swap(r_total[delete_first_monomer_index + monomer_index], r_total[last_polymer_start + monomer_index]);
    }
    
    // 更新cell列表
    if (MN_now > cl_threshold)
    {
        build_cell_list();
    }
    
    // 假设删除，更新状态
    MN_now -= M;
    N_now -= 1;
    
    double W = 1.0;
    
    // 从末端开始回溯，生成候选位置
    for (int monomer_index = M - 1; monomer_index > 0; monomer_index--)
    {
        int current_monomer = MN_now + monomer_index;
        int prev_monomer = MN_now + monomer_index - 1;
        
        // 当前位置
        double r_current[3] = { r_total[current_monomer][0], r_total[current_monomer][1], r_total[current_monomer][2] };
        
        // 生成k_max-1个候选位置（加上当前位置共k_max个）
        std::vector<std::array<double, 3>> r_temp(k_max - 1);
        for (int k = 0; k < k_max - 1; k++)
        {
            r_temp[k] = { r_total[prev_monomer][0], r_total[prev_monomer][1], r_total[prev_monomer][2] };
            add_random_unit_vector(r_temp[k].data(), gen);
        }
        
        // 计算所有候选位置的权重（包括当前位置）
        std::vector<double> Phi_total_try(k_max, 1.0);
        
        // 当前位置的权重
        if (overlap_other_polymer(r_current, delete_polymer_index))
        {
            Phi_total_try[0] = 0.0;
        }
        else
        {
            for (int old_monomer = 0; old_monomer < monomer_index - 1; old_monomer++)
            {
                if (overlap_other_monomer_one(r_total[MN_now + old_monomer], r_current))
                {
                    Phi_total_try[0] = 0.0;
                    break;
                }
            }
            if (Phi_total_try[0] > 0.0)
            {
                Phi_total_try[0] = calculate_boltzmann(r_current);
            }
        }
        
        // 其他候选位置的权重
        for (int k = 0; k < k_max - 1; k++)
        {
            if (overlap_other_polymer(r_temp[k].data(), delete_polymer_index))
            {
                Phi_total_try[k + 1] = 0.0;
                continue;
            }
            
            for (int old_monomer = 0; old_monomer < monomer_index - 1; old_monomer++)
            {
                if (overlap_other_monomer_one(r_total[MN_now + old_monomer], r_temp[k].data()))
                {
                    Phi_total_try[k + 1] = 0.0;
                    break;
                }
            }
            
            if (Phi_total_try[k + 1] > 0.0)
            {
                Phi_total_try[k + 1] = calculate_boltzmann(r_temp[k].data());
            }
        }
        
        // 更新权重W
        W *= (sum_1d(Phi_total_try.data(), k_max) + 0.0) / k_max;
    }
    
    // 计算接受率
    double acc_ratio = 1.0 / W * (N_now + 1) / V * 1.0 / exp_mu_b;
    
    if (acc_ratio > uni_dis(gen))
    {
        acc_delete++;
        // 接受删除，重建cell列表
        build_cell_list();
    }
    else
    {
        // 拒绝删除，恢复状态
        MN_now += M;
        N_now += 1;
        
        // 交换回原来的位置
        for (int monomer_index = 0; monomer_index < M; monomer_index++)
        {
            std::swap(r_total[delete_first_monomer_index + monomer_index], r_total[last_polymer_start + monomer_index]);
        }
        
        // 更新cell列表
        if (MN_now > cl_threshold)
        {
            build_cell_list();
        }
    }
}
