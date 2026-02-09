#include "MuVT_MC_RingPolymer.h"
#include "Utils.h"
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>

double dis_z(double *r1, double *r2, double box_size_xy)
{
    double s = 0;
    s += SQR(r1[2] - r2[2]);
    s += SQR(dis_period(r1[0], r2[0], box_size_xy));
    s += SQR(dis_period(r1[1], r2[1], box_size_xy));
    return sqrt(s);
}

MuVT_MC_RingPolymer::MuVT_MC_RingPolymer(
    std::string configuration, const double mu_b, const int M, int init_N,
    double rho_b, double box_xy, double H, double rcut, int max_N)
    : MuVT_MC_LinearPolymer(configuration, mu_b, M, init_N, rho_b, box_xy, H, rcut, max_N)
{
    std::cout << "ring system is constructed" << std::endl;
}

MuVT_MC_RingPolymer::~MuVT_MC_RingPolymer()
{
    std::cout << "Ring Polymer System destroyed." << std::endl;
}

void MuVT_MC_RingPolymer::init_second()
{
    this->build_topology_map();

    if (this->check_configure_validity())
    {
        std::cout << "Initial configuration is valid" << std::endl;
    }
    else
    {
        std::cout << "error: Initial configuration is not valid, but continuing simulation" << std::endl;
        throw std::runtime_error("Initial configuration is not valid");
    }

    this->print_all_parameters();
}

bool MuVT_MC_RingPolymer::check_configure_validity()
{

    int now_par_index = 0;
    bool ad = false;
    int adj_index_1 = -1;
    int adj_index_2 = -1;

    for (int polymer_index = 0; polymer_index < this->N_now; polymer_index++)
    {
        for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
        {
            now_par_index = monomer_index + polymer_index * M;
            adj_index_1 = polymer_index * M + this->topology_map[monomer_index][0];
            adj_index_2 = polymer_index * M + this->topology_map[monomer_index][1];

            for (int other_mono = now_par_index + 1; other_mono < this->MN_now; other_mono++)
            {
                if (other_mono == adj_index_1 || other_mono == adj_index_2)
                {
                    continue;
                }
                if (this->overlap_other_monomer_one(this->r_total[now_par_index], this->r_total[other_mono]))
                {
                    std::cout << "overlap_index: " << now_par_index << " & " << other_mono << std::endl;
                    std::cout << this->r_total[now_par_index][0] << " " << this->r_total[now_par_index][1] << " " << this->r_total[now_par_index][2] << std::endl;
                    std::cout << this->r_total[other_mono][0] << " " << this->r_total[other_mono][1] << " " << this->r_total[other_mono][2] << std::endl;
                    std::cout << "dis:" << dis_z(this->r_total[now_par_index], this->r_total[other_mono], this->box_xy);
                    return false;
                }
            }
        }
    }

    return true;
}

// 在 RingPolymer 的构造函数中，或者重写 build_topology_map
void MuVT_MC_RingPolymer::build_topology_map()
{
    this->topology_map.resize(this->M);
    for (int i = 0; i < this->M; i++)
    {
        this->topology_map[i].clear();

        // 环状拓扑：邻居总是 (i-1) 和 (i+1)，需要处理周期性
        int prev = (i - 1 + this->M) % this->M;
        int next = (i + 1) % this->M;

        this->topology_map[i].push_back(prev);
        this->topology_map[i].push_back(next);
        // std::cout <<i<<":" <<this->topology_map[i][0] << " "<<this->topology_map[i][1] <<std::endl;
    }
}

void MuVT_MC_RingPolymer::rot_polymer_move(int polymer_index)
{
    std::vector<bool> rot_list(this->M, false);
    // 决定哪些单体需要旋转
    for (int i = 0; i < this->M; i++)
    {
        if (this->uni_dis(this->gen) < this->rot_ratio) // 使用父类的 protected gen
        {
            rot_list[i] = true;
        }
    }

    for (int i = 0; i < this->M; i++)
    {
        if (rot_list[i])
        {
            // 环状聚合物没有端点！所有单体都视为中间单体
            this->rot_mid_move(i, polymer_index);
        }
    }
}

void MuVT_MC_RingPolymer::insert_move(int k_max)
{
    this->num_insert++;

    // 1. 初始化新聚合物的位置
    std::vector<std::array<double, 3>> r_new(this->M);
    double Z_eff = 1.0;

    // Z_eff = z_weight / G/k^M
    // z_weight = \pord_i^M  \sum_k^{k_max} \exp[-\beta u_{ext}(k)] P_{road}(k)
    // G = \pord_i^M  sel_P_{road}(i)

    // 2. 插入每个单体
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &Z_eff, monomer_index, k_max))
        {
            return; // 单体插入失败
        }
    }

    // 3. 计算接受概率
    double acc_ratio = this->exp_mu_b * Z_eff / (this->N_now + 1) * (this->V);

    // 4. 接受或拒绝
    if (acc_ratio > this->uni_dis(this->gen))
    {
        this->acc_insert++;
        // 更新系统状态
        for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
        {
            for (int j = 0; j < 3; j++)
            {
                this->r_total[this->MN_now + monomer_index][j] = r_new[monomer_index][j];
            }
        }
        this->MN_now += this->M;
        this->N_now++;

        // 更新cell列表
        if (this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }

    return;
}

bool MuVT_MC_RingPolymer::insert_one_monomer(std::vector<std::array<double, 3>> &r_new, double *Z_eff, int monomer_index, int k_max)
{
    // monomer_index 正在生长的位置
    // 生成k_max个候选位置
    std::vector<std::array<double, 3>> candidate_positions(k_max);
    std::vector<double> z_weights(k_max, 1.0);

    std::vector<double> P_road(k_max, 1.0);
    int grow_step = this->M - monomer_index;
    if (monomer_index == 0)
    {
        

        for (int k = 0; k < k_max; k++)
        {
            candidate_positions[k][0] = box_size[0] * uni_dis(gen);
            candidate_positions[k][1] = box_size[1] * uni_dis(gen);
            candidate_positions[k][2] = box_size[2] * uni_dis(gen);

            // 检查与其他单体的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_weights[k] = 0.0;
            }
            else
            {
                // 计算外势玻尔兹曼因子
                z_weights[k] = this->Vext_bz(candidate_positions[k].data());
            }
        }
        
    }
    else if (monomer_index == this->M - 1)
    {
          // 使用 Find_Last 方法生成 k_max 个候选位置
        std::vector<double *> candidate_ptrs(k_max);

        for (int i = 0; i < k_max; i++)
        {
            candidate_ptrs[i] = candidate_positions[i].data();
        }
        // 计算从倒数第二个单体到第一个单体的方向向量
        // 使用 Find_Last 方法生成候选位置
        this->Find_Last(r_new[this->M-2].data(), r_new[0].data(), candidate_ptrs.data(), k_max);

        // 检查重叠并计算权重
        for (int k = 0; k < k_max; k++)
        {

            z_weights[k] = this->Vext_bz(candidate_ptrs[k]);
            // 检查与其他聚合物的重叠
            if (this->overlap_all_other_polymer(candidate_ptrs[k]))
            {
                z_weights[k] = 0.0;
                continue;
            }

            // 检查与当前链的重叠
            for (int j = 1; j < this->M-2; j++)
            {
                if (this->overlap_other_monomer_one(r_new[j].data(), candidate_ptrs[k]))
                {
                    z_weights[k] = 0.0;
                    break;
                }
            }


        }

      
    }
    else // middle monomer 
    {
          for (int k = 0; k < k_max; k++)
        {
            candidate_positions[k] = r_new[monomer_index - 1];
            add_random_unit_vector(candidate_positions[k].data(), gen);

            // 检查与其他聚合物的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_weights[k] = 0.0;
                P_road[k] = 0.0;
                continue;
            }

            // 检查与当前链的重叠
            for (int j = 0; j < monomer_index - 1; j++)
            {
                if (this->overlap_other_monomer_one(r_new[j].data(), candidate_positions[k].data()))
                {
                    z_weights[k] = 0.0;
                    P_road[k] = 0.0;
                    break;
                }
            }

            if (z_weights[k] != 0.0)
            {
                // 计算两点之间的距离

                double R = distance_array(r_new[0], candidate_positions[k]);
                
                // 使用打表法计算路径概率 P_road
                P_road[k] = g_p_road_table.get_P_road(R, grow_step);
                // 计算外势玻尔兹曼因子 还有 P_road
                z_weights[k] = this->Vext_bz(candidate_positions[k].data()) * P_road[k];
            }
        }
      
    }

    // 计算总权重
    double sum_z_eff = 0;
    for (int k = 0; k < k_max; k++)
    {
        sum_z_eff += z_weights[k];
    }
    if (sum_z_eff == 0)
    {
        *Z_eff = 0;
        return false; // 所有备选位置都重叠，插入失败
    }
    else
    {
        // 轮盘赌选择
        int select = RouteWheel(z_weights.data(), k_max, gen);

        r_new[monomer_index] = candidate_positions[select];

        *Z_eff *= sum_z_eff / k_max / P_road[select]; // 归一化，除以k_max
        return true;
    }
}

void MuVT_MC_RingPolymer::delete_one_monomer(std::vector<std::array<double, 3>> &r_delete, double *Z_eff, int monomer_index, int k_max)
{
    std::vector<std::array<double, 3>> candidate_positions(k_max);
    std::vector<double> z_eff(k_max, 1.0);
    std::vector<double> P_road(k_max, 1.0); // 把 k_max - 1 设成 老节点的位置
    int grow_step = this->M - monomer_index;

    // 对于最后一个单体，需要考虑与第一个单体的连接
    if (monomer_index == 0) // for first monomer
    {
        
        for (int k = 0; k < k_max-1; k++)
        {
            candidate_positions[k][0] = box_size[0] * uni_dis(gen);
            candidate_positions[k][1] = box_size[1] * uni_dis(gen);
            candidate_positions[k][2] = box_size[2] * uni_dis(gen);

            // 检查与其他单体的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_eff[k] = 0.0;
            }
            else
            {
                // 计算外势玻尔兹曼因子
                z_eff[k] = this->Vext_bz(candidate_positions[k].data());
            }
        }
        z_eff[k_max-1]= this->Vext_bz(r_delete[0].data());
        
    }
    else if (monomer_index == this->M - 1) // for last monomer
    {
        std::vector<double *> candidate_ptrs(k_max-1);

        for (int k = 0; k < k_max - 1; k++)
        {
            candidate_ptrs[k] = candidate_positions[k].data();
        }
        this->Find_Last(r_delete[this->M-2].data(), r_delete[0].data(), candidate_ptrs.data(), k_max-1);
        for( int k = 0; k < k_max - 1; k++)
        {
            z_eff[k] = this->Vext_bz(candidate_positions[k].data());
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_eff[k] = 0.0;
                continue;
            }
            for (int j = 1; j < this->M - 2; j++)
            {
                if (this->overlap_other_monomer_one(r_delete[j].data(), candidate_positions[k].data()))
                {
                    z_eff[k] = 0.0;
                    break;
                }
            }

        }
        z_eff[k_max - 1] = this->Vext_bz(r_delete[monomer_index].data());
    }
    else // for middle monomer
    {
        std::vector<int> cal_index;
        for (int j = 0; j < monomer_index; j++)
        {
            cal_index.push_back(j);
        }
        for (int k = 0; k < k_max - 1; k++)
        {
            candidate_positions[k] = r_delete[monomer_index - 1];
            add_random_unit_vector(candidate_positions[k].data(), gen);
            // this->overlap_insert_polymer(candidate_positions[k], monomer_index + 1, r_delete, cal_index)
            if (this->overlap_insert_polymer(candidate_positions[k], monomer_index - 1, r_delete, cal_index) || this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_eff[k] = 0.0;
                P_road[k] = 0.0;
            }
            else
            {
                double R = distance_array(r_delete[0], candidate_positions[k]);
                P_road[k] = g_p_road_table.get_P_road(R,grow_step);
                z_eff[k] = this->Vext_bz(candidate_positions[k].data()) * P_road[k];
            }
        }


        double R = distance_array(r_delete[0], r_delete[monomer_index]);
        P_road[k_max - 1] = g_p_road_table.get_P_road(R,grow_step);
        z_eff[k_max - 1] = this->Vext_bz(r_delete[monomer_index].data()) * P_road[k_max - 1];
       // std::cout << "R:"<<R << " "<<"try_num:"<<grow_step <<std::endl;
        //std::cout << "P_road:"<<P_road[k_max - 1] << std::endl;
    }

    double sum_z_eff = 0.0;
    for (int k = 0; k < k_max; k++)
    {
        sum_z_eff += z_eff[k];
    }

    *Z_eff *= sum_z_eff / k_max / P_road[k_max - 1];
    return;
}

double MuVT_MC_RingPolymer::get_W_insert_ring(int k_max)
{
    // 1. 初始化新聚合物的位置
    std::vector<std::array<double, 3>> r_new(this->M);
    double Z_eff = 1.0;

    // 2. 插入每个单体，使用环状聚合物的 insert_one_monomer 方法
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &Z_eff, monomer_index, k_max))
        {
            return 0.0; // 单体插入失败，返回 0
        }
    }

    // 3. 返回计算得到的权重
    return Z_eff;
}

double MuVT_MC_RingPolymer::get_G_insert_ring(double insert_z, int k_max) // 不需要指定插入的第一个单体位置，因为环状聚合物没有端点，插入位置对称
{
    std::vector<std::array<double, 3>> r_new(this->M);
    std::vector<int> is_inserted;
    is_inserted.reserve(M); // 不需要申请第一个插入粒子的

    int first_insert_index = 0; // 可以固定为0，因为环状聚合物没有端点，插入位置对称

    // --- 第一步：放置种子节点 ---

    for (int d = 0; d < 2; ++d)
    {
        r_new[first_insert_index][d] = box_size[d] * this->uni_dis(this->gen);
    }
    r_new[first_insert_index][2] = insert_z;
    double Z_eff = 1.0;

    for (int monomer_index = 1; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &Z_eff, monomer_index, k_max))
        {
            return 0.0; // 单体插入失败，返回 0
        }
    }

    return Z_eff;
}

// 删除移动：环状聚合物的删除逻辑
void MuVT_MC_RingPolymer::delete_move(int k_max, int delete_polymer_index)
{
    this->num_delete++;

    int delete_first_monomer_index = delete_polymer_index * this->M;
    int last_polymer_start = this->MN_now - this->M;

    // 2. 交换指针，把要删除的聚合物链放到最后面
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        std::swap(this->r_total[delete_first_monomer_index + monomer_index], this->r_total[last_polymer_start + monomer_index]);
    }

    // 3. 保存要删除的聚合物的位置
    std::vector<std::array<double, 3>> r_delete(this->M);
    for (int i = 0; i < this->M; i++)
    {
        for (int j = 0; j < 3; ++j)
        {
            r_delete[i][j] = this->r_total[last_polymer_start + i][j];
        }
    }
    int original_MN_now = this->MN_now;
    this->MN_now -= this->M;
    this->N_now -= 1;
    // 4. 更新cell列表，因为粒子位置已经改变
    if (this->MN_now > this->cl_threshold)
    {
        this->build_cell_list();
    }

    // 5. 先假设聚合物被删除


    // 6. 使用delete_one_monomer方法计算每个单体的权重
    double Z_eff = 1.0;
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        this->delete_one_monomer(r_delete, &Z_eff, monomer_index, k_max);
    }

    // 7. 计算接受概率
    double acc_ratio = 1.0 / (Z_eff * this->exp_mu_b) * (this->N_now + 1) / this->V;
    // 8. 决定是否接受删除
    if (acc_ratio > uni_dis(this->gen))
    {
        this->acc_delete++;
    }
    else
    {
        this->MN_now += M;
        this->N_now += 1;
        if (this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
}

bool MuVT_MC_RingPolymer::insert_recursive_ring(
    int next_idx, int parent_idx,
    double &Z_eff,
    std::vector<std::array<double, 3>> &r_new,
    std::vector<int> &is_inserted,
    int k_max)
{
    // overlap_test_indices 不是 已插入的单体。 比如：相连的 segment 不需要检测， 比如最后插入的单体，不需要检测头原子

    // 递归出口：环状聚合物终止条件是已插入M个单体
    if (next_idx >= this->M)
    {
        return true;
    }

    std::vector<std::array<double, 3>> candidates(k_max);
    std::vector<double> z_weight(k_max, 0.0);
    std::vector<double> P_road(k_max, 0.0);
    double sum_z_weight = 0.0;
    int grow_step = this->M - next_idx;

    if (next_idx < this->M - 1)
    {
        for (int k = 0; k < k_max; ++k)
        {
            candidates[k] = r_new[parent_idx]; // 基于父节点位置
            add_random_unit_vector(candidates[k].data(), this->gen);
            // 2. 碰撞检测：检查其他链 + 检查本链已插入的部分
            if (this->overlap_all_other_polymer(candidates[k].data()) || this->overlap_insert_polymer(candidates[k], parent_idx, r_new, is_inserted))
            {
                // 为了可读性应该加上下面两行，为了速度注释掉了，因为初始化的时候初始化好了
                // z_weight[k] = 0;
                // P_road[k] = 0.0;
            }
            else
            {
                double R = distance_array(r_new[0], candidates[k]);
                // 使用打表法计算路径概率 P_road
                P_road[k] = g_p_road_table.get_P_road(R, grow_step);
                z_weight[k] = this->Vext_bz(candidates[k].data()) * P_road[k];
                sum_z_weight += z_weight[k];
            }
        }
    }
    else // the last insert monomet
    {

        std::vector<double *> candidate_ptrs(k_max);

        for (int i = 0; i < k_max; i++)
        {
            candidate_ptrs[i] = candidates[i].data();
        }
        // 计算从倒数第二个单体到第一个单体的方向向量
        // 使用 Find_Last 方法生成候选位置
        this->Find_Last(r_new[this->M - 1].data(), r_new[0].data(), candidate_ptrs.data(), k_max);

        for (int k = 0; k < k_max; k++)
        {
            // 最后一个节点 需要单独对 已插入单体的碰撞做检测，因为最后一个节点与第一个节点相连，不能检测第一个节点，但是需要检测其他已插入的单体
            for (int monomer_index = 1; monomer_index < this->M - 1; monomer_index++)
            {
                if (this->overlap_other_monomer_one(r_new[monomer_index].data(), candidates[k].data()))
                {
                    // z_weight[k] = 0;
                    // P_road[k] = 0.0;
                    break;
                }
            }
            if (this->overlap_all_other_polymer(candidate_ptrs[k]))
            {
                // z_weight[k] = 0;
                // P_road[k] = 0.0;
            }
            else
            {
                // 最后一个 不需要计算 P_road
                z_weight[k] = this->Vext_bz(candidate_ptrs[k]);
                sum_z_weight += z_weight[k];
            }
        }
    }

    if (sum_z_weight <= 1e-40)
        return false; // 递归出口

    // 3. 选择并记录位置
    int selected = RouteWheel(z_weight.data(), k_max, this->gen);
    r_new[next_idx] = candidates[selected];
    is_inserted.push_back(next_idx); // 关键：标记为已插入

    // 4. 更新权重
    // Z_eff *=  / static_cast<double>(k_max);

    // 5. 继续向同一方向递归（环状聚合物只向+1方向递归）
    return insert_recursive_ring(next_idx + 1, next_idx, Z_eff, r_new, is_inserted, k_max);
}

// 新的递归插入方法，使用 insert_recursive_ring 方法

void MuVT_MC_RingPolymer::insert_move_recursive_ring(int k_max)
{
    this->num_insert++;

    // 1. 初始化新聚合物的位置
    std::vector<std::array<double, 3>> r_new(this->M);
    std::vector<int> is_inserted;
    is_inserted.reserve(this->M);
    double Z_eff = 1.0;

    // 2. 放置种子节点（第一个单体）在随机位置
    int first_insert_index = 0;
    for (int d = 0; d < 3; ++d)
    {
        r_new[first_insert_index][d] = this->box_size[d] * this->uni_dis(this->gen);
    }

    // 3. 检查种子节点是否与其他聚合物重叠
    if (this->overlap_all_other_polymer(r_new[first_insert_index].data()))
    {
        return; // 种子节点重叠，插入失败
    }

    // 标记种子节点为已插入
    is_inserted.push_back(first_insert_index);

    // 4. 使用 insert_recursive_ring 方法递归插入剩余的单体
    if (!this->insert_recursive_ring(first_insert_index + 1, first_insert_index, Z_eff, r_new, is_inserted, k_max))
    {
        return; // 递归插入失败
    }

    // 5. 计算接受概率
    double acc_ratio = this->exp_mu_b * Z_eff / (this->N_now + 1) * this->V;

    // 6. 决定是否接受插入
    if (acc_ratio > this->uni_dis(this->gen))
    {
        this->acc_insert++;

        // 7. 如果接受，更新系统状态
        for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
        {
            for (int j = 0; j < 3; j++)
            {
                this->r_total[this->MN_now + monomer_index][j] = r_new[monomer_index][j];
            }
        }
        this->MN_now += this->M;
        this->N_now++;

        if (this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }

    return;
}
