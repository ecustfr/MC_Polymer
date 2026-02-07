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

    // Z_eff = Z / G/k^M
    // Z = \pord_i^M  \sum_k^{k_max} \exp[-\beta u_{ext}(k)] P_{road}(k)
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
        // 生成k_max个候选位置
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
    else if (monomer_index < this->M - 1)
    {
        // 从当前单体位置生成k_max个随机方向
        for (int k = 0; k < k_max; k++)
        {
            candidate_positions[k] = r_new[monomer_index - 1];
            add_random_unit_vector(candidate_positions[k].data(), gen);

            // 检查与其他聚合物的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_weights[k] = 0.0;
                continue;
            }

            // 检查与当前链的重叠
            for (int j = 0; j < monomer_index - 1; j++)
            {
                if (this->overlap_other_monomer_one(r_new[j].data(), candidate_positions[k].data()))
                {
                    z_weights[k] = 0.0;
                    break;
                }
            }

            if (z_weights[k] != 0.0)
            {
                // 计算两点之间的距离
                double dx = r_new[0][0] - candidate_positions[k][0];
                double dy = r_new[0][1] - candidate_positions[k][1];
                double dz = r_new[0][2] - candidate_positions[k][2];
                double R = sqrt(dx * dx + dy * dy + dz * dz);
                // 使用打表法计算路径概率 P_road
                P_road[k] = g_p_road_table.get_P_road(R, grow_step);
                // 计算外势玻尔兹曼因子 还有 P_road
                z_weights[k] = this->Vext_bz(candidate_positions[k].data()) * P_road[k];
            }
        }
    }
    else // 最后一个单体，需要闭环
    {
        // 使用 Find_Last 方法生成 k_max 个候选位置
        std::vector<double *> candidate_ptrs(k_max);

        for (int i = 0; i < k_max; i++)
        {
            candidate_ptrs[i] = candidate_positions[i].data();
        }
        // 计算从倒数第二个单体到第一个单体的方向向量
        // 使用 Find_Last 方法生成候选位置
        this->Find_Last(r_new[monomer_index - 1].data(), r_new[0].data(), candidate_ptrs.data(), k_max);

        // 检查重叠并计算权重
        for (int k = 0; k < k_max; k++)
        {

            // 检查与其他聚合物的重叠
            if (this->overlap_all_other_polymer(candidate_ptrs[k]))
            {
                z_weights[k] = 0.0;
                continue;
            }

            // 检查与当前链的重叠
            for (int j = 0; j < monomer_index - 1; j++)
            {
                if (this->overlap_other_monomer_one(r_new[j].data(), candidate_ptrs[k]))
                {
                    z_weights[k] = 0.0;
                    break;
                }
            }

            if (z_weights[k] != 0.0)
            {
                z_weights[k] = this->Vext_bz(candidate_ptrs[k]);
            }
        }

        // 不需要复制，因为 Find_Last 方法已经直接修改了 candidate_positions 的数据
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
    int grow_step = monomer_index + 1;

    int behind_count = std::max(this->M - monomer_index - 2, 0);
    std::vector<int> cal_index(behind_count);
    std::iota(cal_index.begin(), cal_index.end(), monomer_index + 2);

    // 对于最后一个单体，需要考虑与第一个单体的连接

    for (int k = 0; k < k_max - 1; k++)
    {
        candidate_positions[k] = r_delete[monomer_index + 1];
        add_random_unit_vector(candidate_positions[k].data(), gen);
        if (this->overlap_insert_polymer(candidate_positions[k], monomer_index + 1, r_delete, cal_index) || this->overlap_all_other_polymer(candidate_positions[k].data()))
        {
            z_eff[k] = 0.0;
        }
        else
        {
            if (monomer_index == 0)
            {
                z_eff[k] = this->Vext_bz(candidate_positions[k].data());
            }
            else
            {
                double dx = r_delete[M - 1][0] - candidate_positions[k][0];
                double dy = r_delete[M - 1][1] - candidate_positions[k][1];
                double dz = r_delete[M - 1][2] - candidate_positions[k][2];
                double R = sqrt(dx * dx + dy * dy + dz * dz);
                P_road[k] = g_p_road_table.get_P_road(grow_step, R);
                z_eff[k] = this->Vext_bz(candidate_positions[k].data()) * P_road[k];
            }
        }
    }

    if (monomer_index == 0)
    {
        z_eff[k_max - 1] = this->Vext_bz(r_delete[monomer_index].data());
    }
    else
    {
        double dx = r_delete[M - 1][0] - r_delete[monomer_index][0];
        double dy = r_delete[M - 1][1] - r_delete[monomer_index][1];
        double dz = r_delete[M - 1][2] - r_delete[monomer_index][2];
        double R = sqrt(dx * dx + dy * dy + dz * dz);
        P_road[k_max - 1] = g_p_road_table.get_P_road(grow_step, R);
        z_eff[k_max - 1] = this->Vext_bz(r_delete[monomer_index].data()) * P_road[k_max - 1];
    }

    double sum_z_eff = 0.0;
    for (int k = 0; k < k_max; k++)
    {
        sum_z_eff += z_eff[k];
    }

    *Z_eff = sum_z_eff / k_max / P_road[k_max - 1];
    return;
}

double MuVT_MC_RingPolymer::get_W_insert(int k_max)
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

    // 4. 更新cell列表，因为粒子位置已经改变
    if (this->MN_now > this->cl_threshold)
    {
        this->build_cell_list();
    }

    // 5. 先假设聚合物被删除
    int original_MN_now = this->MN_now;
    this->MN_now -= this->M;
    this->N_now -= 1;

    // 6. 使用delete_one_monomer方法计算每个单体的权重
    double Z_eff = 1.0;
    for (int monomer_index = 0; monomer_index < this->M - 1; monomer_index++)
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
