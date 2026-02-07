#include "Bulk_RingPolymer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

Bulk_RingPolymer::Bulk_RingPolymer(std::string configuration, int M, int init_N, 
                                double rho, double box_size[3], double rcut, int max_N)
    : M(M), rho(rho), max_N(max_N)
{
    // 复制盒子大小
    for (int dim = 0; dim < 3; dim++) {
        this->box_size[dim] = box_size[dim];
    }
    this->rcut = rcut;
    
    // 初始化成员变量
    this->easy_cal_init_parameters();
    this->build_memory_r();
    this->read_init_configure(configuration);
    
    // 计算格子数和实际截断半径
    for (int dim = 0; dim < 3; dim++) {
        this->cell_num[dim] = static_cast<int>(floor(this->box_size[dim] / this->rcut));
        this->real_rcut[dim] = this->box_size[dim] / this->cell_num[dim];
    }
    
    this->build_cell_list();
    this->reset_mc_record();
    this->set_sim_parameters(0.1, 0.3, 10); // 默认参数
    
    std::cout << "Bulk Ring Polymer System constructed" << std::endl;
}

Bulk_RingPolymer::~Bulk_RingPolymer() {
    // 释放内存
    for (int i = 0; i < this->max_N * this->M; i++) {
        delete[] this->r_total[i];
        delete[] this->r_init[i];
        delete[] this->r_temp_volumn[i];
    }
    delete[] this->r_total;
    delete[] this->r_init;
    delete[] this->r_temp_volumn;
    
    std::cout << "Bulk Ring Polymer System destroyed." << std::endl;
}

void Bulk_RingPolymer::build_memory_r() {
    // 分配内存
    this->r_total = new double*[this->max_N * this->M];
    this->r_init = new double*[this->max_N * this->M];
    this->r_temp_volumn = new double*[this->max_N * this->M];
    
    for (int i = 0; i < this->max_N * this->M; i++) {
        this->r_total[i] = new double[3];
        this->r_init[i] = new double[3];
        this->r_temp_volumn[i] = new double[3];
    }
}

void Bulk_RingPolymer::easy_cal_init_parameters() {
    // 初始化随机数生成器
    this->seed = std::chrono::system_clock::now().time_since_epoch().count();
    this->gen = std::mt19937(this->seed);
    this->real_dis = std::uniform_real_distribution<double>(-1.0, 1.0);
    this->uni_dis = std::uniform_real_distribution<double>(0.0, 1.0);
    
    // 计算系统状态
    this->N_now = 0; // 初始化为0，read_init_configure会更新
    this->MN_now = 0;
    this->rho_now = 0.0;
    
    std::cout << "Initial parameters calculated." << std::endl;
}

void Bulk_RingPolymer::read_init_configure(std::string fileaddress) {
    std::ifstream infile(fileaddress);
    if (!infile) {
        throw std::runtime_error("Cannot open configuration file: " + fileaddress);
    }
    
    int row = 0;
    std::string line;
    while (std::getline(infile, line) && row < this->max_N * this->M) {
        if (line.empty())
            continue;
        
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
            throw std::runtime_error("Invalid configuration file format.");
        }
        
            // 直接存储位置，不应用周期性边界映射
        this->r_total[row][0] = x;
        this->r_total[row][1] = y;
        this->r_total[row][2] = z;
        
        this->r_init[row][0] = x;
        this->r_init[row][1] = y;
        this->r_init[row][2] = z;
        
        row++;
    }
    
    infile.close();
    
    // 更新系统状态
    this->MN_now = row;
    this->N_now = row / this->M;
    this->rho_now = static_cast<double>(this->MN_now) / (this->box_size[0] * this->box_size[1] * this->box_size[2]);
    
    std::cout << "Initial configuration read. MN_now = " << this->MN_now << ", N_now = " << this->N_now << std::endl;
}

// 边界条件处理方法
double Bulk_RingPolymer::dis_period(double x1, double x2, double box_size) const {
    double dx = fabs(x1 - x2);
    return std::min(dx, box_size - dx);
}

// 计算两点之间的非周期性距离平方
double Bulk_RingPolymer::dis2_no_period(const double* r1, const double* r2) const {
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];
    return dx * dx + dy * dy + dz * dz;
}

// 格子列表相关方法
int Bulk_RingPolymer::get_cell_index(const double* pos) {
    int ix_raw = static_cast<int>(floor(pos[0] / this->real_rcut[0]));
    int iy_raw = static_cast<int>(floor(pos[1] / this->real_rcut[1]));
    int iz_raw = static_cast<int>(floor(pos[2] / this->real_rcut[2]));
    
    int ix = (ix_raw % this->cell_num[0] + this->cell_num[0]) % this->cell_num[0];
    int iy = (iy_raw % this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];
    int iz = (iz_raw % this->cell_num[2] + this->cell_num[2]) % this->cell_num[2];
    
    return ix + iy * this->cell_num[0] + iz * this->cell_num[0] * this->cell_num[1];
}

void Bulk_RingPolymer::build_cell_list() {
    this->cl_total_cells = this->cell_num[0] * this->cell_num[1] * this->cell_num[2];
    
    if (this->cl_head.size() != this->cl_total_cells) {
        this->cl_head.resize(this->cl_total_cells);
    }
    
    if (this->cl_list.size() != this->max_N * this->M) {
        this->cl_list.resize(this->max_N * this->M);
    }
    
    // 清空格子头
    std::fill(this->cl_head.begin(), this->cl_head.end(), -1);
    
    // 填充格子列表
    for (int i = 0; i < this->MN_now; i++) {
        int cell_idx = this->get_cell_index(this->r_total[i]);
        this->cl_list[i] = this->cl_head[cell_idx];
        this->cl_head[cell_idx] = i;
    }
}

// 拓扑构建
void Bulk_RingPolymer::build_topology_map() {
    this->topology_map.resize(this->M);
    for (int i = 0; i < this->M; i++) {
        this->topology_map[i].clear();
        
        // 环状拓扑：每个单体连接到前一个和后一个单体
        int prev = (i - 1 + this->M) % this->M;
        int next = (i + 1) % this->M;
        
        this->topology_map[i].push_back(prev);
        this->topology_map[i].push_back(next);
    }
    std::cout << "Ring topology map constructed." << std::endl;
}

// 初始化第二步
void Bulk_RingPolymer::init_second() {
    this->build_topology_map();
    
    if (this->check_configure_validity()) {
        std::cout << "Initial configuration is valid." << std::endl;
    } else {
        throw std::runtime_error("Initial configuration is invalid.");
    }
    
    this->print_all_parameters();
}

// 设置模拟参数
void Bulk_RingPolymer::set_sim_parameters(double EPS_TRANS, double ROT_RATIO, int K_MAX) {
    this->eps_trans = EPS_TRANS;
    this->rot_ratio = ROT_RATIO;
    this->k_max = K_MAX;
}

// 打印所有参数
void Bulk_RingPolymer::print_all_parameters() {
    std::cout << "------------------ Bulk Ring Polymer Parameters ------------------" << std::endl;
    std::cout << "Polymerization degree M: " << this->M << std::endl;
    std::cout << "Current polymer number N_now: " << this->N_now << std::endl;
    std::cout << "Total monomer number MN_now: " << this->MN_now << std::endl;
    std::cout << "Target density rho: " << this->rho << std::endl;
    std::cout << "Current density rho_now: " << this->rho_now << std::endl;
    std::cout << "Box size: " << this->box_size[0] << "(x), " 
              << this->box_size[1] << "(y), " << this->box_size[2] << "(z)" << std::endl;
    std::cout << "Cutoff rcut: " << this->rcut << std::endl;
    std::cout << "Real rcut: " << this->real_rcut[0] << "(x), " 
              << this->real_rcut[1] << "(y), " << this->real_rcut[2] << "(z)" << std::endl;
    std::cout << "Cell number: " << this->cell_num[0] << "(x), " 
              << this->cell_num[1] << "(y), " << this->cell_num[2] << "(z)" << std::endl;
    std::cout << "eps_trans: " << this->eps_trans << std::endl;
    std::cout << "rot_ratio: " << this->rot_ratio << std::endl;
    std::cout << "k_max: " << this->k_max << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
}

// 重置MC记录
void Bulk_RingPolymer::reset_mc_record() {
    this->acc_trans = 0;
    this->acc_rot = 0;
    this->num_trans = 0;
    this->num_rot = 0;
}

// 结束一个block
void Bulk_RingPolymer::end_block(int block) {
    std::cout << "Block " << block << " completed." << std::endl;
    std::cout << "Translation acceptance ratio: " << static_cast<double>(this->acc_trans) / this->num_trans << std::endl;
    std::cout << "Rotation acceptance ratio: " << static_cast<double>(this->acc_rot) / this->num_rot << std::endl;
    this->reset_mc_record();
}

// 碰撞检测方法
bool Bulk_RingPolymer::overlap_other_monomer_one(const double* r_other, const double* r_try) const {
    double dx = this->dis_period(r_other[0], r_try[0], this->box_size[0]);
    double dy = this->dis_period(r_other[1], r_try[1], this->box_size[1]);
    double dz = this->dis_period(r_other[2], r_try[2], this->box_size[2]);
    
    double dist_sq = dx * dx + dy * dy + dz * dz;
    return dist_sq <= diameter2;
}

bool Bulk_RingPolymer::overlap_all_other_polymer(const double* r_try) {
    // 使用格子列表进行高效碰撞检测
    if (this->MN_now > this->cl_threshold) {
        return this->check_collision_except_polymer(r_try, -1);
    } else {
        // 直接检测所有单体
        for (int i = 0; i < this->MN_now; i++) {
            if (this->overlap_other_monomer_one(this->r_total[i], r_try)) {
                return true;
            }
        }
        return false;
    }
}

bool Bulk_RingPolymer::check_collision_except_polymer(const double* r_try, int ignore_polymer_index) {
    // 计算尝试位置所在的格子
    int cx = static_cast<int>(floor(r_try[0] / this->real_rcut[0]));
    int cy = static_cast<int>(floor(r_try[1] / this->real_rcut[1]));
    int cz = static_cast<int>(floor(r_try[2] / this->real_rcut[2]));
    
    // 遍历相邻格子
    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                // 计算实际格子索引，考虑周期性
                int real_cx = (cx + dx + this->cell_num[0]) % this->cell_num[0];
                int real_cy = (cy + dy + this->cell_num[1]) % this->cell_num[1];
                int real_cz = (cz + dz + this->cell_num[2]) % this->cell_num[2];
                
                int cell_idx = real_cx + real_cy * this->cell_num[0] + real_cz * this->cell_num[0] * this->cell_num[1];
                
                // 遍历格子中的所有单体
                int p_id = this->cl_head[cell_idx];
                while (p_id != -1 && p_id < this->MN_now) {
                    // 检查是否属于要忽略的聚合物
                    if (ignore_polymer_index >= 0) {
                        int polymer_id = p_id / this->M;
                        if (polymer_id == ignore_polymer_index) {
                            p_id = this->cl_list[p_id];
                            continue;
                        }
                    }
                    
                    // 检查碰撞
                    if (this->overlap_other_monomer_one(this->r_total[p_id], r_try)) {
                        return true;
                    }
                    
                    p_id = this->cl_list[p_id];
                }
            }
        }
    }
    return false;
}

bool Bulk_RingPolymer::check_configure_validity() {
    // 检查所有单体对是否重叠
    for (int i = 0; i < this->MN_now - 1; i++) {
        for (int j = i + 1; j < this->MN_now; j++) {
            // 跳过相邻单体
            int polymer_i = i / this->M;
            int polymer_j = j / this->M;
            int local_i = i % this->M;
            int local_j = j % this->M;
            
            bool is_neighbor = false;
            for (int neighbor : this->topology_map[local_i]) {
                if (neighbor == local_j) {
                    is_neighbor = true;
                    break;
                }
            }
            
            if (polymer_i == polymer_j && is_neighbor) {
                continue; // 跳过相邻单体
            }
            
            if (this->overlap_other_monomer_one(this->r_total[i], this->r_total[j])) {
                std::cout << "Overlap detected between monomer " << i << " and " << j << std::endl;
                return false;
            }
        }
    }
    return true;
}

// MC移动方法
void Bulk_RingPolymer::trans_move(int polymer_index) {
    this->num_trans++;
    
    // 生成随机位移
    double delta_move[3];
    for (int j = 0; j < 3; j++) {
        delta_move[j] = this->real_dis(this->gen) * this->box_size[j] * this->eps_trans;
    }
    
    // 检查所有单体是否重叠
    bool overlap = false;
    for (int monomer_index = polymer_index * this->M; monomer_index < (polymer_index + 1) * this->M; monomer_index++) {
        double r_try[3];
        for (int j = 0; j < 3; j++) {
            r_try[j] = this->r_total[monomer_index][j] + delta_move[j];
        }
        
        if (this->check_collision_except_polymer(r_try, polymer_index)) {
            overlap = true;
            break;
        }
    }
    
    // 如果没有重叠，接受移动
    if (!overlap) {
        this->acc_trans++;
        for (int monomer_index = polymer_index * this->M; monomer_index < (polymer_index + 1) * this->M; monomer_index++) {
            for (int j = 0; j < 3; j++) {
                this->r_total[monomer_index][j] += delta_move[j];
            }
        }
        
        // 更新格子列表
        this->build_cell_list();
    }
}

void Bulk_RingPolymer::rot_mid_move(int monomer_index, int polymer_index) {
    this->num_rot++;
    int real_monomer_index = polymer_index * this->M + monomer_index;
    
    // 获取前后单体位置
    double r_prev[3];
    double r_next[3];
    
    int last_monomer = polymer_index * this->M + this->topology_map[monomer_index][0];
    int next_monomer = polymer_index * this->M + this->topology_map[monomer_index][1];
    
    for (int j = 0; j < 3; j++) {
        r_prev[j] = this->r_total[last_monomer][j];
        r_next[j] = this->r_total[next_monomer][j];
    }
    
    // 生成候选位置
    std::vector<std::array<double, 3>> r_temp_new_data(this->k_max);
    std::vector<std::array<double, 3>> r_temp_old_data(this->k_max - 1);
    
    // 创建指针数组
    std::vector<double*> r_temp_new_ptrs(this->k_max);
    std::vector<double*> r_temp_old_ptrs(this->k_max - 1);
    for (int i = 0; i < this->k_max; i++) {
        r_temp_new_ptrs[i] = r_temp_new_data[i].data();
    }
    for (int i = 0; i < this->k_max - 1; i++) {
        r_temp_old_ptrs[i] = r_temp_old_data[i].data();
    }
    
    // 生成候选位置
    this->Find_Last(r_prev, r_next, r_temp_new_ptrs.data(), this->k_max);
    this->Find_Last(r_prev, r_next, r_temp_old_ptrs.data(), this->k_max - 1);
    
    // 合并候选位置：前k_max个为新位置，后k_max-1个为旧位置
    std::vector<std::array<double, 3>> r_temp(2 * this->k_max - 1);
    for (int k = 0; k < this->k_max; k++) {
        r_temp[k] = r_temp_new_data[k];
    }
    for (int k = 0; k < this->k_max - 1; k++) {
        r_temp[k + this->k_max] = r_temp_old_data[k];
    }
    
    // 计算权重
    std::vector<double> w(2 * this->k_max - 1, 1.0);
    for (int k = 0; k < 2 * this->k_max - 1; k++) {
        // 检查碰撞，仅排除当前正在旋转的单体
        if (this->check_collision_except_monomer(r_temp[k].data(), real_monomer_index)) {
            w[k] = 0.0;
        }
    }
    
    // 检查是否所有位置都重叠
    if (std::all_of(w.begin(), w.end(), [](double v) { return v == 0.0; })) {
        return; // 所有候选位置都重叠，旋转失败
    }
    
    // 轮盘赌选择
    double sum_w = 0.0;
    for (double weight : w) {
        sum_w += weight;
    }
    
    double rand_val = this->uni_dis(this->gen) * sum_w;
    double current_sum = 0.0;
    int select_index = -1;
    
    for (int k = 0; k < 2 * this->k_max - 1; k++) {
        current_sum += w[k];
        if (current_sum >= rand_val) {
            select_index = k;
            break;
        }
    }
    
    // 更新位置
    if (select_index < this->k_max) {
        // 接受新位置
        this->acc_rot++;
        for (int j = 0; j < 3; j++) {
            this->r_total[real_monomer_index][j] = r_temp[select_index][j];
        }
        
        // 更新格子列表
        this->build_cell_list();
    }
}

void Bulk_RingPolymer::rot_polymer_move(int polymer_index) {
    // 决定哪些单体需要旋转
    std::vector<bool> rot_list(this->M, false);
    for (int i = 0; i < this->M; i++) {
        if (this->uni_dis(this->gen) < this->rot_ratio) {
            rot_list[i] = true;
        }
    }
    
    // 对选定的单体进行旋转
    for (int i = 0; i < this->M; i++) {
        if (rot_list[i]) {
            // 环状聚合物没有端点，所有单体都是中间单体
            this->rot_mid_move(i, polymer_index);
        }
    }
}

// 单步MC模拟
void Bulk_RingPolymer::mc_one_step() {
    // 对每个聚合物进行平移和旋转
    for (int polymer_index = 0; polymer_index < this->N_now; polymer_index++) {
        this->trans_move(polymer_index);
        this->rot_polymer_move(polymer_index);
    }
}

// 运行多步模拟
void Bulk_RingPolymer::run_simulation(int steps) {
    std::cout << "Running simulation for " << steps << " steps..." << std::endl;
    for (int step = 0; step < steps; step++) {
        this->mc_one_step();
        
        // 每1000步打印一次状态
        if (step % 1000 == 0) {
            std::cout << "Step " << step << ", N_now = " << this->N_now << ", MN_now = " << this->MN_now << std::endl;
        }
    }
    std::cout << "Simulation completed." << std::endl;
}

// 未实现的方法，后续可以扩展
bool Bulk_RingPolymer::check_collision_except_monomer(const double* r_try, int ignore_monomer_index) {
    // 计算尝试位置所在的格子
    int cx = static_cast<int>(floor(r_try[0] / this->real_rcut[0]));
    int cy = static_cast<int>(floor(r_try[1] / this->real_rcut[1]));
    int cz = static_cast<int>(floor(r_try[2] / this->real_rcut[2]));
    
    // === 优化核心：预计算需要忽略的全局 ID 列表 ===
    // 使用一个固定大小的数组存储“黑名单”，因为邻居通常很少
    int ignore_list[10];  // 最多忽略 10 个邻居
    int ignore_count = 0;
    
    if (ignore_monomer_index != -1) {
        // 1. 添加自身到忽略列表
        ignore_list[ignore_count] = ignore_monomer_index;
        ignore_count++;
        
        // 2. 计算所在的聚合物基准 ID 和局部索引
        int local_rank = ignore_monomer_index % this->M;
        int polymer_base = ignore_monomer_index - local_rank;
        
        // 3. 查表获取所有邻居，并添加到忽略列表
        // topology_map[local_rank] 返回的是局部邻居索引
        // 我们加上 polymer_base 转换回全局 ID
        for (int local_neighbor : this->topology_map[local_rank]) {
            ignore_list[ignore_count] = polymer_base + local_neighbor;
            ignore_count++;
        }
    }
    
    // 遍历相邻格子
    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                // 计算实际格子索引，考虑周期性
                int real_cx = (cx + dx + this->cell_num[0]) % this->cell_num[0];
                int real_cy = (cy + dy + this->cell_num[1]) % this->cell_num[1];
                int real_cz = (cz + dz + this->cell_num[2]) % this->cell_num[2];
                
                int cell_idx = real_cx + real_cy * this->cell_num[0] + real_cz * this->cell_num[0] * this->cell_num[1];
                
                // 遍历格子中的所有单体
                int p_id = this->cl_head[cell_idx];
                while (p_id != -1 && p_id < this->MN_now) {
                    // 检查是否在忽略列表中
                    bool is_ignored = false;
                    if (ignore_count > 0) {
                        // 循环检查忽略列表，因为列表很小，速度很快
                        for (int k = 0; k < ignore_count; k++) {
                            if (p_id == ignore_list[k]) {
                                is_ignored = true;
                                break;
                            }
                        }
                    }
                    
                    if (is_ignored) {
                        p_id = this->cl_list[p_id];
                        continue;
                    }
                    
                    // 检查碰撞
                    if (this->overlap_other_monomer_one(this->r_total[p_id], r_try)) {
                        return true;
                    }
                    
                    p_id = this->cl_list[p_id];
                }
            }
        }
    }
    return false;
}

// 辅助函数
void Bulk_RingPolymer::Rot_temp(double* vec, double** r_temp_rot, int k_max) {
    // 计算旋转角度
    double phi = std::acos(std::max(-1.0, std::min(1.0, vec[2])));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    
    if (std::fabs(sin_phi) < 1e-12) {
        sin_phi = 1e-12;
    }
    double theta = std::asin(std::max(-1.0, std::min(1.0, vec[1] / sin_phi)));
    
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    
    // 旋转矩阵
    double R00 = cos_theta * cos_phi;
    double R01 = -sin_theta;
    double R02 = cos_theta * sin_phi;
    
    double R10 = sin_theta * cos_phi;
    double R11 = cos_theta;
    double R12 = sin_theta * sin_phi;
    
    double R20 = -sin_phi;
    double R21 = 0.0;
    double R22 = cos_phi;
    
    // 应用旋转矩阵到所有候选位置
    for (int i = 0; i < k_max; i++) {
        double temp[3] = {r_temp_rot[i][0], r_temp_rot[i][1], r_temp_rot[i][2]};
        r_temp_rot[i][0] = R00 * temp[0] + R01 * temp[1] + R02 * temp[2];
        r_temp_rot[i][1] = R10 * temp[0] + R11 * temp[1] + R12 * temp[2];
        r_temp_rot[i][2] = R20 * temp[0] + R21 * temp[1] + R22 * temp[2];
    }
}

void Bulk_RingPolymer::Find_Last(double* r1, double* r2, double** r_temp_rot, int k_max) {
    // 生成k_max个随机角度
    std::vector<double> temp_theta(k_max);
    for (int i = 0; i < k_max; i++) {
        temp_theta[i] = 2 * M_PI * this->uni_dis(this->gen);
    }
    
    // 计算中间点和距离
    double temp = diameter2 - this->dis2_no_period(r1, r2) / 4.0;
    
    double mid_r[3] = {0.0};
    for (int i = 0; i < 3; i++) {
        mid_r[i] = (r1[i] + r2[i]) / 2;
    }
    
    // 检查是否可以生成候选位置
    if (temp <= 0) {
        // 无法生成候选位置，返回中间点
        for (int k = 0; k < k_max; k++) {
            for (int j = 0; j < 3; j++) {
                r_temp_rot[k][j] = mid_r[j];
            }
        }
        return;
    }
    
    // 生成候选位置
    double R = sqrt(temp);
    for (int i = 0; i < k_max; i++) {
        r_temp_rot[i][0] = R * cos(temp_theta[i]);
        r_temp_rot[i][1] = R * sin(temp_theta[i]);
        r_temp_rot[i][2] = 0.0; // 初始z方向为0
    }
    
    // 计算旋转向量
    double vec[3] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};
    double vec_length = sqrt(this->dis2_no_period(vec, vec));
    
    // 归一化向量
    if (vec[0] >= 0) {
        vec[0] /= vec_length;
        vec[1] /= vec_length;
        vec[2] /= vec_length;
    } else {
        vec[0] /= -vec_length;
        vec[1] /= -vec_length;
        vec[2] /= -vec_length;
    }
    
    // 旋转候选位置
    this->Rot_temp(vec, r_temp_rot, k_max);
    
    // 平移到中间点
    for (int i = 0; i < k_max; i++) {
        for (int j = 0; j < 3; j++) {
            r_temp_rot[i][j] += mid_r[j];
        }
    }
}

// 体积缩放并检查重叠
int Bulk_RingPolymer::change_volumn(int expan_or_shrink, double mc_cond_eps_volumn) {
    // 计算聚合物质心
    std::vector<std::array<double, 3>> r_center(this->N_now);
    for (int i = 0; i < this->N_now; i++) {
        r_center[i][0] = 0.0;
        r_center[i][1] = 0.0;
        r_center[i][2] = 0.0;
        
        for (int par_index = i * this->M; par_index < (i + 1) * this->M; par_index++) {
            r_center[i][0] += this->r_total[par_index][0];
            r_center[i][1] += this->r_total[par_index][1];
            r_center[i][2] += this->r_total[par_index][2];
        }
        
        r_center[i][0] /= this->M;
        r_center[i][1] /= this->M;
        r_center[i][2] /= this->M;
    }
    
    double delta_move[3] = {0.0};
    double temp_eps_volumn = mc_cond_eps_volumn * expan_or_shrink;
    double new_box_size[3];
    for (int dim = 0; dim < 3; dim++) {
        new_box_size[dim] = this->box_size[dim] / (1 + temp_eps_volumn);
    }
    

    // 计算新旧盒子中心
    double old_box_center[3] = {this->box_size[0] / 2, this->box_size[1] / 2, this->box_size[2] / 2};
    double new_box_center[3] = {new_box_size[0] / 2, new_box_size[1] / 2, new_box_size[2] / 2};
    
    // 更新位置
    for (int i = 0; i < this->N_now; i++) {
        // 移动聚合物并映射到新的盒子范围内
        for (int par_index = i * this->M; par_index < (i + 1) * this->M; par_index++) {
            for (int dim = 0; dim < 3; dim++) {
                // 正确的计算逻辑：
                // 1. 计算相对于旧盒子中心的位置
                double rel_pos = this->r_total[par_index][dim] - old_box_center[dim];
                // 2. 缩放这个相对位置
                double scaled_rel_pos = rel_pos / (1 + temp_eps_volumn);
                // 3. 计算相对于新盒子中心的绝对位置
                this->r_temp_volumn[par_index][dim] = new_box_center[dim] + scaled_rel_pos;
                // 4. 将位置映射到新的盒子范围内 [0, new_box_size]
                this->r_temp_volumn[par_index][dim] = fmod(this->r_temp_volumn[par_index][dim], new_box_size[dim]);
                if (this->r_temp_volumn[par_index][dim] < 0) {
                    this->r_temp_volumn[par_index][dim] += new_box_size[dim];
                }
            }
        }
    }
    
    // 基于新的盒子大小计算临时格子参数
    double temp_real_rcut[3];
    int temp_cell_num[3];
    int temp_cl_total_cells;
    
    for (int dim = 0; dim < 3; dim++) {
        temp_cell_num[dim] = static_cast<int>(floor(new_box_size[dim] / this->rcut));
        if (temp_cell_num[dim] < 1) {
            temp_cell_num[dim] = 1;
        }
        temp_real_rcut[dim] = new_box_size[dim] / temp_cell_num[dim];
    }
    
    temp_cl_total_cells = temp_cell_num[0] * temp_cell_num[1] * temp_cell_num[2];
    
    // 构建临时格子列表
    std::vector<int> temp_cl_head(temp_cl_total_cells, -1);
    std::vector<int> temp_cl_list(this->MN_now, -1);
    
    // 临时格子索引计算函数
    auto get_temp_cell_index = [&](const double* pos) -> int {
        int ix_raw = static_cast<int>(floor(pos[0] / temp_real_rcut[0]));
        int iy_raw = static_cast<int>(floor(pos[1] / temp_real_rcut[1]));
        int iz_raw = static_cast<int>(floor(pos[2] / temp_real_rcut[2]));
        
        int ix = (ix_raw % temp_cell_num[0] + temp_cell_num[0]) % temp_cell_num[0];
        int iy = (iy_raw % temp_cell_num[1] + temp_cell_num[1]) % temp_cell_num[1];
        int iz = (iz_raw % temp_cell_num[2] + temp_cell_num[2]) % temp_cell_num[2];
        
        return ix + iy * temp_cell_num[0] + iz * temp_cell_num[0] * temp_cell_num[1];
    };
    
    for (int par_index = 0; par_index < this->MN_now; par_index++) {
        int cell_idx = get_temp_cell_index(this->r_temp_volumn[par_index]);
        
        // 链表插入
        temp_cl_list[par_index] = temp_cl_head[cell_idx];
        temp_cl_head[cell_idx] = par_index;
    }
    
    // 检查是否重叠
    for (int par_index = 0; par_index < this->MN_now; par_index++) {
        // 计算当前单体所在的格子
        int cell_idx = get_temp_cell_index(this->r_temp_volumn[par_index]);
        int cx = cell_idx % temp_cell_num[0];
        int cy = (cell_idx / temp_cell_num[0]) % temp_cell_num[1];
        int cz = cell_idx / (temp_cell_num[0] * temp_cell_num[1]);
        
        // 遍历相邻格子
        for (int dz = -1; dz <= 1; dz++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    // 计算实际格子索引，考虑周期性
                    int real_cx = (cx + dx + temp_cell_num[0]) % temp_cell_num[0];
                    int real_cy = (cy + dy + temp_cell_num[1]) % temp_cell_num[1];
                    int real_cz = (cz + dz + temp_cell_num[2]) % temp_cell_num[2];
                    
                    int neighbor_cell_idx = real_cx + real_cy * temp_cell_num[0] + real_cz * temp_cell_num[0] * temp_cell_num[1];
                    
                    // 遍历格子中的所有单体
                    int p_id = temp_cl_head[neighbor_cell_idx];
                    while (p_id != -1) {
                        // 跳过同一个聚合物的单体
                        if (p_id / this->M == par_index / this->M) {
                            p_id = temp_cl_list[p_id];
                            continue;
                        }
                        
                        // 跳过自己
                        if (p_id == par_index) {
                            p_id = temp_cl_list[p_id];
                            continue;
                        }
                        
                        // 检查碰撞
                        // 位置已经在构建格子列表之前被映射到了新的盒子范围内
                        double dx_dist = this->dis_period(this->r_temp_volumn[par_index][0], this->r_temp_volumn[p_id][0], new_box_size[0]);
                        double dy_dist = this->dis_period(this->r_temp_volumn[par_index][1], this->r_temp_volumn[p_id][1], new_box_size[1]);
                        double dz_dist = this->dis_period(this->r_temp_volumn[par_index][2], this->r_temp_volumn[p_id][2], new_box_size[2]);
                        
                        double temp_dis = dx_dist * dx_dist + dy_dist * dy_dist + dz_dist * dz_dist;
                        if (temp_dis <= diameter2) {
                            return 1; // 重叠了
                        }
                        
                        p_id = temp_cl_list[p_id];
                    }
                }
            }
        }
    }
    
    return 0; // 没有重叠
}

// 插入高分子并计算插入权重
double Bulk_RingPolymer::calculate_insertion_weight(int k_max) {
    // 1. 对第一个粒子使用CBMC采样，生成k_max个候选位置
    std::vector<std::array<double, 3>> r_first_candidates(k_max);
    std::vector<double> weights(k_max, 1.0);
    
    // 生成k_max个候选位置，范围在[0, box_size[i]]
    for (int k = 0; k < k_max; k++) {
        r_first_candidates[k][0] = this->box_size[0] * this->uni_dis(this->gen);
        r_first_candidates[k][1] = this->box_size[1] * this->uni_dis(this->gen);
        r_first_candidates[k][2] = this->box_size[2] * this->uni_dis(this->gen);
        
        // 计算权重：如果重叠则权重为0，否则为1
        if (this->check_collision_except_polymer(r_first_candidates[k].data(), -1)) {
            weights[k] = 0.0;
        }
    }
    
    // 计算总权重
    double sum_weights = 0.0;
    for (int k = 0; k < k_max; k++) {
        sum_weights += weights[k];
    }
    
    if (sum_weights == 0.0) {
        return 0.0; // 所有候选位置都重叠，返回0
    }
    
    // 2. 计算插入权重
    // 对于体相系统，插入权重等于第一个粒子的总权重除以k_max
    // 因为每个候选位置的概率相等，都是1/k_max
    double insertion_weight = sum_weights / k_max;
    
    return insertion_weight;
}