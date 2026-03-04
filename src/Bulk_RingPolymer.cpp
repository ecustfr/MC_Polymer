#include "Bulk_RingPolymer.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <vector>

Bulk_RingPolymer::Bulk_RingPolymer(std::string configuration, int M, int init_N,
                                double rho, double box_size[3], double rcut, int max_N,
                                double mu_b, double rho_b)
    : M(M), rho(rho), max_N(max_N), mu_b(mu_b), rho_b(rho_b),
      exp_mu_b(std::exp(mu_b))
{
    // 复制盒子大小
    for (int dim = 0; dim < 3; dim++) {
        this->box_size[dim] = box_size[dim];
    }
    this->rcut = rcut;
    this->V = this->box_size[0] * this->box_size[1] * this->box_size[2];

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
    this->set_sim_parameters(0.1, 0.3, 5); // 默认参数
    
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

// 计算两点之间的三维周期性距离平方
double Bulk_RingPolymer::dis2_period(const double* r1, const double* r2) const {
    double dx = this->dis_period(r1[0], r2[0], this->box_size[0]);
    double dy = this->dis_period(r1[1], r2[1], this->box_size[1]);
    double dz = this->dis_period(r1[2], r2[2], this->box_size[2]);
    return dx*dx + dy*dy + dz*dz;
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

    // 预分配旋转移动缓存
    int total_candidates = 2 * this->k_max - 1;
    this->rot_mid_cache_data.resize(total_candidates);
    this->rot_mid_cache_ptrs.resize(total_candidates);
    this->rot_mid_cache_weights.resize(total_candidates);

    // 设置指针数组
    for (int i = 0; i < total_candidates; i++) {
        this->rot_mid_cache_ptrs[i] = this->rot_mid_cache_data[i].data();
    }
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
    this->acc_insert = 0;
    this->acc_delete = 0;
    this->num_trans = 0;
    this->num_rot = 0;
    this->num_insert = 0;
    this->num_delete = 0;
}

// 结束一个block
void Bulk_RingPolymer::end_block(int block) {
    std::cout << "Block " << block << " completed." << std::endl;
    std::cout << "Translation acceptance ratio: " << static_cast<double>(this->acc_trans) / this->num_trans << std::endl;
    std::cout << "Rotation acceptance ratio: " << static_cast<double>(this->acc_rot) / this->num_rot << std::endl;
    this->reset_mc_record();
}

// 获取插入接受率
double Bulk_RingPolymer::get_insert_acceptance() const {
    return (this->num_insert > 0) ? static_cast<double>(this->acc_insert) / this->num_insert : 0.0;
}

// 获取删除接受率
double Bulk_RingPolymer::get_delete_acceptance() const {
    return (this->num_delete > 0) ? static_cast<double>(this->acc_delete) / this->num_delete : 0.0;
}

// 获取平移接受率
double Bulk_RingPolymer::get_translation_acceptance() const {
    return (this->num_trans > 0) ? static_cast<double>(this->acc_trans) / this->num_trans : 0.0;
}

// 获取旋转接受率
double Bulk_RingPolymer::get_rotation_acceptance() const {
    return (this->num_rot > 0) ? static_cast<double>(this->acc_rot) / this->num_rot : 0.0;
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
                
                int real_cx = ((cx+dx)%this->cell_num[0] + this->cell_num[0]) % this->cell_num[0]; 
                int real_cy = ((cy+dy)%this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];
                int real_cz = ((cz+dz)%this->cell_num[2] + this->cell_num[2]) % this->cell_num[2];
                
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

    // 距离检查配置（默认关闭以优化性能）
    constexpr bool ENABLE_DISTANCE_CHECK = false;
    constexpr double MIN_ALLOWED = 1.0 - 0.001;
    constexpr double MAX_ALLOWED = 1.0 + 0.001;
    constexpr double MIN_ALLOWED_SQ = MIN_ALLOWED * MIN_ALLOWED;
    constexpr double MAX_ALLOWED_SQ = MAX_ALLOWED * MAX_ALLOWED;

    // 获取前后单体位置
    double r_prev[3];
    double r_next[3];

    int last_monomer = polymer_index * this->M + this->topology_map[monomer_index][0];
    int next_monomer = polymer_index * this->M + this->topology_map[monomer_index][1];

    // 直接复制位置，避免循环开销
    r_prev[0] = this->r_total[last_monomer][0];
    r_prev[1] = this->r_total[last_monomer][1];
    r_prev[2] = this->r_total[last_monomer][2];
    r_next[0] = this->r_total[next_monomer][0];
    r_next[1] = this->r_total[next_monomer][1];
    r_next[2] = this->r_total[next_monomer][2];

    // 使用预分配的缓存
    int total_candidates = 2 * this->k_max - 1;

    // 一次性生成所有候选位置
    if (!this->Find_Last(r_prev, r_next, this->rot_mid_cache_ptrs.data(), total_candidates)) {
        throw std::runtime_error("Find_Last failed in rot_mid_move.");
    }

    // 计算权重 w = exp(-\beta u)，如果不重叠则为1.0（体相系统无外势），重叠则为0
    // 直接使用缓存权重数组
    std::fill(this->rot_mid_cache_weights.begin(), this->rot_mid_cache_weights.end(), 1.0);

    // 检查碰撞并更新权重
    bool all_overlap = true;
    for (int k = 0; k < total_candidates; k++) {
        if (!this->check_collision_except_monomer(this->rot_mid_cache_data[k].data(), real_monomer_index)) {
            all_overlap = false;  // 至少有一个位置不重叠
        } else {
            this->rot_mid_cache_weights[k] = 0.0;
        }
    }

    if (all_overlap) {
        return; // 所有候选位置都重叠，旋转失败
    }

    // 使用 bias 采样判断是否接受（旧权重为1，因为当前是等概率采样）
    // 体相系统无外势，旧位置权重为1.0
    double old_weight = 1.0;
    if (acc_bias_or_not(this->rot_mid_cache_weights.data(), old_weight, this->k_max, this->gen)) {
        this->acc_rot++;
        // 从新位置中选择一个（使用前 k_max 个位置和权重）
        int select_index = RouteWheel(this->rot_mid_cache_weights.data(), this->k_max, this->gen);

        // 直接更新位置
        this->r_total[real_monomer_index][0] = this->rot_mid_cache_data[select_index][0];
        this->r_total[real_monomer_index][1] = this->rot_mid_cache_data[select_index][1];
        this->r_total[real_monomer_index][2] = this->rot_mid_cache_data[select_index][2];

        // 可选的距离检查（默认关闭）
        if (ENABLE_DISTANCE_CHECK) {
            double dx_prev = this->r_total[real_monomer_index][0] - r_prev[0];
            double dy_prev = this->r_total[real_monomer_index][1] - r_prev[1];
            double dz_prev = this->r_total[real_monomer_index][2] - r_prev[2];
            double dist2_prev = dx_prev*dx_prev + dy_prev*dy_prev + dz_prev*dz_prev;

            double dx_next = this->r_total[real_monomer_index][0] - r_next[0];
            double dy_next = this->r_total[real_monomer_index][1] - r_next[1];
            double dz_next = this->r_total[real_monomer_index][2] - r_next[2];
            double dist2_next = dx_next*dx_next + dy_next*dy_next + dz_next*dz_next;

            if (dist2_prev < MIN_ALLOWED_SQ || dist2_prev > MAX_ALLOWED_SQ ||
                dist2_next < MIN_ALLOWED_SQ || dist2_next > MAX_ALLOWED_SQ) {
                std::cerr << "WARNING: Invalid bond length after rotation move." << std::endl;
                std::cerr << "  monomer_index: " << monomer_index << ", polymer_index: " << polymer_index << std::endl;
                std::cerr << "  dist2_prev: " << dist2_prev << " (" << sqrt(dist2_prev) << ")" << std::endl;
                std::cerr << "  dist2_next: " << dist2_next << " (" << sqrt(dist2_next) << ")" << std::endl;
            }
        }

        // 更新格子列表（如果需要）
        if (this->MN_now > this->cl_threshold) {
            this->build_cell_list();
        }
    }
}

void Bulk_RingPolymer::rot_polymer_move(int polymer_index) {
    // 对每个单体独立决定是否旋转，避免额外的vector分配
    for (int i = 0; i < this->M; i++) {
        if (this->uni_dis(this->gen) < this->rot_ratio) {
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
        // this->mc_one_step();
        this->mc_one_step_MuVT();
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
                int real_cx = ((cx+dx)%this->cell_num[0] + this->cell_num[0]) % this->cell_num[0]; 
                int real_cy = ((cy+dy)%this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];
                int real_cz = ((cz+dz)%this->cell_num[2] + this->cell_num[2]) % this->cell_num[2];
                
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

bool Bulk_RingPolymer::Find_Last(double* r1, double* r2, double** r_temp_rot, int k_max) {
    // 计算向量 v = r2 - r1 和中间点 M
    double v[3] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};
    double d2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; // |v|^2

    // 提前检查：如果 d2 > 4.0，则 d > 2.0，无法生成候选位置
    if (d2 > 4.0) {
        double M[3] = {(r1[0] + r2[0]) / 2, (r1[1] + r2[1]) / 2, (r1[2] + r2[2]) / 2};
        for (int i = 0; i < k_max; i++) {
            r_temp_rot[i][0] = M[0];
            r_temp_rot[i][1] = M[1];
            r_temp_rot[i][2] = M[2];
        }
        return false;
    }

    double d = sqrt(d2);
    double M[3] = {(r1[0] + r2[0]) / 2, (r1[1] + r2[1]) / 2, (r1[2] + r2[2]) / 2};

    // 处理 d 过小的情况（r1 和 r2 几乎重合）
    const double EPSILON = 1e-12;
    if (d < EPSILON) {
        // 无法定义正交平面，返回中点
        for (int i = 0; i < k_max; i++) {
            r_temp_rot[i][0] = M[0];
            r_temp_rot[i][1] = M[1];
            r_temp_rot[i][2] = M[2];
        }
        return false;
    }

    // 计算半径 R = sqrt(diameter^2 - (d/2)^2)
    double R_sq = diameter2 - d2 / 4.0;
    double R = sqrt(R_sq);

    // 生成两个正交的单位向量 u, w 垂直于 v
    double u[3], w[3];

    // 方法：选择一个与 v 不正交的临时向量，通过叉积获得正交基
    // 避免数值问题：选择 v 中绝对值最大的分量
    double abs_x = fabs(v[0]), abs_y = fabs(v[1]), abs_z = fabs(v[2]);

    if (abs_x > abs_y) {
        if (abs_x > abs_z) {
            // v 主要沿 x 方向，使用 (0, 1, 0) 作为临时向量
            u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
        } else {
            // v 主要沿 z 方向，使用 (1, 0, 0) 作为临时向量
            u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
        }
    } else {
        if (abs_y > abs_z) {
            // v 主要沿 y 方向，使用 (0, 0, 1) 作为临时向量
            u[0] = 0.0; u[1] = 0.0; u[2] = 1.0;
        } else {
            // v 主要沿 z 方向，使用 (1, 0, 0) 作为临时向量
            u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
        }
    }

    // 计算 w = v × u
    w[0] = v[1]*u[2] - v[2]*u[1];
    w[1] = v[2]*u[0] - v[0]*u[2];
    w[2] = v[0]*u[1] - v[1]*u[0];

    // 归一化 w
    double w_len_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    double w_len = sqrt(w_len_sq);

    if (w_len < 1e-12) {
        // 如果 v 与 u 平行，尝试另一个 u
        u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
        w[0] = v[1]*u[2] - v[2]*u[1];
        w[1] = v[2]*u[0] - v[0]*u[2];
        w[2] = v[0]*u[1] - v[1]*u[0];
        w_len_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
        w_len = sqrt(w_len_sq);
    }

    // 归一化 w：w_unit = w / w_len
    double inv_w_len = 1.0 / w_len;
    w[0] *= inv_w_len; w[1] *= inv_w_len; w[2] *= inv_w_len;

    // 计算 u = w × v (确保正交右手系)
    // 注意：|w × v| = |w| * |v| * sin(90°) = 1 * d * 1 = d
    u[0] = w[1]*v[2] - w[2]*v[1];
    u[1] = w[2]*v[0] - w[0]*v[2];
    u[2] = w[0]*v[1] - w[1]*v[0];

    // 归一化 u：u_unit = u / d（因为 |u| = d）
    double inv_d = 1.0 / d;
    u[0] *= inv_d; u[1] *= inv_d; u[2] *= inv_d;

    // 预计算常数
    const double TWO_PI = 2.0 * M_PI;

    // 生成 k_max 个随机角度的候选位置
    // 使用局部变量减少数组索引开销
    double Mx = M[0], My = M[1], Mz = M[2];
    double ux = u[0], uy = u[1], uz = u[2];
    double wx = w[0], wy = w[1], wz = w[2];

    for (int i = 0; i < k_max; i++) {
        double theta = TWO_PI * this->uni_dis(this->gen);
        double cos_theta, sin_theta;
        // 使用 sincos 同时计算正弦和余弦（如果编译器支持）
        // 否则回退到单独的 cos 和 sin
#ifdef _GNU_SOURCE
        sincos(theta, &sin_theta, &cos_theta);
#else
        cos_theta = cos(theta);
        sin_theta = sin(theta);
#endif

        // 位置公式: B = M + R * (cosθ * u + sinθ * w)
        // 展开计算以减少临时变量
        double term_x = cos_theta * ux + sin_theta * wx;
        double term_y = cos_theta * uy + sin_theta * wy;
        double term_z = cos_theta * uz + sin_theta * wz;

        r_temp_rot[i][0] = Mx + R * term_x;
        r_temp_rot[i][1] = My + R * term_y;
        r_temp_rot[i][2] = Mz + R * term_z;
    }

    // 调试验证：检查第一个候选位置到 r1 和 r2 的距离（仅在调试时启用）
    bool debug_check = false;
    if (debug_check && k_max > 0) {
        double* B = r_temp_rot[0];
        double dx1 = B[0] - r1[0], dy1 = B[1] - r1[1], dz1 = B[2] - r1[2];
        double dx2 = B[0] - r2[0], dy2 = B[1] - r2[1], dz2 = B[2] - r2[2];
        double dist1_sq = dx1*dx1 + dy1*dy1 + dz1*dz1;
        double dist2_sq = dx2*dx2 + dy2*dy2 + dz2*dz2;
        double dist1 = sqrt(dist1_sq);
        double dist2 = sqrt(dist2_sq);

        if (fabs(dist1 - 1.0) > 0.001 || fabs(dist2 - 1.0) > 0.001 || fabs(dist1 - dist2) > 0.001) {
            std::cout << "DEBUG Find_Last validation failed:" << std::endl;
            std::cout << "  d = " << d << ", R = " << R << std::endl;
            std::cout << "  dist to r1 = " << dist1 << ", dist to r2 = " << dist2 << std::endl;
            std::cout << "  difference = " << fabs(dist1 - dist2) << std::endl;
            std::cout << "  v = [" << v[0] << ", " << v[1] << ", " << v[2] << "]" << std::endl;
            std::cout << "  u = [" << ux << ", " << uy << ", " << uz << "]" << std::endl;
            std::cout << "  w = [" << wx << ", " << wy << ", " << wz << "]" << std::endl;
            std::cout << "  M = [" << Mx << ", " << My << ", " << Mz << "]" << std::endl;
            std::cout << "  B = [" << B[0] << ", " << B[1] << ", " << B[2] << "]" << std::endl;
        }
    }

    return true;
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

double Bulk_RingPolymer::get_W_insert_ring(int k_max)
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

double Bulk_RingPolymer::get_W_insert_ring_z(double z, int k_max)
{
    std::vector<std::array<double, 3>> r_new(this->M);
    double Z_eff = 1.0;

    // 放置种子节点（第一个单体）
    for (int d = 0; d < 2; ++d)
    {
        r_new[0][d] = this->box_size[d] * this->uni_dis(this->gen);
    }
    r_new[0][2] = z;

    // 检查种子节点是否与其他聚合物重叠
    if (this->overlap_all_other_polymer(r_new[0].data()))
    {
        return 0.0; // 种子节点重叠，插入失败
    }

    // 插入剩余的单体（从索引1到M-1）
    for (int monomer_index = 1; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &Z_eff, monomer_index, k_max))
        {
            return 0.0; // 单体插入失败，返回 0
        }
    }

    return Z_eff;
}

// ======================================================================
// 巨正则系综模拟方法
// ======================================================================

// 巨正则单步模拟
void Bulk_RingPolymer::mc_one_step_MuVT() {
    this->mc_one_step();
    this->insert_move(this->k_max);
     int delete_N = static_cast<int>(this->N_now * this->uni_dis(this->gen));
    this->delete_move(this->k_max, delete_N);
}

// 巨正则多步模拟
void Bulk_RingPolymer::run_simulation_MuVT(int steps) {
    
    for (int step = 0; step < steps; step++) {
        this->mc_one_step_MuVT();
    }
    
}

// ======================================================================
// 插入移动方法
// ======================================================================

void Bulk_RingPolymer::insert_move(int k_max)
{
    this->num_insert++;

    // 1. 初始化新聚合物的位置
    std::vector<std::array<double, 3>> r_new(this->M);
    double Z_eff = 1.0;

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
        this->rho_now = static_cast<double>(this->MN_now) / (this->box_size[0] * this->box_size[1] * this->box_size[2]);

        // 更新cell列表
        if (this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }

    return;
}

// ======================================================================
// 删除移动方法
// ======================================================================

void Bulk_RingPolymer::delete_move(int k_max, int delete_polymer_index)
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
    this->rho_now = static_cast<double>(this->MN_now) / (this->box_size[0] * this->box_size[1] * this->box_size[2]);
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
    if (acc_ratio > this->uni_dis(this->gen))
    {
        this->acc_delete++;
    }
    else
    {
        this->MN_now += M;
        this->N_now += 1;
        this->rho_now = static_cast<double>(this->MN_now) / (this->box_size[0] * this->box_size[1] * this->box_size[2]);
        if (this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
}

// ======================================================================
// 单体操作方法
// ======================================================================

bool Bulk_RingPolymer::insert_one_monomer(std::vector<std::array<double, 3>> &r_new, double *Z_eff, int monomer_index, int k_max)
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
            candidate_positions[k][0] = this->box_size[0] * this->uni_dis(this->gen);
            candidate_positions[k][1] = this->box_size[1] * this->uni_dis(this->gen);
            candidate_positions[k][2] = this->box_size[2] * this->uni_dis(this->gen);

            // 检查与其他单体的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_weights[k] = 0.0;
            }
            else
            {
                // 体相系统无外势，权重为1
                z_weights[k] = 1.0;
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
        if (!this->Find_Last(r_new[this->M-2].data(), r_new[0].data(), candidate_ptrs.data(), k_max)) {
            // Find_Last失败，无法生成候选位置
            *Z_eff = 0;
            return false;
        }

        // 检查重叠并计算权重
        for (int k = 0; k < k_max; k++)
        {
            // 体相系统无外势，权重为1
            z_weights[k] = 1.0;
            // 检查与其他聚合物的重叠
            if (this->overlap_all_other_polymer(candidate_ptrs[k]))
            {
                z_weights[k] = 0.0;
                continue;
            }

            // 检查与当前链的重叠（跳过第一个单体，因为最后一个单体需要与第一个单体连接）
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
            add_random_unit_vector(candidate_positions[k].data(), this->gen);

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
                // 体相系统无外势，权重等于路径概率
                z_weights[k] = P_road[k];
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
        int select = RouteWheel(z_weights.data(), k_max, this->gen);

        r_new[monomer_index] = candidate_positions[select];

        *Z_eff *= sum_z_eff / k_max / P_road[select]; // 归一化，除以k_max
        return true;
    }
}

void Bulk_RingPolymer::delete_one_monomer(std::vector<std::array<double, 3>> &r_delete, double *Z_eff, int monomer_index, int k_max)
{
    std::vector<std::array<double, 3>> candidate_positions(k_max);
    std::vector<double> z_eff(k_max, 1.0);
    std::vector<double> P_road(k_max, 1.0); // 把 k_max - 1 设成 老节点的位置
    int grow_step = this->M - monomer_index;

    // 对于第一个单体，需要随机生成位置
    if (monomer_index == 0) // for first monomer
    {
        for (int k = 0; k < k_max-1; k++)
        {
            candidate_positions[k][0] = this->box_size[0] * this->uni_dis(this->gen);
            candidate_positions[k][1] = this->box_size[1] * this->uni_dis(this->gen);
            candidate_positions[k][2] = this->box_size[2] * this->uni_dis(this->gen);

            // 检查与其他单体的重叠
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_eff[k] = 0.0;
            }
            else
            {
                // 体相系统无外势，权重为1.0
                z_eff[k] = 1.0;
            }
        }
        // 第k_max-1个位置是原来的位置
        z_eff[k_max-1] = 1.0; // 体相系统无外势

    }
    else if (monomer_index == this->M - 1) // for last monomer
    {
        std::vector<double *> candidate_ptrs(k_max-1);

        for (int k = 0; k < k_max - 1; k++)
        {
            candidate_ptrs[k] = candidate_positions[k].data();
        }
        // 使用Find_Last生成可能的最后一个单体位置（与第一个单体连接）
        if (!this->Find_Last(r_delete[this->M-2].data(), r_delete[0].data(), candidate_ptrs.data(), k_max-1)) {
            // Find_Last失败，将所有候选位置权重设为0
            for (int k = 0; k < k_max - 1; k++) {
                z_eff[k] = 0.0;
            }
            // 原来的位置权重保持为1.0
            z_eff[k_max-1] = 1.0;
        }
        for( int k = 0; k < k_max - 1; k++)
        {
            // 体相系统无外势，权重初始为1.0
            z_eff[k] = 1.0;
            if (this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                z_eff[k] = 0.0;
                continue;
            }
            // 检查与当前链其他单体的重叠（除了第一个和倒数第二个）
            for (int j = 1; j < this->M - 2; j++)
            {
                if (this->overlap_other_monomer_one(r_delete[j].data(), candidate_positions[k].data()))
                {
                    z_eff[k] = 0.0;
                    break;
                }
            }
        }
        z_eff[k_max - 1] = 1.0; // 体相系统无外势
    }
    else // for middle monomer
    {
        // 构建已插入单体的索引列表

        for (int k = 0; k < k_max - 1; k++)
        {
            // 从前一个单体位置开始
            candidate_positions[k] = r_delete[monomer_index - 1];
            // 添加随机单位向量
            add_random_unit_vector(candidate_positions[k].data(), this->gen);

            // 检查碰撞
            bool overlap = false;
            // 检查与当前链的重叠
            for (int j = 0; j < monomer_index - 1; j++)
            {
                if (this->overlap_other_monomer_one(r_delete[j].data(), candidate_positions[k].data()))
                {
                    overlap = true;
                    break;
                }
            }
            // 检查与其他链的重叠
            if (!overlap && this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                overlap = true;
            }

            if (overlap)
            {
                z_eff[k] = 0.0;
                P_road[k] = 0.0;
            }
            else
            {
                // 计算两点之间的距离
                double R = distance_array(r_delete[0], candidate_positions[k]);
                // 使用打表法计算路径概率 P_road
                P_road[k] = g_p_road_table.get_P_road(R, grow_step);
                // 体相系统无外势，权重等于路径概率
                z_eff[k] = P_road[k];
            }
        }

        // 原来的位置
        double R = distance_array(r_delete[0], r_delete[monomer_index]);
        P_road[k_max - 1] = g_p_road_table.get_P_road(R, grow_step);
        z_eff[k_max - 1] = P_road[k_max - 1];
    }

    double sum_z_eff = 0.0;
    for (int k = 0; k < k_max; k++)
    {
        sum_z_eff += z_eff[k];
    }

    *Z_eff *= sum_z_eff / k_max / P_road[k_max - 1];
}