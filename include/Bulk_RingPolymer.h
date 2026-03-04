#define _USE_MATH_DEFINES
#ifndef BULK_RINGPOLYMER_H
#define BULK_RINGPOLYMER_H

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <chrono>

constexpr double diameter = 1.0;
constexpr double diameter2 = 1.0;

class Bulk_RingPolymer {
public:
    // 构造函数和析构函数
    Bulk_RingPolymer(std::string configuration, int M, int init_N,
                    double rho, double box_size[3], double rcut, int max_N,
                    double mu_b = 0.0, double rho_b = 0.0);
    ~Bulk_RingPolymer();
    
    // 初始化和设置
    void init_second();
    void set_sim_parameters(double EPS_TRANS, double ROT_RATIO, int K_MAX);
    void print_all_parameters();
    
    // 模拟核心
    void mc_one_step(); // 单步MC模拟
    void run_simulation(int steps); // 运行多步模拟
    
    // 辅助方法
    void end_block(int block);
    void reset_mc_record();
    
    // Getter方法
    int get_MN_now() const { return MN_now; }
    int get_N_now() const { return N_now; }
    double get_rho_now() const { return rho_now; }
    const double** get_r_total() const { return const_cast<const double**>(r_total); }
    
    // 体积缩放并检查重叠
    int change_volumn(int expan_or_shrink, double mc_cond_eps_volumn);
    
    // 插入高分子并计算插入权重
    double calculate_insertion_weight(int k_max);
    double get_W_insert_ring(int k_max);
    double get_W_insert_ring_z(double z, int k_max);

    // 巨正则模拟方法
    void mc_one_step_MuVT();                    // 巨正则单步模拟
    void run_simulation_MuVT(int steps);        // 巨正则多步模拟

    // 插入/删除移动
    void insert_move(int k_max);
    void delete_move(int k_max, int delete_polymer_index);

    // 单体操作方法
    bool insert_one_monomer(std::vector<std::array<double, 3>> &r_new,
                           double *Z_eff, int monomer_index, int k_max);
    void delete_one_monomer(std::vector<std::array<double, 3>> &r_delete,
                           double *Z_eff, int monomer_index, int k_max);

    // 统计获取方法
    double get_insert_acceptance() const;
    double get_delete_acceptance() const;
    double get_translation_acceptance() const;
    double get_rotation_acceptance() const;
    
private:
    // 系统状态
    int MN_now; // 总单体数
    int N_now; // 当前聚合物数
    double rho_now; // 当前密度
    
    // 系统参数
    const int M; // 每个聚合物的单体数
    const double rho; // 目标密度
    double box_size[3]; // 体相盒子大小 (三个方向都有周期性边界)
    double V; // 系统体积
    double rcut; // 截断半径
    double real_rcut[3]; // 实际截断半径
    int cell_num[3]; // 格子数
    const int max_N; // 最大聚合物数
    const double mu_b; // 化学势
    const double exp_mu_b; // exp(βμ)
    const double rho_b; // 体相参考密度
    
    // 核心数据
    double** r_total; // 所有单体位置
    double** r_init; // 初始单体位置
    double** r_temp_volumn; // 体积缩放临时位置
    
    std::vector<std::vector<int>> topology_map; // 拓扑图
    
    // 格子列表
    std::vector<int> cl_head; // 格子头
    std::vector<int> cl_list; // 格子列表
    int cl_total_cells; // 总格子数
    const int cl_threshold = 0; // 格子法阈值
    
    // MC参数
    double eps_trans; // 平移参数
    double rot_ratio; // 旋转比例
    int k_max; // CBMC尝试次数
    
    // 接受率统计
    double acc_trans; // 平移接受数
    double acc_rot; // 旋转接受数
    double acc_insert; // 插入接受数
    double acc_delete; // 删除接受数
    int num_trans; // 平移尝试次数
    int num_rot; // 旋转尝试次数
    int num_insert; // 插入尝试次数
    int num_delete; // 删除尝试次数
    
    // 随机数生成
    unsigned int seed;
    std::mt19937 gen;
    std::uniform_real_distribution<double> real_dis; // [-1,1]均匀分布
    std::uniform_real_distribution<double> uni_dis; // [0,1]均匀分布

    // 旋转移动缓存（性能优化）
    std::vector<std::array<double, 3>> rot_mid_cache_data;
    std::vector<double*> rot_mid_cache_ptrs;
    std::vector<double> rot_mid_cache_weights;
    
    // 私有方法
    void build_topology_map(); // 构建环状拓扑
    void build_cell_list(); // 构建格子列表
    int get_cell_index(const double* pos); // 获取格子索引
    
    // 初始化辅助
    void read_init_configure(std::string fileaddress);
    void build_memory_r();
    void easy_cal_init_parameters();
    
    // 边界条件处理
    double dis_period(double x1, double x2, double box_size) const; // 周期性距离计算
    
    // 距离计算函数
    double dis2_no_period(const double* r1, const double* r2) const; // 非周期性距离平方
    double dis2_period(const double* r1, const double* r2) const; // 三维周期性距离平方
    
    // 碰撞检测
    bool check_configure_validity();
    bool overlap_other_monomer_one(const double* r_other, const double* r_try) const;
    bool overlap_all_other_polymer(const double* r_try);
    bool check_collision_except_polymer(const double* r_try, int ignore_polymer_index = -1);
    bool check_collision_except_monomer(const double* r_try, int ignore_monomer_index = -1);
    
    // MC移动
    void trans_move(int polymer_index); // 平移聚合物
    void rot_polymer_move(int polymer_index); // 旋转聚合物
    void rot_mid_move(int monomer_index, int polymer_index); // 中间单体旋转
    
    // 辅助函数
    static void Rot_temp(double* vec, double** r_temp_rot, int k_max);
    bool Find_Last(double* r1, double* r2, double** r_temp_rot, int k_max);
    
    // 禁用拷贝构造和赋值
    Bulk_RingPolymer(const Bulk_RingPolymer& other) = delete;
    Bulk_RingPolymer& operator=(const Bulk_RingPolymer& other) = delete;
};

#endif