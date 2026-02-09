#define _USE_MATH_DEFINES
#ifndef MUVT_MC_LINEARPOLYMER_H
#define MUVT_MC_LINEARPOLYMER_H

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <functional>

constexpr double diameter = 1.0;
constexpr double diameter2 = 1.0;

class MuVT_MC_LinearPolymer
{

public:
    // 构造和析构
    MuVT_MC_LinearPolymer(std::string configuration,
        const double mu_b,
        const int M,
        int init_N,
        double rho_b,
        double box_xy,
        double H,
        double rcut,
        int max_N);
    virtual ~MuVT_MC_LinearPolymer();
    
    // 初始化方法
    virtual void init_second();
    
    // 参数设置和打印
    
    
    void set_sim_parameters(double EPS_TRANS, double ROT_RATIO, int K_MAX);// 设置模拟参数
    void print_all_parameters();// 打印所有参数
    
    // 模拟控制方法
    void end_block(int block); // 结束一个 block 的模拟 并重置信息
    void reset_mc_record(); // 重置 mc 接受率等参数
    void mc_one_step_NVT(); // mc 中 NVT 步
    void mc_one_step_MuVT(); // mc 中 MuVT 步
    
    // 插入方法
    void old_insert_move(int k_max);
    virtual void insert_move(int k_max); // 新的插入方法，使用insert_one_monomer
    
    // 删除方法
    void old_delete_move(int k_max,int delete_index);
    virtual void delete_move(int k_max,int delete_index); // 新的删除方法，使用delete_one_monomer

    // 单体操作方法
    virtual bool insert_one_monomer(std::vector<std::array<double, 3>> &r_new, double *W, int monomer_index, int k_max); // 插入单个单体
    virtual void delete_one_monomer(std::vector<std::array<double, 3>> &r_delete, double *W, int monomer_index, int k_max); // 删除单个单体 因为删除总能进行
    
    // 平移和旋转方法
    virtual void trans_move(int polymer_index); // 平移聚合物
    virtual void rot_polymer_move(int polymer_index); // 旋转聚合物
    virtual void rot_mid_move(int monomer_index, int polymer_index); // 转动单个单体
    void rot_end_move(int monomer_index, int polymer_index, int dirct); // 转动末端单体
    
    // 碰撞检测方法
    bool overlap_two_particle(const double *r_try, int monomer_index) ; // 检测两个单体之间有没有重叠
    bool overlap_other_polymer(const double *r_try, int now_polymer_index) ; // 检测一条聚合物与其它聚合物是否重叠
    bool overlap_other_monomer(const double *r_try, int polymer_index, const std::vector<int> &cal_list) ; // 检测单体与本聚合物中的其它单体有无重叠
    bool overlap_insert_polymer(const std::array<double,3> &r_try, int parent_index, const std::vector<std::array<double, 3>> &r_config, const std::vector<int> &cal_list) ; // 检测插入聚合物的碰撞
    bool overlap_other_monomer_one(const double *r_other, const double *r_try) ; // 检测两个单体之间的重叠
    bool overlap_all_other_polymer(const double *r_try) ; // 检测单体与其它所有聚合物有无重叠
    bool check_collision_except_polymer(const double *r_try, int ignore_polymer_index = -1) ; // 检测碰撞，排除指定聚合物
    bool check_collision_except_monomer(const double *r_try, int ignore_monomer_index = -1) ; // 检测碰撞，排除指定单体
    
    // 辅助方法
    void Find_Last(double *r1, double *r2, double **r_temp_rot, int k_max); // 查找最后一个单体位置
    void Generate_Next_Monomer(const double *r_now, std::vector<double> &r_next, int k_max); // 生成下一个单体位置
    virtual bool generate_last_monomer_position(const std::vector<std::array<double, 3>>& r_new, int M, int k_max, std::array<double, 3>& last_position, double& total_weight); // 生成最后一个单体位置（用于环状聚合物闭合）
    // virtual double calculate_insertion_weight(int k_max); // 计算插入聚合物时的权重
    
    // 体积缩放方法
    int get_change_volumn(int expan_or_shrink, double mc_cond_eps_volumn); // 体积缩放并检查重叠
    // 
    virtual double get_W_insert(int k_max); // 计算化学势 
    virtual double get_G_insert(double insert_z, int first_insert_index, int k_max); // 计算插入权重
    virtual bool insert_recursive(int next_idx, int parent_idx, int step, double &total_W, std::vector<std::array<double, 3>> &r_new, std::vector<int> &is_inserted, int k_max); // 递归插入单体 
    virtual double get_Wz_insert(double insert_z ,int first_insert_index, int k_max) ;
    // Getter methods for Python interface
    int get_MN_now() const { return MN_now; }
    int get_N_now() const { return N_now; }
    double get_rho_now() const { return rho_now; }
    int get_M() const { return M; }
    double get_H() const { return H; }
    double get_box_xy() const { return box_xy; }
    int get_max_N() const { return max_N; }
    const double** get_r_total() const { return const_cast<const double**>(r_total); } // Getter for r_total (返回一个指针，用于Python接口）
    
    // Getter methods for acceptance ratios
    double get_trans_acceptance() const;    // 获取平移接受率
    double get_rot_acceptance() const;      // 获取旋转接受率
    double get_insert_acceptance() const;   // 获取插入接受率
    double get_delete_acceptance() const;   // 获取删除接受率
    
    // 随机数设置
    void set_seed(unsigned int seed); // 设置随机数种子

    std::function<double(const double)> V_ext; // 外部势能函数
    std::string V_ext_name; // 外部势能函数名称
    bool is_hs_wall_potential; // 标记是否为 hs_wall 势能
    double Vext_pot(double z) const; // 获取外部势能
    double Vext_bz(const double* pos) const; // 计算外势玻尔兹曼因子
    void set_external_potential(std::function<double(const double)> potential, const std::string& name); // 设置外部势能（lambda 表达式版本）
    void set_external_potential(const std::string& name); // 设置外部势能（势能名称版本）


protected:
    
// --- 模拟状态 (子类需要修改) ---
    int MN_now;
    int N_now;
    double rho_now;
    double N_bulk;

// --- 盒子与系统参数 (子类需要修改) ---
    const int max_N;
    const int M;
    const double mu_b;
    const double exp_mu_b;
    const double rho_b; // 体相对应体积中的聚合物数目

    double V;
    double box_xy;
    double box_size[3];
    double rcut;
    double real_rcut[3];
    double H;
    int cell_num[3];

    // --- 核心数据指针 (子类必须能访问进行修改) ---
    double **r_total;
    double **r_init;
    double **r_insert_temp;
    int **r_total_cell;

    std::vector<std::vector<int>> topology_map;

    std::vector<int> cl_head; // 格子总数
    std::vector<int> cl_list; // max_N * M

    // --- mc move parameters ---
    double eps_trans; // 平移参数 ;
    double rot_ratio; // 单链旋转中 旋转的比例
    int k_max;

    double acc_trans;
    double acc_rot;
    double acc_regrow;
    double acc_insert;
    double acc_delete;

    int num_trans;
    int num_rot;
    int num_regrow;
    int num_insert;
    int num_delete;
    int num_reptate;
    
    // --- 随机数 ---
    unsigned int seed;
    std::mt19937 gen;
    std::uniform_real_distribution<double> real_dis;
    std::uniform_real_distribution<double> uni_dis;

    // --- 外部势能打表 ---
    std::vector<double> z_field_table;
    double table_dz ; // 表步长
    double table_inv_dz;

    
    
//-------------------------------------------------------------------------------------------------------------------------
        
    int cl_total_cells; // 格子总数
    const int cl_threshold = 0; // 阈值，超过该值则使用格子法检测重叠

    virtual void build_topology_map(); // 构建拓扑地图
    

    void build_cell_list(); // 分配内存空间
    
    
    int get_cell_index(const double* pos); // 计算坐标对应的格子ID


    void update_particle_in_cell_list(int particle_id);    
    
    // --- 初始化辅助 --- 
    void read_init_configure(std::string fileaddress); //  
    void build_memory_r(); // 
    
    void easy_cal_init_parameters(); // 
    
    // --- 需要被重写 或 被子类调用的 虚函数 --- 
    
    virtual bool check_configure_validity(); 

    static void Rot_temp(double *vec, double **r_temp_rot, int k_max);

    
private:
    // 对于线性聚合物不需要增添描述结构的数组
    MuVT_MC_LinearPolymer( const MuVT_MC_LinearPolymer &other ) = delete; // 禁用拷贝构造函数
    // 这行代码是什么意思？
    MuVT_MC_LinearPolymer& operator=( const MuVT_MC_LinearPolymer &other ) = delete; // 禁用赋值操作符
};

#endif