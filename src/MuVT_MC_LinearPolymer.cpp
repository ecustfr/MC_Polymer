# include "MuVT_MC_LinearPolymer.h"
# include "Utils.h"
# include <iostream>
# include <string>
# include <algorithm>
# include <stdexcept>   
# include <vector>
# include <fstream>
# include <sstream>
# include <array> 
# include <cmath>
# include <numeric>

MuVT_MC_LinearPolymer::MuVT_MC_LinearPolymer(
    std::string configuration,
    const double mu_b,
    const int M,
    int init_N,
    double rho_b,
    double box_xy,
    double H,
    double rcut,
    int max_N)
    : mu_b(mu_b),
      exp_mu_b(std::exp(mu_b)),
      M(M),
      N_now(init_N),
      rho_b(rho_b),
      rcut(rcut),
      box_xy(box_xy),
      H(H),
      max_N(max_N)
{


    this->easy_cal_init_parameters();

    this->build_memory_r();


    // 阅读初始构型这一步一定要在 build_memory_r 之后
    this->read_init_configure(configuration); 


    for(int dim= 0; dim<3 ;dim++) 
    {
        this->cell_num[dim] = static_cast<int>(floor(this->box_size[dim]/this->rcut));
        this->real_rcut[dim] = this->box_size[dim]/this->cell_num[dim];
    }


    this->build_cell_list();
    
    

    this->reset_mc_record();
    
    this->set_sim_parameters(0.1,0.3,10); // default parameters [eps_trans, rot_ratio , k_max]


}

MuVT_MC_LinearPolymer::~MuVT_MC_LinearPolymer()
{
    for (int i = 0; i < this->max_N * this->M; i++)
    {
        delete[] this->r_total[i];
        delete[] this->r_init[i];
        
    }
    delete[] this->r_total;
    delete[] this->r_init;
    for(int i = 0; i< this->M ; i++)
    {
        delete[] this->r_insert_temp[i];
    }
    delete[] this->r_insert_temp;
}


void MuVT_MC_LinearPolymer::easy_cal_init_parameters()
{
        // --- 随机数 ---
    
    this->seed = std::chrono::system_clock::now().time_since_epoch().count();
    this->gen = std::mt19937(this->seed);
    this->real_dis = std::uniform_real_distribution<double>(-1.0, 1.0);
    this->uni_dis = std::uniform_real_distribution<double>(0.0, 1.0);
    std::cout << "random number generator is initializated" << std::endl;

    this->N_bulk = rho_b * box_xy * box_xy * H/this->M;
    this->MN_now = this->N_now * this->M;
    this->box_size[0] = this->box_xy;
    this->box_size[1] = this->box_xy;
    this->box_size[2] = H;
    this->V = this->box_size[0] * this->box_size[1] * this->box_size[2];
    
    // 初始化 hs_wall 势能标记
    this->is_hs_wall_potential = true;

    std::cout << "easy parameters are calculated" <<std::endl;

}

void MuVT_MC_LinearPolymer::read_init_configure(std::string fileaddress)
{
    std::ifstream infile(fileaddress);
    if (!infile)
    {
        throw std::runtime_error("can't open the file");
    }
    int row = 0;
    std::string line;
    while (std::getline(infile, line) && row < this->N_now*this->M)
    {
        if (line.empty())
            continue;

        std::istringstream iss(line);
        double x, y, z; // 不能动态改变为读两个数
        if (!(iss >> x >> y >> z))
        {
            throw std::runtime_error("文件格式错误");
        }

        // std::cout << x << y << z << std::endl;

        r_total[row][0] = x;
        r_total[row][1] = y;
        r_total[row][2] = z;
        row++;
    }
    infile.close();

    std::cout << "Initial configuration read finished." << std::endl;
}

void MuVT_MC_LinearPolymer::init_second() // 
{
    this->build_topology_map();
    if(this->check_configure_validity())
    {
        std::cout << "Initial configuration is valid" << std::endl;
    }
    else
    {
        // std::cout << "Warning: Initial configuration is not valid, but continuing simulation" << std::endl;
        throw std::runtime_error("Initial configuration is not valid");
    }
   
    this->print_all_parameters();
}

void MuVT_MC_LinearPolymer::build_memory_r()
{
    this->r_total = new double *[this->max_N * this->M];
    this->r_total_cell = new int *[this->max_N * this->M];
    this->r_init = new double *[this->max_N * this->M];

    this->r_insert_temp = new double *[this->M];
    
    for (int i = 0; i < this->max_N * this->M; i++)
    {
        this->r_total[i] = new double[3];
        this->r_total_cell[i] = new int[3];
        this->r_init[i] = new double[3];
        
    }

    for( int i = 0 ;i<this->M ; i++)
    {
        this->r_insert_temp[i] = new double[3];
    }
}
void MuVT_MC_LinearPolymer::print_all_parameters()
{
    std::cout << "------------------ Simulation Parameters ------------------" << std::endl;
    std::cout << "Chemical potential mu_b: " << this->mu_b << std::endl;
    std::cout << "Exp(mu_b): " << this->exp_mu_b << std::endl;
    std::cout << "Polymerization degree M: " << this->M << std::endl;
    std::cout << "Initial polymer number N_now: " << this->N_now << std::endl;
    std::cout << "Polymer N bulk in volume: " << this->N_bulk << std::endl;
    std::cout << "Bulk density rho_b: " << this->rho_b << std::endl;
    std::cout << "Box size in x and y direction: " << this->box_xy << std::endl;
    std::cout << "Box size in z direction H: " << this->H << std::endl;
    std::cout << "Cutoff rcut: " << this->rcut << std::endl;
    std::cout << "Maximum polymer number max_N: " << this->max_N << std::endl;
    std::cout << "Now segment number MN_now: " << this->MN_now << std::endl;   
    std::cout << "real rcut in xyz direction: " << this->real_rcut[0]<<"  "<< this->real_rcut[1]<<" "<<this->real_rcut[2]<<" "<< std::endl;
    std::cout << "Cell number in xyz direction: " << this->cell_num[0]<<"  "<< this->cell_num[1]<<" "<<this->cell_num[2]<<" "<< std::endl;
    std::cout << "eps_trans: " << this->eps_trans << std::endl;
    std::cout << "rot_ratio: " << this->rot_ratio << std::endl;
    std::cout << "k_max: " << this->k_max << std::endl; 
    std::cout << "External potential: " << this->V_ext_name << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

}

// 实现 Vext_pot 方法
double MuVT_MC_LinearPolymer::Vext_pot(double z) const
{
    // 如果是 hs_wall 势能，直接返回 0
    if (this->is_hs_wall_potential) {
        return 0.0;
    }
    
    // 确保 z 在有效范围内
    z = std::max(0.0, std::min(z, this->H));
    
    // 如果打表未初始化，直接使用势能函数计算
    if (this->z_field_table.empty()) {
        return this->V_ext(z);
    }
    
    // 使用打表法计算外部势能
    double scaled_z = z * this->table_inv_dz;
    int idx = static_cast<int>(scaled_z);
    
    // 确保索引在有效范围内
    idx = std::max(0, std::min(idx, static_cast<int>(this->z_field_table.size()) - 2));
    
    // 线性插值
    double frac = scaled_z - idx;
    double V1 = this->z_field_table[idx];
    double V2 = this->z_field_table[idx + 1];
    
    return V1 + frac * (V2 - V1);
}

// 实现 Vext_bz 方法
double MuVT_MC_LinearPolymer::Vext_bz(const double* pos) const
{
    // 如果是 hs_wall 势能，返回 1.0
    if (this->is_hs_wall_potential) {
        return 1.0;
    }
    
    // 计算外势能量
    double z = pos[2];
    double U_ext = this->Vext_pot(z);
    
    // 返回玻尔兹曼因子（假设 beta=1.0）
    return std::exp(-U_ext);
}

// 实现 set_external_potential 方法（lambda 表达式版本）
void MuVT_MC_LinearPolymer::set_external_potential(std::function<double(const double)> potential, const std::string& name)
{
    this->V_ext = potential;
    this->V_ext_name = name;
    
    // 设置 hs_wall 势能标记
    this->is_hs_wall_potential = (name == "hs_wall");
    
    // 如果势能名称是 hs_wall，跳过打表初始化
    if (this->is_hs_wall_potential) {
        std::cout << "External potential set to: " << name << " (skipping table initialization)" << std::endl;
        return;
    }
    
    // 初始化外部势能打表
    if (this->V_ext) {
        // 设置默认表步长
        double dz = 0.01; // 默认步长为 0.01，不使用未初始化的 table_dz
        
        double inv_dz = 1.0 / dz;
        
        // 计算打表点数量
        int num_points = static_cast<int>(std::ceil(this->H / dz)) + 1;
        
        // 确保点数合理
        if (num_points <= 0) {
            num_points = 1001; // 默认点数
        } else if (num_points > 100000) {
            std::cout << "Warning: Too many points for potential table, using 100000 points" << std::endl;
            num_points = 100000;
        }
        
        this->z_field_table.resize(num_points);
        
        // 填充打表数据
        for (int i = 0; i < num_points; i++) {
            double z = i * dz;
            z = std::min(z, this->H); // 确保不超过 H
            this->z_field_table[i] = this->V_ext(z);
        }
        
        // 更新表步长和倒数
        this->table_dz = dz;
        this->table_inv_dz = inv_dz;
        
        std::cout << "External potential table initialized with " << num_points << " points" << std::endl;
    }
    
    std::cout << "External potential set to: " << name << std::endl;
}

// 实现 set_external_potential 方法（势能名称版本）
void MuVT_MC_LinearPolymer::set_external_potential(const std::string& name)
{
    // 这里可以根据势能名称设置不同的势能函数
    // 暂时使用一个默认的势能函数作为示例
    this->V_ext = [](double z) {
        return 0.0; // 默认势能为0
    };
    this->V_ext_name = name;
    
    // 设置 hs_wall 势能标记
    this->is_hs_wall_potential = (name == "hs_wall");
    
    // 如果势能名称是 hs_wall，跳过打表初始化
    if (this->is_hs_wall_potential) {
        std::cout << "External potential set to: " << name << " (skipping table initialization)" << std::endl;
        return;
    }
    
    // 初始化外部势能打表
    if (this->V_ext) {
        // 设置默认表步长
        double dz = 0.01;
        double inv_dz = 1.0 / dz;
        
        // 计算打表点数量
        int num_points = static_cast<int>(std::ceil(this->H / dz)) + 1;
        
        // 确保点数合理
        if (num_points <= 0) {
            num_points = 1001;
        } else if (num_points > 100000) {
            std::cout << "Warning: Too many points for potential table, using 100000 points" << std::endl;
            num_points = 100000;
        }
        
        this->z_field_table.resize(num_points);
        
        // 填充打表数据
        for (int i = 0; i < num_points; i++) {
            double z = i * dz;
            z = std::min(z, this->H);
            this->z_field_table[i] = this->V_ext(z);
        }
        
        // 更新表步长和倒数
        this->table_dz = dz;
        this->table_inv_dz = inv_dz;
        
        std::cout << "External potential table initialized with " << num_points << " points" << std::endl;
    }
    
    std::cout << "External potential set to: " << name << std::endl;
}


// 生成最后一个单体位置（用于环状聚合物闭合）
// 对于线性聚合物，此方法不适用，返回 false
bool MuVT_MC_LinearPolymer::generate_last_monomer_position(const std::vector<std::array<double, 3>>& r_new, int M, int k_max, std::array<double, 3>& last_position, double& total_weight)
{
    return false;
}

void MuVT_MC_LinearPolymer::end_block(int block)
{

    std::cout << "Block:" << block << "is finished." << std::endl;
    if(this->check_configure_validity())
    {
        std::cout << "Configuration is valid" << std::endl;
    }
    std::cout << "Acceptance ratio of translation move: " << static_cast<double>(this->acc_trans) / this->num_trans << std::endl;
    std::cout << "Acceptance ratio of rotation move: " << static_cast<double>(this->acc_rot) / this->num_rot << std::endl;
    std::cout << "Acceptance ratio of insert move:" << static_cast<double>(this->acc_insert)/this->num_insert << std::endl;
    std::cout << "Acceptance ratio of delete move:" << static_cast<double>(this->acc_delete)/this->num_delete << std::endl;
    
    
}

// Getter methods for acceptance ratios
// 计算平移接受率
double MuVT_MC_LinearPolymer::get_trans_acceptance() const
{
    if (this->num_trans == 0) {
        return 0.0;
    }
    return static_cast<double>(this->acc_trans) / this->num_trans;
}

// 计算旋转接受率
double MuVT_MC_LinearPolymer::get_rot_acceptance() const
{
    if (this->num_rot == 0) {
        return 0.0;
    }
    return static_cast<double>(this->acc_rot) / this->num_rot;
}

// 计算插入接受率
double MuVT_MC_LinearPolymer::get_insert_acceptance() const
{
    if (this->num_insert == 0) {
        return 0.0;
    }
    return static_cast<double>(this->acc_insert) / this->num_insert;
}

// 计算删除接受率
double MuVT_MC_LinearPolymer::get_delete_acceptance() const
{
    if (this->num_delete == 0) {
        return 0.0;
    }
    return static_cast<double>(this->acc_delete) / this->num_delete;
}

void MuVT_MC_LinearPolymer::reset_mc_record()
{
    this->acc_trans = 0 ;
    this->acc_rot = 0;
    this->acc_regrow = 0;
    this->acc_insert = 0;
    this->acc_delete = 0;

    this->num_trans = 0;
    this->num_rot = 0;
    this->num_regrow = 0;
    this->num_insert = 0;
    this->num_delete = 0;
}

void MuVT_MC_LinearPolymer::set_sim_parameters(double EPS_TRANS, double ROT_RATIO , int K_MAX)
{
this->eps_trans = EPS_TRANS;
this->rot_ratio = ROT_RATIO;
this->k_max = K_MAX;
}

// 设置随机数种子
void MuVT_MC_LinearPolymer::set_seed(unsigned int seed)
{
    this->seed = seed;
    this->gen = std::mt19937(this->seed);
    this->real_dis = std::uniform_real_distribution<double>(-1.0, 1.0);
    this->uni_dis = std::uniform_real_distribution<double>(0.0, 1.0);
    std::cout << "Random number generator seed set to: " << seed << std::endl;
}

// 删除单个单体的方法

void MuVT_MC_LinearPolymer::Find_Last(double *r1, double *r2, double **r_temp_rot, int k_max)
{
    std::vector<double> temp_theta;
    temp_theta.reserve(k_max);
    for (int i = 0; i < k_max; i++)
    {
        temp_theta.push_back(2 * M_PI * this->uni_dis(this->gen));
    }

    double temp = diameter2 - dis2_no_period(r1, r2) / 4.0;

    double mid_r[3] = {0.0};
    for (int i = 0; i < 3; i++)
    {
        mid_r[i] = (r1[i] + r2[i]) / 2;
    }
    if (temp <= 0)
    {
        for (int k = 0; k < k_max; k++)
        {
            for (int j = 0; j < 3; j++)
            {
                r_temp_rot[k][j] = mid_r[j];
            }
        }
        return;
    }

    double R = sqrt(temp);
    for (int i = 0; i < k_max; i++)
    {
        r_temp_rot[i][0] = R * cos(temp_theta[i]);
        r_temp_rot[i][1] = R * sin(temp_theta[i]);
        r_temp_rot[i][2] = 0.0; // z
    }

    double vec[3] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};
    double vec_length = sqrt(SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]));
    if (vec[0] >= 0)
    {
        vec[0] /= vec_length;
        vec[1] /= vec_length;
        vec[2] /= vec_length;
    }
    else
    {
        vec[0] /= -vec_length;
        vec[1] /= -vec_length;
        vec[2] /= -vec_length;
    }

    Rot_temp(vec, r_temp_rot, k_max);
    for (int i = 0; i < k_max; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            r_temp_rot[i][j] += mid_r[j]; // 平移到中�?
        }
    }
    return;

    // 找到最后一个旋转的单体
}

//-----------------------------------------------------------------------------------------------------------------
//----- 聚合物移动方程 ------
// 对 r_temp_rot 进行旋转，size[r_temp_rot]:k_max * 3
void MuVT_MC_LinearPolymer::Rot_temp(double *vec, double **r_temp_rot, int k_max)
{
    double phi = std::acos(std::max(-1.0, std::min(1.0, vec[2])));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    if (std::fabs(sin_phi) < 1e-12)
    {
        sin_phi = 1e-12;
    }
    double theta = 0.0;
    theta = std::asin(std::max(-1.0, std::min(1.0, vec[1] / sin_phi))); // 计算 theta

    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    double R00 = cos_theta * cos_phi;
    double R01 = -sin_theta;
    double R02 = cos_theta * sin_phi;

    double R10 = sin_theta * cos_phi;
    double R11 = cos_theta;
    double R12 = sin_theta * sin_phi;

    double R20 = -sin_phi;
    double R21 = 0.0;
    double R22 = cos_phi;

    for (int i = 0; i < k_max; i++)
    {
        double temp[3] = {r_temp_rot[i][0], r_temp_rot[i][1], r_temp_rot[i][2]};
        r_temp_rot[i][0] = R00 * temp[0] + R01 * temp[1] + R02 * temp[2];
        r_temp_rot[i][1] = R10 * temp[0] + R11 * temp[1] + R12 * temp[2];
        r_temp_rot[i][2] = R20 * temp[0] + R21 * temp[1] + R22 * temp[2];
    }
}

void MuVT_MC_LinearPolymer::rot_polymer_move(int polymer_index)
{
    std::vector<bool> rot_list(this->M,false);
    // std::cout << "polymer_index: " << polymer_index << std::endl;
    for(int monomer_index =0 ; monomer_index < this->M ; monomer_index++)
    {
        if( this->uni_dis(this->gen) < this->rot_ratio)
        {
            rot_list[monomer_index] = true;
        }
    }
    
    for(int monomer_index = 0; monomer_index < this->M ; monomer_index++)
    {

       if(rot_list[monomer_index])
       {
            if(monomer_index == 0)
            {
               
                this->rot_end_move(monomer_index,polymer_index, 1);
            }
            else if(monomer_index == this->M-1)
            {
                this->rot_end_move(monomer_index,polymer_index,-1);
            }
            else
            {
                this->rot_mid_move(monomer_index,polymer_index);
            }
       }
    }


}

// 向系统中插入一个聚合物链

void MuVT_MC_LinearPolymer::rot_end_move(int monomer_index ,int polymer_index, int direct)
{
    // direct = -1 旋转的尾， direct = 1 旋转的头
    this->num_rot++;
    int real_monomer_index = polymer_index * this->M + monomer_index;
    double r_fixed[3];

    for(int j = 0 ; j < 3 ; j++ )
    {
        r_fixed[j] = this->r_total[real_monomer_index + direct][j];
    }

    std::vector<std::array<double,3>> r_temp(2*this->k_max-1,{r_fixed[0],r_fixed[1],r_fixed[2]}); // 旋转后的位置备选，采用偏执采样方法，前 k_max 个为new，后 k_max-1 个为old
    
    std::vector<double> w(2*this->k_max-1 ,1); // w = exp(-\beta u)

    for(int k = 0 ; k<2*this->k_max-1 ; k++)
    {
        add_random_unit_vector(r_temp[k].data(),this->gen);
    }
    // 待检查
    for(int k = 0 ; k<2*this->k_max-1 ; k++)
    {
        if( this->overlap_all_other_polymer(r_temp[k].data()) )
        {
            w[k] = 0;
            continue;
        }
        for( int other_monomer = 2 ; other_monomer < this->M ; other_monomer++)
        {
            // 尾粒子序号：MN-1
            // 首粒子序号：M(N-1)
            if (direct == -1)
            {
                // std::cout << "monomer_index:"<< polymer_index*this->M+M - other_monomer<< std::endl;
                if( this->overlap_other_monomer_one(this->r_total[(polymer_index+1)*this->M - other_monomer], r_temp[k].data()) ) // 注意代码逻辑
                {
                    w[k] = 0;
                    break;
                }
            }
        
            else
            {
                // std::cout << "monomer_index:"<<(polymer_index)*this->M + other_monomer<< std::endl;
                if( this->overlap_other_monomer_one(this->r_total[(polymer_index)*this->M + other_monomer], r_temp[k].data()) )
                {
                    w[k] = 0;
                    break;
                }
            }
        }
        
        // 如果没有重叠，计算外势玻尔兹曼因子
        if (w[k] != 0)
        {
            w[k] = this->Vext_bz(r_temp[k].data());
        }
    }

    if( std::all_of(w.begin(), w.end(), [](double v) { return v == 0; }) )
    {
        return; // 所有备选位置都重叠，插入失败
    }

    if( acc_bias_or_not(w.data(),this->Vext_bz(this->r_total[real_monomer_index]),this->k_max ,this->gen) )
    {
        this->acc_rot++;
        int select_index = RouteWheel( w.data() , this->k_max, this->gen);
        this->r_total[real_monomer_index][0] = r_temp[select_index][0];
        this->r_total[real_monomer_index][1] = r_temp[select_index][1];
        this->r_total[real_monomer_index][2] = r_temp[select_index][2];

        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }


    
}

// 使用 check_collision_except_monomer 函数重写的 rot_mid_move 函数

void MuVT_MC_LinearPolymer::rot_mid_move(int monomer_index ,int polymer_index)
{
    // std::cout << "rot_mid_move start" << std::endl;
    this->num_rot++;
    int real_monomer_index = polymer_index * this->M + monomer_index;
    
    // std::cout << "rot_mid_move: real_monomer_index = " << real_monomer_index << std::endl;
    
    // 获取前一个和后一个单体的位置
    double r_prev[3];
    double r_next[3];
    int last_monomer = polymer_index*this->M + this->topology_map[monomer_index][0];
    int next_monomer = polymer_index*this->M + this->topology_map[monomer_index][1];

    int now_monomer = polymer_index*this->M + monomer_index;
    for(int j = 0; j < 3; j++)
    {
        r_prev[j] =  this->r_total[last_monomer][j] ; // this->r_total[real_monomer_index - 1][j];
        r_next[j] = this->r_total[next_monomer][j];
    }

    
    // std::cout << "rot_mid_move: got r_prev and r_next" << std::endl;
    
    // 使用 vector 存储数据（更安全，自动管理内存）
    std::vector<std::array<double,3>> r_temp_new_data(this->k_max);
    std::vector<std::array<double,3>> r_temp_old_data(this->k_max-1);
    
    // 创建指针数组用于 Find_Last（需要 double** 接口）
    std::vector<double*> r_temp_new_ptrs(this->k_max);
    std::vector<double*> r_temp_old_ptrs(this->k_max-1);
    for(int i = 0; i < this->k_max; i++)
    {
        r_temp_new_ptrs[i] = r_temp_new_data[i].data();
    }
    for(int i = 0; i < this->k_max-1; i++)
    {
        r_temp_old_ptrs[i] = r_temp_old_data[i].data();
    }
    
    //std::cout << "rot_mid_move: created vectors" << std::endl;
    
    // 生成 k_max 个新位置（连接前一个和后一个单体）
    Find_Last(r_prev, r_next, r_temp_new_ptrs.data(), this->k_max);
    
    // 生成 k_max-1 个旧位置（也从相同的前后位置生成，用于bias采样）
    Find_Last(r_prev, r_next, r_temp_old_ptrs.data(), this->k_max-1);
    
    // std::cout << "rot_mid_move: generated positions" << std::endl;
    
    // 将所有候选位置合并到一个向量中：前 k_max 个为新位置，后 k_max-1 个为旧位置
    std::vector<std::array<double,3>> r_temp(2*this->k_max-1);
    for(int k = 0; k < this->k_max; k++)
    {
        r_temp[k] = r_temp_new_data[k];
    }
    for(int k = 0; k < this->k_max-1; k++)
    {
        r_temp[k + this->k_max] = r_temp_old_data[k];
    }
    
    // std::cout << "rot_mid_move: merged positions" << std::endl;
    
    // 计算权重 w = exp(-\beta u)，如果不重叠则为玻尔兹曼因子，重叠则为0
    std::vector<double> w(2*this->k_max-1, 1.0);
    
    for(int k = 0; k < 2*this->k_max-1; k++)
    {
        if(this->check_collision_except_monomer(r_temp[k].data(), real_monomer_index))
        {
            w[k] = 0;
        }
        else
        {
            // 计算外势玻尔兹曼因子
            w[k] = this->Vext_bz(this->r_total[now_monomer]);
        }
    }
    
    // std::cout << "rot_mid_move: calculated weights" << std::endl;
    
    // 检查是否所有位置都重叠
    if(std::all_of(w.begin(), w.end(), [](double v) { return v == 0.0; }))
    {
        return; // 所有备选位置都重叠，旋转失败
    }
    
    // std::cout << "rot_mid_move: not all positions overlap" << std::endl;
    
    // 使用 bias 采样判断是否接受（旧权重为1，因为当前是等概率采样）
    if(acc_bias_or_not(w.data(), this->Vext_bz(this->r_total[now_monomer]), this->k_max, this->gen))
    {
        this->acc_rot++;
        // 从新位置中选择一个（使用前 k_max 个位置和权重）
        int select_index = RouteWheel(w.data(), this->k_max, this->gen);
        this->r_total[real_monomer_index][0] = r_temp[select_index][0];
        this->r_total[real_monomer_index][1] = r_temp[select_index][1];
        this->r_total[real_monomer_index][2] = r_temp[select_index][2];

        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
    
    // std::cout << "rot_mid_move end" << std::endl;
}

void MuVT_MC_LinearPolymer::old_insert_move(int k_max)
{
    this->num_insert++ ;
    this->r_insert_temp[0][0] = this->box_size[0]*this->uni_dis(this->gen);
    this->r_insert_temp[0][1] = this->box_size[1]*this->uni_dis(this->gen);
    this->r_insert_temp[0][2] = this->box_size[2]*this->uni_dis(this->gen);

    if(this->overlap_all_other_polymer(this->r_insert_temp[0]))
    {
        return ;
    }

    double W = 1;
    double r_now[3] = {this->r_insert_temp[0][0], this->r_insert_temp[0][1], this->r_insert_temp[0][2]};

    std::vector<double> Phi_HS_try(k_max,1.0);// exp(-\beta U_{hs})

    // std::vector<double> Phi_ext_try(k_max,1.0); 
    
    std::vector<std::array<double,3>>  r_temp(k_max,std::array<double,3>{r_now[0], r_now[1], r_now[2]}); //插入链的单体位置备选
    
    int select_index = 0;

    for(int monomer_index = 1; monomer_index<this->M; monomer_index++)
    {
        std::fill(Phi_HS_try.begin(), Phi_HS_try.end(), 1.0);
        // 生成 k_max 个随机方向
        for(int k = 0 ; k< k_max ; k++)
        {
            add_random_unit_vector(r_temp[k].data(), this->gen);
        }

        // 判断k_max 个备选位置是否与其它聚合物重叠
        // Phi_try = 1.0; 不重叠
        // Phi_try = 0.0; 重叠
        for(int k=0 ; k<k_max ;k++)
        {
            if( this->overlap_all_other_polymer(r_temp[k].data()) )
            {
                Phi_HS_try[k] =0.0;
                continue;
            }
            for( int old_monomer = 0 ; old_monomer < monomer_index -1; old_monomer++)
            {
                // std::cout << "old_monomer: " << old_monomer << std::endl;
                if( this->overlap_other_monomer_one(this->r_insert_temp[old_monomer], r_temp[k].data()) )
                {                   
                    
                    Phi_HS_try[k] = 0.0;
                    break;
                }
            }
        }

        if( std::all_of(Phi_HS_try.begin(), Phi_HS_try.end(), [](double v) { return v == 0.0; }) )
        {
            return; // 所有备选位置都重叠，插入失败
        }

        W = W*sum_1d(Phi_HS_try.data(), k_max)/k_max;
        
        select_index = RouteWheel(Phi_HS_try.data(), k_max, this->gen);
        this->r_insert_temp[monomer_index][0] = r_temp[select_index][0];
        this->r_insert_temp[monomer_index][1] = r_temp[select_index][1];
        this->r_insert_temp[monomer_index][2] = r_temp[select_index][2]; 

        // 重置 r_temp
        std::fill( r_temp.begin(), r_temp.end(), std::array<double,3> {r_temp[select_index][0],r_temp[select_index][1], r_temp[select_index][2]} );
            
    }


    // 所有单体都可以被插入
    double acc_ratio = this->exp_mu_b*W/(this->N_now+1)*(this->V);
    if (acc_ratio  >this->uni_dis(this->gen))
    {
        this->acc_insert++;
        for(int monomer_index =0 ; monomer_index < this->M ; monomer_index++)
        {
            for(int j =0 ; j<3 ; j++)
            {
                this->r_total[this->MN_now + monomer_index][j] = this->r_insert_temp[monomer_index][j];
            }   
        }
        this->MN_now += this->M;
        this->N_now++ ;

        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
    
    return;

    // 插入一个聚合物链
    /*
    this->num_insert++ ;
    
    // 为第一个粒子生成k_max个候选位置
    std::vector<std::array<double,3>> r_temp_first(k_max);
    std::vector<double> Phi_HS_try_first(k_max, 1.0);
    
    for(int k = 0; k < k_max; k++)
    {
        r_temp_first[k][0] = this->box_size[0] * this->uni_dis(this->gen);
        r_temp_first[k][1] = this->box_size[1] * this->uni_dis(this->gen);
        r_temp_first[k][2] = this->box_size[2] * this->uni_dis(this->gen);
        
        // 检查与其他聚合物的重叠
        if(this->overlap_all_other_polymer(r_temp_first[k].data()))
        {
            Phi_HS_try_first[k] = 0.0;
        }
    }
    
    // 检查是否所有位置都重叠
    if(std::all_of(Phi_HS_try_first.begin(), Phi_HS_try_first.end(), [](double v) { return v == 0.0; }))
    {
        return; // 所有备选位置都重叠，插入失败
    }
    
    // 计算权重并选择一个位置
    double sum_phi_first = sum_1d(Phi_HS_try_first.data(), k_max);
    double W = sum_phi_first / k_max; // 归一化，除以k_max
    int select_index = RouteWheel(Phi_HS_try_first.data(), k_max, this->gen);
    
    // 更新第一个粒子的位置
    this->r_insert_temp[0][0] = r_temp_first[select_index][0];
    this->r_insert_temp[0][1] = r_temp_first[select_index][1];
    this->r_insert_temp[0][2] = r_temp_first[select_index][2];

    double r_now[3] = {this->r_insert_temp[0][0], this->r_insert_temp[0][1], this->r_insert_temp[0][2]};

    std::vector<double> Phi_HS_try(k_max,1.0);// exp(-\beta U_{hs})

    // std::vector<double> Phi_ext_try(k_max,1.0); 
    
    std::vector<std::array<double,3>>  r_temp(k_max,std::array<double,3>{r_now[0], r_now[1], r_now[2]}); //插入链的单体位置备选

    for(int monomer_index = 1; monomer_index<this->M; monomer_index++)
    {
        std::fill(Phi_HS_try.begin(), Phi_HS_try.end(), 1.0);
        // 生成 k_max 个随机方向
        for(int k = 0 ; k< k_max ; k++)
        {
            add_random_unit_vector(r_temp[k].data(), this->gen);
        }

        // 判断k_max 个备选位置是否与其它聚合物重叠
        // Phi_try = 1.0; 不重叠
        // Phi_try = 0.0; 重叠
        for(int k=0 ; k<k_max ;k++)
        {
            if( this->overlap_all_other_polymer(r_temp[k].data()) )
            {
                Phi_HS_try[k] =0.0;
                continue;
            }
            for( int old_monomer = 0 ; old_monomer < monomer_index -1; old_monomer++)
            {
                // std::cout << "old_monomer: " << old_monomer << std::endl;
                if( this->overlap_other_monomer_one(this->r_insert_temp[old_monomer], r_temp[k].data()) )
                {                   
                    
                    Phi_HS_try[k] = 0.0;
                    break;
                }
            }
        }

        if( std::all_of(Phi_HS_try.begin(), Phi_HS_try.end(), [](double v) { return v == 0.0; }) )
        {
            return; // 所有备选位置都重叠，插入失败
        }

        W = W*sum_1d(Phi_HS_try.data(), k_max)/k_max; // 归一化，除以k_max
        
        select_index = RouteWheel(Phi_HS_try.data(), k_max, this->gen);
        this->r_insert_temp[monomer_index][0] = r_temp[select_index][0];
        this->r_insert_temp[monomer_index][1] = r_temp[select_index][1];
        this->r_insert_temp[monomer_index][2] = r_temp[select_index][2]; 

        // 重置 r_temp
        std::fill( r_temp.begin(), r_temp.end(), std::array<double,3> {r_temp[select_index][0],r_temp[select_index][1], r_temp[select_index][2]} );
            
    }


    // 所有单体都可以被插入
    double acc_ratio = this->exp_mu_b*W/(this->N_now+1)*(this->V);
    if (acc_ratio  >this->uni_dis(this->gen))
    {
        this->acc_insert++;
        for(int monomer_index =0 ; monomer_index < this->M ; monomer_index++)
        {
            for(int j =0 ; j<3 ; j++)
            {
                this->r_total[this->MN_now + monomer_index][j] = this->r_insert_temp[monomer_index][j];
            }   
        }
        this->MN_now += this->M;
        this->N_now++ ;

        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
    
    return;

    // 插入一个聚合物链
    */
}

// 插入单个单体
bool MuVT_MC_LinearPolymer::insert_one_monomer( std::vector<std::array<double, 3>> &r_new, double *W ,int monomer_index , int k_max )
{
    // 生成k_max个候选位置
    std::vector<std::array<double, 3>> candidate_positions(k_max);
    std::vector<double> weights(k_max, 1.0);
    
    
    if(monomer_index == 0)
    {
        // 简化第一个粒子的位置生成：直接在盒子内生成一个随机位置
        r_new[0][0] = box_size[0] * this->uni_dis(this->gen);
        r_new[0][1] = box_size[1] * this->uni_dis(this->gen);
        r_new[0][2] = box_size[2] * this->uni_dis(this->gen);
        

        // 检查与其他单体的重叠
        if ( this->overlap_all_other_polymer(r_new[0].data())) 
        {
            *W = 0.0;
            return false;
        }
        else
        {
            // 计算外势玻尔兹曼因子
            
            *W *=  this->Vext_bz(r_new[0].data()); 
            return true;
        }
        
    }
    else
    {
        // 创建需要碰撞检测的粒子
        std::vector<int> cal_index (monomer_index-1) ;
        std::iota(cal_index.begin(),cal_index.end(),0);


        for(int k=0;k<k_max;k++)
        {
           // 生成随机单位向量并添加到前一个单体位置
           candidate_positions[k] = r_new[monomer_index-1];
           add_random_unit_vector(candidate_positions[k].data(), this->gen);
            
           // 检测 插入本链的所有粒子，
           if(this->overlap_insert_polymer(candidate_positions[k], monomer_index - 1, r_new , cal_index) || this->overlap_all_other_polymer(candidate_positions[k].data()) )
           {
               weights[k] = 0.0;
           }
           else
        {
            // 计算外势玻尔兹曼因子
            weights[k] = this->Vext_bz(candidate_positions[k].data());
        }

        }

    }

    double sum_w = 0;
    for (int k = 0; k < k_max; k++) {
        sum_w += weights[k];
    }

    if( sum_w==0 )
    {
        *W = 0;
        return false; // 所有备选位置都重叠，插入失败
    }
    else
    {
        int select = RouteWheel(weights.data(), k_max, this->gen);
        r_new[monomer_index] = candidate_positions[select];
        *W *= sum_w / k_max; // 归一化，除以k_max，与insert_move方法一致
        return true;
    }

}


void MuVT_MC_LinearPolymer::insert_move(int k_max)
{
    this->num_insert++;
    
    // 1. 初始化新聚合物的位置
    std::vector<std::array<double, 3>> r_new(this->M);
    double W = 1.0;
    
    
    // 3. 插入剩余的单体
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &W, monomer_index, k_max))
        {
            return; // 单体插入失败
        }
    }
    
    // 4. 所有单体都可以被插入
    double acc_ratio = this->exp_mu_b * W / (this->N_now + 1) * (this->V);
    
    if (acc_ratio > this->uni_dis(this->gen))
    {
        this->acc_insert++;
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


void MuVT_MC_LinearPolymer::delete_one_monomer(std::vector<std::array<double, 3>> &r_delete, double *W, int monomer_index, int k_max)
{
    std::vector<std::array<double,3>>  candidate_positions(k_max-1);
    std::vector<double> weights(k_max-1, 1.0);
    
    if(monomer_index == 0)
    { 
        // 简化第一个粒子的权重计算：直接使用当前位置的外势玻尔兹曼因子
        // 不需要生成多个候选位置，只使用当前位置
        *W *= this->Vext_bz(r_delete[monomer_index].data());
        return ;
        // 只需要计算当前位置的权重
    }
    else
    {
        // int behind_count = std::max(, 0);
        std::vector<int> cal_index(monomer_index - 1);
        std::iota(cal_index.begin(), cal_index.end(), 0);

        for(int k=0;k<k_max-1;k++)
        {
            candidate_positions[k] = r_delete[monomer_index-1];
            add_random_unit_vector(candidate_positions[k].data(), this->gen);

            if(this->overlap_insert_polymer(candidate_positions[k], monomer_index - 1, r_delete , cal_index)||this->overlap_all_other_polymer(candidate_positions[k].data()))
            {
                weights[k] = 0.0;
            }
            else
            {
                weights[k] = this->Vext_bz(candidate_positions[k].data());
            }
        
        }
    }

    double sum_w = 0;
    for(int k = 0; k < k_max - 1; k++)
    {
        sum_w += weights[k];
    }
    
    sum_w = sum_w +  this->Vext_bz(r_delete[monomer_index].data());
    
    *W *=sum_w/k_max; 
    return;

}

// 新的删除方法，使用delete_one_monomer
void MuVT_MC_LinearPolymer::delete_move(int k_max, int delete_index)
{
    this->num_delete++;
    
    int delete_first_monomer_index = delete_index * this->M;
    int last_polymer_start = this->MN_now - this->M;
    
    // 交换指针，把要删除的聚合物链放到最后面
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++) {
        std::swap(this->r_total[delete_first_monomer_index + monomer_index], this->r_total[last_polymer_start + monomer_index]);
    }
    
    std::vector<std::array<double, 3>> r_delete(this->M);
    for (int i = 0; i < this->M; i++) 
    {
        for (int j = 0; j < 3; ++j) {
            r_delete[i][j] = this->r_total[last_polymer_start + i][j];
        }
    }

    
    // 先假设聚合物被删除
    this->MN_now -= this->M;
    this->N_now -= 1;
    
        // 更新cell列表，因为粒子位置已经改变

    if (this->MN_now > this->cl_threshold) {
        this->build_cell_list();
    }
    double W = 1.0;

    
    // 使用delete_one_monomer方法计算每个单体的权重
    // 最后一次删除 没办法删除 所以不计算多个权重
    for(int monomer_index = 0; monomer_index < this->M; monomer_index++) 
    {
        this->delete_one_monomer(r_delete, &W, monomer_index, k_max);
    }

    double acc_ratio = 1/W * (this->N_now+1)/(this->V) * 1/this->exp_mu_b;
   
    if(acc_ratio > this->uni_dis(this->gen))
    {
        this->acc_delete++;
        // 接受删除，重建cell列表
        // this->build_cell_list();
    }
    else
    {
        // 拒绝删除，恢复MN_now和N_now
        this->MN_now += this->M;
        this->N_now += 1;
    
        
        // 更新cell列表，恢复正确的粒子位置
        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
}


void MuVT_MC_LinearPolymer::old_delete_move(int k_max, int delete_polymer_index)
{
        int delete_first_monomer_index = (delete_polymer_index) * this->M;
    int last_polymer_start = this->MN_now - this->M;

    this->num_delete++;

    // 交换指针，把要删除的聚合物链放到最后面
    for(int monomer_index = 0 ; monomer_index < this->M ; monomer_index++)
    {
        std::swap( this->r_total[delete_first_monomer_index + monomer_index] , this->r_total[last_polymer_start + monomer_index] );
    }
    
    // 更新cell列表，因为粒子位置已经改变
    if(this->MN_now > this->cl_threshold)
    {
        this->build_cell_list();
    }

    // 先假设聚合物被删除
    this->MN_now -= this->M;
    this->N_now -= 1;
     
    double W = 1.0;
    std::vector<std::array<double,3>> r_temp(k_max-1); // 被删除单体的位置备选

    std::vector<double> Phi_HS_try(k_max-1,1.0);

    std::fill(r_temp.begin(), r_temp.end(), std::array<double,3>{this->r_total[this->MN_now][0], this->r_total[this->MN_now][1], this->r_total[this->MN_now][2]}); 

    for(int monomer_index = 1; monomer_index < this->M; monomer_index++)
    {
        std::fill(Phi_HS_try.begin(), Phi_HS_try.end(), 1.0);
        for(int k=0 ; k<k_max-1 ; k++)
        {
            add_random_unit_vector(r_temp[k].data(),this->gen);
        }

        for(int k=0 ; k<k_max-1 ; k++)
        {
            if( this->overlap_all_other_polymer(r_temp[k].data()) ) // 把聚合物链放到最后就是为了方便这一步
            {
                Phi_HS_try[k] =0.0;
                continue;
            }
            for( int old_monomer = 0 ; old_monomer <monomer_index-1  ; old_monomer++)
            {
                if( this->overlap_other_monomer_one(this->r_total[this->MN_now + old_monomer], r_temp[k].data()) )
                {
                    Phi_HS_try[k] = 0.0;
                    break;
                }
            }
        }

        W = W*(sum_1d(Phi_HS_try.data(),k_max-1)+1)/k_max;
        std::fill( r_temp.begin(), r_temp.end(), std::array<double,3> 
        {this->r_total[this->MN_now+monomer_index][0],this->r_total[this->MN_now+monomer_index][1], this->r_total[this->MN_now+monomer_index][2]} );

    }
    double acc_ratio = 1/W * (this->N_now+1)/(this->V) * 1/this->exp_mu_b;
    if(acc_ratio > uni_dis(gen))
    {
        this->acc_delete++;
        // 接受删除，重建cell列表
        this->build_cell_list();
    }
    else
    {
        // 拒绝删除，恢复MN_now和N_now
        this->MN_now += this->M;
        this->N_now += 1;
        
        // 将聚合物链交换回原来的位置
        for(int monomer_index = 0 ; monomer_index < this->M ; monomer_index++)
        {
            std::swap( this->r_total[delete_first_monomer_index + monomer_index] , this->r_total[last_polymer_start + monomer_index] );
        }
        
        // 更新cell列表，恢复正确的粒子位置
        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
    /*
    int delete_first_monomer_index = (delete_polymer_index) * this->M;
    int last_polymer_start = this->MN_now - this->M;

    this->num_delete++;

    // 交换指针，把要删除的聚合物链放到最后面
    for(int monomer_index = 0 ; monomer_index < this->M ; monomer_index++)
    {
        std::swap( this->r_total[delete_first_monomer_index + monomer_index] , this->r_total[last_polymer_start + monomer_index] );
    }
    
    // 更新cell列表，因为粒子位置已经改变
    if(this->MN_now > this->cl_threshold)
    {
        this->build_cell_list();
    }

    // 先假设聚合物被删除
    this->MN_now -= this->M;
    this->N_now -= 1;
     
    double W = 1.0;
    std::vector<std::array<double,3>> r_temp(k_max-1); // 被删除单体的位置备选

    std::vector<double> Phi_HS_try(k_max-1,1.0);

    std::fill(r_temp.begin(), r_temp.end(), std::array<double,3>{this->r_total[this->MN_now][0], this->r_total[this->MN_now][1], this->r_total[this->MN_now][2]}); 

    for(int monomer_index = 1; monomer_index < this->M; monomer_index++)
    {
        std::fill(Phi_HS_try.begin(), Phi_HS_try.end(), 1.0);
        for(int k=0 ; k<k_max-1 ; k++)
        {
            add_random_unit_vector(r_temp[k].data(),this->gen);
        }

        for(int k=0 ; k<k_max-1 ; k++)
        {
            if( this->overlap_all_other_polymer(r_temp[k].data()) ) // 把聚合物链放到最后就是为了方便这一步
            {
                Phi_HS_try[k] =0.0;
                continue;
            }
            for( int old_monomer = 0 ; old_monomer <monomer_index-1  ; old_monomer++)
            {
                if( this->overlap_other_monomer_one(this->r_total[this->MN_now + old_monomer], r_temp[k].data()) )
                {
                    Phi_HS_try[k] = 0.0;
                    break;
                }
            }
        }

        W = W*(sum_1d(Phi_HS_try.data(),k_max-1)+1)/k_max;
        std::fill( r_temp.begin(), r_temp.end(), std::array<double,3> 
        {this->r_total[this->MN_now+monomer_index][0],this->r_total[this->MN_now+monomer_index][1], this->r_total[this->MN_now+monomer_index][2]} );

    }
    double acc_ratio = 1/W * (this->N_now+1)/(this->V) * 1/this->exp_mu_b;
    if(acc_ratio > uni_dis(gen))
    {
        this->acc_delete++;
        // 接受删除，重建cell列表
        this->build_cell_list();
    }
    else
    {
        // 拒绝删除，恢复MN_now和N_now
        this->MN_now += this->M;
        this->N_now += 1;
        
        // 将聚合物链交换回原来的位置
        for(int monomer_index = 0 ; monomer_index < this->M ; monomer_index++)
        {
            std::swap( this->r_total[delete_first_monomer_index + monomer_index] , this->r_total[last_polymer_start + monomer_index] );
        }
        
        // 更新cell列表，恢复正确的粒子位置
        if(this->MN_now > this->cl_threshold)
        {
            this->build_cell_list();
        }
    }
    */
}



void MuVT_MC_LinearPolymer::trans_move(int polymer_index)
{
    this->num_trans++;
    double delta_move[3];
    double r_try[3];
    // int r_try_cell[3];
    double temp_dis = 0.0;

    for (int j = 0; j < 3; j++) // 生成 聚合物 随机 位移
    {
        delta_move[j] = this->real_dis(this->gen) * this->box_size[j] * this->eps_trans;
    }

    // 检查重叠
    for (int monomer_index = polymer_index * this->M; monomer_index < (polymer_index + 1) * this->M; monomer_index++)
    {
        for (int j = 0; j < 3; j++)
        {
            r_try[j] = this->r_total[monomer_index][j] + delta_move[j];
            // r_try_cell[j] = static_cast<int>(floor(r_try[j] / this->rcut)) % this->cell_num[j];
        }
        if (this->overlap_other_polymer(r_try, polymer_index))
        {
            return;
        }
    }
    
    // 计算外势能量差
    double delta_U_ext = 0.0;
    for (int monomer_index = polymer_index * this->M; monomer_index < (polymer_index + 1) * this->M; monomer_index++)
    {
        // 计算移动后的位置
        double r_new[3];
        for (int j = 0; j < 3; j++)
        {
            r_new[j] = this->r_total[monomer_index][j] + delta_move[j];
        }
        
        // 计算外势能量差
        delta_U_ext += this->Vext_pot(r_new[2]) - this->Vext_pot(this->r_total[monomer_index][2]);
    }
    
    // 计算玻尔兹曼因子
    double boltz_factor = std::exp(-delta_U_ext);
    
    // 使用玻尔兹曼因子决定是否接受移动
    if (boltz_factor > this->uni_dis(this->gen))
    {
        // std::cout <<delta_move[0] << " " << delta_move[1] << " " << delta_move[2] << std::endl;
        this->acc_trans++;
        for (int monomer_index = polymer_index * this->M; monomer_index < (polymer_index + 1) * this->M; monomer_index++)
        {
            for (int j = 0; j < 3; j++)
            {
                this->r_total[monomer_index][j] += delta_move[j];
            }
        }
        // 

        if(this->MN_now > this-> cl_threshold)
        {
            this->build_cell_list();
        }
    }
}
//-------------------------------------------------------------------------------------------------------------------------
//------判断重叠类函数 ------------------------------------------------------------------------------------------------------
bool MuVT_MC_LinearPolymer::overlap_two_particle(double *r_try, int monomer_index)
{
    int temp_dis_int;
    double temp_dis=0 ;

    if(r_try[2] > this->H-0.5 || r_try[2] < 0.5)
    {
        return true;
    }
    
    
    for(int j = 0 ;j<2;j++)
    {
        temp_dis += SQR(dis_period(this->r_total[monomer_index][j],r_try[j], this->box_size[j]));

        if(temp_dis > diameter2)
        {
            return false;
        }
    }
    temp_dis += SQR(this->r_total[monomer_index][2]-r_try[2]);

    if(temp_dis <= diameter2)
    {
       // std::cout << "overlap distance2: " << temp_dis << std::endl;
        return true;
    }
    else
    {
    return false;
    }
    
    
}

// 判断插入聚合物时是否与其他单体重叠
bool MuVT_MC_LinearPolymer::overlap_insert_polymer(std::array<double,3> &r_try,int parent_index,std::vector< std::array<double, 3> > &r_config , std::vector<int> &cal_list )
{
    for(const auto& index : cal_list)
    {
        // 跳过父节点
        if(index == parent_index)
        {
            continue;
        }
        
        // 检查是否与其他单体重叠
        if(this->overlap_other_monomer_one(r_try.data(), r_config[index].data() ))
        {
            return true;
        }
    }
    return false;
}

bool MuVT_MC_LinearPolymer::overlap_other_polymer(double *r_try, int now_polymer_index)
{
    // 判断当前 聚合物与 其它聚合物有无重叠
    if(this->MN_now > this-> cl_threshold)
    {
        return this->check_collision_except_polymer(r_try, now_polymer_index);
    }
    else
    {
        for(int polymer_index = 0 ; polymer_index< this->N_now;polymer_index++)
        {
            if(polymer_index == now_polymer_index)
            {
                continue;
            }
            else
            {
                for(int par_index = (polymer_index)*this->M ; par_index< (polymer_index+1)*this->M ; par_index++)
                {
                    if(overlap_two_particle(r_try , par_index))
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
bool MuVT_MC_LinearPolymer::overlap_other_monomer(double *r_try, int polymer_index , std::vector<int> &cal_list)
{
    for(const auto& index : cal_list)
    {
        if(this->overlap_two_particle(r_try, polymer_index*this->M+index))
        {
            return true;
        }
    }
    return false;
}
bool MuVT_MC_LinearPolymer::overlap_other_monomer_one(const double *r_try, const double *r_other)
{

    if(r_try[2] > this->H - 0.5 || r_try[2] < 0.5)
    {
        return true;
    }
    double temp_dis = 0;
    for(int j = 0; j < 2; j++)
    {
        temp_dis += SQR(dis_period(r_other[j], r_try[j], this->box_size[j]));
        if(temp_dis > diameter2)
        {
            return false;
        }
    }
    temp_dis += SQR(r_other[2]-r_try[2]);
    if(temp_dis <= diameter2)
    {
        return true;
    }
    else
    {
        return false;
    }


}


bool MuVT_MC_LinearPolymer::overlap_all_other_polymer(const double *r_try)
{
    if(this->MN_now > this-> cl_threshold)
    {
        return this->check_collision_except_polymer(r_try,-1);
    }
    else
    {
        for(int par_index = 0; par_index < this->MN_now ;par_index++)
        {
            // if(this->overlap_two_particle(r_try , par_index))
            if(this->overlap_other_monomer_one(this->r_total[par_index], r_try))
            {
                return true;
            }
            
        }
    }


    return false;
}
bool MuVT_MC_LinearPolymer::check_configure_validity()
{
    for(int par_index = 0; par_index < this->MN_now-1; par_index++)
    {
        if(par_index+1%this->M ==0) // 末端单体没有连接下一个单体
        {
            if(this->overlap_two_particle(this->r_total[par_index], par_index+1))
            {
                return false;
            }
        }

        for(int other_par_index = par_index+2; other_par_index < this->MN_now; other_par_index++)
        {
            if(this->overlap_two_particle(this->r_total[par_index], other_par_index))
            {
                std::cout << this->r_total[par_index][0] << " " << this->r_total[par_index][1] << " " << this->r_total[par_index][2] <<std::endl;
                std::cout << this->r_total[other_par_index][0] << " " << this->r_total[other_par_index][1] << " " << this->r_total[other_par_index][2] <<std::endl;

                // 打印cell_num lattice信息
                std::cout << "=== Cell Lattice Information ===" << std::endl;
                std::cout << "cell_num: " << this->cell_num[0] << "(x) " << this->cell_num[1] << "(y) " << this->cell_num[2] << "(z)" << std::endl;
                std::cout << "real_rcut: " << this->real_rcut[0] << "(x), " << this->real_rcut[1] << "(y), " << this->real_rcut[2] << "(z)" << std::endl;
                
                // 计算并打印两个粒子所在的cell索引
                int cell_idx1 = this->get_cell_index(this->r_total[par_index]);
                int cell_idx2 = this->get_cell_index(this->r_total[other_par_index]);
                std::cout << "Monomer " << par_index << " cell index: " << cell_idx1 << std::endl;
                std::cout << "Monomer " << other_par_index << " cell index: " << cell_idx2 << std::endl;
                
                // 计算并打印cell坐标
                int ix1 = cell_idx1 % this->cell_num[0];
                int iy1 = (cell_idx1 / this->cell_num[0]) % this->cell_num[1];
                int iz1 = cell_idx1 / (this->cell_num[0] * this->cell_num[1]);
                std::cout << "Monomer " << par_index << " cell coordinates: " << ix1 << ", " << iy1 << ", " << iz1 << std::endl;
                
                int ix2 = cell_idx2 % this->cell_num[0];
                int iy2 = (cell_idx2 / this->cell_num[0]) % this->cell_num[1];
                int iz2 = cell_idx2 / (this->cell_num[0] * this->cell_num[1]);
                std::cout << "Monomer " << other_par_index << " cell coordinates: " << ix2 << ", " << iy2 << ", " << iz2 << std::endl;
                
                std::cout << "overlap between monomer " << par_index << " and monomer " << other_par_index << std::endl;
                return false;
            }
        }
    }

    return true;
}

int MuVT_MC_LinearPolymer::get_cell_index(const double* r)
{

    int ix_raw = static_cast<int>( floor( r[0] / this->real_rcut[0] )) ;
    int iy_raw = static_cast<int>( floor( r[1] / this->real_rcut[1] )) ; 

    int ix = (ix_raw % this->cell_num[0] + this->cell_num[0]) % this->cell_num[0];
    int iy = (iy_raw % this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];

    
    int iz = static_cast<int>(floor(r[2] / real_rcut[2])); // z 不会出界的

    // if (iz < 0 || iz >= this->cell_num[2]) std::cout << "z方向越界" << std::endl; // z方向越界
    
    return ix + iy * this->cell_num[0] + iz * this->cell_num[0] * this->cell_num[1];
}
//---------------------------------------------------------------------------------------------------------------------------

void MuVT_MC_LinearPolymer::build_cell_list()
{
    // 初始化大小


    this->cl_total_cells = this->cell_num[0] * this->cell_num[1] * this->cell_num[2];
    if(this->cl_head.size() != this->cl_total_cells)
    {
        this->cl_head.resize(this->cl_total_cells); // 给链表头分配空间
    }

    if(this->cl_list.size()!= this->max_N*this->M)
    {
        this->cl_list.resize(this->max_N*this->M); // 给链表分配空间
    }

    // 1. 清空Head
    std::fill(this->cl_head.begin(),this->cl_head.end(),-1); // 很重要
    
    // 2. 填充链表
    for(int i= 0 ; i< this->MN_now ; i++)
    {
        int cell_idx = this->get_cell_index(this->r_total[i]);

        this->cl_list[i] = this->cl_head[cell_idx];
        this->cl_head[cell_idx] = i;
    }

}

void MuVT_MC_LinearPolymer::build_topology_map()
{
    this->topology_map.resize(this->M);
    for(int i = 0; i < this->M;i++)
    {
        // 前一个邻居
        if(i > 0)
        {
            this->topology_map[i].push_back(i-1);
        }
        // 后一个邻居
        // his->topology_map[i].push_back(i);

        if(i < this->M-1)
        {
            this->topology_map[i].push_back(i+1);
        }

    }
    std::cout << "topology map is OK" <<std::endl;
}

bool MuVT_MC_LinearPolymer::check_collision_except_polymer( const double* r_try, int ignore_polymer_index ) // 使用cell list 进行碰撞检测
{
    if(r_try[2] > this->H - 0.5 || r_try[2] < 0.5) return true;

    int cx = static_cast<int> (floor(r_try[0] / this->real_rcut[0]));
    int cy = static_cast<int> (floor(r_try[1] / this->real_rcut[1]));
    int cz = static_cast<int> (floor(r_try[2] / this->real_rcut[2]));

    /*
    cx = (cx % this->cell_num[0] + this->cell_num[0]) % this->cell_num[0];
    cy = (cy % this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];
    cz = (cz % this->cell_num[2] + this->cell_num[2]) % this->cell_num[2];
    */

    // if(cz < 0 || cz >= this->cell_num[2]) std::cout << "z方向越界" << std::endl; // z方向越界
    int z_min = std::max(0,cz-1);
    int z_max = std::min(cz+1,this->cell_num[2]-1);
    for(int z = z_min;z<=z_max ; z++ )
    {
        for(int y = cy-1; y<= cy+1 ; y++)
        {
            for(int x = cx-1 ; x <= cx+1 ; x++)
            {
                int real_x = (x%this->cell_num[0] + this->cell_num[0]) % this->cell_num[0]; 
                int real_y = (y%this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];

                int cell_idx = real_x + real_y*this->cell_num[0] + z*this->cell_num[0]*this->cell_num[1]; 
                 
                // 遍历该格子的链表
                int p_id = this->cl_head[cell_idx] ; 
                // std::cout << 0 << p_id << std::endl;
                
                while( p_id != -1 && p_id < this->MN_now)
                {
                    // 逻辑：如果这个粒子属于我们应该忽略的聚合物，就跳过检测
                    if(ignore_polymer_index != -1 && (p_id / this->M) == ignore_polymer_index)
                    {
                        // std::cout << 1 << p_id << std::endl;
                        p_id = this->cl_list[p_id]; // 得到下个粒子的索引
                        // std::cout << 2 <<p_id << std::endl;
                        continue;
                    }

                    // 计算距离
                    // 可以展开优化速度
                    
                    if(this->overlap_other_monomer_one(this->r_total[p_id], const_cast<double*>(r_try)))
                    {
                        return true; // 碰撞发生
                    }

                    // std::cout << 3  <<p_id << std::endl;
                    p_id = this->cl_list[p_id]; // 下一个粒子索引
                    // std::cout << 4  <<p_id << std::endl;
                }

            }
        }
    }
    return false; // 无碰撞
}
// 数据处理类，只处理和记录数据 不保存数据
bool MuVT_MC_LinearPolymer::check_collision_except_monomer(const double* r_try, int ignore_monomer_index)
{
    // ... 边界检查和 Cell 计算 (保持不变) ...

    // === 优化核心：预计算需要忽略的全局 ID 列表 ===
    // 我们不需要在 while 循环里做复杂的 prev/next 逻辑判断
    // 而是使用一个固定大小的数组存储“黑名单”
    
    // 为了性能，用一个固定小数组代替 vector，因为邻居通常很少
    if(r_try[2] > this->H - 0.5 || r_try[2] < 0.5) return true;


    int cx = static_cast<int> (floor(r_try[0] / this->real_rcut[0]));
    int cy = static_cast<int> (floor(r_try[1] / this->real_rcut[1]));
    int cz = static_cast<int> (floor(r_try[2] / this->real_rcut[2]));

    /*
    cx = (cx % this->cell_num[0] + this->cell_num[0]) % this->cell_num[0];
    cy = (cy % this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];
    cz = (cz % this->cell_num[2] + this->cell_num[2]) % this->cell_num[2];
    */
    int ignore_list[10];  // 最多忽略 10 个邻居
    int ignore_count = 0;

    if (ignore_monomer_index != -1)
    {
        // 1. 添加自身
        ignore_list[ignore_count] = ignore_monomer_index;
        ignore_count++;

        // 2. 计算所在的聚合物基准 ID (Polymer Base ID)
        int local_rank = ignore_monomer_index % this->M;
        int polymer_base = ignore_monomer_index - local_rank;

        // 3. 查表！获取所有邻居
        // topology_map[local_rank] 返回的是局部邻居 (如 0, 2)
        // 我们加上 polymer_base 转换回全局 ID
        for (int local_neighbor : this->topology_map[local_rank]) // 遍历
        {
            ignore_list[ignore_count] = polymer_base + local_neighbor;
            ignore_count++;
        }
    }

    // ... Cell 循环 (z, y, x) 保持不变 ...
    int z_min = std::max(0,cz-1);
    int z_max = std::min(this->cell_num[2]-1 , cz+1);
    for (int z = z_min; z <= z_max; z++)
    {

        for(int y = cy-1; y<= cy+1 ; y++)
        {
            for(int x = cx-1 ; x <= cx+1 ; x++)
            {
                int real_x = (x%this->cell_num[0] + this->cell_num[0]) % this->cell_num[0]; 
                int real_y = (y%this->cell_num[1] + this->cell_num[1]) % this->cell_num[1];

                int cell_idx = real_x + real_y*this->cell_num[0] + z*this->cell_num[0]*this->cell_num[1]; 
                 
                // 遍历该格子的链表
                int p_id = this->cl_head[cell_idx] ; 
                
                while( p_id != -1 && p_id < this->MN_now)
                {
                    // --- 极速黑名单检查 ---
                    // 检查 p_id 是否在 ignore_list 中
                    bool is_ignored = false;
                    
                    if (ignore_count > 0)
                    {
                        // 循环展开非常快，因为 ignore_count 通常只有 2 或 3
                        for (int k = 0; k < ignore_count; k++)
                        {
                            if (p_id == ignore_list[k])
                            {
                                is_ignored = true;
                                break;
                            }
                        }
                    }

                    if (is_ignored)
                    {
                        p_id = this->cl_list[p_id];
                        continue;
                    }

                    // --- 距离检测 ---
                    // 调用你的单体距离检查函数
                    // 注意：r_total[p_id] 是已存在的粒子，r_try 是尝试位置
                    if (this->overlap_other_monomer_one(this->r_total[p_id], const_cast<double*>(r_try)))
                    {
                        return true; // 碰撞发生
                    }

                    // 移动到链表下一个
                    p_id = this->cl_list[p_id];
                }
            }
        }
    }
    return false;// 安全，无碰撞
}

void MuVT_MC_LinearPolymer::mc_one_step_NVT()
{
    int now_N = this->N_now;
    for(int polymer_index =0 ; polymer_index < now_N ; polymer_index++)
    {
        this->trans_move(polymer_index);
        this->rot_polymer_move(polymer_index);
        
    } 
};

void MuVT_MC_LinearPolymer::mc_one_step_MuVT()
{
    this->insert_move(this->k_max);
    //this->old_insert_move(this->k_max);

    int delete_N = static_cast<int>(this->N_now*this->uni_dis(this->gen));

    //this->old_delete_move(this->k_max,delete_N);

    this->delete_move(this->k_max, delete_N);
}


// 体积缩放并检查重叠
int MuVT_MC_LinearPolymer::get_change_volumn(int expan_or_shrink, double mc_cond_eps_volumn)
{
    // 计算聚合物质心
    std::vector<std::array<double, 3>> r_center(N_now);
    for (int i = 0; i < N_now; i++) {
        r_center[i][0] = 0.0;
        r_center[i][1] = 0.0;
        r_center[i][2] = 0.0;
        
        for (int par_index = i * M; par_index < (i + 1) * M; par_index++) {
            r_center[i][0] += r_total[par_index][0];
            r_center[i][1] += r_total[par_index][1];
            r_center[i][2] += r_total[par_index][2];
        }
        
        r_center[i][0] /= M;
        r_center[i][1] /= M;
        r_center[i][2] /= M;
    }
    
    double delta_move[3] = {0.0};
    double temp_eps_volumn = mc_cond_eps_volumn * expan_or_shrink;
    double new_box_size[3];
    for (int dim = 0; dim < 3; dim++) {
        new_box_size[dim] = box_size[dim] / (1 + temp_eps_volumn);
    }
    
    // 计算盒子中心
    double box_center[3] = {box_size[0] / 2, box_size[1] / 2, box_size[2] / 2};
    
    // 分配临时内存
    double** r_total_volumn = new double*[MN_now];
    int** r_total_cell_volumn = new int*[MN_now];
    for (int i = 0; i < MN_now; i++) {
        r_total_volumn[i] = new double[3];
        r_total_cell_volumn[i] = new int[3];
    }
    
    // 更新位置和格子索引
    for (int i = 0; i < N_now; i++) {
        // 计算移动量
        for (int dim = 0; dim < 3; dim++) {
            delta_move[dim] = 1/(1+temp_eps_volumn) * (r_center[i][dim] - box_center[dim]) + box_center[dim] - r_center[i][dim];
        }
        
        // 移动聚合物到盒子中心
        for (int par_index = i * M; par_index < (i + 1) * M; par_index++) {
            for (int dim = 0; dim < 3; dim++) {
                r_total_volumn[par_index][dim] = r_total[par_index][dim] + delta_move[dim];
                r_total_cell_volumn[par_index][dim] = static_cast<int>(floor(r_total_volumn[par_index][dim] / real_rcut[dim])) % cell_num[dim];
                if (r_total_cell_volumn[par_index][dim] < 0) {
                    r_total_cell_volumn[par_index][dim] += cell_num[dim];
                }
            }
        }
    }
    
    // 构建连接链表
    std::vector<std::vector<int>> cell_list(cl_total_cells);
    for (int par_index = 0; par_index < MN_now; par_index++) {
        int cell_idx = get_cell_index(r_total_volumn[par_index]);
        cell_list[cell_idx].push_back(par_index);
    }
    
    // 检查是否重叠
    double temp_dis = 0.0;
    int next_polymer_fir_monomer = 0;
    int temp_dis_int = 0;
    
    for (int par_index = 0; par_index < MN_now - M; par_index++) {
        next_polymer_fir_monomer = (par_index / M + 1) * M;
        
        for (int other_par_index = next_polymer_fir_monomer; other_par_index < MN_now; other_par_index++) {
            temp_dis = 0.0;
            bool check_distance = true;
            
            // 检查格子是否相邻
            for (int dim = 0; dim < 3; dim++) {
                temp_dis_int = (r_total_cell_volumn[par_index][dim] - r_total_cell_volumn[other_par_index][dim]) % cell_num[dim];
                if (temp_dis_int < 0) {
                    temp_dis_int += cell_num[dim];
                }
                
                if (!(temp_dis_int == 0 || temp_dis_int == 1 || temp_dis_int == cell_num[dim] - 1)) {
                    check_distance = false;
                    break;
                }
            }
            
            if (!check_distance) {
                continue;
            }
            
            // 计算周期性距离
            for (int dim = 0; dim < 3; dim++) {
                double dx = fabs(r_total_volumn[par_index][dim] - r_total_volumn[other_par_index][dim]);
                dx = std::min(dx, new_box_size[dim] - dx);
                temp_dis += dx * dx;
                
                if (temp_dis > diameter2) {
                    break;
                }
            }
            
            // 检查是否重叠
            if (temp_dis <= diameter2) {
                // 释放内存
                for (int i = 0; i < MN_now; i++) {
                    delete[] r_total_volumn[i];
                    delete[] r_total_cell_volumn[i];
                }
                delete[] r_total_volumn;
                delete[] r_total_cell_volumn;
                
                return 1; // 重叠了
            }
        }
    }
    
    // 释放内存
    for (int i = 0; i < MN_now; i++) {
        delete[] r_total_volumn[i];
        delete[] r_total_cell_volumn[i];
    }
    delete[] r_total_volumn;
    delete[] r_total_cell_volumn;
    
    return 0; // 没有重叠
}


double MuVT_MC_LinearPolymer::get_W_insert(int k_max) 
{
    double W_insert = 1.0;
    std::vector<std::array<double, 3>> r_new(this->M);
    for (int monomer_index = 0; monomer_index < this->M; monomer_index++)
    {
        if (!this->insert_one_monomer(r_new, &W_insert, monomer_index, k_max))
        {
            W_insert = 0;
            return W_insert; // 单体插入失败
        }
    }
    return W_insert;
}

double MuVT_MC_LinearPolymer::get_G_insert(double insert_z ,int first_insert_index, int k_max) 
{
    
    std::vector<std::array<double,3>> r_new(this->M);
    std::vector<int> is_inserted;
    is_inserted.reserve(M-1);// 不需要申请第一个插入粒子的
    double W = 1.0;
    

    // --- 第一步：放置种子节点 ---
    for (int d = 0; d < 2; ++d)
    { 
        r_new[first_insert_index][d] = box_size[d] * this->uni_dis(this->gen);
    }
    r_new[first_insert_index][2] = insert_z;
    
    
    // --- 第二步：向左递归 (step = -1) ---
    if (!insert_recursive(first_insert_index - 1, first_insert_index, -1, W, r_new, is_inserted, k_max)) {
        W = 0.0;
        return W;
    }

    // --- 第三步：向右递归 (step = 1) ---
    if (!insert_recursive(first_insert_index + 1, first_insert_index, 1, W, r_new, is_inserted, k_max)) {
        W = 0.0;
        return W;
    }

    return W;
}

double MuVT_MC_LinearPolymer::get_Wz_insert(double insert_z ,int first_insert_index, int k_max) 
{
    double Wz = 1.0;
    std::vector<std::array<double,3>> r_new(this->M);
    std::vector<int> is_inserted;
    is_inserted.reserve(M);// 不需要申请第一个插入粒子的
    is_inserted.push_back(first_insert_index); // 标记第一个粒子已经插入


    // --- 第一步：放置种子节点 ---
    for (int d = 0; d < 2; ++d)
    { 
        r_new[first_insert_index][d] = box_size[d] * this->uni_dis(this->gen);
    }
    r_new[first_insert_index][2] = insert_z;
    
    if(this->overlap_all_other_polymer(r_new[first_insert_index].data()))
    {
        Wz = 0.0;
        return Wz;
    }
    else
    {
        Wz = this->Vext_bz(r_new[first_insert_index].data());
    }
    // --- 第二步：向左递归 (step = -1) ---
    if (!insert_recursive(first_insert_index - 1, first_insert_index, -1, Wz, r_new, is_inserted, k_max)) {
        Wz = 0.0;
        return Wz;
    }

    // --- 第三步：向右递归 (step = 1) ---
    if (!insert_recursive(first_insert_index + 1, first_insert_index, 1, Wz, r_new, is_inserted, k_max)) {
        Wz = 0.0;
        return Wz;
    }

    return Wz;
}

bool MuVT_MC_LinearPolymer::insert_recursive(
    int next_idx, int parent_idx, int step, 
    double &total_W, 
    std::vector<std::array<double, 3>> &r_new, 
    std::vector<int> &is_inserted,
    int k_max)
{
    // 递归出口
    if (next_idx < 0 || next_idx >= this->M) {return true;}

    std::vector<std::array<double, 3>> candidates(k_max);
    std::vector<double> w_trial(k_max, 0.0);
    double sum_w = 0.0;

    // 1. 尝试 k_max 个候选位置
    for (int k = 0; k < k_max; ++k) {
        candidates[k] = r_new[parent_idx]; // 基于父节点位置
        add_random_unit_vector(candidates[k].data(), this->gen);

        // 2. 碰撞检测：检查其他链 + 检查本链已插入的部分
        if (this->overlap_all_other_polymer(candidates[k].data()) || this->overlap_insert_polymer(candidates[k], parent_idx, r_new, is_inserted)) 
        {
            w_trial[k] = 0; // 其实是多余的，但是为了保证
        }
        else
        {
            w_trial[k] = this->Vext_bz(candidates[k].data());
            sum_w += w_trial[k];
        }
    }

    if (sum_w <= 1e-20) return false;// 递归出口

    // 3. 选择并记录位置
    int selected = RouteWheel(w_trial.data(), k_max, this->gen);
    r_new[next_idx] = candidates[selected];
    is_inserted.push_back(next_idx); // 关键：标记为已插入
    
    // 4. 更新权重
    total_W *= (sum_w / static_cast<double>(k_max));

    // 5. 继续向同一方向递归
    return insert_recursive(next_idx + step, next_idx, step, total_W, r_new, is_inserted, k_max);
}

// 将打表数据输出到文件
void MuVT_MC_LinearPolymer::export_potential_table(const std::string& filename) const
{
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl;
        return;
    }

    // 写入表头
    outfile << "# z \t V_ext(z)" << std::endl;
    outfile << "# Potential table for " << V_ext_name << std::endl;
    outfile << "# Table step dz: " << table_dz << std::endl;
    outfile << "# Number of points: " << z_field_table.size() << std::endl;

    // 写入数据
    for (size_t i = 0; i < z_field_table.size(); ++i) {
        double z = i * table_dz;
        z = std::min(z, H); // 确保不超过 H
        outfile << z << "\t" << z_field_table[i] << std::endl;
    }

    outfile.close();
    std::cout << "Potential table exported to " << filename << std::endl;
}

// 返回打表数据给 Python 接口
std::vector<double> MuVT_MC_LinearPolymer::get_potential_table() const
{
    return z_field_table;
}

// 返回打表对应的 z 坐标
std::vector<double> MuVT_MC_LinearPolymer::get_potential_table_z() const
{
    std::vector<double> z_values(z_field_table.size());
    for (size_t i = 0; i < z_field_table.size(); ++i) {
        double z = i * table_dz;
        z_values[i] = std::min(z, H); // 确保不超过 H
    }
    return z_values;
}