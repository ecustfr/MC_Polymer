#define _USE_MATH_DEFINES
#ifndef MUVT_MC_LINEARPOLYMER_EXT_H
#define MUVT_MC_LINEARPOLYMER_EXT_H

#include "MuVT_MC_LinearPolymer.h"
#include <vector>
#include <string>

class MuVT_MC_LinearPolymer_Ext : public MuVT_MC_LinearPolymer
{
public:
    // 构造和析构
    MuVT_MC_LinearPolymer_Ext(std::string configuration,
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
        double C);
    virtual ~MuVT_MC_LinearPolymer_Ext();

    // 势场计算方法
    double calculate_potential(const double* pos) const;
    double calculate_boltzmann(const double* pos) const;

    // 重写MC move方法
    virtual void trans_move(int polymer_index) override;
    virtual void rot_mid_move(int monomer_index, int polymer_index) override;
    void rot_end_move(int monomer_index, int polymer_index, int dirct); // 注意：父类中不是virtual
    virtual void insert_move(int k_max) override;
    virtual void delete_move(int k_max, int delete_polymer_index) override;

protected:
    // 势场参数
    std::vector<double> An;           // 正弦分量振幅
    std::vector<double> phi_n;        // 正弦分量相位
    std::vector<std::vector<double>> Vlin_par;  // 线性段参数 [2][4]
    std::vector<std::vector<double>> x_tar;     // 线性段区间 [2][4]
    double C;                         // 振幅缩放因子
    const double kB = 1.0;            // 玻尔兹曼常数
    const double beta = 1.0;          // 1/(kB*T)
};

#endif 