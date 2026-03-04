#include "Bulk_RingPolymer.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdlib> // for atof
#include <iomanip> // for std::setprecision
#include <sstream> // for std::stringstream
#include <limits>   // for std::numeric_limits

int main(int argc, char* argv[]) {
    // 测试Bulk_RingPolymer类，统计平衡时系统的密度
    std::cout << "=== 测试 Bulk_RingPolymer 密度计算 ===" << std::endl;

    // 检查命令行参数


    // 从命令行读取化学势
    double mu_b = std::atof(argv[1]);

    // 设置参数（根据用户要求）
    std::string configuration = "D:\\Desktop\\NewMuVTPolymer\\input\\init_config_Z\\TrivialN64M8H20.00rho0.10.dat";
    int M = 8; // 每个聚合物的单体数
    int init_N = 64; // 初始聚合物数 (64/8=8)
    double rho = 0.1; // 目标密度
    double box_size[3] = {20.0, 20.0, 20.0}; // 体相盒子大小
    
    /*
    if (mu_b <1.0)
    {
        box_size[0] = 40.0;
        box_size[1] = 40.0;
        box_size[2] = 40.0;
    }
    */
    
    double rcut = 1.0; // 截断半径
    int max_N = 1000; // 最大聚合物数
    double rho_b = 0.1; // 体相参考密度

    std::cout << "参数设置:" << std::endl;
    std::cout << "  输入文件: " << configuration << std::endl;
    std::cout << "  单体数/聚合物 (M): " << M << std::endl;
    std::cout << "  初始聚合物数: " << init_N << std::endl;
    std::cout << "  目标密度: " << rho << std::endl;
    std::cout << "  盒子大小: [" << box_size[0] << ", " << box_size[1] << ", " << box_size[2] << "]" << std::endl;
    std::cout << "  截断半径: " << rcut << std::endl;
    std::cout << "  最大聚合物数: " << max_N << std::endl;
    std::cout << "  化学势 (mu_b): " << mu_b << std::endl;
    std::cout << "  初始密度 (rho_b): " << rho_b << std::endl;

    try {
        // 创建Bulk_RingPolymer对象
        std::cout << "\n创建 Bulk_RingPolymer 对象..." << std::endl;
        Bulk_RingPolymer bulk_ring(configuration, M, init_N, rho, box_size, rcut, max_N, mu_b, rho_b);

        // 初始化第二步
        std::cout << "执行第二步初始化..." << std::endl;
        bulk_ring.init_second();

        // 设置模拟参数
        double EPS_TRANS = 0.02;
        double ROT_RATIO = 0.3;
        int K_MAX = 10;
        bulk_ring.set_sim_parameters(EPS_TRANS, ROT_RATIO, K_MAX);

        // 打印参数
        bulk_ring.print_all_parameters();

        // 弛豫阶段
        int equilibration_steps = 5000;
        std::cout << "\n开始弛豫阶段 (" << equilibration_steps << " 步)..." << std::endl;
        for (int step = 0; step < equilibration_steps; step++) {
            if (step % 1000 == 0) {
                std::cout << "  弛豫步数: " << step << "/" << equilibration_steps << std::endl;
            }
            bulk_ring.run_simulation_MuVT(1); // 运行1步MuVT模拟
        }

        // 采样阶段（分块进行，每块结束后调整平移步长）
        int production_steps = 40000;
        int blocks = 4; // 将总步数分为4个块
        int steps_per_block = production_steps / blocks;
        int sample_interval = 5;

        double target_acceptance = 0.5; // 目标平移接受率
        // 使用已有的 EPS_TRANS, ROT_RATIO, K_MAX（已在弛豫阶段设置）
        // 在块间调整 EPS_TRANS
        bulk_ring.set_sim_parameters(EPS_TRANS, ROT_RATIO, K_MAX);

        // 存储每个块的平均值
        std::vector<double> block_avg_density;
        std::vector<double> block_avg_polymer_count;
        std::vector<double> block_avg_W;
        std::vector<double> block_avg_mu_ex; // 每个块的过剩化学势 μ_ex = -log(<W>)
        std::vector<double> block_avg_mu_id; // 每个块的理想化学势 μ_id = log(ρ/M)
        std::vector<int> block_valid_W_counts; // 每个块中W>0的样本数

        // 存储所有样本（用于整体平均值，可选）
        std::vector<double> all_density_samples;
        std::vector<double> all_polymer_count_samples;
        std::vector<double> all_insertion_weight_samples;

        std::cout << "\n开始采样阶段 (" << production_steps << " 步，分为 " << blocks << " 个块)..." << std::endl;

        for (int block = 0; block < blocks; block++) {
            std::cout << "\n=== 块 " << block << " (" << steps_per_block << " 步) ===" << std::endl;

            // 重置当前块的统计计数器（可选，但为了清晰起见）
            bulk_ring.reset_mc_record();

            // 当前块的样本
            std::vector<double> block_density_samples;
            std::vector<double> block_polymer_count_samples;
            std::vector<double> block_insertion_weight_samples;

            for (int step = 0; step < steps_per_block; step++) {
                if (step % 2000 == 0) {
                    std::cout << "  块内步数: " << step << "/" << steps_per_block << std::endl;
                }
                bulk_ring.run_simulation_MuVT(1); // 运行1步MuVT模拟

                if (step % sample_interval == 0) {
                    // 记录当前密度和聚合物数
                    double density = bulk_ring.get_rho_now();
                    int polymer_count = bulk_ring.get_N_now();
                    block_density_samples.push_back(density);
                    block_polymer_count_samples.push_back(polymer_count);
                    all_density_samples.push_back(density);
                    all_polymer_count_samples.push_back(polymer_count);

                    // 计算插入权重 W
                    double W = bulk_ring.get_W_insert_ring(K_MAX);
                    block_insertion_weight_samples.push_back(W);
                    all_insertion_weight_samples.push_back(W);
                }
            }

            // 结束当前块，打印接受率
            bulk_ring.end_block(block);

            // 计算当前块的平均密度和平均聚合物数
            if (!block_density_samples.empty()) {
                double block_avg_dens = std::accumulate(block_density_samples.begin(), block_density_samples.end(), 0.0) / block_density_samples.size();
                double block_avg_poly = std::accumulate(block_polymer_count_samples.begin(), block_polymer_count_samples.end(), 0.0) / block_polymer_count_samples.size();

                // 计算当前块的平均插入权重（排除零值）
                double sum_W = 0.0;
                int count_nonzero_W = 0;
                for (double W : block_insertion_weight_samples) {
                    sum_W += W;
                    count_nonzero_W++;
                }
                double block_avg_W_val = (count_nonzero_W > 0) ? sum_W / count_nonzero_W : 0.0;

                block_avg_density.push_back(block_avg_dens);
                block_avg_polymer_count.push_back(block_avg_poly);
                block_avg_W.push_back(block_avg_W_val);
                block_valid_W_counts.push_back(count_nonzero_W);

                // 计算当前块的过剩化学势 μ_ex = -log(<W>)
                double block_mu_ex = (block_avg_W_val > 0.0) ? -std::log(block_avg_W_val) : std::numeric_limits<double>::quiet_NaN();
                block_avg_mu_ex.push_back(block_mu_ex);

                // 计算当前块的理想化学势 μ_id = log(ρ/M)
                double block_mu_id = (block_avg_dens > 0.0) ? std::log(block_avg_dens / M) : std::numeric_limits<double>::quiet_NaN();
                block_avg_mu_id.push_back(block_mu_id);

                std::cout << "  块平均密度: " << block_avg_dens << std::endl;
                std::cout << "  块平均聚合物数: " << block_avg_poly << std::endl;
                if (block_avg_W_val > 0.0) {
                    std::cout << "  块平均插入权重 <W>: " << block_avg_W_val << " (有效样本: " << count_nonzero_W << ")" << std::endl;
                } else {
                    std::cout << "  块平均插入权重 <W>: 0 (无有效样本)" << std::endl;
                }
            }

            // 根据平移接受率调整平移步长 EPS_TRANS
            double trans_accept = bulk_ring.get_translation_acceptance();
            if (trans_accept > 0.0) {
                // 调整因子：目标接受率/实际接受率，但限制调整幅度
                double adjust_factor = trans_accept / target_acceptance;
                // 限制调整因子在合理范围内（例如0.8~1.2）
                if (adjust_factor < 0.2) adjust_factor = 0.2;
                if (adjust_factor > 3.0) adjust_factor = 3.0;
                EPS_TRANS *= adjust_factor;
                // 限制 EPS_TRANS 在合理范围内
                if (EPS_TRANS < 0.001) EPS_TRANS = 0.001;
                if (EPS_TRANS > 0.5) EPS_TRANS = 0.5;

                bulk_ring.set_sim_parameters(EPS_TRANS, ROT_RATIO, K_MAX);
                std::cout << "  平移接受率: " << trans_accept << ", 调整 EPS_TRANS 为: " << EPS_TRANS << std::endl;
            } else {
                std::cout << "  平移接受率为0，不调整 EPS_TRANS" << std::endl;
            }
        }

        // ===== 整体统计结果（基于所有样本）=====
        if (!all_density_samples.empty()) {
            // 整体平均值（所有样本）
            double avg_density = std::accumulate(all_density_samples.begin(), all_density_samples.end(), 0.0) / all_density_samples.size();
            double avg_polymer_count = std::accumulate(all_polymer_count_samples.begin(), all_polymer_count_samples.end(), 0.0) / all_polymer_count_samples.size();

            // 整体插入权重的平均值（排除零值）
            double sum_W = 0.0;
            int count_nonzero_W = 0;
            for (double W : all_insertion_weight_samples) {
                sum_W += W;
                count_nonzero_W++;
            }
            double avg_W = (count_nonzero_W > 0) ? sum_W / count_nonzero_W : 0.0;

            // 整体样本的标准差（单样本波动）
            double sum_sq_diff = 0.0;
            for (double d : all_density_samples) {
                sum_sq_diff += (d - avg_density) * (d - avg_density);
            }
            double std_density = std::sqrt(sum_sq_diff / all_density_samples.size());

            // ===== 块平均统计 =====
            int num_blocks = block_avg_density.size();
            if (num_blocks > 0) {
                // 计算块平均值的平均值
                double block_mean_density = std::accumulate(block_avg_density.begin(), block_avg_density.end(), 0.0) / num_blocks;
                double block_mean_W = std::accumulate(block_avg_W.begin(), block_avg_W.end(), 0.0) / num_blocks;

                // 计算块平均值的标准差（用于误差估计）
                double block_sum_sq_diff_dens = 0.0;
                double block_sum_sq_diff_W = 0.0;
                for (int i = 0; i < num_blocks; i++) {
                    block_sum_sq_diff_dens += (block_avg_density[i] - block_mean_density) * (block_avg_density[i] - block_mean_density);
                    block_sum_sq_diff_W += (block_avg_W[i] - block_mean_W) * (block_avg_W[i] - block_mean_W);
                }
                double block_std_density = std::sqrt(block_sum_sq_diff_dens / num_blocks);
                double block_std_W = std::sqrt(block_sum_sq_diff_W / num_blocks);
                // 标准误差 = 标准差 / sqrt(块数)
                double block_se_density = block_std_density / std::sqrt(num_blocks);
                double block_se_W = block_std_W / std::sqrt(num_blocks);

                // 计算过剩化学势 μ_ex 的块统计（只考虑有效值）
                std::vector<double> valid_mu_ex;
                for (double mu_ex : block_avg_mu_ex) {
                    if (!std::isnan(mu_ex)) {
                        valid_mu_ex.push_back(mu_ex);
                    }
                }
                int num_valid_mu_ex = valid_mu_ex.size();
                double block_mean_mu_ex = 0.0;
                double block_std_mu_ex = 0.0;
                double block_se_mu_ex = 0.0;

                if (num_valid_mu_ex > 0) {
                    block_mean_mu_ex = std::accumulate(valid_mu_ex.begin(), valid_mu_ex.end(), 0.0) / num_valid_mu_ex;

                    double block_sum_sq_diff_mu_ex = 0.0;
                    for (double mu_ex : valid_mu_ex) {
                        block_sum_sq_diff_mu_ex += (mu_ex - block_mean_mu_ex) * (mu_ex - block_mean_mu_ex);
                    }
                    block_std_mu_ex = std::sqrt(block_sum_sq_diff_mu_ex / num_valid_mu_ex);
                    block_se_mu_ex = block_std_mu_ex / std::sqrt(num_valid_mu_ex);
                }

                // 计算理想化学势 μ_id 的块统计（只考虑有效值）
                std::vector<double> valid_mu_id;
                for (double mu_id : block_avg_mu_id) {
                    if (!std::isnan(mu_id)) {
                        valid_mu_id.push_back(mu_id);
                    }
                }
                int num_valid_mu_id = valid_mu_id.size();
                double block_mean_mu_id = 0.0;
                double block_std_mu_id = 0.0;
                double block_se_mu_id = 0.0;

                if (num_valid_mu_id > 0) {
                    block_mean_mu_id = std::accumulate(valid_mu_id.begin(), valid_mu_id.end(), 0.0) / num_valid_mu_id;

                    double block_sum_sq_diff_mu_id = 0.0;
                    for (double mu_id : valid_mu_id) {
                        block_sum_sq_diff_mu_id += (mu_id - block_mean_mu_id) * (mu_id - block_mean_mu_id);
                    }
                    block_std_mu_id = std::sqrt(block_sum_sq_diff_mu_id / num_valid_mu_id);
                    block_se_mu_id = block_std_mu_id / std::sqrt(num_valid_mu_id);
                }

                std::cout << "\n===== 块平均统计结果 =====" << std::endl;
                std::cout << "块数: " << num_blocks << std::endl;
                std::cout << "密度块平均值: " << block_mean_density << " ± " << block_se_density << " (标准误差)" << std::endl;
                std::cout << "密度块间标准差: " << block_std_density << std::endl;
                if (block_mean_W > 0.0) {
                    std::cout << "插入权重块平均值 <W>: " << block_mean_W << " ± " << block_se_W << " (标准误差)" << std::endl;
                    std::cout << "插入权重块间标准差: " << block_std_W << std::endl;
                } else {
                    std::cout << "插入权重块平均值 <W>: 0 (无有效样本)" << std::endl;
                }
                if (num_valid_mu_ex > 0) {
                    std::cout << "过剩化学势块平均值 μ_ex: " << block_mean_mu_ex << " ± " << block_se_mu_ex << " (标准误差)" << std::endl;
                    std::cout << "过剩化学势块间标准差: " << block_std_mu_ex << std::endl;
                } else {
                    std::cout << "过剩化学势块平均值 μ_ex: 无有效数据 (所有块的 <W> ≤ 0)" << std::endl;
                }
                if (num_valid_mu_id > 0) {
                    std::cout << "理想化学势块平均值 μ_id: " << block_mean_mu_id << " ± " << block_se_mu_id << " (标准误差)" << std::endl;
                    std::cout << "理想化学势块间标准差: " << block_std_mu_id << std::endl;
                } else {
                    std::cout << "理想化学势块平均值 μ_id: 无有效数据 (所有块的密度 ≤ 0)" << std::endl;
                }

                // 计算化学势（基于整体平均值或块平均值？这里使用整体平均值）
                if (avg_W > 0.0) {
                    double mu_calc = -std::log(avg_W) + std::log(avg_density / M);
                    std::cout << "\n===== 化学势计算 =====" << std::endl;
                    std::cout << "整体平均插入权重 <W>: " << avg_W << std::endl;
                    std::cout << "整体有效采样点数 (W>0): " << count_nonzero_W << std::endl;
                    std::cout << "输入化学势 mu_input: " << mu_b << std::endl;
                    std::cout << "计算化学势 mu_calc: " << mu_calc << std::endl;
                    std::cout << "差值 (mu_calc - mu_input): " << (mu_calc - mu_b) << std::endl;
                } else {
                    std::cout << "\n警告: 未获得有效的插入权重 (W>0)" << std::endl;
                }

                // 保存结果到文件（文件名包含化学势）
                std::stringstream filename_stream;
                filename_stream << "bulk_density_results_mu_" << std::fixed << std::setprecision(3) << mu_b << ".txt";
                std::string filename = filename_stream.str();

                std::ofstream outfile(filename);
                if (outfile) {
                    outfile << "# mu_b: " << mu_b << std::endl;
                    outfile << "# ave_density: " << avg_density << std::endl;
                    outfile << "# std_density: " << std_density << std::endl;
                    outfile << "# sample_times: " << all_density_samples.size() << std::endl;
                    outfile << "# block_count: " << block_avg_density.size() << std::endl;
                    outfile << "# block_mean_density: " << block_mean_density << std::endl;
                    outfile << "# block_se_density: " << block_se_density << std::endl;
                    outfile << "# block_std_density: " << block_std_density << std::endl;
                    if (avg_W > 0.0) {
                        double mu_calc = -std::log(avg_W) + std::log(avg_density / M);
                        outfile << "# <W>: " << avg_W << std::endl;
                        outfile << "# mu_ex: " << -std::log(avg_W) << std::endl;
                        outfile << "# mu_id: " << std::log(avg_density / M) << std::endl;
                        outfile << "# mu_calc: " << mu_calc << std::endl;
                        outfile << "# (mu_calc - mu_input): " << (mu_calc - mu_b) << std::endl;
                        outfile << "# block_mean_W: " << block_mean_W << std::endl;
                        outfile << "# block_se_W: " << block_se_W << std::endl;
                        outfile << "# block_std_W: " << block_std_W << std::endl;
                        if (num_valid_mu_ex > 0) {
                            outfile << "# block_mean_mu_ex: " << block_mean_mu_ex << std::endl;
                            outfile << "# block_se_mu_ex: " << block_se_mu_ex << std::endl;
                            outfile << "# block_std_mu_ex: " << block_std_mu_ex << std::endl;
                        }
                        // 输出每个块的 mu_ex 值
                        outfile << "# block_mu_ex_values:";
                        for (double mu_ex_val : block_avg_mu_ex) {
                            if (std::isnan(mu_ex_val)) {
                                outfile << " nan";
                            } else {
                                outfile << " " << mu_ex_val;
                            }
                        }
                        outfile << std::endl;
                        if (num_valid_mu_id > 0) {
                            outfile << "# block_mean_mu_id: " << block_mean_mu_id << std::endl;
                            outfile << "# block_se_mu_id: " << block_se_mu_id << std::endl;
                            outfile << "# block_std_mu_id: " << block_std_mu_id << std::endl;
                        }
                        // 输出每个块的 mu_id 值
                        outfile << "# block_mu_id_values:";
                        for (double mu_id_val : block_avg_mu_id) {
                            if (std::isnan(mu_id_val)) {
                                outfile << " nan";
                            } else {
                                outfile << " " << mu_id_val;
                            }
                        }
                        outfile << std::endl;
                    }
                    outfile << "# final_EPS_TRANS: " << EPS_TRANS << std::endl;
                    outfile << "# final_trans_acceptance: " << bulk_ring.get_translation_acceptance() << std::endl;
                    outfile << "# final_rotation_acceptance: " << bulk_ring.get_rotation_acceptance() << std::endl;
                    outfile << "# final_insertion_acceptance: " << bulk_ring.get_insert_acceptance() << std::endl;
                    outfile << "# final_deletion_acceptance: " << bulk_ring.get_delete_acceptance() << std::endl;
                    outfile.close();
                    std::cout << "结果已保存到 " << filename << std::endl;
                }
            } // 结束 if (num_blocks > 0)

            std::cout << "\n===== 整体样本统计结果 =====" << std::endl;
            std::cout << "总采样点数: " << all_density_samples.size() << std::endl;
            std::cout << "整体平均单体密度: " << avg_density << std::endl;
            std::cout << "整体样本标准差: " << std_density << std::endl;
            std::cout << "相对波动: " << (std_density / avg_density * 100) << "%" << std::endl;
            std::cout << "整体平均聚合物数: " << avg_polymer_count << std::endl;
            std::cout << "当前聚合物数: " << bulk_ring.get_N_now() << std::endl;
            std::cout << "当前单体数: " << bulk_ring.get_MN_now() << std::endl;
        } else {
            std::cout << "\nWarning: No sample data collected." << std::endl;
        }

        // 保存最终系统状态
        /*
        std::stringstream config_filename_stream;
        config_filename_stream << "final_config_mu_" << std::fixed << std::setprecision(3) << mu_b << ".dat";
        std::string config_filename = config_filename_stream.str();

        std::ofstream config_file(config_filename);
        if (config_file) {
            const double** r_total = bulk_ring.get_r_total();
            int MN_now = bulk_ring.get_MN_now();
            
            config_file << std::fixed << std::setprecision(8);
            for (int i = 0; i < MN_now; i++) {
                config_file << r_total[i][0] << " " << r_total[i][1] << " " << r_total[i][2] << std::endl;
            }
            config_file.close();
            std::cout << "最终系统状态已保存到 " << config_filename << std::endl;
            std::cout << "总单体数: " << MN_now << std::endl;
            std::cout << "总聚合物数: " << bulk_ring.get_N_now() << std::endl;
            std::cout << "最终密度: " << bulk_ring.get_rho_now() << std::endl;
        } else {
            std::cout << "警告: 无法保存最终系统状态到文件 " << config_filename << std::endl;
        }

        std::cout << "\n测试完成!" << std::endl;
        */

    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}