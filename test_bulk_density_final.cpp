#include "Bulk_RingPolymer.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdlib> // for atof
#include <iomanip> // for std::setprecision
#include <sstream> // for std::stringstream

int main(int argc, char* argv[]) {
    // 测试Bulk_RingPolymer类，统计平衡时系统的密度
    std::cout << "=== 测试 Bulk_RingPolymer 密度计算 ===" << std::endl;

    // 检查命令行参数
    if (argc < 2) {
        std::cerr << "用法: " << argv[0] << " <化学势>" << std::endl;
        std::cerr << "示例: " << argv[0] << " 1.23" << std::endl;
        std::cerr << "注意: 化学势为必选参数" << std::endl;
        return 1;
    }

    // 从命令行读取化学势
    double mu_b = std::atof(argv[1]);

    // 设置参数（根据用户要求）
    std::string configuration = "input/init_config_Z/TrivialN64M8H20.00rho0.10.dat";
    int M = 8; // 每个聚合物的单体数
    int init_N = 64; // 初始聚合物数 (64/8=8)
    double rho = 0.1; // 目标密度
    double box_size[3] = {20.0, 20.0, 20.0}; // 体相盒子大小
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
    std::cout << "  体相参考密度 (rho_b): " << rho_b << std::endl;

    try {
        // 创建Bulk_RingPolymer对象
        std::cout << "\n创建 Bulk_RingPolymer 对象..." << std::endl;
        Bulk_RingPolymer bulk_ring(configuration, M, init_N, rho, box_size, rcut, max_N, mu_b, rho_b);

        // 初始化第二步
        std::cout << "执行第二步初始化..." << std::endl;
        bulk_ring.init_second();

        // 设置模拟参数
        double EPS_TRANS = 0.1;
        double ROT_RATIO = 0.3;
        int K_MAX = 5;
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

        // 采样阶段
        int production_steps = 20000;
        int sample_interval = 10;
        std::vector<double> density_samples;
        std::vector<double> polymer_count_samples;
        std::vector<double> insertion_weight_samples;

        std::cout << "\n开始采样阶段 (" << production_steps << " 步)..." << std::endl;
        for (int step = 0; step < production_steps; step++) {
            if (step % 2000 == 0) {
                std::cout << "  采样步数: " << step << "/" << production_steps << std::endl;
            }
            bulk_ring.run_simulation_MuVT(1); // 运行1步MuVT模拟

            if (step % sample_interval == 0) {
                // 记录当前密度和聚合物数
                double density = bulk_ring.get_rho_now();
                int polymer_count = bulk_ring.get_N_now();
                density_samples.push_back(density);
                polymer_count_samples.push_back(polymer_count);

                // 计算插入权重 W
                double W = bulk_ring.get_W_insert_ring(K_MAX);
                insertion_weight_samples.push_back(W);
            }
        }

        // 结束block
        bulk_ring.end_block(0);

        // 统计结果
        if (!density_samples.empty()) {
            // 计算平均值
            double avg_density = std::accumulate(density_samples.begin(), density_samples.end(), 0.0) / density_samples.size();
            double avg_polymer_count = std::accumulate(polymer_count_samples.begin(), polymer_count_samples.end(), 0.0) / polymer_count_samples.size();

            // 计算插入权重的平均值（排除零值）
            double sum_W = 0.0;
            int count_nonzero_W = 0;
            for (double W : insertion_weight_samples) {
                    sum_W += W;
                    count_nonzero_W++;
            }
            double avg_W = (count_nonzero_W > 0) ? sum_W / count_nonzero_W : 0.0;

            // 计算标准差
            double sum_sq_diff = 0.0;
            for (double d : density_samples) {
                sum_sq_diff += (d - avg_density) * (d - avg_density);
            }
            double std_density = std::sqrt(sum_sq_diff / density_samples.size());

            std::cout << "\n===== 密度统计结果 =====" << std::endl;
            std::cout << "采样点数: " << density_samples.size() << std::endl;
            std::cout << "平均单体密度: " << avg_density << std::endl;
            std::cout << "标准差: " << std_density << std::endl;
            std::cout << "相对波动: " << (std_density / avg_density * 100) << "%" << std::endl;
            std::cout << "平均聚合物数: " << avg_polymer_count << std::endl;
            std::cout << "当前聚合物数: " << bulk_ring.get_N_now() << std::endl;
            std::cout << "当前单体数: " << bulk_ring.get_MN_now() << std::endl;

            // 计算化学势
            if (avg_W > 0.0) {
                double mu_calc = -std::log(avg_W) + std::log(avg_density / M);
                std::cout << "\n===== 化学势计算 =====" << std::endl;
                std::cout << "平均插入权重 <W>: " << avg_W << std::endl;
                std::cout << "有效采样点数 (W>0): " << count_nonzero_W << std::endl;
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
                outfile << "# Bulk_RingPolymer 密度采样数据" << std::endl;
                outfile << "# 化学势 (mu_b): " << mu_b << std::endl;
                outfile << "# 平均密度: " << avg_density << std::endl;
                outfile << "# 标准差: " << std_density << std::endl;
                outfile << "# 采样点数: " << density_samples.size() << std::endl;
                if (avg_W > 0.0) {
                    double mu_calc = -std::log(avg_W) + std::log(rho / M);
                    outfile << "# 平均插入权重 <W>: " << avg_W << std::endl;
                    outfile << "# 有效W采样点数: " << count_nonzero_W << std::endl;
                    outfile << "# 计算化学势 mu_calc: " << mu_calc << std::endl;
                    outfile << "# 差值 (mu_calc - mu_input): " << (mu_calc - mu_b) << std::endl;
                }
                outfile << "# 步数 密度 聚合物数 插入权重" << std::endl;
                for (size_t i = 0; i < density_samples.size(); i++) {
                    outfile << i * sample_interval << " " << density_samples[i] << " " << polymer_count_samples[i] << " " << insertion_weight_samples[i] << std::endl;
                }
                outfile.close();
                std::cout << "结果已保存到 " << filename << std::endl;
            }
        } else {
            std::cout << "未采集到密度数据。" << std::endl;
        }

        // 保存最终系统状态
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

    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}