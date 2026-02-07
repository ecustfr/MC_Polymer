#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

// === 参数设置 (请根据你的模拟参数修改) ===
const double BOX_Z = 10.0;   // 盒子 Z 方向高度 (H)
const double MARGIN = 0.5;   // 边界预留距离 (即 0.5)

int main() {
    // 1. 设置文件路径
    std::string filename = "mc_sim_Linear_H10.0_mub2.75_M6_trajectory.xyz"; 
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }

    // 2. 计算有效范围 [0.5, H-0.5]
    double min_valid_z = MARGIN;
    double max_valid_z = BOX_Z - MARGIN;

    std::cout << "Starting Boundary Check..." << std::endl;
    std::cout << "Box Height (H): " << BOX_Z << std::endl;
    std::cout << "Valid Z Range: [" << min_valid_z << ", " << max_valid_z << "]" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    std::string line;
    int current_step = -1;
    int particle_index = 0;
    long long total_violations = 0;
    bool reading_particles = false;

    while (std::getline(file, line)) {
        // 去除首尾空白
        size_t first = line.find_first_not_of(" \t\r\n");
        if (std::string::npos == first) continue; // 空行
        line = line.substr(first, line.find_last_not_of(" \t\r\n") - first + 1);

        // 检测 Step 标题行
        if (line.find("Step:") != std::string::npos) {
            std::stringstream ss(line);
            std::string temp;
            ss >> temp >> current_step; // 读取 Step 数值
            
            reading_particles = true;
            particle_index = 0; // 重置粒子索引
        }
        else if (reading_particles) {
            // 读取坐标行 (x, y, z)
            std::stringstream ss(line);
            double x, y, z;
            
            // 尝试读取三个浮点数
            if (ss >> x >> y >> z) {
                // === 核心检查逻辑 ===
                bool violation = false;
                std::string type = "";

                if (z < min_valid_z) {
                    violation = true;
                    type = "LOWER BOUND (z < 0.5)";
                } 
                else if (z > max_valid_z) {
                    violation = true;
                    type = "UPPER BOUND (z > H-0.5)";
                }

                if (violation) {
                    std::cout << "[Step " << current_step << "] "
                              << "Particle " << particle_index << " OUT OF BOUNDS! "
                              << "z = " << z << " -> " << type << std::endl;
                    total_violations++;
                }

                particle_index++;
            }
        }
    }

    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Check Finished." << std::endl;
    if (total_violations == 0) {
        std::cout << "✅ PERFECT! No particles crossed the Z-boundaries." << std::endl;
    } else {
        std::cout << "❌ FOUND " << total_violations << " VIOLATIONS." << std::endl;
    }

    return 0;
}