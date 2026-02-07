#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

// === 参数设置 ===
const double BOX_X = 19.595917942265423; // 盒子 X 方向长度 (请根据您的模拟设置修改)
const double BOX_Y = 19.595917942265423; // 盒子 Y 方向长度 (请根据您的模拟设置修改)
const double BOX_Z = 10.0; // 盒子 Z 方向长度 (请根据您的模拟设置修改)
const double DIAMETER = 1.0; // 粒子直径 (即重叠判断的阈值)
const double DIAMETER_SQ = DIAMETER * DIAMETER;

struct Particle {
    double x, y, z;
};

// 计算考虑周期性边界条件的距离平方
double distance_sq_pbc(const Particle& p1, const Particle& p2) {
    double dx = std::abs(p1.x - p2.x);
    double dy = std::abs(p1.y - p2.y);
    double dz = std::abs(p1.z - p2.z);

    if (dx > BOX_X * 0.5) dx = BOX_X - dx;
    if (dy > BOX_Y * 0.5) dy = BOX_Y - dy;

    return dx * dx + dy * dy + dz * dz;
}

int main() {
    std::string filename = "mc_sim_Linear_H10.0_mub2.75_M6_trajectory.xyz"; // 您的轨迹文件名
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }

    std::string line;
    int current_step = -1;
    int num_polymers = 0;
    std::vector<Particle> particles;
    bool reading_particles = false;

    std::cout << "Starting overlap check..." << std::endl;
    std::cout << "Box Size: " << BOX_X << " x " << BOX_Y << " x " << BOX_Z << std::endl;
    std::cout << "Overlap Cutoff (Diameter): " << DIAMETER << std::endl;

    while (std::getline(file, line)) {
        // 移除行首尾空格（可选）
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        // 检查是否是 Step 标题行
        if (line.find("Step:") != std::string::npos) {
            // 如果之前有读取粒子数据，先处理上一帧
            if (!particles.empty()) {
                bool overlap_found = false;
                for (size_t i = 0; i < particles.size(); ++i) {
                    for (size_t j = i + 1; j < particles.size(); ++j) {
                        // 【注意】这里假设相邻粒子（同一聚合物内的键连粒子）允许重叠或距离很近
                        // 如果您的模型是切线硬球链，相邻粒子距离固定为 1.0，不算重叠。
                        // 如果需要排除相邻粒子（例如 i 和 i+1），请取消下面这行的注释：
                        // if (j == i + 1 && (i + 1) % 6 != 0) continue; // 假设链长为 6，需根据实际 M 修改
                        
                        // 或者更通用的做法：只检查距离是否小于 DIAMETER * 0.99 以容忍浮点误差
                        // 这里我们使用严格的 DIAMETER_SQ，如果相邻粒子刚好接触，可能会误报。
                        // 建议将阈值设得稍微小一点，比如 0.95 * DIAMETER
                        
                        double dist2 = distance_sq_pbc(particles[i], particles[j]);
                        if (dist2 < 0.99 * DIAMETER_SQ) { // 使用 0.99 避免浮点误差误报接触的粒子
                             std::cout << "[Step " << current_step << "] Overlap found between particle " 
                                       << i << " and " << j << ", dist^2: " << dist2 << std::endl;
                             overlap_found = true;
                             // 找到一个重叠就跳出当前两层循环（如果只想知道这一帧有没有重叠）
                             // goto next_frame; 
                        }
                    }
                }
                if (!overlap_found) {
                    // std::cout << "[Step " << current_step << "] No overlap." << std::endl;
                }
            }

            // next_frame:; 

            // 解析新的 Step 信息
            particles.clear();
            std::stringstream ss(line);
            std::string temp;
            ss >> temp >> current_step; // 读取 "Step:" 和 数字
            
            // 跳过中间可能的逗号或文本，找到 "Polymers:"
            size_t poly_pos = line.find("Polymers:");
            if (poly_pos != std::string::npos) {
                std::stringstream ss2(line.substr(poly_pos + 9)); // 9 是 "Polymers:" 的长度
                ss2 >> num_polymers;
            }
            
            // std::cout << "Processing Step: " << current_step << ", Polymers: " << num_polymers << std::endl;
            reading_particles = true;
        } 
        else if (reading_particles) {
            // 读取坐标行
            std::stringstream ss(line);
            double x, y, z;
            if (ss >> x >> y >> z) {
                particles.push_back({x, y, z});
            }
        }
    }

    // 处理文件最后一帧（因为循环结束后最后一帧还没处理）
    if (!particles.empty()) {
        bool overlap_found = false;
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                 double dist2 = distance_sq_pbc(particles[i], particles[j]);
                 if (dist2 < 0.99 * DIAMETER_SQ) {
                      std::cout << "[Step " << current_step << "] Overlap found between particle " 
                                << i << " and " << j << ", dist^2: " << dist2 << std::endl;
                      overlap_found = true;
                 }
            }
        }
        if (!overlap_found) {
            // std::cout << "[Step " << current_step << "] No overlap." << std::endl;
        }
    }

    std::cout << "Check finished." << std::endl;
    return 0;
}