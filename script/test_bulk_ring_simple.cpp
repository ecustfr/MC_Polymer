#include "Bulk_RingPolymer.h"
#include <iostream>
#include <fstream>
#include <cmath>

// 创建一个简单的测试配置文件
void create_test_config_file(const std::string& filename, int N, int M, double box_size) {
    std::ofstream outfile(filename);
    if (!outfile) {
        throw std::runtime_error("Cannot create test configuration file.");
    }
    
    const double spacing = 1.5; // 单体间距，大于直径1.0
    
    // 生成N个环状聚合物，每个聚合物有M个单体
    for (int polymer = 0; polymer < N; polymer++) {
        // 为每个聚合物生成M个不重叠的单体位置
        for (int monomer = 0; monomer < M; monomer++) {
            // 生成不重叠的位置
            double x = polymer * spacing * 3.0 + monomer * spacing;
            double y = polymer * spacing * 3.0 + monomer * spacing;
            double z = polymer * spacing * 3.0;
            
            // 确保位置在盒子范围内
            x = x < box_size ? x : x - box_size;
            y = y < box_size ? y : y - box_size;
            z = z < box_size ? z : z - box_size;
            
            outfile << x << " " << y << " " << z << std::endl;
        }
    }
    
    outfile.close();
    std::cout << "Test configuration file created: " << filename << std::endl;
}

int main() {
    // 测试Bulk_RingPolymer类
    std::cout << "Testing Bulk_RingPolymer class with simple configuration..." << std::endl;
    
    // 创建测试配置文件
    std::string configuration = "test_bulk_ring_config.dat";
    int N = 2; // 2个聚合物
    int M = 4; // 每个聚合物4个单体
    double box_size = 10.0; // 盒子大小
    create_test_config_file(configuration, N, M, box_size);
    
    // 设置参数
    double rho = 0.1; // 目标密度
    double box_size_3d[3] = {box_size, box_size, box_size}; // 体相盒子大小
    double rcut = 1.0; // 截断半径
    int max_N = 100; // 最大聚合物数
    
    try {
        // 创建Bulk_RingPolymer对象
        Bulk_RingPolymer bulk_ring(configuration, M, N, rho, box_size_3d, rcut, max_N);
        
        // 初始化第二步
        bulk_ring.init_second();
        
        // 设置模拟参数
        bulk_ring.set_sim_parameters(0.1, 0.3, 10);
        
        // 打印参数
        bulk_ring.print_all_parameters();
        
        // 运行100步模拟
        bulk_ring.run_simulation(100);
        
        // 结束block
        bulk_ring.end_block(0);
        
        std::cout << "Test completed successfully!" << std::endl;
        
        // 清理测试文件
        std::remove(configuration.c_str());
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}