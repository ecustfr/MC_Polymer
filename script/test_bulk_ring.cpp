#include "Bulk_RingPolymer.h"
#include <iostream>

int main() {
    // 测试Bulk_RingPolymer类
    std::cout << "Testing Bulk_RingPolymer class..." << std::endl;
    
    // 设置参数
    std::string configuration = "input/init_config_Z/TrefoilN64M40H20.00rho0.10.dat"; // 使用现有文件进行测试
    int M = 40; // 每个聚合物的单体数
    int init_N = 1; // 初始聚合物数
    double rho = 0.1; // 目标密度
    double box_size[3] = {20.0, 20.0, 20.0}; // 体相盒子大小
    double rcut = 1.0; // 截断半径
    int max_N = 100; // 最大聚合物数
    
    try {
        // 创建Bulk_RingPolymer对象
        Bulk_RingPolymer bulk_ring(configuration, M, init_N, rho, box_size, rcut, max_N);
        
        // 初始化第二步
        bulk_ring.init_second();
        
        // 设置模拟参数
        bulk_ring.set_sim_parameters(0.1, 0.3, 10);
        
        // 打印参数
        bulk_ring.print_all_parameters();
        
        // 运行1000步模拟
        bulk_ring.run_simulation(1000);
        
        // 结束block
        bulk_ring.end_block(0);
        
        std::cout << "Test completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}