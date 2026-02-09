#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试外势打表功能
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# 添加项目根目录到 Python 路径
# 导入 pymcpolymer 模块
try:
    import pymcpolymer
    print("pymcpolymer imported successfully!")
except Exception as e:
    print(f"Error importing pymcpolymer: {e}")
    import traceback
    traceback.print_exc()

# 测试函数
def test_potential_table():
    """测试外势打表功能"""
    print("=== 测试外势打表功能 ===")
    
    # 创建一个简单的线性聚合物模拟对象
    # 使用一个不存在的配置文件，因为我们只测试外势打表功能
    config_file = "input/init_config_Z/LinearN64M10H10rhob0.1.dat"
    
    try:
        # 从 JSON 文件中读取模拟参数
        import json
        config_file_json = "config_0000_M6_Linear_H6.0_mu-1.53.json"
        with open(config_file_json, 'r') as f:
            config = json.load(f)
        
        input_params = config['input_params']
        
        # 创建模拟对象
        sim = pymcpolymer.MuVT_MC_LinearPolymer(
            input_params['configuration'],
            mu_b=input_params['mu_b'],
            M=input_params['M'],
            init_N=1,  # 使用小的 init_N 以便快速初始化
            rho_b=input_params['rho_b'],
            box_xy=input_params['box_xy'],
            H=input_params['H'],
            rcut=input_params['rcut'],
            max_N=10  # 使用小的 max_N 以便快速初始化
        )
        
        # 初始化系统
        print("初始化系统...")
        sim.init_second()
        
        # 从 JSON 文件中读取 Vext_params 参数
        import json
        config_file = "config_0000_M6_Linear_H6.0_mu-1.53.json"
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        Vext_params = config['Vext_params']
        An = Vext_params['An']
        phi_n = Vext_params['phi_n']
        Vlin_par = Vext_params['Vlin_par']
        x_tar = Vext_params['x_tar']
        C = Vext_params['C']
        H = config['input_params']['H']
        
        # 创建 custom 势函数
        print("设置 custom 外部势能...")
        from external_potential import create_potential_function
        custom_potential = create_potential_function(An, phi_n, H, Vlin_par, x_tar, C)
        
        # 设置外部势能
        sim.set_external_potential(custom_potential, "custom")
        
        # 测试 1: 导出打表数据到文件
        print("\n测试 1: 导出打表数据到文件")
        output_file = "potential_table_test.txt"
        sim.export_potential_table(output_file)
        print(f"打表数据已导出到: {output_file}")
        
        # 测试 2: 通过 Python 接口获取打表数据
        print("\n测试 2: 通过 Python 接口获取打表数据")
        pot_values = sim.get_potential_table()
        z_values = sim.get_potential_table_z()
        print(f"打表数据点数量: {len(pot_values)}")
        print(f"z 范围: {z_values[0]:.2f} 到 {z_values[-1]:.2f}")
        print(f"势能范围: {min(pot_values):.2f} 到 {max(pot_values):.2f}")
        
        # 测试 3: 验证打表数据是否正确
        print("\n测试 3: 验证打表数据是否正确")
        # 计算理论值并与打表值比较
        max_error = 0.0
        for z, pot in zip(z_values, pot_values):
            theory_pot = custom_potential(z)
            error = abs(pot - theory_pot)
            max_error = max(max_error, error)
        print(f"最大误差: {max_error:.6f}")
        if max_error < 1e-10:
            print("打表数据正确！")
        else:
            print("打表数据可能存在问题！")
        
        # 测试 4: 绘制势能分布图
        print("\n测试 4: 绘制势能分布图")
        plt.figure(figsize=(10, 6))
        plt.plot(z_values, pot_values, 'b-', label='Table values')
        
        
        plt.xlabel('z')
        plt.ylabel('V_ext(z)')
        plt.title('External Potential Distribution')
        plt.ylim(-10, 10)
        plt.legend()
        plt.grid(True)
        plt.savefig('potential_distribution.png')
        print("势能分布图已保存为: potential_distribution.png")
        
        # 清理临时文件
        if os.path.exists(output_file):
            os.remove(output_file)
        
        print("\n=== 测试完成 ===")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_potential_table()
