# -*- coding: utf-8 -*-
"""
System analyzer modules for Monte Carlo simulation results.

This module provides two classes for analyzing simulation data using custom functions:
1. PolymerAnalyzer: For calculating per-polymer averages (like radius of gyration)
2. DistributionAnalyzer: For calculating system-wide distributions (like density profiles)
"""

import numpy as np
from typing import Any, Callable,Union

class GlobalPropertyAnalyzer:
    """
    负责统计系统的全局性质。
    支持标量（float）或固定形状的数组（ndarray）。
    """
    def __init__(self, name: str = "Property"):
        self.name = name
        self.acc_val = 0.0          # 累积值
        self.total_samples = 0      # 采样计数
        self._is_initialized = False

    def reset(self):
        """重置所有统计数据"""
        self.acc_val = 0.0
        self.total_samples = 0
        self._is_initialized = False

    def accumulate(self, sim: Any, func: Callable):
        """
        通过传入的函数采样并累积数据。
        """
        val = func(sim)
        
        # 自动处理数组初始化：如果是第一次采样且结果是数组
        if not self._is_initialized:
            if isinstance(val, np.ndarray):
                self.acc_val = np.zeros_like(val, dtype=float)
            self._is_initialized = True

        self.acc_val += val
        self.total_samples += 1

    @property
    def average(self) -> Union[float, np.ndarray]:
        """
        获取采样平均值: $\langle A \rangle = \frac{1}{N} \sum_{i=1}^{N} A_i$
        """
        if self.total_samples == 0:
            return 0.0
        return self.acc_val / self.total_samples

    def __repr__(self):
        return f"<GlobalPropertyAnalyzer(name='{self.name}', samples={self.total_samples})>"


class PolymerAnalyzer:
    """
    Class for calculating per-polymer averages using custom functions.
    """
    
    def __init__(self, name,M: int, N: int, trace_or_not: bool = True, max_samples: int = 1000):
        """
        Initialize the PolymerAnalyzer.
        
        Args:
            M: Number of monomers per polymer
            N: Maximum number of polymers
            max_samples: Maximum number of samples to store
            trace_or_not: Whether to enable sampling of calculated variables
        """
        self.name = name
        self.M = M
        self.N = N
        self.max_samples = max_samples
        
        self.accumulated_values = np.zeros(N)
        self.sample_count = 0
        self.trace_or_not = trace_or_not
        self.trace_count = 0
        
        if self.trace_or_not:
            self.calculated_variables = np.zeros((max_samples, N))
        else:
            self.calculated_variables = np.zeros((1, N))
   
    
    def reset(self):
        """Reset all accumulated data and sample count."""
        self.accumulated_values.fill(0.0)
        self.calculated_variables.fill(0.0)
        self.sample_count = 0
        self.trace_count = 0
    
    def calculate(self, sim: Any, calculation_func: Callable[[np.ndarray], float]) -> np.ndarray:
        """
        Calculate values for each polymer using a custom function.
        
        Args:
            sim: Simulation object with r_total (positions) and get_N_now() methods
            calculation_func: Function that takes monomer coordinates (shape: (M, 3)) and returns a float
            
        Returns:
            values: Array of values for each polymer (shape: (N_now,))
        """
        n_now = sim.get_N_now()
        result = np.zeros(n_now)
        r_total = sim.r_total
        
        for polymer_idx in range(n_now):
            start_idx = polymer_idx * self.M
            end_idx = start_idx + self.M
            r_polymer = r_total[start_idx:end_idx]
            result[polymer_idx] = calculation_func(r_polymer)
        
        return result
    
    def accumulate(self, sim: Any, calculation_func: Callable[[np.ndarray], float]):
        """
        Accumulate per-polymer values using a custom function.
        
        Args:
            sim: Simulation object with r_total and get_N_now() methods
            calculation_func: Function that takes monomer coordinates (shape: (M, 3)) and returns a float
        """
        current_values = self.calculate(sim, calculation_func)
        min_size = min(len(self.accumulated_values), len(current_values))
        
        # Store currently calculated variables
        if self.trace_or_not:
            if self.trace_count < self.max_samples:
                self.calculated_variables[self.trace_count, :min_size] = current_values[:min_size]
                self.trace_count += 1
            else:
                raise ValueError("Exceeded maximum number of samples")
    
        self.accumulated_values[:min_size] += current_values[:min_size]
        self.sample_count += 1
    
    def get_average_values(self) -> np.ndarray:
        """
        Get the average of accumulated per-polymer values.
        
        Returns:
            avg_values: Array representing the average values for each polymer (shape: (N,))
        """
        if self.sample_count == 0:
            return np.zeros(self.N)
        return self.accumulated_values / self.sample_count
    
    def get_calculated_variables(self, use_all_samples: bool = True) -> np.ndarray:
        """
        Get the stored calculated variables.
        
        Args:
            use_all_samples: If True, returns all stored samples; if False, returns only valid samples
            
        Returns:
            variables: Array of stored variables (shape: (max_samples, N) or (trace_count, N))
        """
        if use_all_samples or self.trace_count >= self.max_samples:
            return self.calculated_variables
        return self.calculated_variables[:self.trace_count]
        

class DistributionAnalyzer:
    """
    一个更健壮的高分子系统分布分析器
    """
    def __init__(self, dz: float, bins: int, name):
        self.dz = dz
        self.bins = bins
        self.name = name 
        # 初始化存储空间
        self.acc_val = np.zeros(bins, dtype=float)  # 累积值
        self.acc_cnt = np.zeros(bins, dtype=float)  # 针对每个 bin 的计数
        self.total_samples = 0        

    def reset(self):
        """彻底重置所有数据"""
        self.acc_val.fill(0.0)
        self.acc_cnt.fill(0.0)
        self.total_samples = 0

    def accumulate(self, sim: Any, func: Callable[[Any, float, int], np.ndarray]):
        """
        通用的累积方法。支持函数返回 1D (仅数据) 或 2D (计数, 数据) 结果。
        """
        profile = func(sim, self.dz, self.bins)
        
        if profile.ndim == 1:
            # 如果只返回了数据，我们认为全系统均匀计数 +1
            self.acc_val += profile
            self.acc_cnt += 1
        elif profile.ndim == 2 and profile.shape[1] == 2:
            # 如果返回的是 [bins, 2]，第一列为计数，第二列为数据
            self.acc_cnt += profile[:, 0]
            self.acc_val += profile[:, 1]
        else:
            raise ValueError(f"不支持的返回形状: {profile.shape}。请返回 (bins,) 或 (bins, 2)")
            
        self.total_samples += 1

    @property
    def average(self) -> np.ndarray:
        """
        获取平均分布。使用 property 伪装成属性，调用更方便。
        """
        return np.divide(
            self.acc_val, 
            self.acc_cnt, 
            out=np.zeros_like(self.acc_val), 
            where=self.acc_cnt != 0
        )



def cal_polymer_rg(r_polymer: np.ndarray) -> float:
    """
    计算单条聚合物的回旋半径
    """
    r_cm = np.mean(r_polymer, axis=0)
    sq_distances = np.sum((r_polymer - r_cm) ** 2, axis=1)
    rg_sq = np.mean(sq_distances)
    return np.sqrt(rg_sq)

def cal_density_profile(sim: Any, dz: float, n_bins: int) -> np.ndarray:
    """
    计算系统的密度分布
    """
    profile = np.zeros((n_bins, 2))
    
    # 2. 获取所有粒子的坐标
    # 假设 r 的形状是 (N_total, 3)
    r = sim.r_total
    M = sim.get_M()
    
    # 3. 核心：使用切片步长提取所有链的第 monomer_index 个单体
    # r[start:end:step] -> 从 monomer_index 开始，每隔 M 个取一个，取第 2 列(z)
    z_coords = r[:, 2]
    
    # 4. 统计直方图
    # range 必须从 0 到 n_bins * dz，确保跟 Analyzer 对齐
    counts, _ = np.histogram(z_coords, bins=n_bins, range=(0, n_bins * dz))
    
    # 5. 填充结果
    bin_volume = dz * sim.get_box_xy() * sim.get_box_xy()
    profile[:, 1] = counts.astype(float)/bin_volume  # 存储分布数值
    profile[:, 0] = 1.0                   # 这一帧的采样权重记为 1
    
    return profile



def cal_mono_density_profile(sim: Any, dz: float, n_bins: int, monomer_index: int) -> np.ndarray:
    """
    计算系统中特定单体的密度分布。数据结构假设为：[链1单体0, 1...M-1, 链2单体0, 1...M-1, ...]
    
    Args:
        sim: 模拟对象，需提供 get_r_total() 返回所有粒子的坐标 (N_total, 3)
        dz: z轴分段宽度
        n_bins: 总分段数
        monomer_index: 目标单体在链内的索引 (0 到 M-1)
        M: 每条高分子链包含的单体总数
    """
    # 1. 初始化 (n_bins, 2) 数组
    profile = np.zeros((n_bins, 2))
    
    # 2. 获取所有粒子的坐标
    # 假设 r 的形状是 (N_total, 3)
    r = sim.r_total
    M = sim.get_M()
    
    # 3. 核心：使用切片步长提取所有链的第 monomer_index 个单体
    # r[start:end:step] -> 从 monomer_index 开始，每隔 M 个取一个，取第 2 列(z)
    z_coords = r[monomer_index::M, 2]
    
    # 4. 统计直方图
    # range 必须从 0 到 n_bins * dz，确保跟 Analyzer 对齐
    counts, _ = np.histogram(z_coords, bins=n_bins, range=(0, n_bins * dz))
    
    # 5. 填充结果
    bin_volume = dz * sim.get_box_xy() * sim.get_box_xy()
    profile[:, 1] = counts.astype(float)/bin_volume  # 存储分布数值
    profile[:, 0] = 1.0                   # 这一帧的采样权重记为 1
    
    return profile

def cal_G_profile(sim:Any,dz:float,n_bins:int,monomer_index:int, k_max = 10 , insert_time=20) -> np.ndarray:
    """
    计算系统的G值分布
    """

    profile = np.zeros((n_bins,2)) # 第一行存次数
    h = sim.get_H()
    # insert_time = 20
    
    z = np.random.rand(insert_time) * (h-1)+0.5
    for pos in z:
        Npos = int(pos/dz)
        g_value = sim.get_G_insert(pos, monomer_index, k_max)
        profile[Npos,1] += g_value
        profile[Npos,0] += 1


    return profile

def cal_Wz_profile(sim:Any,dz:float,n_bins:int,monomer_index:int, k_max = 10 , insert_time=20) -> np.ndarray:
    """
    计算系统的Wz值分布
    """

    profile = np.zeros((n_bins,2)) # 第一行存次数
    h = sim.get_H()
    # insert_time = 20
    
    z = np.random.rand(insert_time) * (h-1)+0.5
    for pos in z:
        Npos = int(pos/dz)
        g_value = sim.get_Wz_insert(pos, monomer_index, k_max)
        profile[Npos,1] += g_value
        profile[Npos,0] += 1


    return profile

def cal_W(sim:Any,k_max:int):
    return sim.get_W_insert(k_max)

    
    
if __name__ == "__main__":
    pass