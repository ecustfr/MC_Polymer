# -*- coding: utf-8 -*-
"""
Radius of Gyration calculation module for Monte Carlo simulation results.

This module provides classes and functions for calculating and analyzing
radius of gyration values from Monte Carlo simulation data.
"""

import numpy as np
from typing import Any, List


class RadiusOfGyrationCalculator:
    """
    Class for calculating and analyzing radius of gyration values.
    
    Attributes:
        accumulated_rg: Accumulated radius of gyration values per polymer
        sample_count: Number of samples taken
        n_polymers: Number of polymers in the system
        m_per_polymer: Monomers per polymer
    """
    
    def __init__(self, n_polymers: int, m_per_polymer: int):
        """
        Initialize the RadiusOfGyrationCalculator.
        
        Args:
            n_polymers: Number of polymers in the system
            m_per_polymer: Number of monomers per polymer
        """
        self.n_polymers = n_polymers
        self.m_per_polymer = m_per_polymer
        self.accumulated_rg = np.zeros(n_polymers)
        self.sample_count = 0
    
    def reset(self):
        """Reset the accumulated RG values and sample count."""
        self.accumulated_rg.fill(0.0)
        self.sample_count = 0
    
    def calculate_rg_single_polymer(self, r_polymer: np.ndarray) -> float:
        """
        Calculate radius of gyration for a single polymer.
        
        Args:
            r_polymer: Array of monomer positions with shape (M, 3)
            
        Returns:
            rg: Radius of gyration
        """
        # Calculate center of mass
        r_cm = np.mean(r_polymer, axis=0)
        
        # Calculate squared distances from center of mass
        sq_distances = np.sum((r_polymer - r_cm) ** 2, axis=1)
        
        # Calculate radius of gyration squared and take square root
        rg_sq = np.mean(sq_distances)
        return np.sqrt(rg_sq)
    
    def calculate_rg_all_polymers(self, sim: Any) -> np.ndarray:
        """
        Calculate radius of gyration for all polymers in the simulation.
        
        Args:
            sim: Simulation object with r_total (positions), get_N_now(), and get_M() methods
            
        Returns:
            rg_values: Array of radius of gyration values for each polymer
        """
        # 使用 getter 方法获取 N_now 和 M
        n_now = sim.get_N_now()
        m = sim.get_M()
        
        rg_values = np.zeros(n_now)
        
        for polymer_index in range(n_now):
            # Extract monomer positions for current polymer
            start_idx = polymer_index * m
            end_idx = start_idx + m
            r_polymer = sim.r_total[start_idx:end_idx]
            
            # Calculate RG for this polymer
            rg_values[polymer_index] = self.calculate_rg_single_polymer(r_polymer)
        
        return rg_values
    
    def accumulate_rg(self, sim: Any):
        """
        Accumulate RG values from current simulation state.
        
        Args:
            sim: Simulation object with r_total, N_now, and M attributes
        """
        current_rg = self.calculate_rg_all_polymers(sim)
        
        # Ensure we don't exceed the initialized size
        min_size = min(len(self.accumulated_rg), len(current_rg))
        self.accumulated_rg[:min_size] += current_rg[:min_size]
        self.sample_count += 1
    
    def get_average_rg(self) -> np.ndarray:
        """
        Get the average RG values over all accumulated samples.
        
        Returns:
            average_rg: Average radius of gyration values per polymer
        """
        if self.sample_count == 0:
            return np.zeros(self.n_polymers)
        return self.accumulated_rg / self.sample_count
    
    def get_overall_average_rg(self) -> float:
        """
        Get the overall average RG across all polymers and samples.
        
        Returns:
            overall_average_rg: Overall average radius of gyration
        """
        if self.sample_count == 0:
            return 0.0
        total_average = np.mean(self.accumulated_rg)
        return total_average / self.sample_count
    
    def save_rg_values(self, file_path: str):
        """
        Save the average RG values to a file.
        
        Args:
            file_path: Path to output file
        """
        average_rg = self.get_average_rg()
        try:
            np.savetxt(file_path, average_rg, delimiter=' ')
        except IOError as e:
            print(f"Error saving RG values to {file_path}: {e}")
    
    def calculate_rg_distribution(self, bins: int = 50, range_min: float = 0.0, range_max: float = 20.0) -> tuple:
        """
        Calculate the distribution of RG values.
        
        Args:
            bins: Number of bins in the histogram
            range_min: Minimum RG value for the histogram
            range_max: Maximum RG value for the histogram
            
        Returns:
            hist: Histogram of RG values
            bin_edges: Bin edges for the histogram
        """
        average_rg = self.get_average_rg()
        hist, bin_edges = np.histogram(average_rg, bins=bins, range=(range_min, range_max), density=True)
        return hist, bin_edges