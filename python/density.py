# -*- coding: utf-8 -*-
"""
Density calculation module for Monte Carlo simulation results.

This module provides classes and functions for calculating and analyzing
density profiles from Monte Carlo simulation data.
"""

import numpy as np
from typing import Any, List, Optional


class DensityProfileCalculator:
    """
    Class for calculating and analyzing density profiles from simulation data.
    
    Attributes:
        dz: Bin width along z-axis
        n_bins: Number of bins in the profile
        rho_profile: Current density profile
        sample_count: Number of samples taken
    """
    
    def __init__(self, dz: float, n_bins: int):
        """
        Initialize the DensityProfileCalculator.
        
        Args:
            dz: Bin width along z-axis
            n_bins: Number of bins in the profile
        """
        self.dz = dz
        self.n_bins = n_bins
        self.rho_profile = np.zeros(n_bins)
        self.sample_count = 0
    
    def reset(self):
        """Reset the density profile and sample count."""
        self.rho_profile.fill(0.0)
        self.sample_count = 0
    
    def calculate_profile(self, sim: Any) -> np.ndarray:
        """
        Calculate density profile for a single simulation state.
        
        Args:
            sim: Simulation object with r_total (positions) attribute
            
        Returns:
            profile: Density profile for this simulation state
        """
        profile = np.zeros(self.n_bins)
        
        # Extract z-coordinates and calculate bin indices
        # 注意：请确保 sim.r_total 只包含有效的粒子（例如已切片 sim.r_total[:sim.N_now]）
        z_coords = sim.r_total[:, 2]
        
        # 强制转换为整数索引
        bin_indices = (z_coords / self.dz).astype(int)
        
        # 直接使用 np.add.at 进行累加
        # 1. 如果 bin_indices 中有任何值超出 [0, n_bins-1]，这里会直接抛出 IndexError
        # 2. 它能正确处理多个粒子落在同一个 bin 里的情况（向量化操作，比 for 循环快）
        np.add.at(profile, bin_indices, 1)
        
        return profile
    
    """
    def calculate_profile(self, sim: Any) -> np.ndarray:

        profile = np.zeros(self.n_bins)
        
        # Extract z-coordinates and calculate bin indices
        z_coords = sim.r_total[:, 2]
        bin_indices = (z_coords / self.dz).astype(int)
        
        # Ensure bin indices are within bounds and count monomers
        valid_indices = bin_indices[(bin_indices >= 0) & (bin_indices < self.n_bins)]
        for idx in valid_indices:
            profile[idx] += 1
        
        return profile
    """
    
    def accumulate_profile(self, sim: Any):
        """
        Accumulate density profile from current simulation state.
        
        Args:
            sim: Simulation object with r_total (positions) attribute
        """
        current_profile = self.calculate_profile(sim)
        self.rho_profile += current_profile
        self.sample_count += 1
    
    def get_average_profile(self) -> np.ndarray:
        """
        Get the average density profile over all accumulated samples.
        
        Returns:
            average_profile: Average density profile
        """
        if self.sample_count == 0:
            return np.zeros(self.n_bins)
        return self.rho_profile / self.sample_count
    
    def normalize_profile(self, volume_per_bin: float) -> np.ndarray:
        """
        Normalize the density profile by volume per bin.
        
        Args:
            volume_per_bin: Volume of each bin
            
        Returns:
            normalized_profile: Normalized density profile
        """
        average_profile = self.get_average_profile()
        return average_profile / volume_per_bin
    
    def save_profile(self, file_path: str, normalized: bool = False, volume_per_bin: Optional[float] = None):
        """
        Save the density profile to a file.
        
        Args:
            file_path: Path to output file
            normalized: Whether to save normalized profile
            volume_per_bin: Volume per bin (required if normalized=True)
        """
        if normalized and volume_per_bin is None:
            raise ValueError("volume_per_bin must be provided for normalized profile")
        
        if normalized:
            profile = self.normalize_profile(volume_per_bin)
        else:
            profile = self.get_average_profile()
        
        try:
            np.savetxt(file_path, profile)
        except IOError as e:
            print(f"Error saving density profile to {file_path}: {e}")
    
    def get_z_coordinates(self) -> np.ndarray:
        """
        Get the z-coordinates corresponding to the bin centers.
        
        Returns:
            z_coords: Array of z-coordinates for each bin center
        """
        return np.arange(self.n_bins) * self.dz + self.dz / 2
