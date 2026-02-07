# -*- coding: utf-8 -*-
"""
Data analysis module for Monte Carlo simulation results.

This module provides classes and functions for analyzing Monte Carlo simulation
results, including radius of gyration calculation, density profile analysis,
and data recording.
"""

import numpy as np
import os
from typing import List, Dict, Tuple, Any, Optional


class SimulationData:
    """
    Container class for simulation data storage and analysis.
    
    Attributes:
        block: Current simulation block number
        sample_address: File path for sample data output
        times: Number of samples taken
    """
    
    def __init__(self, block: int, sample_address: str):
        """
        Initialize the SimulationData object.
        
        Args:
            block: Current simulation block number
            sample_address: File path for sample data output
        """
        self.block = block
        self.sample_address = sample_address
        self.times = 0
        self.sample_file = None
        self._open_sample_file()
    
    def _open_sample_file(self):
        """Open the sample file for writing."""
        try:
            self.sample_file = open(self.sample_address, 'w')
        except IOError as e:
            print(f"Error opening file {self.sample_address}: {e}")
    
    def close(self):
        """Close the sample file if it's open."""
        if self.sample_file and not self.sample_file.closed:
            self.sample_file.close()
    
    def __del__(self):
        """Destructor to ensure file is closed."""
        self.close()
    
    def begin_block(self):
        """Write block beginning information to file and console."""
        print(f"--------------------------------block: {self.block}------------------------------------------------")
        if self.sample_file and not self.sample_file.closed:
            self.sample_file.write(f"block:{self.block} begin\n")
    
    def end_block(self, new_block: int):
        """Write block ending information to file and console, then update block number."""
        print("---------------------------------------------end-----------------------------------------------------------")
        if self.sample_file and not self.sample_file.closed:
            self.sample_file.write(f"block:{self.block}end\n")
        self.block = new_block


def record_simulation_trace(sim: Any, time: int, sample_file: Optional[Any] = None, flush: bool = True):
    """
    Record simulation trace data to file.
    
    Args:
        sim: Simulation object containing r_total, MN_now attributes
        time: Current simulation time step
        sample_file: File object to write to, or None to use print
        flush: Whether to flush the file buffer immediately
    """
    data_str = f"{time} {sim.r_total.shape} {sim.MN_now} 3\n"
    if sample_file:
        sample_file.write(data_str)
        if flush:
            sample_file.flush()
    else:
        print(data_str, end="")


def calculate_density_profile(sim: Any, dz: float, n_bins: int) -> np.ndarray:
    """
    Calculate density profile along z-axis.
    
    Args:
        sim: Simulation object with r_total (positions) and MN_now (number of monomers)
        dz: Bin width along z-axis
        n_bins: Number of bins in the profile
        
    Returns:
        rho_profile: Density profile array
    """
    rho_profile = np.zeros(n_bins)
    
    # Calculate bin indices for each monomer's z-coordinate
    z_coords = sim.r_total[:, 2]  # Get all z-coordinates
    bin_indices = (z_coords / dz).astype(int)
    
    # Ensure bin indices are within bounds
    valid_indices = bin_indices[bin_indices < n_bins]
    
    # Count monomers in each bin
    for idx in valid_indices:
        rho_profile[idx] += 1
    
    return rho_profile


def calculate_radius_of_gyration(r_polymer: np.ndarray) -> float:
    """
    Calculate radius of gyration for a single polymer.
    
    Args:
        r_polymer: Array of monomer positions with shape (M, 3)
        
    Returns:
        rg: Radius of gyration
    """
    M = r_polymer.shape[0]
    if M == 0:
        return 0.0
    
    # Calculate center of mass
    r_cm = np.mean(r_polymer, axis=0)
    
    # Calculate squared distances from center of mass
    sq_distances = np.sum((r_polymer - r_cm) ** 2, axis=1)
    
    # Calculate radius of gyration
    rg_sq = np.mean(sq_distances)
    rg = np.sqrt(rg_sq)
    
    return rg


def calculate_rg_for_all_polymers(sim: Any) -> np.ndarray:
    """
    Calculate radius of gyration for all polymers in the simulation.
    
    Args:
        sim: Simulation object with r_total, get_N_now(), and get_M() methods
        
    Returns:
        rg_values: Array of radius of gyration values for each polymer
    """
    # 使用 getter 方法获取 N_now 和 M
    n_now = sim.get_N_now()
    m = sim.get_M()
    
    rg_values = np.zeros(n_now)
    
    for polymer_index in range(n_now):
        # Extract positions for current polymer
        start_idx = polymer_index * m
        end_idx = start_idx + m
        r_polymer = sim.r_total[start_idx:end_idx]
        
        # Calculate radius of gyration
        rg_values[polymer_index] = calculate_radius_of_gyration(r_polymer)
    
    return rg_values


def calculate_center_z(r_polymer: np.ndarray) -> float:
    """
    Calculate average z-coordinate (center of mass in z-direction) for a polymer.
    
    Args:
        r_polymer: Array of monomer positions with shape (M, 3)
        
    Returns:
        center_z: Average z-coordinate
    """
    if r_polymer.shape[0] == 0:
        return 0.0
    return np.mean(r_polymer[:, 2])


def calculate_rc_z_for_all_polymers(sim: Any) -> np.ndarray:
    """
    Calculate average z-coordinate for all polymers in the simulation.
    
    Args:
        sim: Simulation object with r_total, get_N_now(), and get_M() methods
        
    Returns:
        z_centers: Array of average z-coordinates for each polymer
    """
    # 使用 getter 方法获取 N_now 和 M
    n_now = sim.get_N_now()
    m = sim.get_M()
    
    z_centers = np.zeros(n_now)
    
    for polymer_index in range(n_now):
        # Extract positions for current polymer
        start_idx = polymer_index * m
        end_idx = start_idx + m
        r_polymer = sim.r_total[start_idx:end_idx]
        
        # Calculate center z-coordinate
        z_centers[polymer_index] = calculate_center_z(r_polymer)
    
    return z_centers


def save_data_to_file(file_path: str, data: Any, mode: str = 'w', delimiter: str = ' '):
    """
    Save data to file using numpy's savetxt function.
    
    Args:
        file_path: Path to output file
        data: Data to save (numpy array)
        mode: File mode ('w' for write, 'a' for append)
        delimiter: Delimiter for data columns
    """
    try:
        np.savetxt(file_path, data, delimiter=delimiter)
    except IOError as e:
        print(f"Error saving data to {file_path}: {e}")


def load_data_from_file(file_path: str, delimiter: str = ' ') -> np.ndarray:
    """
    Load data from file using numpy's loadtxt function.
    
    Args:
        file_path: Path to input file
        delimiter: Delimiter for data columns
        
    Returns:
        data: Loaded data as numpy array
    """
    try:
        return np.loadtxt(file_path, delimiter=delimiter)
    except IOError as e:
        print(f"Error loading data from {file_path}: {e}")
        return np.array([])
