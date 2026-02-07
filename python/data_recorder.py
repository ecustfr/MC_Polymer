# -*- coding: utf-8 -*-
"""
Data Recorder module for Monte Carlo simulation results.

This module provides a class for saving data to files.
"""

import numpy as np


class DataRecorder:
    """
    Class for saving data to files.
    
    This class provides a simple interface for saving arrays to files.
    """
    
    def __init__(self):
        """
        Initialize the DataRecorder.
        """
        pass
    
    def save(self, file_path: str, data: np.ndarray, format: str = "txt"):
        """
        Save data to a file.
        
        Args:
            file_path: Path to the output file
            data: Array data to save
            format: File format (default: "txt")
        """
        format = format.lower()
        
        if format == "txt":
            self.save_txt(file_path, data)
        elif format == "npy":
            self.save_npy(file_path, data)
        else:
            raise ValueError(f"Unsupported file format: {format}")
    
    def save_txt(self, file_path: str, data: np.ndarray):
        """
        Save data to a text file.
        
        Args:
            file_path: Path to the output text file
            data: Array data to save
        """
        try:
            np.savetxt(file_path, data, delimiter=' ')
            print(f"Data saved to {file_path}")
        except IOError as e:
            print(f"Error saving data to {file_path}: {e}")
    
    def save_npy(self, file_path: str, data: np.ndarray):
        """
        Save data to a numpy binary file.
        
        Args:
            file_path: Path to the output numpy file
            data: Array data to save
        """
        try:
            np.save(file_path, data)
            print(f"Data saved to {file_path}")
        except IOError as e:
            print(f"Error saving data to {file_path}: {e}")
    
    def save_multiple_arrays(self, file_path: str, arrays: list, names: list = None):
        """
        Save multiple arrays to a single file.
        
        Args:
            file_path: Path to the output file
            arrays: List of arrays to save
            names: List of array names (optional)
        """
        try:
            if names is None:
                names = [f"array_{i}" for i in range(len(arrays))]
            
            with open(file_path, 'w') as f:
                for name, array in zip(names, arrays):
                    f.write(f"# {name}\n")
                    np.savetxt(f, array, delimiter=' ')
                    f.write("\n")
            print(f"Multiple arrays saved to {file_path}")
        except IOError as e:
            print(f"Error saving multiple arrays to {file_path}: {e}")