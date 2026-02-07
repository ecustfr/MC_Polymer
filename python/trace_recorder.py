# -*- coding: utf-8 -*-
"""
Trace Recorder module for Monte Carlo simulation results.

This module provides a class for writing trajectory files.
"""

import numpy as np


class TraceRecorder:
    """
    Class for writing trajectory files.
    
    This class only handles writing trajectory files, not reading.
    """
    
    def __init__(self, file_path: str, format: str = "xyz"):
        """
        Initialize the TraceRecorder.
        
        Args:
            file_path: Path to the output trajectory file
            format: Trajectory file format (default: "xyz")
        """
        self.file_path = file_path
        self.format = format.lower()
        self.file = None
        self._open_file()
    
    def _open_file(self):
        """Open the trajectory file for writing."""
        try:
            self.file = open(self.file_path, 'w')
        except IOError as e:
            raise IOError(f"Failed to open trajectory file {self.file_path}: {e}")
    
    def write_frame(self, positions: np.ndarray, n_polymers: int, step: int):
        """
        Write a frame to the trajectory file.
        
        Args:
            positions: Array of particle positions (shape: (N, 3) where N is the number of particles)
            n_polymers: Number of polymers in the frame
            step: Simulation step number
        """
        if self.file is None:
            raise ValueError("Trajectory file is not open")
        
        if self.format == "xyz":
            self._write_xyz_frame(positions, n_polymers, step)
        else:
            raise ValueError(f"Unsupported trajectory format: {self.format}")
    
    def _write_xyz_frame(self, positions: np.ndarray, n_polymers: int, step: int):
        """
        Write a frame in XYZ format.
        
        Args:
            positions: Array of particle positions (shape: (N, 3))
            n_polymers: Number of polymers in the frame
            step: Simulation step number
        """
        n_particles = positions.shape[0]
        
        # Write header
        self.file.write(f"{n_particles}\n")
        self.file.write(f"Step: {step}, Polymers: {n_polymers}\n")
        
        # Write particle positions
        for i, pos in enumerate(positions):
            # Only write coordinates, no particle type included
            self.file.write(f"{pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    
    def close(self):
        """
        Close the trajectory file.
        """
        if self.file is not None:
            self.file.close()
            self.file = None
    
    def __enter__(self):
        """
        Enter context manager.
        """
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context manager and close the file.
        """
        self.close()