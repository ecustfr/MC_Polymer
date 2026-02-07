#!/usr/bin/env python3
"""
External potential calculation module

This module provides functions for calculating and generating external potentials
for polymer simulations. All potential calculations use the [0, H] range for z-coordinates.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def calculate_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar):
    """
    Calculate external potential at position z

    Parameters:
    -----------
    z : float
        Position along z-axis (should be in [0, box_size_z])
    An : list of float
        Amplitudes of sine components
    phi_n : list of float
        Phases of sine components
    box_size_z : float
        Box size in z-direction
    Vlin_par : list of list of float
        Parameters for linear segments (2x4 matrix)
    x_tar : list of list of float
        Target intervals for linear segments (2x4 matrix)

    Returns:
    --------
    float
        External potential value at position z
    """
    # Calculate Vext1: sum of sine functions
    Vext1 = 0.0
    for n in range(4):
        Vext1 += An[n] * np.sin(2 * np.pi * (n+1) * z / box_size_z + phi_n[n])

    # Calculate Vext2: sum of linear functions
    Vext2 = 0.0
    for n in range(4):
        x1, x2 = x_tar[0][n], x_tar[1][n]
        if x1 <= z <= x2:
            V_temp = Vlin_par[0][n] + (Vlin_par[1][n] - Vlin_par[0][n]) * (z - x1) / (x2 - x1)
            Vext2 += V_temp

    # Total potential
    V = Vext1 + Vext2

    # Boundary condition (using [0, box_size_z] range)
    if z < 0.5 or z > box_size_z - 0.5:
        V = 1e20

    return V

def create_potential_function(An, phi_n, box_size_z, Vlin_par, x_tar, C=1.0):
    """
    Create potential function with given parameters

    Parameters:
    -----------
    An : list of float
        Amplitudes of sine components
    phi_n : list of float
        Phases of sine components
    box_size_z : float
        Box size in z-direction
    Vlin_par : list of list of float
        Parameters for linear segments (2x4 matrix)
    x_tar : list of list of float
        Target intervals for linear segments (2x4 matrix)
    C : float, optional
        Amplitude scaling factor

    Returns:
    --------
    function
        Potential function that takes z as input
    """
    def potential(z):
        return C * calculate_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar)

    return potential

def generate_vext_params(H, seed=None):
    """
    Generate external potential parameters

    Parameters:
    -----------
    H : float
        Box size in z-direction
    seed : int, optional
        Random seed for reproducibility

    Returns:
    --------
    dict
        Dictionary of external potential parameters
    """
    rng = np.random.default_rng(seed)  # Using independent random stream
    
    # Ensure correct sorting
    x_tar = np.sort(rng.random((2, 4)) * H, axis=0)
    
    return {
        "An": rng.normal(0, np.sqrt(2.5), 4).tolist(),
        "phi_n": (2 * np.pi * rng.random(4)).tolist(),
        "Vlin_par": rng.normal(0, 2, (2, 4)).tolist(),
        "x_tar": x_tar.tolist(),
        "C": 0.25
    }

def plot_potential(An, phi_n, box_size_z, Vlin_par, x_tar, C=1.0, output_file=None):
    """
    Plot external potential distribution

    Parameters:
    -----------
    An : list of float
        Amplitudes of sine components
    phi_n : list of float
        Phases of sine components
    box_size_z : float
        Box size in z-direction
    Vlin_par : list of list of float
        Parameters for linear segments (2x4 matrix)
    x_tar : list of list of float
        Target intervals for linear segments (2x4 matrix)
    C : float, optional
        Amplitude scaling factor
    output_file : str, optional
        Output file path for saving the plot

    Returns:
    --------
    None
    """
    # Generate z values in [0, box_size_z]
    z_min = 0
    z_max = box_size_z
    z_values = np.linspace(z_min, z_max, 1000)
    
    # Calculate potential values
    potential_values = []
    for z in z_values:
        V = C * calculate_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar)
        # Cap potential at a reasonable value for plotting
        V = min(V, 50.0)
        potential_values.append(V)
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(z_values, potential_values, 'b-', linewidth=2)
    plt.title('External Potential Distribution')
    plt.xlabel('z Position')
    plt.ylabel('Potential Energy')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(z_min, z_max)
    plt.ylim(-10,10)
    
    # Add plot info
    plt.text(0.02, 0.95, f'Box size (H): {box_size_z}', transform=plt.gca().transAxes, 
             bbox=dict(facecolor='white', alpha=0.8))
    plt.text(0.02, 0.90, f'External potential: custom', transform=plt.gca().transAxes, 
             bbox=dict(facecolor='white', alpha=0.8))
    
    # Save plot if output file is specified
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Potential plot saved to: {output_file}")
    
    # Show plot
    plt.show()
