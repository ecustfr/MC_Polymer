#!/usr/bin/env python3
"""
Python version of Vext_z.m for external potential calculation
"""

import numpy as np
from pymcpolymer import MuVT_MC_LinearPolymer


def calculate_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar):
    """
    Calculate external potential at position z
    
    Parameters:
    -----------
    z : float
        Position along z-axis
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
    
    # Boundary condition
    half_box = box_size_z / 2
    if z < -half_box + 0.5 or z > half_box - 0.5:
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


def example_usage():
    """
    Example usage of the potential function with C++ interface
    """
    # Example parameters (similar to MATLAB code)
    np.random.seed(42)  # For reproducibility
    
    # Parameters
    phi_n = 2 * np.pi * np.random.rand(4)
    An = np.random.normal(0, np.sqrt(2.5), 4)
    Vlin_par = np.random.normal(0, 2, (2, 4))
    box_size_z = 10.0  # Example box size
    x_tar = np.sort(np.random.rand(2, 4), axis=0) * box_size_z - 0.5 * box_size_z
    C = 1.0
    
    # Create potential function
    potential_func = create_potential_function(An, phi_n, box_size_z, Vlin_par, x_tar, C)
    
    # Example: Test the potential function
    test_z = 0.0
    print(f"Potential at z={test_z}: {potential_func(test_z)}")
    
    # Example: How to use with C++ interface
    print("\nTo use with C++ interface:")
    print("1. Create a polymer system")
    print("   polymer = pymcpolymer.MuVT_MC_LinearPolymer(...)")
    print("2. Set the external potential")
    print("   polymer.set_external_potential(potential_func, 'custom_potential')")
    print("3. Run simulations as usual")


if __name__ == "__main__":
    example_usage()
