#!/usr/bin/env python3
"""
External potential calculation module

This module provides functions for calculating and generating external potentials
for polymer simulations. All potential calculations use the [0, H] range for z-coordinates.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def calculate_custom_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar):
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
        return C * calculate_custom_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar)

    return potential

def create_potential_function_from_params(params, box_size_z):
    """
    Create potential function from parameters dictionary

    Parameters:
    -----------
    params : dict
        Dictionary containing potential parameters, must include 'potential_type'
    box_size_z : float
        Box size in z-direction

    Returns:
    --------
    function
        Potential function that takes z as input
    """
    potential_type = params.get('potential_type', 'custom')

    if potential_type == 'step':
        # Step potential parameters
        boundaries = params['boundaries']
        potentials = params['potentials']
        C = params.get('C', 1.0)
        return create_step_potential_function(boundaries, potentials, box_size_z, C)
    elif potential_type == 'custom':
        # Custom potential parameters (sine + linear)
        An = params['An']
        phi_n = params['phi_n']
        Vlin_par = params['Vlin_par']
        x_tar = params['x_tar']
        C = params.get('C', 1.0)
        return create_potential_function(An, phi_n, box_size_z, Vlin_par, x_tar, C)
    else:
        raise ValueError(f"Unknown potential type: {potential_type}")

def generate_vext_params(H, potential_type="custom", n_steps=3, seed=None, **kwargs):
    """
    Generate external potential parameters

    Parameters:
    -----------
    H : float
        Box size in z-direction
    potential_type : str, optional
        Type of potential: "custom" (sine+linear) or "step"
    n_steps : int, optional
        Number of steps for step potential (if potential_type="step")
    seed : int, optional
        Random seed for reproducibility
    **kwargs : dict
        Additional parameters passed to specific potential generators
        For step potential: potential_mean, potential_std

    Returns:
    --------
    dict
        Dictionary of external potential parameters
    """
    if potential_type == "step":
        # Extract step-specific parameters from kwargs
        potential_mean = kwargs.get('potential_mean', 0.0)
        potential_std = kwargs.get('potential_std', 2.0)
        C = kwargs.get('C', 0.1)  # Default C value for step potential
        return generate_step_vext_params(
            H,
            n_steps=n_steps,
            seed=seed,
            potential_mean=potential_mean,
            potential_std=potential_std,
            C=C
        )
    elif potential_type == "custom":
        rng = np.random.default_rng(seed)  # Using independent random stream

        # Ensure correct sorting
        x_tar = np.sort(rng.random((2, 4)) * H, axis=0)

        return {
            "potential_type": "custom",
            "An": rng.normal(0, np.sqrt(2.5), 4).tolist(),
            "phi_n": (2 * np.pi * rng.random(4)).tolist(),
            "Vlin_par": rng.normal(0, 2, (2, 4)).tolist(),
            "x_tar": x_tar.tolist()
        }
    else:
        raise ValueError(f"Unknown potential type: {potential_type}")

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
        V = C * calculate_custom_Vext(z, An, phi_n, box_size_z, Vlin_par, x_tar)
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
    plt.ylim(-10*C-0.5,10*C+0.5)
    
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

def calculate_step_Vext(z, boundaries, potentials, box_size_z):
    """
    Calculate step potential at position z

    Parameters:
    -----------
    z : float
        Position along z-axis (should be in [0, box_size_z])
    boundaries : list of float
        Boundary points of step intervals (must include 0 and box_size_z)
    potentials : list of float
        Potential values for each interval (length must be len(boundaries)-1)
    box_size_z : float
        Box size in z-direction

    Returns:
    --------
    float
        Step potential value at position z
    """
    # Hard wall boundary condition (same as existing system)
    if z < 0.5 or z > box_size_z - 0.5:
        return 1e20

    # Validate input parameters
    if len(boundaries) < 2:
        raise ValueError("boundaries must have at least 2 elements (0 and H)")
    if len(potentials) != len(boundaries) - 1:
        raise ValueError(f"potentials length ({len(potentials)}) must be boundaries length ({len(boundaries)}) - 1")
    if boundaries[0] != 0.0:
        raise ValueError(f"First boundary must be 0.0, got {boundaries[0]}")
    if boundaries[-1] != box_size_z:
        raise ValueError(f"Last boundary must be box_size_z ({box_size_z}), got {boundaries[-1]}")

    # Find which interval z belongs to using binary search
    # boundaries is sorted: [0, x1, x2, ..., xn, H]
    left, right = 0, len(boundaries) - 1
    while left <= right:
        mid = (left + right) // 2
        if mid < len(boundaries) - 1 and boundaries[mid] <= z < boundaries[mid + 1]:
            return potentials[mid]
        elif z < boundaries[mid]:
            right = mid - 1
        else:
            left = mid + 1

    # Handle edge case: z == boundaries[-1] (z == H)
    if z == boundaries[-1]:
        return potentials[-1]

    # This should not happen if boundaries cover [0, H]
    raise ValueError(f"Position z={z} not found in boundaries {boundaries}")

def create_step_potential_function(boundaries, potentials, box_size_z, C=1.0):
    """
    Create step potential function with given parameters

    Parameters:
    -----------
    boundaries : list of float
        Boundary points of step intervals (must include 0 and box_size_z)
    potentials : list of float
        Potential values for each interval (length must be len(boundaries)-1)
    box_size_z : float
        Box size in z-direction
    C : float, optional
        Amplitude scaling factor

    Returns:
    --------
    function
        Step potential function that takes z as input
    """
    # Validate parameters
    if len(boundaries) < 2:
        raise ValueError("boundaries must have at least 2 elements")
    if len(potentials) != len(boundaries) - 1:
        raise ValueError(f"potentials length ({len(potentials)}) must be boundaries length ({len(boundaries)}) - 1")
    if not all(boundaries[i] < boundaries[i+1] for i in range(len(boundaries)-1)):
        raise ValueError("boundaries must be strictly increasing")
    if boundaries[0] != 0.0:
        raise ValueError(f"First boundary must be 0.0, got {boundaries[0]}")
    if boundaries[-1] != box_size_z:
        raise ValueError(f"Last boundary must be box_size_z ({box_size_z}), got {boundaries[-1]}")

    def potential(z):
        return C * calculate_step_Vext(z, boundaries, potentials, box_size_z)

    return potential

def generate_step_vext_params(H, n_steps=3, seed=None, potential_mean=0.0, potential_std=2.0, C=0.1):
    """
    Generate step potential parameters

    Parameters:
    -----------
    H : float
        Box size in z-direction
    n_steps : int, optional
        Number of steps (intervals), default 3
    seed : int, optional
        Random seed for reproducibility
    potential_mean : float, optional
        Mean value for potential heights, default 0.0
    potential_std : float, optional
        Standard deviation for potential heights, default 2.0
    C : float, optional
        Amplitude scaling factor, default 0.1

    Returns:
    --------
    dict
        Dictionary of step potential parameters
    """
    rng = np.random.default_rng(seed)

    # Generate random boundaries (excluding 0 and H)
    if n_steps < 1:
        raise ValueError("n_steps must be at least 1")

    # Generate n_steps-1 random points between 0 and H
    if n_steps > 1:
        internal_points = np.sort(rng.random(n_steps - 1) * H)
    else:
        internal_points = np.array([])

    # Create boundaries array: [0, internal_points..., H]
    boundaries = np.concatenate([[0.0], internal_points, [H]])

    # Generate random potentials for each interval
    potentials = rng.normal(potential_mean, potential_std, n_steps).tolist()

    return {
        "potential_type": "step",
        "boundaries": boundaries.tolist(),
        "potentials": potentials,
        "C": C  # Use the provided scaling factor
    }

def plot_step_potential(boundaries, potentials, box_size_z, C=1.0, output_file=None):
    """
    Plot step potential distribution

    Parameters:
    -----------
    boundaries : list of float
        Boundary points of step intervals (must include 0 and box_size_z)
    potentials : list of float
        Potential values for each interval (length must be len(boundaries)-1)
    box_size_z : float
        Box size in z-direction
    C : float, optional
        Amplitude scaling factor
    output_file : str, optional
        Output file path for saving the plot

    Returns:
    --------
    None
    """
    # Validate parameters
    if len(boundaries) < 2:
        raise ValueError("boundaries must have at least 2 elements")
    if len(potentials) != len(boundaries) - 1:
        raise ValueError(f"potentials length ({len(potentials)}) must be boundaries length ({len(boundaries)}) - 1")

    # Generate z values in [0, box_size_z]
    z_min = 0
    z_max = box_size_z
    z_values = np.linspace(z_min, z_max, 1000)

    # Calculate potential values
    potential_values = []
    for z in z_values:
        try:
            V = C * calculate_step_Vext(z, boundaries, potentials, box_size_z)
            # Cap potential at a reasonable value for plotting
            if V > 50.0:
                V = 50.0
            potential_values.append(V)
        except ValueError:
            potential_values.append(50.0)  # For plotting purposes

    # Create step plot
    plt.figure(figsize=(10, 6))

    # For step plot, we need to create step-shaped data
    step_z = []
    step_V = []
    for i in range(len(potentials)):
        # Add boundary points
        step_z.extend([boundaries[i], boundaries[i+1]])
        step_V.extend([potentials[i] * C, potentials[i] * C])

    plt.step(step_z, step_V, 'b-', linewidth=2, where='post')

    # Add hard wall regions shading
    plt.axvspan(0, 0.5, alpha=0.2, color='red', label='Hard wall (z<0.5)')
    plt.axvspan(box_size_z - 0.5, box_size_z, alpha=0.2, color='red', label='Hard wall (z>H-0.5)')

    plt.title('Step Potential Distribution')
    plt.xlabel('z Position')
    plt.ylabel('Potential Energy')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(z_min, z_max)
    plt.ylim(min(potentials)*C-0.5, max(potentials)*C+0.5)

    # Add plot info
    plt.text(0.02, 0.95, f'Box size (H): {box_size_z}', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8))
    plt.text(0.02, 0.90, f'Potential type: step ({len(potentials)} steps)', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8))
    plt.text(0.02, 0.85, f'Scaling factor C: {C}', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8))

    # Add legend
    plt.legend(loc='upper right')

    # Save plot if output file is specified
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Step potential plot saved to: {output_file}")

    # Show plot
    plt.show()


def _create_and_set_potential(Vext_params, box_size_z, mc_sys, external_potential):
    """
    Create potential function from parameters and set it to the MC system.

    Parameters:
    -----------
    Vext_params : dict
        Potential parameters dictionary
    box_size_z : float
        Box size in z-direction
    mc_sys : object
        Monte Carlo simulation system object
    external_potential : str
        External potential type ('custom')

    Returns:
    --------
    function
        Created potential function
    """
    # Create potential function using the unified factory function
    potential_func = create_potential_function_from_params(Vext_params, box_size_z)

    # Set external potential
    mc_sys.set_external_potential(potential_func, f'{external_potential}_potential')
    print(f"{external_potential.capitalize()} external potential set successfully")

    return potential_func


def _save_potential_values_to_file(potential_func, box_size_z, output_dir, external_potential, Vext_params):
    """
    Calculate potential values and save them to a file.

    Parameters:
    -----------
    potential_func : function
        Potential function that takes z as input
    box_size_z : float
        Box size in z-direction
    output_dir : str
        Output directory for saving the file
    external_potential : str
        External potential type ('custom')
    Vext_params : dict
        Potential parameters dictionary for the file header

    Returns:
    --------
    None
    """
    dz = 0.05
    z_values = np.arange(0, box_size_z + dz/2, dz)
    potential_values = []
    for z in z_values:
        V = potential_func(z)
        potential_values.append(V)

    # Save to output directory
    output_file = os.path.join(output_dir, f'{external_potential}_potential_values.txt')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"# {external_potential.capitalize()} External Potential Values\n")
        f.write(f"# Box size H = {box_size_z}\n")
        f.write(f"# Potential type: {external_potential}\n")
        f.write(f"# Parameters: {Vext_params}\n")
        f.write("# z\tpotential\n")
        for z, V in zip(z_values, potential_values):
            f.write(f"{z:.6f}\t{V:.6f}\n")

    print(f"{external_potential.capitalize()} potential values saved to: {output_file}")
    print(f"Number of points: {len(z_values)}, dz={dz}")


def _generate_potential_plot_file(Vext_params, box_size_z, output_dir, external_potential):
    """
    Generate and save a plot of the potential.

    Parameters:
    -----------
    Vext_params : dict
        Potential parameters dictionary
    box_size_z : float
        Box size in z-direction
    output_dir : str
        Output directory for saving the plot
    external_potential : str
        External potential type ('custom')

    Returns:
    --------
    None
    """
    plot_file = os.path.join(output_dir, f'{external_potential}_potential_plot.png')
    potential_type = Vext_params.get('potential_type', 'custom')

    if potential_type == 'custom':
        # Custom potential plot
        An = Vext_params['An']
        phi_n = Vext_params['phi_n']
        Vlin_par = Vext_params['Vlin_par']
        x_tar = Vext_params['x_tar']
        C = Vext_params.get('C', 1.0)
        plot_potential(An, phi_n, box_size_z, Vlin_par, x_tar, C=C, output_file=plot_file)
    elif potential_type == 'step':
        # Step potential plot
        boundaries = Vext_params['boundaries']
        potentials = Vext_params['potentials']
        C = Vext_params.get('C', 1.0)
        plot_step_potential(boundaries, potentials, box_size_z, C=C, output_file=plot_file)
    else:
        raise ValueError(f"Unknown potential type in Vext_params: {potential_type}")

    print(f"{external_potential.capitalize()} potential plot saved to: {plot_file}")


def setup_external_potential(mc_sys, input_params, config, output_dir):
    """
    Set up external potential for simulation based on configuration.

    Parameters:
    -----------
    mc_sys : object
        Monte Carlo simulation system object (MuVT_MC_RingPolymer or MuVT_MC_LinearPolymer)
    input_params : dict
        Input parameters dictionary containing 'external_potential' and 'H'
    config : dict
        Full configuration dictionary containing Vext_params if external_potential is 'custom'
    output_dir : str
        Output directory for saving potential values and plots

    Returns:
    --------
    None

    Notes:
    ------
    - Supports 'custom' and 'None' external_potential types from input_params
    - For 'custom' type, Vext_params must contain 'potential_type' ('custom' or 'step')
    - For 'custom' potential_type: requires An, phi_n, Vlin_par, x_tar parameters
    - For 'step' potential_type: requires boundaries and potentials parameters
    """
    # Set external potential if specified
    external_potential = input_params['external_potential']
    print(f"External potential type: {external_potential}")

    if external_potential == 'custom':
        Vext_params = config['Vext_params']
        box_size_z = input_params['H']  # Use H as box size in z-direction

        # Create and set the potential function, then save values to file
        _save_potential_values_to_file(
            _create_and_set_potential(Vext_params, box_size_z, mc_sys, external_potential),
            box_size_z, output_dir, external_potential, Vext_params
        )

        # Generate and save potential plot
        _generate_potential_plot_file(Vext_params, box_size_z, output_dir, external_potential)

    elif external_potential == 'None':
        mc_sys.set_external_potential("hs_wall")
        print("No external potential will be used")
