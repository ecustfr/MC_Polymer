#!/usr/bin/env python3
"""
Matlab interface module for ring polymer theory calculations.

Provides simple functions to call Tian_slit_ringpoly_HS_Vext.m from Python.
"""

import json
import numpy as np
from pathlib import Path
import warnings

try:
    import matlab.engine
    MATLAB_ENGINE_AVAILABLE = True
except ImportError:
    MATLAB_ENGINE_AVAILABLE = False


class MatlabTheoryCalculator:
    """
    Calculator for ring polymer theory using Matlab function.

    This class provides methods to:
    1. Load configuration from JSON files
    2. Call Tian_slit_ringpoly_HS_Vext.m Matlab function
    3. Return results as numpy arrays
    """

    def __init__(self, matlab_script_path=None, use_engine=None):
        """
        Initialize the calculator.

        Parameters:
        -----------
        matlab_script_path : str or Path, optional
            Path to Tian_slit_ringpoly_HS_Vext.m
            If None, searches in current directory
        use_engine : bool, optional
            Force use of MATLAB Engine API (True) or subprocess (False)
            If None, uses engine if available
        """
        if matlab_script_path is None:
            matlab_script_path = Path(__file__).parent / "Tian_slit_ringpoly_HS_Vext.m"
        self.matlab_script_path = Path(matlab_script_path)

        if not self.matlab_script_path.exists():
            raise FileNotFoundError(f"Matlab script not found: {self.matlab_script_path}")

        self.use_engine = use_engine if use_engine is not None else MATLAB_ENGINE_AVAILABLE
        self.eng = None

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - cleanup."""
        self.close()

    def close(self):
        """Close MATLAB engine if open."""
        if self.eng is not None:
            self.eng.quit()
            self.eng = None

    def load_config(self, config_path):
        """
        Load configuration from JSON file.

        Parameters:
        -----------
        config_path : str or Path
            Path to JSON configuration file

        Returns:
        --------
        dict
            Configuration dictionary
        """
        config_path = Path(config_path)
        with open(config_path, 'r') as f:
            config = json.load(f)
        return config

    def extract_params(self, config):
        """
        Extract Matlab function parameters from configuration.

        Parameters:
        -----------
        config : dict
            Configuration dictionary

        Returns:
        --------
        tuple
            (rhob, H, m, dz, Vext_stru)
        """
        input_params = config["input_params"]
        sim_params = config["simulation_params"]
        vext_params = config["Vext_params"]

        rhob = input_params["rho_b"]
        H = input_params["H"]
        m = input_params["M"]
        dz = sim_params["dz"]

        Vext_stru = {
            "phi_n": vext_params["phi_n"],
            "An": vext_params["An"],
            "Vlin_par": vext_params["Vlin_par"],
            "x_tar": vext_params["x_tar"],
            "C": vext_params.get("C", 1.0)
        }

        return rhob, H, m, dz, Vext_stru

    def _run_via_engine(self, rhob, H, m, dz, Vext_stru):
        """Run via MATLAB Engine API."""
        if self.eng is None:
            self.eng = matlab.engine.start_matlab()
            self.eng.addpath(str(self.matlab_script_path.parent), nargout=0)

        # Convert parameters to MATLAB types
        rhob_matlab = float(rhob)
        H_matlab = float(H)
        m_matlab = float(m)
        dz_matlab = float(dz)

        # Create MATLAB struct
        Vext_stru_matlab = self.eng.struct()
        for key, value in Vext_stru.items():
            if isinstance(value, list):
                if key in ["Vlin_par", "x_tar"]:
                    value_matlab = matlab.double(value)
                else:
                    value_matlab = matlab.double(value)
            else:
                value_matlab = float(value)
            Vext_stru_matlab[key] = value_matlab

        # Call function
        rho_profile_matlab, G_matlab = self.eng.Tian_slit_ringpoly_HS_Vext(
            rhob_matlab, H_matlab, m_matlab, dz_matlab, Vext_stru_matlab,
            nargout=2
        )

        # Convert to numpy
        rho_profile = np.array(rho_profile_matlab).flatten()
        G = np.array(G_matlab).flatten()

        return rho_profile, G

    def _run_via_subprocess(self, rhob, H, m, dz, Vext_stru):
        """Run via MATLAB subprocess."""
        import tempfile
        import subprocess
        import scipy.io

        temp_dir = tempfile.mkdtemp()
        temp_script = Path(temp_dir) / "run_theory.m"

        # Create MATLAB code
        matlab_code = f"""
addpath('{self.matlab_script_path.parent}');

rhob = {rhob};
H = {H};
m = {m};
dz = {dz};

phi_n = {Vext_stru['phi_n']};
An = {Vext_stru['An']};
Vlin_par = {Vext_stru['Vlin_par']};
x_tar = {Vext_stru['x_tar']};
C = {Vext_stru.get('C', 1.0)};

Vext_stru = struct();
Vext_stru.phi_n = phi_n;
Vext_stru.An = An;
Vext_stru.Vlin_par = Vlin_par;
Vext_stru.x_tar = x_tar;
Vext_stru.C = C;

[rho_profile, G] = Tian_slit_ringpoly_HS_Vext(rhob, H, m, dz, Vext_stru);

save('{temp_dir}/results.mat', 'rho_profile', 'G');
exit;
"""

        with open(temp_script, 'w') as f:
            f.write(matlab_code)

        # Run MATLAB
        cmd = ['matlab', '-batch', f"run('{temp_script}')"]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            raise RuntimeError(f"MATLAB failed: {result.stderr}")

        # Load results
        results = scipy.io.loadmat(f'{temp_dir}/results.mat')
        rho_profile = results['rho_profile'].flatten()
        G = results['G'].flatten()

        # Cleanup
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)

        return rho_profile, G

    def calculate(self, config_path=None, config_dict=None):
        """
        Calculate density profile from configuration.

        Parameters:
        -----------
        config_path : str or Path, optional
            Path to JSON configuration file
        config_dict : dict, optional
            Configuration dictionary directly

        Returns:
        --------
        tuple
            (z, rho_profile, G) as numpy arrays

        Notes:
        ------
        Must provide either config_path or config_dict.
        """
        if config_path is not None:
            config = self.load_config(config_path)
        elif config_dict is not None:
            config = config_dict
        else:
            raise ValueError("Must provide either config_path or config_dict")

        # Extract parameters
        rhob, H, m, dz, Vext_stru = self.extract_params(config)

        # Run calculation
        if self.use_engine:
            rho_profile, G = self._run_via_engine(rhob, H, m, dz, Vext_stru)
        else:
            rho_profile, G = self._run_via_subprocess(rhob, H, m, dz, Vext_stru)

        # Generate z values
        z = np.arange(0, H + dz/2, dz)

        # Ensure same length
        min_len = min(len(z), len(rho_profile), len(G))
        z = z[:min_len]
        rho_profile = rho_profile[:min_len]
        G = G[:min_len]

        return z, rho_profile, G, config


def run_theory_from_json(config_path, output_dir=None, use_engine=None):
    """
    Simple function to run theory calculation from JSON file.

    Parameters:
    -----------
    config_path : str or Path
        Path to JSON configuration file
    output_dir : str or Path, optional
        Output directory for results
    use_engine : bool, optional
        Force use of MATLAB Engine API

    Returns:
    --------
    tuple
        (z, rho_profile, G, config)
    """
    calculator = MatlabTheoryCalculator(use_engine=use_engine)

    try:
        z, rho_profile, G, config = calculator.calculate(config_path=config_path)

        # Save results if output_dir provided
        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Save text file
            output_file = output_dir / "theory_results.txt"
            with open(output_file, 'w') as f:
                f.write("# Theory results from Tian_slit_ringpoly_HS_Vext.m\n")
                f.write(f"# H = {config['input_params']['H']}, dz = {config['simulation_params']['dz']}\n")
                f.write(f"# rhob = {config['input_params']['rho_b']}, M = {config['input_params']['M']}\n")
                f.write("# z\trho\tG\n")
                for zi, rhoi, Gi in zip(z, rho_profile, G):
                    f.write(f"{zi:.6f}\t{rhoi:.6f}\t{Gi:.6f}\n")

            # Save numpy arrays
            np.save(output_dir / "z.npy", z)
            np.save(output_dir / "rho_profile.npy", rho_profile)
            np.save(output_dir / "G.npy", G)

            print(f"Results saved to: {output_dir}")

        return z, rho_profile, G, config

    finally:
        calculator.close()


def compare_with_simulation(theory_results, simulation_results, output_dir=None):
    """
    Compare theory results with simulation results.

    Parameters:
    -----------
    theory_results : tuple
        (z_theory, rho_theory, G_theory) from calculate()
    simulation_results : tuple or dict
        Either (z_sim, rho_sim) or dict with 'z' and 'rho_profile' keys
    output_dir : str or Path, optional
        Output directory for comparison plots

    Returns:
    --------
    dict
        Comparison metrics
    """
    z_theory, rho_theory, G_theory = theory_results[:3]

    # Extract simulation data
    if isinstance(simulation_results, tuple):
        z_sim, rho_sim = simulation_results[:2]
    else:  # dict
        z_sim = simulation_results.get('z')
        rho_sim = simulation_results.get('rho_profile')

    if z_sim is None or rho_sim is None:
        raise ValueError("Simulation results must contain z and rho_profile")

    # Interpolate to common grid if needed
    if not np.array_equal(z_theory, z_sim):
        from scipy import interpolate
        interp_func = interpolate.interp1d(z_sim, rho_sim, bounds_error=False, fill_value="extrapolate")
        rho_sim_interp = interp_func(z_theory)
    else:
        rho_sim_interp = rho_sim

    # Calculate metrics
    mse = np.mean((rho_theory - rho_sim_interp)**2)
    mae = np.mean(np.abs(rho_theory - rho_sim_interp))
    max_error = np.max(np.abs(rho_theory - rho_sim_interp))
    correlation = np.corrcoef(rho_theory, rho_sim_interp)[0, 1]

    metrics = {
        'mse': mse,
        'mae': mae,
        'max_error': max_error,
        'correlation': correlation,
        'theory_mean': np.mean(rho_theory),
        'simulation_mean': np.mean(rho_sim_interp)
    }

    # Plot comparison if output_dir provided
    if output_dir is not None:
        try:
            import matplotlib.pyplot as plt

            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            plt.figure(figsize=(12, 8))

            # Plot density profiles
            plt.subplot(2, 2, 1)
            plt.plot(z_theory, rho_theory, 'b-', linewidth=2, label='Theory')
            plt.plot(z_sim, rho_sim, 'r--', linewidth=2, label='Simulation')
            plt.xlabel('z')
            plt.ylabel('Density')
            plt.title('Density Profile Comparison')
            plt.legend()
            plt.grid(True, alpha=0.3)

            # Plot difference
            plt.subplot(2, 2, 2)
            plt.plot(z_theory, rho_theory - rho_sim_interp, 'g-', linewidth=2)
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
            plt.xlabel('z')
            plt.ylabel('Difference (Theory - Sim)')
            plt.title(f'Difference (MSE={mse:.6f})')
            plt.grid(True, alpha=0.3)

            # Plot G
            plt.subplot(2, 2, 3)
            plt.plot(z_theory, G_theory, 'purple-', linewidth=2)
            plt.xlabel('z')
            plt.ylabel('G')
            plt.title('Propagator G from Theory')
            plt.grid(True, alpha=0.3)

            # Plot correlation
            plt.subplot(2, 2, 4)
            plt.scatter(rho_theory, rho_sim_interp, alpha=0.6)
            min_val = min(rho_theory.min(), rho_sim_interp.min())
            max_val = max(rho_theory.max(), rho_sim_interp.max())
            plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
            plt.xlabel('Theory Density')
            plt.ylabel('Simulation Density')
            plt.title(f'Correlation = {correlation:.4f}')
            plt.grid(True, alpha=0.3)
            plt.axis('equal')

            plt.tight_layout()
            plt.savefig(output_dir / 'theory_simulation_comparison.png', dpi=150)

            # Save metrics
            with open(output_dir / 'comparison_metrics.json', 'w') as f:
                json.dump(metrics, f, indent=4)

            print(f"Comparison plot saved to: {output_dir}/theory_simulation_comparison.png")

        except ImportError:
            print("Matplotlib not available. Skipping plots.")

    return metrics


# Example usage
if __name__ == "__main__":
    # Example 1: Simple usage
    config_file = "input/Ring_configs/config_0000_M8_Trivial_H8.0_mu0.50.json"

    try:
        # Create calculator
        with MatlabTheoryCalculator() as calculator:
            # Load and calculate
            z, rho_profile, G, config = calculator.calculate(config_path=config_file)

            print(f"Calculation completed:")
            print(f"  Config: {config_file}")
            print(f"  z range: {z[0]:.2f} to {z[-1]:.2f} ({len(z)} points)")
            print(f"  Density range: {rho_profile.min():.4f} to {rho_profile.max():.4f}")

            # Save results
            output_dir = Path("theory_results")
            output_dir.mkdir(exist_ok=True)

            np.save(output_dir / "z.npy", z)
            np.save(output_dir / "rho_profile.npy", rho_profile)
            np.save(output_dir / "G.npy", G)

            print(f"Results saved to {output_dir}/")

    except Exception as e:
        print(f"Error: {e}")