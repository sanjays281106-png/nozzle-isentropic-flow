"""
Analysis of Isentropic Flow Through a Converging-Diverging Nozzle Using Python
------------------------------------------------------------------------------
Quasi-one-dimensional, steady, isentropic flow of a calorically perfect gas
through a de Laval nozzle. Solves for Mach number, static temperature,
pressure, density, speed of sound, velocity and mass flow rate at 400
stations along the nozzle axis.

Run:  python nozzle_flow.py
Outputs (written to ./outputs/):
    nozzle_flow_results1.png   - six-panel plot of A, M, p, T, rho, mdot vs x
    program_output.txt         - printed summary values
    nozzle_flow_data.csv       - tabulated results at every station
"""

import os
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")  # works without a display; remove for interactive use
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# 1. Gas / flow properties
# ---------------------------------------------------------------------
gamma = 1.69   # specific heat ratio (value used in the case study)
R = 267.0      # specific gas constant, J/(kg-K)
T0 = 690.0     # stagnation temperature, K
p0 = 6.9e5     # stagnation pressure, Pa


# ---------------------------------------------------------------------
# 2. Nozzle geometry definition
# ---------------------------------------------------------------------
def nozzle_area(x, A_throat=1.5, A_inlet=6.90, A_exit=4.5, x_throat=0.0,
                x_start=-1.0, x_end=1.0):
    A = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi <= x_throat:
            # converging section: A_inlet -> A_throat
            frac = (xi - x_start) / (x_throat - x_start)
            A[i] = A_throat + (A_inlet - A_throat) * 0.5 * (1 + np.cos(np.pi * frac))
        else:
            # diverging section: A_throat -> A_exit
            frac = (xi - x_throat) / (x_end - x_throat)
            A[i] = A_throat + (A_exit - A_throat) * 0.5 * (1 - np.cos(np.pi * frac))
    return A


# ---------------------------------------------------------------------
# 3. Area-Mach relation (from energy + isentropic relations)
# ---------------------------------------------------------------------
def area_ratio(M, gamma):
    """Returns A/A* for a given Mach number M."""
    term = (2.0 / (gamma + 1.0)) * (1.0 + (gamma - 1.0) / 2.0 * M**2)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / M) * term**exponent


def solve_mach(A_over_Astar, gamma, supersonic):
    """
    Inverts the area-Mach relation for M, given A/A* and a branch choice.
    Uses Brent's method (robust bracketed root finder).
    """
    f = lambda M: area_ratio(M, gamma) - A_over_Astar
    if supersonic:
        return brentq(f, 1.0 + 1e-8, 15.0)
    else:
        return brentq(f, 1e-4, 1.0 - 1e-8)


# ---------------------------------------------------------------------
# 4. Build the solution along the nozzle axis
# ---------------------------------------------------------------------
def solve_nozzle(x, A, gamma, R, T0, p0):
    A_star = A.min()  # throat area = minimum area
    i_throat = np.argmin(A)
    M = np.zeros_like(x)
    for i in range(len(x)):
        ratio = A[i] / A_star
        if ratio <= 1.0 + 1e-6:
            M[i] = 1.0
        else:
            supersonic_branch = i > i_throat
            M[i] = solve_mach(ratio, gamma, supersonic=supersonic_branch)

    # Isentropic relations (from energy equation)
    T = T0 / (1.0 + (gamma - 1.0) / 2.0 * M**2)
    p = p0 / (1.0 + (gamma - 1.0) / 2.0 * M**2)**(gamma / (gamma - 1.0))
    rho = p / (R * T)
    a = np.sqrt(gamma * R * T)  # local speed of sound
    V = M * a
    return M, T, p, rho, V, A_star


# ---------------------------------------------------------------------
# 5. Run
# ---------------------------------------------------------------------
if __name__ == "__main__":
    out_dir = "outputs"
    os.makedirs(out_dir, exist_ok=True)

    x = np.linspace(-1.0, 1.0, 400)
    A = nozzle_area(x)
    M, T, p, rho, V, A_star = solve_nozzle(x, A, gamma, R, T0, p0)

    # Mass flow conservation check
    mdot = rho * A * V

    lines = [
        f"Throat area A* = {A_star:.4f}",
        f"Mass flow rate: min = {mdot.min():.6f} kg/s, "
        f"max = {mdot.max():.6f} kg/s (should be ~constant)",
        f"Exit Mach number: {M[-1]:.4f}",
        f"Exit pressure: {p[-1]:.2f} Pa (p/p0 = {p[-1]/p0:.4f})",
        f"Exit temperature: {T[-1]:.2f} K",
        f"Exit velocity: {V[-1]:.2f} m/s",
    ]
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(out_dir, "program_output.txt"), "w") as fh:
        fh.write(summary + "\n")

    # Save station-by-station results as CSV
    data = np.column_stack([x, A, M, p, T, rho, V, mdot])
    np.savetxt(os.path.join(out_dir, "nozzle_flow_data.csv"), data, delimiter=",",
               header="x_m,A_m2,Mach,p_Pa,T_K,rho_kg_m3,V_m_s,mdot_kg_s",
               comments="")

    # -------------------------------------------------------------
    # 6. Plots
    # -------------------------------------------------------------
    fig, axs = plt.subplots(3, 2, figsize=(11, 9))

    axs[0, 0].plot(x, A, 'k')
    axs[0, 0].set_title("Nozzle Area A(x)")
    axs[0, 0].set_xlabel("x")
    axs[0, 0].set_ylabel("A")

    axs[0, 1].plot(x, M, 'b')
    axs[0, 1].axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    axs[0, 1].set_title("Mach Number M(x)")
    axs[0, 1].set_xlabel("x")
    axs[0, 1].set_ylabel("M")

    axs[1, 0].plot(x, p / 1e5, 'r')
    axs[1, 0].set_title("Static Pressure p(x)")
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("p [bar]")

    axs[1, 1].plot(x, T, 'g')
    axs[1, 1].set_title("Static Temperature T(x)")
    axs[1, 1].set_xlabel("x")
    axs[1, 1].set_ylabel("T [K]")

    axs[2, 0].plot(x, rho, 'm')
    axs[2, 0].set_title("Density rho(x)")
    axs[2, 0].set_xlabel("x")
    axs[2, 0].set_ylabel("rho [kg/m^3]")

    axs[2, 1].plot(x, mdot, 'c')
    axs[2, 1].set_title("Mass Flow Rate (conservation check)")
    axs[2, 1].set_xlabel("x")
    axs[2, 1].set_ylabel("mdot [kg/s]")

    plt.tight_layout()
    fig_path = os.path.join(out_dir, "nozzle_flow_results1.png")
    plt.savefig(fig_path, dpi=150)
    print(f"\nSaved plot to {fig_path}")
