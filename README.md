# Analysis of Isentropic Flow Through a Converging–Diverging Nozzle Using Python

A case study that computes quasi-one-dimensional, steady, isentropic flow of a
calorically perfect gas through a converging–diverging (de Laval) nozzle.

## Objectives
- Define a smooth cosine-blend nozzle area distribution A(x) with a single throat.
- Invert the area–Mach relation numerically with Brent's method, choosing the subsonic branch before the throat and the supersonic branch after it.
- Compute static temperature, pressure, density, speed of sound and velocity from the isentropic relations.
- Verify the solution through mass-flow conservation (ρAV = constant), the closed-form choked mass flow and the energy equation.
- Plot A, M, p, T, ρ and ṁ along the nozzle axis.

## Project structure
```
nozzle_project/
├── nozzle_flow.py            # main program
├── requirements.txt          # Python dependencies
├── README.md
└── outputs/
    ├── nozzle_flow_results1.png   # six-panel visualization
    ├── program_output.txt         # printed summary values
    └── nozzle_flow_data.csv       # results at all 400 stations
```
**Dataset:** none. All inputs are numerical constants set at the top of `nozzle_flow.py`.

## Inputs (base case)
| Parameter | Symbol | Value |
|---|---|---|
| Specific heat ratio | γ | 1.69 |
| Gas constant | R | 267 J/(kg·K) |
| Stagnation temperature | T0 | 690 K |
| Stagnation pressure | p0 | 6.9 bar |
| Inlet / throat / exit area | A_in / A_t / A_e | 6.9 / 1.5 / 4.5 m² |
| Axial range | x | −1 to +1 m (400 stations) |

## Requirements
- Python 3.8+
- NumPy, SciPy, Matplotlib (see `requirements.txt`)

## How to run
```bash
pip install -r requirements.txt
python nozzle_flow.py
```
Results are written to the `outputs/` folder and printed to the console. The script uses the non-interactive `Agg` backend so it works on any machine; to display the figure in a window, delete the `matplotlib.use("Agg")` line and add `plt.show()` at the end.

To try other cases, edit `gamma`, `R`, `T0`, `p0` at the top of the file, or the area arguments of `nozzle_area()`.

## Functions
| Function | Purpose |
|---|---|
| `nozzle_area` | Cosine-blend area A(x) (converging half, then diverging half) |
| `area_ratio` | A/A* as a function of Mach number and γ |
| `solve_mach` | Inverts the area–Mach relation with `scipy.optimize.brentq` |
| `solve_nozzle` | Finds the throat, solves M at every station, returns M, T, p, ρ, V, A* |

## Expected results (base case)
| Quantity | Value |
|---|---|
| Throat area A* | 1.5000 m² |
| Mass flow rate (constant) | 1759.1457 kg/s |
| Exit Mach number | 3.0342 |
| Exit pressure | 20815.98 Pa (p/p0 = 0.0302) |
| Exit temperature | 165.22 K |
| Exit velocity | 828.47 m/s |

## Results plot
![Nozzle flow results](outputs/nozzle_flow_results1.png)

## Assumptions and limitations
- 1-D, steady, adiabatic, reversible flow; constant γ; no shocks, friction, heat transfer or real-gas effects.
- Area values are interpreted as m² and x as metres.
- The throat is taken as the minimum of the sampled area; because 400 stations do not include x = 0 exactly, A* = 1.50005 and ṁ is about 0.003 % above the theoretical value. An odd number of stations (e.g. 401) samples the throat exactly.
- The solution is valid only if the back pressure equals the exit pressure (fully expanded design condition).
- γ = 1.69 and R = 267 J/(kg·K) are not the values for air (γ = 1.4, R = 287 J/(kg·K)); results correspond to the values entered.

## References
- J. D. Anderson, *Modern Compressible Flow*.
- G. P. Sutton and O. Biblarz, *Rocket Propulsion Elements*.
