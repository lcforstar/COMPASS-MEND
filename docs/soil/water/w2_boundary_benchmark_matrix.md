# W2 water-boundary benchmark matrix

Date: 2026-07-11

## Common numerical baseline

The first comparison tier uses the HYDRUS-1D nonlinear acceptance rules:

- `MaxIt = 20`;
- `TolTh = 1e-4` for the maximum change in water content at unsaturated nodes;
- `TolH = 0.1 cm` for the maximum change in pressure head at saturated nodes;
- modified-Picard mixed-form storage;
- no Anderson acceleration;
- all prescribed boundary values are audited at every retry;
- accepted-step water-balance residual must be at most `0.05 mm`.

ClimaLand Newton runs remain a second reference. A case is not accepted merely
because one model converges: profiles, boundary fluxes, storage, and cumulative
balance are compared on a common grid and time axis.

## Case matrix

| ID | Capability | Reference case | Boundary definition | Comparison target | Status |
|---|---|---|---|---|---|
| W2-01 | Constant top head/state plus free drainage | ClimaLand `experiments/standalone/Soil/richards_comparison.jl`; HYDRUS mirror; Bonan clay artifact | Top `theta=0.494`, equivalent `h=-1.428821 cm`; bottom free drainage | Fresh COMPASS, ClimaLand, HYDRUS profiles and storage | Fresh three-model run complete. COMPASS-HYDRUS final RMSE `0.001442`; COMPASS balance residual `0.002107 mm`. Current ClimaLand's own Bonan artifact test fails at RMSE `0.027383`, so the artifact is historical evidence, not a golden gate. |
| W2-02 | Closed column / zero flux | ClimaLand `docs/src/tutorials/standalone/Soil/richards_equation.jl` | Zero flux at top and bottom | Exact constant storage plus internal redistribution | Next implementation target. An internal mixed-form sign test exists, but the sourced full-column run still needs a shared fixture. |
| W2-03 | Constant top and bottom flux | ClimaLand `experiments/standalone/Soil/water_conservation.jl` | Prescribed `flux_in` and `flux_out` | Exact integral storage identity and common-grid profiles | Source locked; mirror fixtures pending. Includes negative top flux for controlled evaporation. |
| W2-04 | Constant head at both ends / saturated Darcy flow | Analytical Darcy steady state plus a one-material HYDRUS mirror | Fixed top and bottom pressure heads | Flux equals analytical Darcy value; linear total-head profile | Fixture pending. This case also tests whether positive saturated pressure head is preserved by the solver. |
| W2-05 | Upward infiltration / capillary rise | HYDRUS `DATA/Examples/Direct/UPINFIL` | Zero-flux top; fixed-tension bottom; `0-10800 s` | HYDRUS nodal profiles, cumulative bottom inflow, storage | Official local case locked; COMPASS and ClimaLand mirrors pending. |
| W2-06 | Horizontal infiltration | HYDRUS `DATA/Examples/Inverse/HORIZINF` (George Vachaud data) | Horizontal column (`CosAlpha=0`) with prescribed inlet state | Observations and HYDRUS profiles; gravity term must vanish | Official local case locked; use forward parameters after separating inverse calibration inputs. |
| W2-07 | Drainage with fixed lower head | HYDRUS `DATA/Examples/Direct/1DRAINAG` | Zero-flux top; fixed-head bottom; `0-100 h` | HYDRUS drainage and profile series | Official local case locked; mirror fixtures pending. |
| W2-08 | Variable head/flux series | HYDRUS `DATA/Examples/Inverse/MULTSTEP` and atmospheric examples | Piecewise boundary schedule | Endpoint head interpolation and interval-integrated flux invariance under retries | Source locked; COMPASS already has CSV materialization tests, but a full three-model trajectory is pending. |
| W2-09 | Controlled evaporation | HYDRUS `DATA/Examples/Inverse/EVAPOR` (Ole Wendroth data) | Upward surface removal, represented in COMPASS by negative top flux | Cumulative evaporation, profiles, storage | Source locked; forward-only fixture pending. |
| W2-10 | Layered interface infiltration | Huang et al. SV62 data and ClimaLand `docs/src/tutorials/standalone/Soil/layered_soil.jl` | Fixed wet top; free drainage bottom; sampled hydraulic parameters; `0-60 min` | Observed profiles, ClimaLand, HYDRUS reference, COMPASS | The sampled-depth interpretation is now part of canonical initialization: `[soil_input].vertical_parameter_mode="sample_linear"` maps primitive `soil.csv` parameters to element centers with held endpoints, while `horizon_bounds` retains legacy tiling. The benchmark adapter no longer overwrites soil parameters. A fresh 60 min run accepted all steps, closed the COMPASS ledger to `6.32e-5 mm`, and has COMPASS-ClimaLand profile RMSE `0.01140`, `0.01114`, and `0.03021 m3/m3` at 8, 24, and 60 min. The remaining mismatch is a spatial/interface-flux discretization audit, not an input-tiling, convergence-tolerance, or mass-loss defect. Local HYDRUS 4.17 still rejects 16 materials; the bundled reference is retained. |
| W2-11 | Experimental ponded infiltration | Touma and Vauclin (1986), reproduced in Zeng (2013), *Coupled Dynamics in Soil*, Sect. 4.4.1 | 100 cm column; initial equilibrium over a water table 120 cm below the surface; top head `+2.3 cm`; open-column bottom head `-5.8 cm` | Observed water-content profiles and the published open-column solution; wetting front near 50 cm at 0.4 h | Source locked from the local book. The open column is the water-only target because the book reports agreement with the Richards approximation. The bounded-column variant is reserved for coupled air-water validation. |

## Execution order

1. W2-02 closed column: isolate conservation with no boundary input.
2. W2-04 saturated Darcy: pin fixed-head and saturated-head handling.
3. W2-03 constant flux: pin the public positive-into-domain sign and storage.
4. W2-05 upward infiltration and W2-07 drainage: cover opposite vertical flow directions.
5. W2-06 horizontal infiltration: isolate removal of gravity.
6. W2-08 and W2-09: exercise time-series interpolation and outward flux.
7. W2-11 TV86 open-column ponding: compare against an experimental strong
   wetting front after the fixed-head and homogeneous convergence gates pass.
8. W2-10 layered SV62: return after homogeneous discretization and convergence gates are stable.

No vegetation, heat, vapor, solute, or root process is enabled in these W2
fixtures. Natural atmospheric switching is validated after the laboratory
boundary matrix passes.
