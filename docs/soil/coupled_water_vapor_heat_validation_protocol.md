# Coupled water-vapor-air-heat validation protocol

Date: 2026-07-11  
Status: preregistered design; do not tune formulas against the final curves

## 1. Purpose

This protocol separates constitutive derivatives, conservation equations,
boundary handling, and nonlinear coupling. A full natural simulation is not a
diagnostic test because compensating coefficient errors can produce a plausible
trajectory.

No coupled-physics implementation is accepted unless it passes the lower tiers
before the integrated cases.

## 2. Acceptance order

1. symbolic identities and finite-difference coefficient checks;
2. zero-dimensional closed-cell conservation checks;
3. one-equation column tests;
4. two-equation coupling tests;
5. three-equation water-vapor-air-heat tests;
6. freeze-thaw tests;
7. published experiments and field trajectories.

Every case records the source equation, units, sign convention, active
switches, boundary fluxes, storage change, nonlinear residual, retry history,
and accepted/rejected state behavior.

## 3. Coefficient-level gates

| ID | Test | Perturbation | Required invariant |
|---|---|---|---|
| C-01 | Signed VG tangent | Central difference of `theta(h)` at unsaturated negative heads | `dtheta/dh > 0`; analytic and finite-difference values agree away from branch points. |
| C-02 | Strict-SFCC identities | Independent perturbations of `h` and `T` in `theta_u=F(h)` and `theta_l=F(h+h_frz)` | Each public derivative matches its declared mathematical field; no crossed naming is allowed. |
| C-03 | Phase-volume closure | Perturb `h` and `T` in a frozen cell | Derivatives satisfy `theta_u=theta_l+(rho_i/rho_l)theta_i` and `theta_g=porosity-theta_l-theta_i`. |
| C-04 | Vapor phase-branch derivative | Perturb the free/lower/upper Kelvin regions and the frozen ice-equilibrium branch | Analytic derivatives match each active piecewise function; the IAPWS 230 K pressure matches Table 3 and one-sided checks are used at phase kinks. |
| C-05 | Dry-air density | Perturb `h`, `T`, and `Pg` independently | `rho_a` derivatives match A22-A23 and remain consistent with the active vapor branch. |
| C-06 | Water storage Jacobian | Finite difference of conserved `Mw` | `Chh` and `ChT` equal derivatives of total water mass, unfrozen and frozen. |
| C-07 | Dry-air storage Jacobian | Finite difference of conserved `Ma` | `Cah`, `CaT`, and `Caa` equal derivatives of total dry-air mass. |
| C-08 | Energy storage Jacobian | Finite difference of conserved `E` | `CTh`, `CTT`, and `CTa` equal derivatives of one total-energy function. |
| C-09 | Clapeyron continuity | Approach `Tcrit(h)` from both sides | `h+h_frz`, `theta_l`, and conserved energy are continuous; branch derivative behavior is documented. |
| C-10 | Unit conversion | Evaluate the same physical state in SI and internal centimetre/gram units | Coefficients transform exactly; no empirical factor 2, 10, or 1e4 remains unexplained. |

Coefficient gates use several representative states, not random-only testing:

- wet unfrozen, but not saturated;
- dry unfrozen above the HR floor;
- saturated HR clamp;
- dry HR floor;
- just above and below `Tcrit`;
- partially frozen with both liquid and ice nonzero;
- almost fully frozen but above residual liquid content.

## 4. Closed-cell conservation gates

These cases have no spatial gradient and no boundary or root flux.

| ID | Active physics | Required result |
|---|---|---|
| Z-01 | Water only | Total water is unchanged to roundoff. |
| Z-02 | Water plus vapor equilibrium | Phase partition may change after a temperature perturbation, but total water is unchanged. |
| Z-03 | Dry air | Total free plus dissolved dry air is unchanged. |
| Z-04 | Sensible heat | Temperature is constant without a source. |
| Z-05 | Freeze-thaw energy | Liquid and ice may repartition, but total water and total energy are unchanged. |

These tests must evaluate conserved densities directly, not only assembled
linear residuals.

## 5. Single-equation column gates

| ID | Equation isolated | Boundary and initial condition | Reference or invariant |
|---|---|---|---|
| P-01 | Isothermal liquid water | Closed top and bottom; nonuniform `h` | Monotone redistribution; exact water conservation. |
| P-02 | Saturated liquid water | Fixed head at both ends | Linear total-head profile and analytical Darcy flux. |
| P-03 | Pure heat conduction | Fixed temperatures or insulated closed column | Analytical slab solution; exact energy accounting. |
| P-04 | Vapor diffusion by `h` | Fixed uniform `T`, prescribed head gradient, liquid conductivity disabled | Flux direction follows `-rho_v,h grad(h)` and integrated vapor mass closes. |
| P-05 | Thermal vapor diffusion | Fixed uniform `h`, prescribed temperature gradient, liquid thermal flow disabled | Flux follows the active `rho_v,T` and enhancement branch. |
| P-06 | Dry-air pressure diffusion | Uniform `h,T`, prescribed `Pg` gradient | Ideal-gas diffusion direction and dry-air mass closure. |

## 6. Coupling gates

| ID | Active equations | Discriminating condition | Required result |
|---|---|---|---|
| K-01 | Water and heat, no vapor/air | Closed freezing column | Phase redistribution cannot create water; latent heat has the opposite sign to melting/freezing. |
| K-02 | Water and vapor, fixed temperature and pressure | Dry head gradient | Liquid plus vapor mass closes even when their fluxes oppose. |
| K-03 | Water and dry air, isothermal | Ponded infiltration with bounded versus vented gas | Open case tends to Richards; bounded case slows from gas compression without mass leakage. |
| K-04 | Heat and vapor, fixed gas pressure | Temperature pulse in dry soil | Vapor latent transport changes the heat front while total energy closes. |
| K-05 | Water, vapor, and heat | Milly strong-coupling problem from Zeng Sect. 4.4.2 | Profiles and flux decomposition match the published benchmark within preregistered tolerances. |

## 7. Freeze-thaw gates

| ID | Case | Source | Required comparison |
|---|---|---|---|
| FT-01 | SFCC curve and derivative | Yu 2018 A1-A5 | `theta_u`, `theta_l`, `theta_i`, and all phase derivatives against direct formulas. |
| FT-02 | Zero-flux freezing cell | Conserved-density derivation | Total water and total energy close through the transition. |
| FT-03 | One-dimensional freezing front | Yu 2018 Sect. 3.3 flux analysis | Liquid and thermal-vapor fluxes move toward the freezing front; ice impedance suppresses liquid transport in the frozen zone. |
| FT-04 | Thawing reversal | Same column with reversed thermal forcing | Latent heat and phase rates reverse without hysteretic numerical leakage unless a physical hysteresis model is explicitly active. |
| FT-05 | Ice-vapor closure under freezing | IAPWS R14-08(2011) Eq. (6), plus the active SFCC branch | Frozen `rho_v` follows ice sublimation equilibrium, `rho_v,h=0` within the branch, and no false storage Jacobian is introduced at HR or phase switches. |

## 8. Published integrated cases

| ID | Case | Model role | Gate |
|---|---|---|---|
| R-01 | Touma and Vauclin open column | Water-only Richards reference | Reuse W2-11 after W2 homogeneous gates pass. |
| R-02 | Touma and Vauclin bounded column | Coupled water-air reference | Match infiltration retardation, pressure response, profiles, and separate water/air balances. |
| R-03 | Milly heat-moisture benchmark | Coupled water-vapor-heat reference | Match published profiles on a common time/depth grid. |
| R-04 | Yu 2018 Maqu freeze-thaw period | Full frozen-process reference | Compare temperature, liquid water, ice, latent heat, and flux-component directions; do not calibrate before conservation gates pass. |
| R-05 | HYDRUS freeze-thaw mirror | Independent numerical reference | Use identical retention, thermal, grid, time, and boundary inputs; compare conserved storage and fronts before field calibration. |

## 9. Numerical and transactional requirements

Each PDE case is repeated with at least two grids and two maximum time steps.
An accepted result must satisfy all of the following:

- conservation residual decreases under refinement;
- the physical front or profile converges under refinement;
- rejected trials do not change accepted state, time, boundary accumulators, or
  diagnostics;
- boundary series are materialized at trial time and integrated over the
  accepted interval;
- convergence is checked in physically scaled state increments and conserved
  residuals, not only the Thomas-solver algebraic residual;
- no result is accepted because it matches a legacy curve while violating a
  conserved-density gate.

## 10. Implementation increments after approval

1. Introduce explicit phase-derivative auxiliaries without changing any
   solver result.
2. Add conserved water, dry-air, and energy density evaluators plus
   coefficient finite-difference tests.
3. Correct the VG sign and SFCC derivative identity at the producer, retaining
   a branch-explicit legacy adapter for comparison.
4. Make vapor density and its derivatives use the same piecewise branch.
5. Generate water and dry-air storage coefficients from conserved densities.
6. Generate energy storage coefficients from total energy and remove legacy
   frozen overrides only after coefficient and closed-cell gates pass.
7. Wire the approved gas and heat solves into the transactional point-column
   loop with accepted/rejected state tests.
8. Run the published integrated cases, then natural field trajectories.

No increment may alter the accepted Srt ownership or the W1 boundary request
contract.
