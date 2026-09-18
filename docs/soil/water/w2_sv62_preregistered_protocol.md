# W2 SV62 preregistered comparison protocol

Date frozen: 2026-07-11

## Question

Can the current COMPASS water-only Richards path reproduce the infiltration
front and storage response of the sourced SV62 layered-soil experiment under
the same material parameters and static boundary conditions as ClimaLand and a
fresh HYDRUS-1D run?

## Frozen physical setup

- vertical one-dimensional column, depth 110 cm;
- 16 stepwise hydraulic horizons from the ClimaLand artifact;
- van Genuchten-Mualem retention and conductivity, `l = 0.5`;
- wetting alpha equal to twice the artifact drying alpha;
- measured initial water-content profile;
- constant surface pressure head of +5 cm;
- free-drainage lower boundary;
- duration 60 min;
- requested comparison times 0, 8, 16, 24, 32, 40, and 60 min;
- no heat, vapor, solute, root uptake, or plant processes.

Each model keeps its native spatial discretization and nonlinear solver. All
profiles are interpolated only during analysis onto a common 0-100 cm,
1-cm comparison grid; no model output is used to initialize another model.

HYDRUS uses its standard modified-Picard controls (`TolTh = 1e-4`,
`TolH = 0.1 cm`, maximum 20 iterations) after a pre-output pilot with
`TolTh = 1e-6` and `TolH = 1e-3 cm` entered persistent minimum-step reduction.
The failed pilot produced no profiles and did not inform any physical setting
or comparison threshold. Its project uses seconds because that convention is
exercised by the installed official `UPINFIL` case; Ksat and print times are
converted exactly from cm/min and minutes.

COMPASS uses its current production head-increment tolerance of 50 cm. A
pre-output pilot at `1e-5 cm` rejected the first 30 s outer step after nine
binary retries down to 0.1171875 s; the last 80-iteration residual was
11.060253 cm and the accepted state correctly remained unchanged. The pilot
therefore diagnoses the already scheduled W3 convergence-contract gap and did
not inform any physical input or comparison threshold.

## Confirmatory checks

1. Input identity: horizon depths, theta_r, theta_s, wetting alpha, n, Ksat,
   top head, bottom policy, and duration must match after unit conversion.
2. Boundary identity: every accepted COMPASS trial must request fixed head at
   the top and free drainage at the bottom; HYDRUS and ClimaLand must be
   configured with the same policies.
3. Initial-profile control: each model's time-zero profile is plotted and its
   common-grid RMSE from the measured initial profile is reported.
4. Cross-model profile agreement: at every requested time, report pairwise
   RMSE and maximum absolute difference in volumetric water content.
5. Observation agreement: report RMSE and mean error against the measured
   profiles at the 11 measurement depths. These metrics are descriptive because
   the field experiment is not an exact manufactured solution.
6. COMPASS transaction check: every outer step must accept; rejected trials
   may be retried but cannot advance time or leak state.
7. COMPASS conservation check: report initial/final storage and the integrated
   actual top and bottom boundary fluxes using the solver's signed actual-flux
   records.

## Acceptance rule for this W2 case

The case passes boundary correctness only if checks 1, 2, and 6 are exact and
the COMPASS storage residual computed from actual boundary fluxes is no larger
than 0.05 mm over 60 min. Numerical parity is reported, not forced: COMPASS is
flagged for W3 investigation when its RMSE from the fresh HYDRUS run exceeds
the ClimaLand-to-HYDRUS RMSE by more than 0.01 m3/m3 at any requested time.

The historical HYDRUS CSV bundled with ClimaLand is shown separately. It is a
provenance control rather than an acceptance oracle because the file's surface
saturated water content does not equal the current tutorial parameter array.

## No post-hoc changes

The physical parameters, comparison depths/times, error metrics, and thresholds
above are frozen before examining new three-model outputs. Any later change is
recorded as a new protocol version rather than overwriting this one.
