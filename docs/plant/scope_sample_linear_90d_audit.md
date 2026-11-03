# 90-day sample-linear soil and SCOPE output audit

## Run metadata

- Site: Changwu
- Period: 90 days, 2160 hourly steps
- Soil input mapping: `sample_linear`
- Coupling: SCOPE plus water-only soil solve; soil heat, dry-air, nitrogen, and freezing solves disabled
- Standard NetCDF: `C:\Users\zjzx\AppData\Local\Temp\scope_soil_sample_linear_90d\nc`
- Captured SCOPE diagnostics: `docs/figures/scope_diag_sample_linear_90d.mat`

## Data quality

All 2160 soil steps completed and all captured SCOPE fields used in this audit were finite. The
standard `runOneStep` NetCDF path still does not write the SCOPE radiation, energy, carbon, canopy,
or hydraulic groups; these fields were therefore captured directly from `workspace.scope_step` by
`tools/captureScopeWindow.m`. The standard `soil_moisture` diagnostic is also still all missing, so
liquid water content was reconstructed from matric head and the initialized van Genuchten inputs.

## Soil-water result

- Reconstructed liquid water content: 0.1802 to 0.3091 m3 m-3.
- Written versus independently reconstructed liquid storage: maximum difference 0.0303 mm; this is
  the expected small vapor-equivalent contribution in the written storage definition.
- Column storage: 2420.654 mm initially and 2346.350 mm finally.
- Maximum absolute hourly water-balance residual: 9.55e-5 mm h-1.
- Ninety-day totals: throughfall 42.850 mm, soil evaporation 88.318 mm, transpiration 27.341 mm,
  and bottom water output 1.501 mm.
- Relative to the old horizon mapping, the interpolated profile removes the artificial parameter
  steps. Final storage is 2.911 mm lower; evaporation is 1.367 mm higher; transpiration changes by
  only -0.050 mm. The water solve remains stable and conservative.

## SCOPE result

Most ecosystem-scale flux ranges are physically plausible for this winter, low-LAI run:

- Net radiation: -90.39 to 491.62 W m-2.
- Latent heat: -8.61 to 290.84 W m-2.
- Sensible heat: -92.32 to 145.13 W m-2.
- Ground heat: -18.82 to 126.15 W m-2.
- Ninety-day GPP: 71.97 g C m-2 at constant LAI = 0.293.
- Plant hydraulic supply ratio remains 1; mean leaf water potential ranges from -1.29 to -0.157 in
  the model's hydraulic output units, so the failures are not caused by soil-water supply stress.

The full SCOPE result is nevertheless not acceptable:

- 11 hourly steps reached the 200-iteration energy-balance limit.
- The first failure is day 33.333 (KT=801); failures occur at local early-morning model hours 07-08.
- Maximum local sunlit-leaf energy residual is 1166.79 W m-2.
- Maximum ecosystem closure residual `Rn - LE - H - G` is 56.88 W m-2.
- A sunlit-leaf angular class exceeds 70 degC in 24 steps and reaches 99.998 degC; shaded-leaf and
  soil-skin maxima remain 20.33 and 31.52 degC, respectively.

## Located cause

The failure is independent of the soil interpolation mode. The solar zenith angle is staged in
`scope_step.angles.tts`, but `vegetation.energy.ebalClm45` checks `meteo.tts`. Because that field is absent,
the code records a 90-degree solar altitude for all 2160 hours and disables the intended low-sun
guard. At the failing dawn steps, effective sunlit LAI lies just above the separate 0.05 cutoff;
individual direct-beam leaf-angle classes then receive roughly 2500-4000 W m-2, producing the
near-100 degC leaf temperatures and non-convergence.

## Recommendation

Keep `sample_linear` for the soil inputs. Fix the SCOPE geometry connection by supplying the real
`angles.tts` to the low-sun check (or explicitly adding it to the ebal input contract), then rerun the
discriminating KT=801 case before another 90-day soak. Also restore `diagnostics.writeScopeStructs`
to the `runOneStep` accepted-step diagnostic path so future runs write the registered SCOPE groups.

## Fix verification

The solar-geometry wiring was corrected on 2026-07-12 in both the new-chain and legacy non-profile
and profile solvers. A complete 90-day regression produced zero failed energy-balance steps. Maximum
iteration count fell from 201 to 12 and maximum local residual from 1166.79 to 4.99928 W m-2. The
resolved solar altitude spans 0 to 58.95 degrees instead of remaining fixed at 90 degrees. Three
converged dawn steps at solar altitude 5.13-5.85 degrees still contain individual sunlit angular bins
at 76.64-80.82 degC; this is a separate threshold-policy issue at the current 5-degree low-sun cutoff,
not the former missing-angle wiring failure.
