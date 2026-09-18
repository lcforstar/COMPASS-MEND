# Heat boundary contract

This package separates a heat boundary into four objects. The separation is
important because a requested condition is not always the equation ultimately
imposed, and an imposed temperature is not an actual energy flux.

1. `normalizeSpec` converts site configuration into an immutable `Spec`.
2. `materializeRequest` freezes the forcing over one exact trial interval and
   produces a `Request`.
3. `solveThermalTrial` converts the request into the matrix row actually used
   by the finite-element solve and returns `Applied`.
4. `buildResult` combines the request, applied row kind, and accepted boundary
   reaction into `Actual`.

The public `[boundary]` keys are `heat_top_kind`, `heat_bottom_kind`,
`heat_top_value`, and `heat_bottom_value`. Legacy `NBCT` integers are generated
only by `init.buildBoundary` for the retained matrix adapter.

The public flux convention is positive into the soil column. Therefore the
legacy reactions are converted as

```text
q_U,top^in    = -QET
q_U,bottom^in =  QEB
```

Both fluxes are total energy fluxes in `J cm^-2 s^-1`. They are not conductive
temperature fluxes alone. For a prescribed temperature, the solver first
imposes the Dirichlet correction row

```text
-a_h delta_h + a_U delta_U - a_P delta_Pg = T_bc - T_iter
```

and then recovers the total energy reaction from the complete coupled row. For
a prescribed energy flux, the frozen request is inserted directly into the
Neumann row. A rejected trial never enters the accepted boundary accumulator.

Supported policies are:

| Side | Request kinds | Applied equation |
|---|---|---|
| top | `atmospheric`, `surface_energy_response`, `surface_temperature`, `specified_heat_flux` | natural matched flux, affine response, fixed temperature, or fixed total-energy flux |
| bottom | `seasonal_bottom_temperature`, `fixed_temperature`, `specified_heat_flux`, `free_drainage_internal_energy` | fixed temperature, fixed total-energy flux, or liquid internal-energy export paired with bottom free drainage |

`surface_energy_response` is an experimental atmospheric-surface contract for
bare, partial-snow, and full-snow states. It applies

```text
q_top^in(T0) = q_ref + dq_dT (T0 - T_ref),  dq_dT <= 0
```

inside every complete `h/U/P_g` Picard map. Bare radiation, sensible heat and
latent demand are linearized at the SCOPE surface solution; covered ground
adds snow-to-soil conduction from prognostic snow enthalpy. The dry-water
active set replaces potential latent cooling with accepted evaporation before
the same heat row is assembled. Every boundary switch clears acceleration
history and rebuilds water, heat, and dry-air rows at the same trial time.

SCOPE's empirical force-restore `Gtot` is a predictor/diagnostic only. It is
not the flux prescribed to the soil. The accepted boundary summary stores the
time-weighted soil-surface temperature used to evaluate the response because
adaptive substeps make the final endpoint temperature insufficient for
reconstructing the outer-step mean flux.

`surface_temperature` is an explicit diagnostic/legacy operator-split contract
and remains a frozen Dirichlet condition, even when an affine response exists.
`atmospheric` is the production natural-surface policy. The outer surface-soil
coupler matches the surface tangent to the upper-20-cm dynamic heat response,
supplies its total-energy flux, and materialization converts
`W m^-2 -> J cm^-2 s^-1` without imposing the diagnostic skin temperature.

`free_drainage_internal_energy` is valid only with
`water_bottom_kind = "free_drainage"`. In each coupled Picard map it applies

```text
q_U,bottom^in = q_w,bottom^in rho_l c_l (T_bottom - 0 degC)
```

as a Robin row, where both inward fluxes are negative during drainage. The
water flux and bottom temperature come from the same `h/U/P_g` map and are
refreshed on the next Picard iteration. Soil wetting energy is not exported as
liquid enthalpy. This follows ClimaLand `EnergyWaterFreeDrainage` in
`D:\OneDrive\应用\Model\ClimaLand\src\standalone\Soil\boundary_conditions.jl`,
lines 585-615, translated to COMPASS's 0 degC liquid reference and CGS units.
