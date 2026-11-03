# COMPASS partial-snow surface implementation and validation contract

Updated: 2026-07-19

This document defines the implemented CLM5/ClimaLand-inspired surface split
over one shared one-dimensional COMPASS soil column. It is an engineering
reduction, not a claim of numerical identity with CLM5 or ClimaLand.

## 1. Domain ownership

| Domain | Prognostic authority | Surface responsibility |
|---|---|---|
| atmosphere above canopy | hourly forcing | incoming radiation, precipitation, air state |
| canopy | SCOPE/plant/canopy water | radiation, photosynthesis, leaf exchange, interception |
| land surface | `+surface`, `+snow` | bare/snow flux split, snow enthalpy, ponding/runoff, boundary contracts |
| vadose soil | `+soil` | coupled `h/U/P_g`, phase change, root sinks |
| groundwater/deep boundary | typed bottom requests | drainage, pressure and thermal export/boundary |

`runOneStep` is the accepted-step orchestrator. Every surface and soil
candidate is constructed from the same accepted state and committed or rolled
back together.

## 2. Radiation and canopy contract

There is one canopy RTMo/RTMt solve. Effective lower-boundary reflectance is

```text
r_eff = (1-f_snow) r_bare + f_snow r_snow.
```

Thermal emission uses fourth-power averaging, not a linear temperature mean:

```text
eps_eff Trad^4 = (1-f) eps_bare T_bare^4 + f eps_snow T_snow^4.
```

Snow therefore feeds back to canopy radiation, PAR/aPAR, photosynthesis, SIF,
and longwave exchange through the common radiative solve. PAR/aPAR remains one
canopy field and is never multiplied by `(1-f)` or `f` again.

After the common solve, ground absorption is diagnosed separately for bare
and snow optical properties. Every public tile contribution is already on a
grid-area basis and must satisfy `X_ground=X_bare+X_snow`.

## 3. Surface temperatures and turbulent exchange

`T_bare` is the SCOPE bare-ground skin solution. `T_snow` is diagnosed from
prognostic snow internal energy. `T_rad` is only the effective radiation
interface temperature. These three quantities are never aliases.

Bare and snow pathways separately calculate:

```text
H_bare, LE_bare, qv_bare
H_snow, LE_snow, qv_snow.
```

Atmosphere-facing `H/LE/qv` are positive upward. Bare evaporation demand is
sent to the atmospheric soil-water active set; bare dew becomes liquid supply.
Snow sublimation/deposition changes only the snow mass ledger. Canopy wet
evaporation and transpiration remain separate canopy/plant pathways.

## 4. Soil-top total-energy response

The shared soil column never receives two tile temperatures. It receives one
inward-positive affine total-energy response:

```text
q_top(T0) = q_ref + dq_dT (T0-T_ref), dq_dT <= 0.
```

Bare external exchange contributes

```text
F_bare(T0) = Rn_bare(T0) - H_bare(T0) - LE_bare,demand(T0).
```

Snow contributes internal snow-to-soil conduction

```text
G_snow_to_soil(T0) = C_snow (T_snow-T0).
```

The grid response is the once-area-weighted sum plus any conservative
terminal snow-energy redistribution. SCOPE's force-restore `G` is used only in
the SCOPE skin-temperature solve and diagnostic comparison; it is absent from
the soil boundary response.

The response is inserted into the U-primary top row as a Robin condition. The
full `h/U/P_g` Picard system therefore sees the same surface-temperature
feedback without a separate surface-soil outer iteration.

## 5. Water active set and latent coupling

The atmospheric water boundary chooses among supply-limited, finite storage,
capacity/runoff, and residual-water dry states. It never switches solely on
`P>Ksat`; wet-side feasibility is diagnosed from the coupled soil response.

When the dry state limits actual evaporation below the SCOPE demand, the
unrealized latent cooling and its temperature tangent are removed from the
heat request inside the same Picard mapping. Liquid water, vapor and dry-air
equations remain active. Any branch change discards the rejected `h/U/P_g`
trial, clears acceleration history and reassembles all three equations at the
same substep time.

## 6. Snow transaction

Bulk snow owns total SWE `S`, liquid water `S_l`, and internal energy `U_s`.
Rain/snow partition uses air temperature and relative humidity. Liquid holding
is pore-limited; excess drains. Frozen depression water transfers to snow only
on outer-step acceptance.

The staged snow update uses the reference snow-to-soil flux. After the soil
hour accepts, the snow enthalpy is reconciled with the time-weighted accepted
soil-surface temperature. If this changes snow-domain existence, phase
resolution, or cover materially, the complete SCOPE/surface/soil candidate
graph is rebuilt once from the original accepted state. A second event request
rejects the outer step; no unlimited outer iteration is allowed.

Terminal snow removal is energy based. A solid trace exits only if its latent
energy is below the configured coupled energy resolution. Remaining liquid
water and enthalpy are delivered once; there is no simple SWE cutoff.

## 7. Code flow

```text
stageCanopyInputs
  -> stageCanopyEnergy
       -> snow surface predictor
       -> RTMo/RTMt
       -> ebal bare/snow radiation and turbulent responses
  -> stagePostEbal
       -> advanceCanopyWaterSnow
            -> partitionPrecipitation
            -> advanceBulkSnow
            -> computeGroundHeatResponse
            -> buildSoilTopHeatResponse
  -> materializeSoilBoundaryTrial
       -> water/heat/dry-air typed requests
  -> attemptSolveSoilEquations
       -> water active-set reaction
       -> realized-evaporation heat correction
       -> coupled h/U/P_g Picard map
  -> advanceSoilEquationOuterStep
       -> accepted substep means and conservation ledgers
  -> reconcileSnowGroundHeat
       -> optional one-time complete restage
  -> commitPending + commitOuterStepLedgers
```

## 8. Public output rules

- Names state direction, area basis, units, and whether a value is demand,
  accepted reaction, empirical diagnostic, or internal transfer.
- `soil_top_energy_response_evaluated_flux` is the accepted response;
  `soil_top_energy_flux_total_actual` is the U-primary reaction.
- `snow_to_soil_internal_energy_flux` is internal to the combined control
  volume and must cancel in the snow+soil total ledger.
- `scope_ground_heat_flux_empirical` is never presented as the actual soil
  flux.
- Ambiguous old public names such as potential/actual infiltration and generic
  `Gtot` are not part of the authoritative NetCDF contract.

## 9. Acceptance tests

The implementation is accepted only if all of the following remain green:

1. zero-snow `snow.enabled=true/false` complete-path parity;
2. exact bare+snow tile aggregation identities;
3. evaporation, dew, sublimation and deposition mass/sign identities;
4. cold full-snow end-to-end closure;
5. warm complete-meltout water and enthalpy single delivery;
6. failed soil solve rolls back snow, surface, SCOPE history and `h/U/T/P_g`;
7. wet and dry boundary active sets remain converged and transactional;
8. water, soil energy, dry-air and combined snow+soil ledgers meet tolerances;
9. no NaN/Inf and no nonphysical temperature/state excursions in short and
   annual runs;
10. Changwu completes `8760/8760`, followed by comparison with preserved
    pre-snow output, legacy output and measured soil moisture.

## 10. Explicit non-goals

- no second canopy biochemical/radiative column;
- no unbounded SCOPE-surface-soil outer Picard;
- no Nitsche/mixed finite-element rewrite;
- no rainfall-event breakpoint scheduler in this change;
- no watershed routing or MACRO/LAST preferential-flow implementation;
- no hydraulic calibration before the thermal/surface baseline is accepted.
