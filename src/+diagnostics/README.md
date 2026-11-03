# +diagnostics

COMPASS diagnostic output system (ClimaLand-style registry, CF-Convention NetCDF).
Replaces the legacy "103 places scattered Sim_output.X(...) = ..." anti-pattern.

## Why

Old layout (规范 §7.2 触雷反例):

```matlab
% Scattered across src/COMPASS.m 103 lines:
Sim_output.SoilMoisture(out_idx, :) = SoilVariables.Theta_LL(...);
Sim_output.NitrificationFlux(out_idx, 1) = NitrogenFlux.cvSnit(KT);
Sim_output.fwet_sun(out_idx) = ...;
% Then writeCompassOutputs.m loops over 99 fields and dumps each as CSV.
```

Problems:
- No central place to see "what variables exist + what units"
- New field needs editing 3 files
- Units / metadata lives in side `.md` files that drift from code
- No NetCDF, no CF Convention → not paper / data-archive ready

## New design

3-phase lifecycle:

```matlab
% (1) REGISTER (before main loop)
reg = diagnostics.Registry(tS, NL, NN);
diagnostics.registerAllDiagnostics(reg, cfg);

% (2) WRITE (inside main loop, replaces every Sim_output.X = ...)
reg.write('soil_moisture', SoilVariables.Theta_LL, KT);
reg.write('nitrification_flux', NitrogenFlux.cvSnit(KT), KT);

% (3) FLUSH (after main loop)
diagnostics.flushToNetCDF(reg, output_dir, run_metadata);
diagnostics.flushToCsv(reg, output_dir);   % compatibility for legacy plots
```

## Output layout

```
data/output_<tag>/
  metadata.yaml                  cross-group summary (human-readable)
  soil_water.nc                  group 1
  soil_thermal.nc                group 2
  soil_freezing.nc               group 3 (cfg.lFreezing=1)
  soil_air.nc                    group 4 (cfg.lSoilairefc=1)
  soil_nitrogen.nc               group 5 (cfg.lNefc=1)
  plant_canopy.nc                group 6
  plant_hydraulics.nc            group 7
  surface_water.nc               group 8
  evapotranspiration.nc          group 9
  root_uptake.nc                 group 10
  dual_porosity.nc               group 11 (cfg.lDualPor=1)
  conservation.nc                group 12 (cross-cutting residuals)
  numerics.nc                    group 13
  radiation.nc                   (stub, future)
  solver_internal.nc             (stub, future)
  Sim_output_*.csv               legacy compatibility (CSV mirror, kill in v34)
```

## Per-variable metadata (CF-1.10)

Each `register(...)` call sets:

| field | example |
|---|---|
| `name` | `soil_moisture` (snake_case, NetCDF var name) |
| `long_name` | `Volumetric soil moisture (liquid phase)` |
| `standard_name` | `mass_content_of_water_in_soil_layer` (CF table) |
| `units` | `m3 m-3` (UDUNITS-2) |
| `dims` | `{'time', 'node'}` |
| `group` | `soil_water` → goes into `soil_water.nc` |
| `source` | `SoilVariables.Theta_LL` (in-memory origin) |
| `source_code` | `src/+soil/+water/calculateMatricCoefficients.m` |
| `transform` | `edge_to_node` / `none` / `flip_node` / `scalar` |
| `comment` | free-text caveat |
| `fill_value` | `-9999` |

## Adding a new variable

1. Choose a group (or add a new `register*.m` for a new module)
2. Add one `reg.register(struct(...))` block
3. In the main loop, one line: `reg.write('your_name', value, KT)`
4. (Optional) add legacy CSV mapping in `flushToCsv.m::legacyCaseFor`

No edits to 3 files. No metadata drift.

## Python read-side

```python
import xarray as xr
ds = xr.open_dataset('data/output_xxx/soil_water.nc')
ds['soil_moisture'].sel(node=10).plot()
ds['soil_moisture'].plot.contourf(x='time', y='depth_node')
```

## See also

- `docs/netcdf_output_design_decisions.md` — 8 design decisions
- `docs/output_catalog.md` — full variable catalog (legacy + new)
- 规范 §7 — diagnostic container conventions
- 规范 §12 — refactor priority #2
