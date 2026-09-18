# Water boundary input in `site.toml`

All user-facing water-boundary settings stay in the single `[boundary]` table.
Each endpoint has one `kind`; the kind selects the physical rule, data source,
CSV column, and interpolation policy.

## Default Changwu boundary

```toml
[boundary]
water_boundary_file = "water_boundary.csv"
water_top_kind       = "atmospheric"
water_bottom_kind    = "free_drainage"
```

This is physically identical to the historical `water_top=3` and
`water_bottom=3` setup. All prescribed values live in the CSV, not site.toml.

## Public modes

| Position | `kind` | Input | Interpolation |
|---|---|---|---|
| top | `atmospheric` | Throughfall and evaporation demand | natural surface switching |
| top/bottom | `constant_head` | CSV `<position>_value` [cm], unchanged | automatic `constant` |
| top/bottom | `variable_head` | CSV `<position>_value` [cm] | automatic `linear` |
| top/bottom | `constant_flux` | CSV `<position>_value` [cm/s], unchanged | automatic `constant` |
| top/bottom | `variable_flux` | CSV `<position>_value` [cm/s] | automatic `previous` |
| top/bottom | `zero_flux` | none | none |
| bottom | `free_drainage` | none | gravity drainage |

Users do not configure `source`, `series_column`, or `interpolation`.

## Select a field or laboratory site

The repository-root launcher reads the named folder under `data/sites` and
dispatches from its explicit `[simulation]` block:

```matlab
model = run_site("changwu");
model = run_site("lab_water");
model = run_site("lab_water", struct('stop_after_initialization', true));
```

`point_context="field"` uses the natural forcing/canopy/plant loop.
`point_context="laboratory"` uses the standalone Richards loop. The minimal
`data/sites/lab_water` example contains only `site.toml`, `soil.csv`,
`soil_initial_profile.csv`, and `water_boundary.csv`; it intentionally has no meteorological,
plant, root, vegetation, or SCOPE input. Do not infer the context by deleting
files: the explicit context and module switches are the only selectors.

## Prescribed boundaries and CSV

Constant and variable prescribed boundaries use the same CSV. For a constant
mode, repeat the same value for the entire corresponding column; validation
fails if that column changes. Variable modes may change by row.

The public flux sign is always positive into the soil domain:

- positive top flux means downward entry;
- negative top flux means upward removal;
- positive bottom flux means upward recharge;
- negative bottom flux means downward drainage.

Legacy solver signs remain confined to the internal adapter.

The historical top-boundary `BCh/QMT` convention was upward-positive. Do not
copy that sign directly into the current CSV. An old prescribed top evaporation
rate `+E` must be entered as `top_value=-E`; the solver adapter converts it back
to the legacy internal sign and converts the actual result to the public sign
again. The bottom legacy flux already follows the current into-domain-positive
direction.

For a controlled laboratory evaporation experiment, select a prescribed flux
and write a negative top value. For example:

```toml
[boundary]
water_top_kind = "constant_flux"
water_bottom_kind = "zero_flux"
```

```csv
time_s,top_value,bottom_value
0,-1.0e-6,0
86400,-1.0e-6,0
```

This requests upward removal at `1.0e-6 cm/s`. It is a prescribed experimental
flux, not an atmospheric evaporation demand. Natural evaporation remains a
nonnegative demand produced by the energy-balance chain and is converted to an
actual signed soil flux by the natural surface-boundary transaction.

`water_boundary_file` may be site-relative or absolute. The CSV has one common
time column and exactly two boundary-value columns:

```csv
time_s,top_value,bottom_value
0,0,0
3600,1.0e-6,-100
7200,0,-80
```

`time_s` is elapsed time from `[time].start` and must cover the full simulation.
The selected kind assigns units independently at each endpoint:

- `constant_head`: corresponding value is [cm] and every row must be equal;
- `variable_head`: corresponding value is [cm] and is linearly interpolated at
  the trial endpoint;
- `constant_flux`: corresponding value is [cm/s] and every row must be equal;
- `variable_flux`: corresponding value is [cm/s] and is piecewise constant on
  `[time_s(i), time_s(i+1))`; the scheduler integrates the exact trial interval.

The CSV keyword `previous` is internal policy metadata and is unrelated to the
solver state container `workspace.step.previous`.

## Historical migration

| Historical input | Current `[boundary]` mode |
|---|---|
| `water_top=1`, `h_top=x` | `water_top_kind="constant_head"`; repeat `x` in CSV `top_value` |
| `water_top=2`, constant flux `x` | `water_top_kind="constant_flux"`; repeat `x` in CSV `top_value` |
| `water_top=3` | `water_top_kind="atmospheric"` |
| `water_top=4` | `water_top_kind="zero_flux"` |
| bottom fixed head | `water_bottom_kind="constant_head"` |
| bottom fixed flux | `water_bottom_kind="constant_flux"` |
| `water_bottom=3` | `water_bottom_kind="free_drainage"` |
| `water_bottom=4` | `water_bottom_kind="zero_flux"` |

The removed split `source/series_column/interpolation` schema fails with an
explicit migration error. The loader converts the selected public mode into the
canonical internal policy and transitional Richards aliases.

`h_sur_cm` remains an observation. Irrigation remains a management input and
must not be duplicated as a prescribed top flux in the same scenario.
