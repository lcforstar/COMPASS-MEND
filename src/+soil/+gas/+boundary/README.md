# Dry-air boundary contract

This package gives the dry-air equation the same four-layer boundary contract
as water and heat:

1. `normalizeSpec` creates an immutable `Spec` from configuration.
2. `materializeRequest` freezes pressure or mass-flux drivers for one trial.
3. `solveDryAirTrial` applies the corresponding pressure or mass-flux row and
   returns `Applied`.
4. `buildResult` records the request, applied row kind, and accepted `Actual`
   boundary reaction.

The public `[boundary]` keys are `dry_air_top_kind`, `dry_air_bottom_kind`,
`dry_air_top_value`, and `dry_air_bottom_value`. Legacy `NBCP` integers are
generated only by `init.buildBoundary` for the retained matrix adapter.

The conserved endpoint storage is

```text
B_a = rho_a (theta_g + H_c theta_l)                 [g cm^-3]
```

and the accepted column ledger is

```text
M_a = sum_e dz_e (B_a,e,1 + B_a,e,2) / 2            [g cm^-2]
R_a = M_a^(n+1) - M_a^n
      - dt (q_a,top^in + q_a,bottom^in)              [g cm^-2].
```

The public sign convention is positive into the soil column. Legacy reactions
are converted only at the package boundary:

```text
q_a,top^in    = -QAT
q_a,bottom^in =  QAB
```

For prescribed pressure, the actual mass flux must be recovered from the full
coupled boundary row after `h`, diagnostic `T`, and `P_g` have reached the
accepted state. A pressure-only reaction is incomplete because dry-air storage
and transport depend on all three coordinates. The production path therefore
uses `calculateCoupledBoundaryFlux`; `calculateBoundaryFlux` remains a legacy
single-coordinate compatibility helper.

Supported policies are:

| Side | Request kinds | Applied equation |
|---|---|---|
| top | `fixed_pressure`, `specified_mass_flux`, `atmospheric_pressure` | fixed pressure or fixed dry-air mass flux |
| bottom | `hydrostatic_pressure`, `specified_mass_flux`, `zero_mass_flux` | fixed pressure or fixed dry-air mass flux |

The pressure supplied by users and forcing is in Pa. Conversion to the internal
`g cm^-1 s^-2` coordinate is confined to the solver adapter. A rejected trial
does not update either the boundary transaction or the column mass ledger.
