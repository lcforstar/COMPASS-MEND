# Coupled water-vapor-air-heat equation audit

Date: 2026-07-11  
Status: derivation and code audit; no coupled-physics change is authorized by
  this document; implementation approval granted 2026-07-11 for the separate
  conservative coefficient path

## 1. Decision and scope

The coupled water-vapor-air-heat implementation must not be selected by
matching a legacy trajectory. The target equations are first derived from the
conserved water mass, dry-air mass, and total energy. Code coefficients are
then classified against those derivatives. Only after discriminating tests
identify the failing term may the final kernel be changed.

The accepted current-step root-water wiring is outside this audit and must not
  be redesigned. The temporary `abs(DTheta_LLh)` conversion has been removed
  from the preserved legacy assemblers. The new transactional path uses named
  phase partials and a separate conservative assembler.

Primary sources:

- Zeng, Y. (2013), *Coupled Dynamics in Soil*, Sect. 4.2-4.4, especially
  Eqs. 4.4, 4.10, 4.13, 4.15, 4.17, 4.57-4.65,
  DOI: 10.1007/978-3-642-34073-4.
- Yu, L., Zeng, Y., Wen, J., and Su, Z. (2018), "Liquid-Vapor-Air Flow in
  the Frozen Soil," Eqs. 1-3 and Appendix Eqs. A1-A5 and A21-A23,
  DOI: 10.1029/2018JD028502.

## 2. State convention

The independent variables are signed liquid pressure head `h`, soil
temperature `T`, and mixed pore-gas pressure `Pg`.

- `h < 0` in unsaturated soil, as defined by Zeng Eq. 4.4.
- `theta_l` is liquid volume per bulk volume.
- `theta_i` is ice volume per bulk volume.
- `theta_g = porosity - theta_l - theta_i` is gas-filled pore volume.
- `theta_u = theta_l + (rho_i/rho_l) theta_i` is total water expressed as
  liquid-water-equivalent volume.
- Unfrozen `rho_v` follows the clipped Kelvin liquid-water branch; when
  `hh_frez < 0`, `rho_v` follows the IAPWS ice-Ih sublimation curve.
- `rho_a = Pg/(Ra Tk) - rho_v Rv/Ra` follows Dalton's and ideal-gas laws.

The names `Theta_V` and `Theta_g` in the current code both denote gas-filled
pore volume, not a separate volumetric amount of pure vapor.

## 3. Conserved densities and exact storage derivatives

### 3.1 Water mass

Normalize total water mass by liquid density:

```text
Mw = theta_u + (rho_v/rho_l) theta_g
```

For any primary variable `x` in `{h,T}`:

```text
theta_g,x = -(rho_l/rho_i) theta_u,x
            + (rho_l/rho_i - 1) theta_l,x

Mw,x = theta_u,x
       + (rho_v/rho_l) theta_g,x
       + (theta_g/rho_l) rho_v,x
```

This reduces in unfrozen soil to the published coefficient

```text
Mw,x = (1-rho_v/rho_l) theta_l,x
       + (theta_g/rho_l) rho_v,x .
```

The frozen form must not be inferred from a field name. It requires both the
total-water and liquid-water derivatives because ice and liquid occupy
different volumes.

### 3.2 Dry-air mass

The conserved dry-air density, including Henry-law dissolved air, is

```text
Ma = rho_a (theta_g + Hc theta_l) .
```

For `x` in `{h,T}`:

```text
Ma,x = rho_a,x (theta_g + Hc theta_l)
       + rho_a [-(rho_l/rho_i) theta_u,x
                + (Hc + rho_l/rho_i - 1) theta_l,x] .
```

For gas pressure, phase fractions are held fixed:

```text
Ma,Pg = rho_a,Pg (theta_g + Hc theta_l) .
```

The familiar unfrozen term `rho_a (Hc-1) theta_l,x` is therefore not a valid
three-phase derivative.

### 3.3 Total energy

With reference temperature zero degrees Celsius, define phase energy
densities per phase volume:

```text
e_l = rho_l (c_l T - W)
e_v = rho_v (c_v T + L0) = rho_v (c_l T + L(T))
e_a = rho_a c_a T
e_i = rho_i (c_i T - Lf) .
```

The conserved energy is

```text
E = E_solid + e_l theta_l + e_v theta_g
    + e_a theta_g + e_i theta_i .
```

For `x=h`, the exact storage coefficient is

```text
E,h = (e_l-e_v-e_a) theta_l,h
      + (e_i-e_v-e_a) theta_i,h
      + theta_g [(c_v T+L0) rho_v,h + c_a T rho_a,h] .
```

For `x=T`, add the bulk sensible heat capacity to the same chain-rule form:

```text
E,T = C_bulk
      + (e_l-e_v-e_a) theta_l,T
      + (e_i-e_v-e_a) theta_i,T
      + theta_g [(c_v T+L0) rho_v,T + c_a T rho_a,T] .
```

For `x=Pg`:

```text
E,Pg = theta_g c_a T rho_a,Pg .
```

These expressions provide one auditable source for `CTh`, `CTT`, and `CTa`.
Separate hand-written frozen and unfrozen energy formulas are unnecessary once
all phase derivatives are explicit.

## 4. SFCC derivative identities

Yu et al. Eqs. A1-A2 specify

```text
theta_u = F(h)
theta_l = F(h + h_frz) .
```

Therefore, within a smooth active branch:

```text
theta_u,h = F'(h)
theta_l,h = F'(h+h_frz)
theta_l,T = F'(h+h_frz) h_frz,T .
```

The printed Eq. A3 and its stated Heaviside convention are not physically
self-consistent: they activate on the warm side of `Tcrit`. A continuous
Clapeyron closure is obtained from

```text
h_eq(T) = Lf (T-T0)/(g T0)
Tcrit(h) = T0 + g h T0/Lf
h_frz = H(Tcrit-T) [h_eq(T)-h] .
```

Then `h+h_frz` is continuous at `Tcrit` and equals the equilibrium freezing
head below the transition. Unit conversion to centimetres contributes the
existing factor `1e4` when `Lf` is in J/kg and `g` is in cm/s2.

Before 2026-07-12, `computeFreezingPoint` used this continuous difference form
but divided the A4 depression by an additional factor of two. The current
implementation removes that factor and satisfies exact Clapeyron closure on
the smooth, unclamped frozen branch. The analytical and finite-difference
regression is `tests/test_computeFreezingPoint_oracle.m`.

## 5. Current-code mapping

| Mathematical quantity | Current producer or consumer | Audit result |
|---|---|---|
| `theta_u=F(h)` | `calculateTheta_UU.m`, `computeEndpointHydraulics.m` | Formula path agrees with A1 in the strict SFCC branch. |
| `theta_l=F(h+h_frz)` | `calculateTheta_LL.m`, `computeEndpointHydraulics.m` | Formula path agrees with A2 in the strict SFCC branch. |
| `theta_u,h=F'(h)` | Produced under the name `DTheta_LLh` in strict SFCC | Confirmed semantic cross-over. The derivative routine also returns the negative derivative with respect to `abs(h)`, not the positive derivative with respect to signed `h`. |
| `theta_l,h=F'(h+h_frz)` | Produced under the name `DTheta_UUh` in strict SFCC | Confirmed semantic cross-over and the same sign error. |
| `rho_v,h`, `rho_v,T` | `soil.gas.Density_V` | Resolved 2026-07-13: derivatives now follow the active clipped-Kelvin or IAPWS ice-equilibrium branch; the old -20 degC density plateau was removed. |
| `rho_a` derivatives | `soil.gas.Density_DA` | Algebra agrees with A22-A23 if the vapor derivatives supplied to it are correct. |
| Unfrozen `Chh`, `ChT` | `soil.water.calculateMatricCoefficients` | Algebra agrees with Zeng 4.57 after correcting the signed-head derivative. The local `abs` is only a compatibility conversion. |
| Frozen water storage | Same function | Dominant `theta_u,h` is obtained accidentally through the crossed field meaning. The gas-volume correction still assumes liquid and ice have the same volume, and the public names are unsafe. |
| Frozen dry-air storage | `soil.gas.calculateDryAirParameters` | Confirmed incomplete three-phase derivative. It uses actual three-phase pore volume but differentiates it with the unfrozen `(Hc-1) DTheta_LLh` expression. |
| Unfrozen energy storage | `soil.thermal.calculateEnergyParameters` | Algebra agrees with Zeng 4.64. `c_l*T+L(T)` is exactly `c_v*T+L0`. The public upstream line-continuation loss has already been repaired in COMPASS. |
| Frozen energy storage | Strict and sinusoidal overrides in `calculateEnergyParameters` | Not derivable from one conserved energy density as written. It depends on crossed derivative names, drops all phase-composition terms from `CTh`, and uses separate density factors in the latent branch. Classified as a high-risk unresolved defect pending coefficient-level tests. |
| Time-old water state | `solveWaterTrial` maps accepted `Y.h/Y.T` to legacy `h/T`; trial values map to `hh/TT` | Correct for a substep transaction. `workspace.step.previous` is the outer A3 extrapolation anchor and must not replace the accepted substep-start state globally. |
| Time-old gas and heat state | Legacy gas/thermal solvers and B3b adapter | The algebra uses accepted substep-start `h,T,Pg` against the current trial. Current `runOneStep` does not yet execute a coupled gas/thermal solve; B2/B3 are coefficient/read-back adapters. |

## 6. Severity-ranked findings

### F1. Signed hydraulic derivative is wrong at the producer

Severity: critical for nonlinear Richards behavior.  
Confidence: confirmed analytically.

The VG functions use `abs(h)`, but `calculateDTheta` omits the derivative of
`abs(h)` for `h<0`. The returned tangent is negative although
`d theta/d h` must be positive. A local matrix `abs` does not repair dry-air,
energy, diagnostics, secants, or future consumers.

### F2. Strict-SFCC derivative identities are crossed

Severity: critical architectural ambiguity.  
Confidence: confirmed from A1-A2 and both scalar/vectorized producers.

`DTheta_LLh` evaluates `F'(h)` and `DTheta_UUh` evaluates
`F'(h+h_frz)`. Some legacy consumers appear to rely on the actual crossed
meaning. A global search-and-replace would therefore introduce new errors.

### F3. Frozen dry-air storage is not the derivative of dry-air mass

Severity: high when air flow and freezing are both active.  
Confidence: confirmed by chain rule.

The current storage volume includes ice, but `Cah` and `CaT` omit the ice-volume
derivative and do not distinguish total-water from liquid-water derivatives.

### F4. Frozen energy coefficients are not closed from total energy

Severity: potentially critical for freeze-front temperature and conservation.  
Confidence: high for structural inconsistency; term-by-term numerical impact
still requires tests.

The strict branch discards liquid/ice composition contributions to `CTh` and
reconstructs `CTT` with legacy special terms. The result cannot be shown equal
to the exact derivatives in Sect. 3.3 using the actual definitions of
`theta_u`, `theta_l`, and `theta_i`.

### F5. Resolved: vapor-density branches and derivatives described different functions

Original severity: high near saturation, at the dry cutoff, and below -20 degC.  
Status: resolved in `soil.gas.Density_V` on 2026-07-13.

The lower and upper relative-humidity plateaus now have zero Kelvin activity,
the unfrozen temperature derivative differentiates the returned liquid-water
curve, and frozen nodes use the IAPWS ice-Ih sublimation curve with its analytic
temperature derivative and zero local matric-head derivative. Temperatures
outside [-100, 150] degC are rejected instead of silently projected.

### F6. Resolved: the freezing-point depression contained an unexplained factor 1/2

Original severity: high for transition timing and frozen-zone extent.  
Status: resolved in `computeFreezingPoint` on 2026-07-12; full coupled-run impact remains to be quantified.

The continuous difference form of `h_frz` is preferable to the internally
inconsistent printed A3, but `Tcrit` must still satisfy the Clapeyron equality.
The removed factor two prevented that equality. The corrected smooth branch now
has `dh_frz/dh=-1` and `d(h+h_frz)/dh=0`; numerical clamp branches retain zero
derivatives by convention.

### F7. Current coupled wiring is incomplete in the new transactional loop

Severity: blocking for declaring a final coupled model, not a defect in the
water-only W2 kernel.  
Confidence: confirmed by call-graph inspection.

`runOneStep` executes water-only attempts. B2/B3 compute and audit legacy
thermal/gas coefficients, but the new loop does not yet solve and transactionally
commit gas pressure and soil temperature. Coupled coefficient repair and I6
wiring must be separate increments.

## 7. Required final interface

The final auxiliary state must expose mathematical meanings directly:

```text
theta_liquid
theta_ice
theta_water_equiv
theta_gas

dtheta_liquid_dh
dtheta_liquid_dT
dtheta_water_equiv_dh
dtheta_water_equiv_dT
dtheta_ice_dh
dtheta_ice_dT
dtheta_gas_dh
dtheta_gas_dT

drho_v_dh
drho_v_dT
drho_a_dh
drho_a_dT
drho_a_dPg
```

Water, dry-air, and energy storage coefficients must be generated from the
three conserved-density functions in Sect. 3. Legacy names may exist only in a
temporary adapter and must assert the active branch they represent.

## 8. Approval boundary

This audit authorizes no change to the governing equations. Before editing the
coupled kernels, the following must be presented for approval:

1. coefficient-level tests that fail for exactly F2-F6;
2. the proposed conserved-density helper signatures and units;
3. a term-by-term old/new coefficient table for unfrozen and frozen states;
4. a staged compatibility plan that leaves W1 root uptake and W2 boundary
   contracts unchanged.
