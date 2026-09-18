# W2 van Genuchten Water-Capacity Sign Audit

Date: 2026-07-11  
Status: superseded implementation boundary; derivation remains valid

## Decision

For the pressure-head convention used by STEMMUS and COMPASS, unsaturated
pressure head is negative and becomes less negative as soil wets. Therefore
the Richards storage capacity

```text
C(h) = d theta / d h
```

must be non-negative on the monotone van Genuchten retention curve. The
negative analytical value returned by the public STEMMUS implementation is
not a separate physical convention supported by its governing equations. It
is an inherited analytical-derivative sign defect.

The final W2 boundary preserves the legacy `calculateMatricCoefficients*`
functions and their inherited `DTheta_LLh` convention unchanged. The new
transactional water chain now uses named positive phase derivatives and
`calculateWaterEquationCoefficients`; it does not convert the legacy field
with `abs` and does not overwrite it with a time-layer secant.

## Governing-equation check

The STEMMUS technical report defines, for `h < 0`,

```text
theta(h) = theta_r + (theta_s - theta_r)
           * (1 + |alpha h|^n)^(-m).
```

Differentiating with the chain rule gives

```text
d theta / d h = -(theta_s - theta_r) * m * n * alpha
                * |alpha h|^(n-1)
                * (1 + |alpha h|^n)^(-m-1)
                * sign(h).
```

Because `h < 0`, `sign(h) = -1`, so `d theta / d h > 0`.

The public formula is exactly the derivative with respect to the positive
suction magnitude `psi = |h|`, for which `d theta / d psi < 0`. It is then
stored and used as though it were the derivative with respect to the signed
pressure head `h`. The missing conversion is

```text
d theta / d h = (d theta / d psi) * (d |h| / d h)
              = (d theta / d psi) * sign(h).
```

For unsaturated `h < 0`, this missing factor reverses the sign.

## Public-source trace

The current official source and the local public checkout agree:

- `calculateDTheta.m`, cases 1 and 5, return the negative analytical value.
- `calculateMatricCoefficients.m` inserts `DTheta_LLh` into `Chh` without a
  compensating sign change.
- `assembleCoefficientMatrices.m` adds `Chh` directly to `C1`.
- `calculateTimeDerivatives.m` uses `C1 / dt` on the left-hand matrix and the
  old-head term on the right-hand side, again without a compensating sign.

The comment referring to STEMMUS_SCOPE issue 181 does not justify the sign.
Item 7 of that issue asks whether `theta_s` or `theta_m` should be used. The
maintainer response concerns extremely dry conditions and numerical
convergence; it does not discuss the missing `sign(h)` factor.

The same formula already appears in the original `yijianzeng/STEMMUS`
repository. Its only `CondL_h.m` history entry is the initial 2020 public
upload with commit message `SourceCode`; no public commit explanation for the
sign was found. That original file explicitly describes `h` as suction head
with a negative sign, and its initialization creates negative `h` values from
water content. Thus the suction-magnitude interpretation above is a
mathematical reconstruction from the source, not a claim about undocumented
author intent.

Official references:

- Technical report: <https://ris.utwente.nl/ws/files/5182404/STEMMUS_Report.pdf>
- Original derivation: Yijian Zeng (2013), *Coupled Dynamics in Soil*,
  Sect. 4.2-4.4, DOI `10.1007/978-3-642-34073-4`; local copy:
  `D:/OneDrive/文档/最近读文献/STEMMUS/试验观测与案例/Coupled Dynamics in Soil.pdf`.
- Current source: <https://github.com/EcoExtreML/STEMMUS_SCOPE/blob/main/src/%2Bconductivity/%2BhydraulicConductivity/calculateDTheta.m>
- Issue 181, item 7: <https://github.com/EcoExtreML/STEMMUS_SCOPE/issues/181#issuecomment-1658570621>
- Original implementation: <https://github.com/yijianzeng/STEMMUS/blob/master/CondL_h.m>

The book removes the remaining sign-convention ambiguity:

- Eq. 4.4 defines the modeled head from pore-water minus pore-air pressure,
  divided by the specific weight of water; it is negative in unsaturated soil.
- Eq. 4.17 applies the van Genuchten curve on the explicit `h < 0` branch.
- Eq. 4.57 expands storage using `partial theta_L / partial h`.
- Eqs. 4.48-4.53 place that storage coefficient on the same positive mass-
  matrix side used by the MATLAB `h_MAT` and `h_EQ` implementation.

Therefore, the negative analytical source value cannot be reconciled by
claiming that the book solves for a positive suction variable.

## Why the public model can still run

The negative analytical branch is not the only derivative used by the public
solver. For the common VG path, its control flow is:

1. When the new and previous heads are almost equal, use the negative
   analytical tangent.
2. Once a Picard iterate changes head sufficiently, use a secant derivative,
   such as `(theta_new - theta_old) / (h_new - h_old)`.
3. The monotone retention curve makes that secant derivative positive.
4. With thermal head correction enabled, `calculateMatricCoefficients.m` can
   overwrite the analytical derivative with another positive secant estimate.

Consequently, a normal field run may receive a negative initial tangent and
then recover through positive secants, small forcing increments, damping, or
time-step reduction. This explains successful public simulations; it does not
make the analytical negative sign correct. A strong laboratory wetting front
is a discriminating case because the wrong first linearization can destabilize
the iteration before the secant path rescues it.

## COMPASS migration finding

Before the vectorized endpoint-kernel migration, COMPASS contained the full
chain-rule expression with `sign(h)`, which was positive for unsaturated
pressure head. Commit `53a6a191` replaced it with the public negative formula
for source parity. The later endpoint kernel then evaluated that negative
analytical derivative on every refresh and no longer reproduced the public
solver's tangent-to-secant switching. This made the latent public defect more
persistent in the new water-only chain.

## Numerical evidence

At the Bonan clay initial head `h = -566.010378836 cm`:

```text
centered finite-difference d theta / d h  > 0
independent analytical VG derivative     = +8.62807340181e-05 1/cm
inherited public DTheta_LLh              = -8.62807340181e-05 1/cm
```

With the negative capacity inserted into the Richards matrix, the unaccelerated
Picard iteration oscillates and diverges. With the positive capacity, the same
physics, grid, boundary conditions, and time step converge monotonically. The
completed 1000-step Bonan run then matches HYDRUS at final-time profile RMSE
`0.0014423935` in volumetric water content, with water-balance residual
`0.0021066 mm`.

## Superseded compatibility boundary

The temporary implementation was:

```text
shared DTheta_LLh auxiliary: retain inherited signed value
Richards Chh matrix input:    use abs(DTheta_LLh)
```

That temporary `abs` conversion has now been removed from the preserved legacy
assemblers. The implemented boundary is:

```text
legacy calculateMatricCoefficients*: unchanged signed behavior
new phase auxiliary interface:       named positive physical partials
new Richards Chh/ChT:                 derivative of conserved water mass
new mixed storage:                    liquid + ice equivalent + vapor
```

The complete derivation and freeze-thaw extension are in
`docs/water_heat_freeze_conservative_derivation.md`.

## Remaining work

Before full coupled capability is declared, audit every non-Richards consumer
of `DTheta_LLh` and `DTheta_UUh` independently. In particular:

- heat-capacity and latent-heat coefficients;
- vapor and dry-air storage coefficients;
- freezing/SFCC analytical derivatives;
- thermal-head correction and its secant overwrite;
- previous-step versus current-iterate derivative ownership.

No global sign replacement is authorized by this W2 result.
