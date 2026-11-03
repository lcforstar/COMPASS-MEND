# W2 Bonan clay three-model protocol

Date frozen: 2026-07-11

## Purpose

This case supplies the fresh COMPASS-ClimaLand-HYDRUS comparison that the
16-material SV62 case cannot provide with the locally installed HYDRUS 4.17
binary. It is the homogeneous clay case already used by ClimaLand as a
regression against Bonan's independent MATLAB result.

## Frozen setup

- 150 cm homogeneous vertical column;
- 150 one-centimetre elements in COMPASS and ClimaLand, 151 one-centimetre
  nodes in HYDRUS;
- theta_s=0.495, theta_r=0.124, alpha=0.026 1/cm, n=1.43, Mualem l=0.5;
- Ksat=0.0443 cm/hour;
- uniform initial theta=0.24;
- top theta=0.494, represented in head-based models by
  h=-1.4288210312335 cm from the same VG curve;
- free-drainage bottom;
- duration 1,000,000 s;
- profiles at 0 and every 100,000 s;
- water only: no heat, vapor, solute, roots, or plants.

ClimaLand uses its source regression controls: ARS111, 1000 s step, Newton
maximum 50, and maximum-error tolerance 1e-8. HYDRUS uses its standard
modified-Picard controls: initial 10 s step, 0.01-1000 s range, maximum 20,
TolTh=1e-4 and TolH=0.1 cm. COMPASS uses its current production head-increment
tolerance of 50 cm and transactional adaptive substeps.

## Frozen metrics and gates

1. Exact input and boundary-policy audit.
2. Every requested output time must exist and all values must be finite.
3. Report pairwise RMSE and maximum absolute theta difference on the common
   one-centimetre grid for all eleven times.
4. At final time report RMSE from the preserved Bonan MATLAB profile.
5. Report COMPASS actual-boundary storage residual; boundary correctness
   requires absolute residual <=0.05 mm.
6. Flag COMPASS for W3 when its HYDRUS RMSE exceeds the ClimaLand-HYDRUS RMSE
   by more than 0.01 m3/m3 at any requested time.

No physical parameter, metric, or threshold is changed after model outputs are
examined. A failed run remains part of the result rather than being hidden by
calibration.

## Recorded outcomes

The preregistered COMPASS production head-form run remains preserved as a
failure: outer step 1 was rejected after the minimum trial step reached
3.90625 s. Reducing the minimum to 0.01 s and disabling Anderson did not rescue
that formulation.

The diagnosis found an inherited negative signed VG derivative convention in
the endpoint refresh, while the Richards matrix requires positive water
capacity. The shared endpoint auxiliary, retention curve, conductivity curve,
and coupled golden behavior were left unchanged; only matrix assembly converts
the signed auxiliary with `C=abs(DTheta_LLh)`. This made the Picard map monotone. A separate
same-standard run then used mixed-form storage, no Anderson, `MaxIt=20`,
`TolTh=1e-4`, and `TolH=0.1 cm`. It completed all 1,000 steps without rejection
and closed cumulative water balance to 0.002107 mm.

The current ClimaLand source test against its preserved Bonan artifact was run
directly and failed at RMSE 0.0273825355, exactly matching the cached fresh
ClimaLand runner. The Bonan artifact is therefore reported but is not used as a
golden acceptance gate for this source snapshot.
