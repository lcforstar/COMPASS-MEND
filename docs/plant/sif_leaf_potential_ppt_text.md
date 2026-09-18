# SIF-Based Leaf Water Potential Prediction

## Physically Interpretable Formula

The normalized fluorescence-yield model was written as:

```text
Ψleaf = Ψref
      + KFY (SIFyield - SIFyield,ref)
      - KD  (VPD - VPDref)
      - KT  (Tcave - Tcave,ref)
```

where:

```text
SIFyield = SIF / aPAR
```

## Reference State

The reference state was directly calculated from the annual valid daytime samples, not fitted:

```text
Ψref          = mean(Ψleaf)
SIFyield,ref  = mean(SIF / aPAR)
VPDref        = mean(VPD)
Tcave,ref     = mean(Tcave)
```

For this simulation:

```text
Ψref         = -0.9149 MPa
SIFyield,ref = 0.001729
VPDref       = 0.8576 kPa
Tcave,ref    = 20.20 °C
```

## Coefficient Calculation

```text
KFY = 630.42 MPa
KD  = 0.1661 MPa kPa^-1
KT  = 0.06154 MPa °C^-1
```

These coefficients were not assigned arbitrarily. They can be estimated from the annual valid simulation samples as sensitivity coefficients around the reference state:

```text
KFY =  ∂Ψleaf / ∂SIFyield
KD  = -∂Ψleaf / ∂VPD
KT  = -∂Ψleaf / ∂Tcave
```

In practice, the coefficients were calculated by fitting the centered sensitivity equation:

```text
Ψleaf - Ψref
    = KFY (SIFyield - SIFyield,ref)
    - KD  (VPD - VPDref)
    - KT  (Tcave - Tcave,ref)
```

This form is physically more meaningful than fitting a raw intercept, because all predictors are expressed as deviations from the annual reference state.

After calculating the reference state and sensitivity coefficients, the prediction equation becomes:

```text
Ψleaf = -0.9149
      + 630.42 (SIFyield - 0.001729)
      - 0.1661 (VPD - 0.8576)
      - 0.06154 (Tcave - 20.20)
```

## Units

```text
Ψleaf    : MPa
SIF      : W m^-2
aPAR     : W m^-2
SIFyield : dimensionless
VPD      : kPa
Tcave    : °C
KFY      : MPa
KD       : MPa kPa^-1
KT       : MPa °C^-1
```

## Physical Meaning

`Ψref` is the reference leaf water potential under the average annual daytime state. It is not an artificial intercept at zero radiation, zero VPD, and zero temperature.

`KFY` is the fluorescence-yield hydraulic sensitivity coefficient. Because `SIFyield = SIF/aPAR` is dimensionless, the unit of `KFY` is MPa. A positive `KFY` means that higher normalized fluorescence yield corresponds to less negative leaf water potential, indicating better plant water status.

`KD` is the atmospheric drought drawdown coefficient. It quantifies how strongly atmospheric water demand lowers leaf water potential. Here, when VPD increases by 1 kPa relative to the reference state, `Ψleaf` decreases by about 0.166 MPa.

`KT` is the canopy thermal stress coefficient. It quantifies the decrease in leaf water potential associated with canopy warming. Here, when canopy temperature increases by 1 °C relative to the reference state, `Ψleaf` decreases by about 0.0615 MPa.

## Prediction Validation

After substituting the calculated coefficients back into the equation, the predicted leaf water potential was compared with the simulated annual `Ψleaf`:

```text
Model: Ψleaf = Ψref + KFY(SIFyield - SIFyield,ref)
              - KD(VPD - VPDref)
              - KT(Tcave - Tcave,ref)

R²              = 0.687
RMSE            = 0.431 MPa
Blocked CV R²   = 0.653
Blocked CV RMSE = 0.454 MPa
```

This indicates that the coefficients are calculable from the model outputs, and after being inserted into the centered physical equation, the formula can reproduce the annual variation of leaf water potential with moderate skill.

Therefore, the result supports the following interpretation:

```text
SIFyield, VPD, and Tcave can be used to construct a physically interpretable
remote-sensing proxy for Ψleaf. The calculated sensitivity coefficients have
clear units and signs, and the resulting equation gives a reasonable prediction
of simulated leaf water potential.
```

For comparison, a less directly interpretable but more predictive remote-practical model using raw SIF and absorbed radiation was:

```text
Ψleaf = a + b SIF + c aPAR + d VPD + e Tcave
```

with:

```text
R²              = 0.790
RMSE            = 0.353 MPa
Blocked CV R²   = 0.763
Blocked CV RMSE = 0.375 MPa
```

Therefore, the SIF-yield model is preferred for physical interpretation, while the raw-SIF + aPAR model is preferred for prediction accuracy.

## Suggested PPT Conclusion

```text
The reference state and sensitivity coefficients were calculated from annual valid simulation samples. After substituting these coefficients into the centered equation, SIFyield, VPD, and canopy temperature predicted Ψleaf with R² = 0.687 and blocked cross-validation RMSE = 0.454 MPa.

This suggests that normalized SIF can provide a physically interpretable remote-sensing proxy for leaf water potential when combined with atmospheric drought demand and canopy thermal stress.
```
