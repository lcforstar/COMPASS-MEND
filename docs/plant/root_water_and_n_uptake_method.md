# Root Water And Nitrogen Uptake Method

## Scope

This document explains the current root water uptake and root nitrogen uptake implementation in the MATLAB workflow, with one-to-one correspondence to the code.

It covers:

- plant hydraulic soil-root conductance `gsr`
- node-based root water uptake `RWU`
- layer-endpoint root uptake array `Srt`
- passive NH4 / NO3 uptake
- demand-driven active NO3 uptake
- call order inside the main workflow

Main code locations:

- [src/+plant/+PlantHydraulic/Plantrootconductivity.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/Plantrootconductivity.m:1)
- [src/+plant/+PlantHydraulic/calSPAC.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/calSPAC.m:1)
- [src/+plant/+PlantHydraulic/solveEbalHydraulics.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/solveEbalHydraulics.m:1)
- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:1)
- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:577)
- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:1)
- [src/N_PARM.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/N_PARM.m:1)

## Variable Conventions

### Grid conventions

- `NL`: number of soil layers
- `NN = NL + 1`: number of soil nodes
- `DeltZ(NLx1)`: layer thickness in cm

### Root uptake conventions

- `gsr(NNx1)`: nodal soil-root conductance, unit `[m s^-1 MPa^-1]`
- `RWU(NNx1)`: nodal root water uptake sink/source strength, unit `[1/s]`
- `Srt(NLx2)`: layer-endpoint root water uptake array

The mapping from node uptake to layer-endpoint uptake is:

```text
Srt(ML,1) = RWU(ML)
Srt(ML,2) = RWU(ML+1)
```

implemented in [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:592):

```matlab
Srt = [RWU(1:NL), RWU(2:NN)];
```

### Root nitrogen uptake conventions

- `SinkS(NLx2x2)`: root solute sink terms
- `SinkS(:,:,1)`: NH4 uptake sink
- `SinkS(:,:,2)`: NO3 uptake sink
- `SPUptake(1)`: depth-integrated passive NH4 uptake
- `SPUptake(2)`: depth-integrated passive NO3 uptake
- `SAUptakeP`: potential active NO3 uptake still required after passive uptake

## Runtime Call Order

The runtime chain is:

1. `AMainloop` calls `ebal_clm45`
2. `ebal_clm45` calls `solveEbalHydraulics`
3. `solveEbalHydraulics` calls `Plant_hydraulics`
4. `Plant_hydraulics` calls `calSPAC`
5. `calSPAC` calls `Plantrootconductivity`
6. `Plantrootconductivity` computes nodal `gsr`
7. `ebal_clm45` computes nodal `RWU`
8. `AMainloop` maps `RWU(NNx1)` to `Srt(NLx2)`
9. `AMainloop` calls `applyRootUptake`
10. `applyRootUptake` computes `SinkS`
11. `N_sub -> N_PARM` integrates `SinkS` into nitrogen mass balance

The main call sites are:

- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:577): call to `ebal_clm45`
- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:592): `RWU -> Srt`
- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:735): call to `applyRootUptake`

## Part 1: Soil-Root Conductance `gsr`

### Function

`Plantrootconductivity` computes the root-weighted soil potential `psis` and the nodal soil-root conductance `gsr`.

See [src/+plant/+PlantHydraulic/Plantrootconductivity.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/Plantrootconductivity.m:1).

### Inputs

- `Khh(NLx2)`: soil hydraulic conductivity at layer endpoints, unit `cm/s`
- `RAI`: root area index
- `hh(NNx1)`: nodal soil pressure head, unit `cmH2O`
- `DeltZ(NLx1)`: layer thickness

### RAI endpoint expansion

If `RAI` is supplied as `NLx1`, it is first converted to nodal support and then mapped back to `NLx2` endpoints:

```matlab
RAI_node(1) = 0.5 * RAI_layer(1);
RAI_node(MN) = 0.5 * (RAI_layer(MN - 1) + RAI_layer(MN));
RAI_node(NN) = 0.5 * RAI_layer(NL);
RAI = [RAI_node(1:NL), RAI_node(2:NN)];
```

Code location:

- [src/+plant/+PlantHydraulic/Plantrootconductivity.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/Plantrootconductivity.m:56)

This prevents a layer-total `RAI` profile from being counted twice at the two endpoints.

### Conductance formula

For each layer endpoint:

```text
gsr_endpoint = K_eff * sqrt(RAI) / (Zrj * pi * rho_w * g) * 10^6
```

where:

- `K_eff = K_m_s * f_contact`
- `K_m_s = Khh / 100` converts `cm/s` to `m/s`
- `f_contact` is the root-soil contact reduction factor
- `Zrj = DeltZ / 100` converts `cm` to `m`
- `10^6` converts conductance from per Pa to per MPa

Implemented in:

- [src/+plant/+PlantHydraulic/Plantrootconductivity.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/Plantrootconductivity.m:111)

```matlab
K_m_s = Khh(ML, ND) / 100;
f_contact = rootSoilContactFactor(hhMPA(MN), rootContact);
K_eff_m_s = K_m_s * f_contact;
layer_endpoint_gsr = K_eff_m_s * sqrt(RAI(ML, ND)) / depth_factor * pa_per_mpa;
gsr(MN) = gsr(MN) + layer_endpoint_gsr;
```

### Root-zone weighted soil potential

After all nodal conductances are assembled:

```text
psis = sum(gsr_i * psi_soil_i) / sum(gsr_i)
```

Implemented in:

- [src/+plant/+PlantHydraulic/Plantrootconductivity.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/Plantrootconductivity.m:143)

## Part 2: Plant Hydraulic Solver

### Function role

`solveEbalHydraulics` packages the current soil state and transpiration demand and calls the plant hydraulic solver.

Key output fields:

- `hydraulic.gsr`
- `hydraulic.psis`
- `hydraulic.psir`
- `hydraulic.Tp_t`
- `hydraulic.supply_ratio`

See:

- [src/+plant/+PlantHydraulic/solveEbalHydraulics.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/solveEbalHydraulics.m:55)

### Shared-supply logic

`calSPAC` solves the actual transpiration supply under one shared soil-root-stem hydraulic path for sunlit and shaded leaves.

Demand conversion:

```text
demand_sun_m_s = Tp_sun / 100
demand_sha_m_s = Tp_sha / 100
demand_total_m_s = demand_sun_m_s + demand_sha_m_s
```

Actual root potential:

```text
psir = psis - actual_total_m_s / sumgsr
```

Implemented in:

- [src/+plant/+PlantHydraulic/calSPAC.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+PlantHydraulic/calSPAC.m:82)

## Part 3: Node-Based Root Water Uptake `RWU`

### Function

`ebal_clm45` converts nodal conductance and potential gradients into nodal root water uptake.

See:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:878)

### Potential gradient

```text
grad_node = psi_soil_node - psi_root
```

Code:

```matlab
grad_node = hhMPA - psir_sc;
```

Location:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:891)

### Node tributary length

Because `RWU` is now nodal, each node uses its own tributary control length:

```text
DeltZ_node(1)   = DeltZ(1)/2
DeltZ_node(MN)  = (DeltZ(MN-1) + DeltZ(MN))/2
DeltZ_node(NN)  = DeltZ(NL)/2
```

Code:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:952)

### Passive redistribution branch

The nodal hydraulic redistribution flux is:

```text
flux_hr_i = gsr_i / (rectifier_i * dry_penalty_i) * grad_i
```

Then convert flux to nodal volumetric sink/source strength:

```text
RWU_hr_raw_i = flux_hr_i / DeltZ_node_i
```

Implemented in:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:961)

### Active transpiration branch

Only positive soil-to-root flow contributes:

```text
flux_active_i = gsr_i * max(grad_i, 0)
RWU_active_capacity_i = flux_active_i / DeltZ_node_i
```

Then scale the active capacity so that the integrated active uptake equals the actual transpiration supply:

```text
target_uptake = Tp_t / 100
RWU_active = RWU_active_capacity * target_uptake / sum(RWU_active_capacity * DeltZ_node)
```

Implemented in:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:991)

### Day/night logic

Nighttime:

- no net plant uptake
- redistribution is allowed
- mean nodal redistribution is removed so the net integral is zero

Code:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:981)

Daytime:

- active uptake satisfies `Tp_t`
- redistribution is computed separately
- redistribution is further suppressed by `daytime_hr_suppression`

Code:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:988)

### Unit conversion

Before return:

```text
RWU = RWU * 100
```

because:

```text
[m/s/cm] * 100 = [cm/s/cm] = [1/s]
```

Code:

- [src/ebal_clm45.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/ebal_clm45.m:1021)

## Part 4: Mapping `RWU` To `Srt`

After `ebal_clm45` returns `RWU(NNx1)`, `AMainloop` maps it to `Srt(NLx2)`:

```matlab
Srt = [RWU(1:NL), RWU(2:NN)];
```

This means:

- left column: top node of each layer
- right column: bottom node of each layer

Code:

- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:592)

This `Srt(NLx2)` format is then used consistently by the nitrogen uptake routine and by output writing:

- [src/AMainloop.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/AMainloop.m:336)

## Part 5: Root Nitrogen Uptake

### Function

`applyRootUptake` computes passive NH4 / NO3 uptake and optional active NO3 uptake.

See:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:1)

### Input compatibility

The routine accepts:

- `NNx1` nodal `RWU`
- `NLx1` layer uptake
- `NLx2` endpoint `Srt`

and internally converts to:

```text
Srt_nodes(NLx2)
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:30)

### Passive uptake

Current passive uptake uses pure mass-flow logic:

```text
SinkS_NH4(ML,ND) = Srt(ML,ND) * max(C_NH4(MN), 0)
SinkS_NO3(ML,ND) = Srt(ML,ND) * max(C_NO3(MN), 0)
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:100)

This means:

- passive NH4 and passive NO3 are computed separately
- no empirical concentration cap is applied
- only non-negative concentrations are allowed

### Passive uptake integration

Depth-integrated passive uptake is:

```text
SPUptake(1) = sum over layers [ mean(endpoint NH4 sink) * DeltZ ]
SPUptake(2) = sum over layers [ mean(endpoint NO3 sink) * DeltZ ]
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:105)

### Demand-driven active uptake

Active uptake is now demand-driven rather than driven by the legacy `SPOT * wt` proxy.

The active uptake trigger is:

```text
SAUptakeP = max(N_demand - passiveTotal, 0)
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:121)

If:

- `lActRSU ~= 1`, or
- `N_demand <= 0`

the routine returns after passive uptake.

### Active uptake allocation weight

The active uptake spatial weight follows the already allocated water uptake profile:

```text
weightDenom = sum over layers [ mean(Srt endpoints) * DeltZ ]
uptakeWeight(ML,ND) = Srt(ML,ND) * DeltZ(ML) / weightDenom
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:83)

This means active uptake is not redistributed by a separate empirical root profile. It follows the current root water uptake allocation.

### Active NO3 uptake formula

Only NO3 has an active uptake term.

The local concentration response is Michaelis-Menten type:

```text
f_NO3 = cc / (rKM + cc)
cc = max(C_NO3 - Nmin, 0)
```

Then:

```text
AUptakeA(ML,ND) = f_NO3 * uptakeWeight(ML,ND) * SAUptakeP
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:133)

### Compensation factor

The uncompensated active uptake sum is:

```text
SAUptake = sum over layers [ mean(AUptakeA endpoints) * DeltZ ]
```

Then:

```text
pai = SAUptake / SAUptakeP
Compen = max(pai, paic) when 0 < pai < paic
Compen = pai otherwise
```

with `paic = 1` in the current code.

Final compensated active uptake added into nitrate sink:

```text
AUptakeA = f_NO3 * uptakeWeight * SAUptakeP / Compen
SinkS(:,:,2) = SinkS(:,:,2) + AUptakeA
```

Code:

- [src/+plant/+root/applyRootUptake.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/+plant/+root/applyRootUptake.m:159)

### Why only NO3 is active

The current implementation keeps:

- passive NH4 uptake
- passive NO3 uptake
- active NO3 uptake only

This matches the legacy `Evap_Cal_FAO` logic and avoids adding an unsupported NH4 active-uptake parameterization before a full plant N-demand model exists.

## Part 6: Nitrogen Mass Balance Coupling

`SinkS` is passed into the nitrogen solver and integrated to obtain bulk uptake fluxes:

```matlab
SinkS1_avg = (SinkS(:,1,1) + SinkS(:,2,1)) * 0.5;
SinkS2_avg = (SinkS(:,1,2) + SinkS(:,2,2)) * 0.5;
cvSinkS(KT)  = DeltZ' * SinkS1_avg;
cvSinkS1(KT) = DeltZ' * SinkS2_avg;
```

Code:

- [src/N_PARM.m](/f:/STEMMUS_Zhujizixing%20-%E8%8B%B9%E6%9E%9C17%20-%20%E5%8F%98%E6%AD%A5%E9%95%BF/src/N_PARM.m:135)

This means:

- `cvSinkS` is depth-integrated NH4 uptake
- `cvSinkS1` is depth-integrated NO3 uptake

## Part 7: Current Design Decisions

### Decision 1: `RWU` is node-based

Reason:

- plant hydraulic conductance `gsr` is assembled at nodes
- pressure head `hh` is node-based
- using `RWU(NNx1)` keeps the hydraulic calculation internally consistent

### Decision 2: `Srt` is `NLx2`

Reason:

- solute uptake and old root-uptake code use layer-endpoint form
- `Srt(NLx2)` is the natural bridge between node-based hydraulics and layer-based nitrogen transport

### Decision 3: active uptake is demand-driven

Reason:

- the old `SPOT * wt` term mixed supply-side and demand-side meaning
- the new interface leaves a clean hook for future plant N-demand models

### Decision 4: passive uptake has no empirical concentration cap

Reason:

- passive uptake is treated as `Srt * C`
- this keeps the physical meaning simple
- only non-negative concentration clipping remains

## Part 8: Remaining Interface Hook

The current `applyRootUptake` interface already reserves:

```matlab
N_demand
```

Accepted forms:

- `[]`
- numeric scalar
- struct with field `.total`

Current behavior:

- if `N_demand` is empty, active uptake is off
- therefore the present model is effectively passive-only unless an external plant N-demand module is added

## Part 9: Legacy Relation To `Evap_Cal_FAO`

The current implementation retains several ideas from `Evap_Cal_FAO`:

- passive NH4 / NO3 uptake from `Srt`
- Michaelis-Menten type active NO3 response
- compensation-style scaling through `pai`

But it changes several core points:

- no legacy empirical `bx` root profile inside solute uptake
- no `SPOT * wt` demand proxy
- no passive concentration cap
- explicit `N_demand` interface
- node-based hydraulic `RWU` with `NLx2` `Srt` bridge

## Summary

The present workflow is:

1. build nodal soil-root conductance `gsr`
2. solve plant hydraulic supply
3. compute nodal `RWU`
4. map `RWU` to layer-endpoint `Srt`
5. compute passive NH4 / NO3 uptake
6. if plant N demand is supplied, compute active NO3 uptake
7. pass `SinkS` into the nitrogen transport and balance solver

This is the current root water and root nitrogen uptake implementation used by the new `AMainloop + ebal_clm45 + applyRootUptake` path.
