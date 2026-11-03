function [stateNext, flux] = advanceBulkSnow(state, parameters, drivers, dt_s)
%ADVANCEBULKSNOW  Advance and thermodynamically close one bulk-snow step.
%
%   This MATLAB port preserves ClimaLand v1.0 SnowModel's prognostic states
%   S/S_l/U, auxiliary-update order, runoff and precipitation energy, snow
%   cover weighting, water/energy clipping, explicit phase-change flux, and
%   liquid-state clipping. State is stored in COMPASS units (cm, J cm-2), but
%   all ClimaLand equations are evaluated in their native SI units.  Because
%   COMPASS accepts one outer-step transaction after this explicit update,
%   the updated liquid/ice partition is projected at fixed S and U onto the
%   bulk phase-equilibrium manifold.  This prevents a long Euler step from
%   accepting the impossible combination of finite ice and T > freezing.
%
%   Driver conventions for this adapter are positive downward/input:
%     rainfall_cm_s, snowfall_cm_s, snow_to_atmosphere_vapor_cm_s,
%     surface_closure_flux_W_m2, snow_to_soil_heat_flux_W_m2.
%   frozen_pond_transfer_cm is added as solid precipitation over this step.
%
%   Source: ClimaLand.jl v1.0, src/standalone/Snow/{Snow.jl,
%   snow_parameterizations.jl,boundary_fluxes.jl}, Apache-2.0,
%   Copyright 2022 California Institute of Technology.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-07-17 Jizixing Zhu -- accept separate snow-tile and
%     bare-soil surface closures from the SCOPE ground-flux partition.
%   Change: 2026-07-17 Jizixing Zhu -- conservatively drain an ice-free bulk-snow
%     state so liquid meltwater cannot persist as a summertime snow tile.
%   Change: 2026-07-17 Jizixing Zhu -- allow the current-step predicted
%     surface cover to weight SCOPE heat/vapor exchange while storage still
%     advances exactly once from the accepted snow state.
%   Change: 2026-07-18 Jizixing Zhu -- conservatively project the accepted
%     liquid/ice partition onto the phase-equilibrium manifold at fixed SWE
%     and internal energy, and include that conversion in the liquid ledger.
%   Change: 2026-07-18 Jizixing Zhu -- drain every ice-free remnant together with
%     its retained enthalpy instead of requiring numerically zero energy.
%   Change: 2026-07-18 Jizixing Zhu -- drain a near-liquid tail when its
%     remaining solid latent energy is below the coupled energy resolution.
%   Change: 2026-07-18 Jizixing Zhu -- expose bare and snow contributions
%     to the one aggregate energy flux entering the shared soil column.
%   Change: 2026-07-19 Jizixing Zhu -- expose the signed snow external-energy
%     terms required by the combined surface-snow-soil conservation ledger.

    assert(dt_s > 0 && isfinite(dt_s), ...
        'advanceBulkSnow:badDt', 'dt_s must be positive and finite.');
    old = surface.snow.diagnoseState(state, parameters);
    S = old.SWE_cm / 100;
    Sl = old.liquid_water_cm / 100;
    U = old.energy_J_cm2 * 1e4;
    rain = scalarDriver(drivers, 'rainfall_cm_s', 0) / 100;
    snowfall = scalarDriver(drivers, 'snowfall_cm_s', 0) / 100;
    vapor = scalarDriver(drivers, 'snow_to_atmosphere_vapor_cm_s', 0) / 100;
    surfaceClosure = scalarDriver(drivers, 'surface_closure_flux_W_m2', 0);
    bareSurfaceClosure = scalarDriver( ...
        drivers, 'bare_surface_closure_flux_W_m2', surfaceClosure);
    snowToSoil = scalarDriver(drivers, 'snow_to_soil_heat_flux_W_m2', 0);
    airT_C = scalarDriver(drivers, 'air_temperature_C', 0);
    frozenPondM = max(scalarDriver( ...
        drivers, 'frozen_pond_transfer_cm', 0), 0) / 100;

    % ClimaLand precipitation drivers are negative toward the surface.
    P_snow = -snowfall - frozenPondM / dt_s;
    P_liq = -rain;
    scf = old.cover_fraction;
    if isfield(drivers, 'surface_cover_fraction') ...
            && ~isempty(drivers.surface_cover_fraction)
        scf = double(drivers.surface_cover_fraction);
        assert(isscalar(scf) && isfinite(scf) && scf >= 0 && scf <= 1, ...
            'advanceBulkSnow:badSurfaceCover', ...
            'surface_cover_fraction must be a finite scalar in [0,1].');
    end
    rhoL = parameters.rho_liquid_kg_m3;
    cpL = parameters.specific_heat_liquid_J_kg_K;
    TairK = airT_C + parameters.freezing_temperature_K;

    qmax = maximumLiquidFraction(old, parameters);
    excessLiquid = max(Sl - qmax * max(S, 0), 0);
    tau = max(dt_s, old.depth_m ...
        / parameters.wet_snow_hydraulic_conductivity_m_s);
    waterRunoff = -excessLiquid / tau; % [m/s], negative leaves snow downward
    liquidEnergyDensity = rhoL * cpL ...
        * (old.temperature_K - parameters.reference_temperature_K);
    energyRunoff = waterRunoff * liquidEnergyDensity; % [W/m2]

    totalWater = P_snow + (P_liq + vapor - waterRunoff) * scf;
    liquidWater = (P_liq + vapor * old.liquid_fraction ...
        - waterRunoff) * scf;
    fallingSnowEnergy = -rhoL * parameters.latent_heat_fusion_J_kg * P_snow;
    fallingRainEnergy = rhoL * cpL ...
        * (TairK - parameters.reference_temperature_K) * P_liq;
    % ClimaLand fluxes are positive toward the atmosphere. SCOPE's closure
    % flux and the snow-to-soil flux are positive downward.
    totalEnergy = fallingSnowEnergy + ( ...
        -surfaceClosure - energyRunoff + snowToSoil ...
        + fallingRainEnergy) * scf;

    appliedWater = clipWaterFlux(S, totalWater, dt_s);
    appliedEnergy = clipEnergyFlux(U, S, totalEnergy, totalWater, dt_s);
    phaseFlux = phaseChangeFlux( ...
        U, S, old.liquid_fraction, appliedEnergy, parameters, dt_s);
    liquidWater = liquidWater + phaseFlux * scf;
    appliedLiquid = clipLiquidFlux( ...
        Sl, S, liquidWater, appliedWater, dt_s);

    Snext = S - appliedWater * dt_s;
    Slnext = Sl - appliedLiquid * dt_s;
    Unext = U - appliedEnergy * dt_s;
    scale = max([1, abs(S), abs(Snext)]);
    if abs(Snext) <= 32 * eps(scale), Snext = 0; end
    if abs(Slnext) <= 32 * eps(scale), Slnext = 0; end
    if abs(Unext) <= 32 * eps(max([1, abs(U), abs(Unext)])), Unext = 0; end
    Snext = max(Snext, 0);
    Slnext = min(max(Slnext, 0), Snext);
    Unext = min(Unext, 0);
    % ClimaLand normally integrates the snow ODE with internal solver steps.
    % A single COMPASS outer-hour Euler update can otherwise leave S_l lagging
    % U, e.g. retaining ice after U already implies T > T_freeze.  Restore the
    % algebraic phase closure without changing total SWE or snow energy.  The
    % adjustment is an internal melt/freeze flux and therefore belongs in the
    % liquid-state ledger, not the total-water or energy ledgers.
    [Slnext, equilibriumLiquidFlux] = projectPhaseEquilibrium( ...
        Snext, Slnext, Unext, parameters, dt_s);
    appliedLiquid = appliedLiquid + equilibriumLiquidFlux;
    % The ClimaLand flux kernel asymptotically drains a fully liquid remnant
    % because its runoff is area-weighted by snow-cover fraction.  In the
    % COMPASS adapter that remnant would otherwise remain a snow tile and
    % switch the entire soil column to a snow heat-flux boundary throughout
    % summer. Once the solid phase is exhausted, route the remaining liquid
    % and retained enthalpy to the soil-facing interface in this same
    % conservative transaction. This is a phase-domain closure, not a SWE
    % cutoff: any energetically resolved ice keeps the prognostic snowpack,
    % while an ice-free or sub-resolution solid tail belongs to the
    % surface/soil liquid domain even at nonzero retained enthalpy.  The
    % criterion is solid latent energy, not total SWE, so a large liquid store
    % cannot remain a snow tile solely because of an unresolvable ice trace.
    terminalDrainage = 0; % [m water over the grid cell]
    terminalEnergy = 0; % [J m-2], carried with terminal liquid drainage
    terminalPhaseAdjustment = 0; % [m s-1], negative creates terminal liquid
    solidDepth = max(Snext - Slnext, 0); % [m water equivalent]
    solidEnergy = rhoL * parameters.latent_heat_fusion_J_kg ...
        * solidDepth; % [J m-2]
    solidResolution = solidPhaseEnergyResolution(parameters) * 1e4; % [J m-2]
    solidRoundoff = 32 * eps(max([1, abs(Snext), abs(Slnext)]));
    unresolvedSolid = solidDepth <= solidRoundoff ...
        || solidEnergy <= solidResolution;
    if Snext > 0 && unresolvedSolid
        % Reclassify the unresolved final solid mass as liquid before the
        % complete store is drained.  This closes the liquid-phase ledger;
        % total water and U are unchanged, and the original U is carried by
        % terminalEnergy into the receiving surface/soil interface.
        terminalPhaseAdjustment = -solidDepth / dt_s;
        appliedLiquid = appliedLiquid + terminalPhaseAdjustment;
        terminalDrainage = Snext;
        terminalEnergy = Unext;
        Snext = 0;
        Slnext = 0;
        Unext = 0;
    end
    stateNext = struct('SWE', Snext * 100, ...
        'liquid_water', Slnext * 100, 'energy', Unext / 1e4);
    next = surface.snow.diagnoseState(stateNext, parameters);

    excessWater = totalWater - appliedWater;
    excessEnergy = totalEnergy - appliedEnergy;
    liquidSupply = -waterRunoff * scf + rain * (1 - scf) ...
        + max(-excessWater, 0) + terminalDrainage / dt_s;
    extraEvap = max(excessWater, 0);
    bareSoilHeatContribution = (1 - scf) * bareSurfaceClosure;
    snowSoilHeatContribution = scf * snowToSoil ...
        - excessEnergy + terminalEnergy / dt_s;
    soilHeat = bareSoilHeatContribution + snowSoilHeatContribution;
    waterResidualM = (Snext - S) + appliedWater * dt_s ...
        + terminalDrainage;
    liquidResidualM = (Slnext - Sl) + appliedLiquid * dt_s ...
        + terminalDrainage;
    energyResidual = (Unext - U) + appliedEnergy * dt_s + terminalEnergy;

    flux = struct( ...
        'throughfall_liquid_cm_s', liquidSupply * 100, ...
        'additional_soil_evaporation_cm_s', extraEvap * 100, ...
        'soil_evaporation_fraction', 1 - scf, ...
        'snowfall_cm_s', snowfall * 100, ...
        'rainfall_cm_s', rain * 100, ...
        'rain_on_snow_cm_s', rain * scf * 100, ...
        'snow_vapor_flux_cm_s', vapor * scf * 100, ...
        'snow_runoff_cm_s', ...
            (-waterRunoff * scf + terminalDrainage / dt_s) * 100, ...
        'complete_meltout_cm_s', terminalDrainage / dt_s * 100, ...
        'complete_meltout_energy_W_m2', terminalEnergy / dt_s, ...
        'terminal_meltout_triggered', terminalDrainage > 0, ...
        'terminal_solid_energy_J_cm2', solidEnergy / 1e4, ...
        'solid_phase_energy_resolution_J_cm2', solidResolution / 1e4, ...
        'frozen_pond_transfer_cm', frozenPondM * 100, ...
        'total_water_flux_m_s', totalWater, ...
        'applied_water_flux_m_s', appliedWater, ...
        'liquid_water_flux_m_s', appliedLiquid, ...
        'total_energy_flux_W_m2', totalEnergy, ...
        'applied_energy_flux_W_m2', appliedEnergy, ...
        'surface_cover_fraction', scf, ...
        'snow_surface_external_energy_flux_W_m2', surfaceClosure, ...
        'falling_snow_energy_flux_W_m2', fallingSnowEnergy, ...
        'falling_rain_energy_flux_W_m2', fallingRainEnergy, ...
        'snow_liquid_runoff_energy_flux_W_m2', energyRunoff, ...
        'excess_energy_flux_W_m2', excessEnergy, ...
        'phase_change_flux_m_s', phaseFlux, ...
        'phase_equilibrium_adjustment_m_s', equilibriumLiquidFlux, ...
        'terminal_phase_adjustment_m_s', terminalPhaseAdjustment, ...
        'snow_soil_heat_flux_W_m2', snowToSoil, ...
        'bare_soil_heat_contribution_W_m2', bareSoilHeatContribution, ...
        'snow_soil_heat_contribution_W_m2', snowSoilHeatContribution, ...
        'soil_heat_flux_W_m2', soilHeat, ...
        'bare_ground_downward_soil_heat_flux_contribution_W_m2', ...
            bareSoilHeatContribution, ...
        'snow_downward_soil_energy_flux_contribution_W_m2', ...
            snowSoilHeatContribution, ...
        'soil_top_downward_energy_flux_W_m2', soilHeat, ...
        'excess_heat_to_soil_W_m2', -excessEnergy, ...
        'water_balance_residual_cm', waterResidualM * 100, ...
        'liquid_balance_residual_cm', liquidResidualM * 100, ...
        'energy_balance_residual_J_cm2', energyResidual / 1e4, ...
        'state', next);
end

function resolution = solidPhaseEnergyResolution(parameters)
% Read the nonnegative latent-energy resolution for solid snow [J cm-2].
% Editor: Jizixing Zhu; Date: 2026-07-18.
    resolution = 1e-2;
    if isfield(parameters, 'solid_phase_energy_resolution_J_cm2') ...
            && ~isempty(parameters.solid_phase_energy_resolution_J_cm2)
        resolution = double( ...
            parameters.solid_phase_energy_resolution_J_cm2);
    end
    assert(isscalar(resolution) && isfinite(resolution) && resolution >= 0, ...
        'advanceBulkSnow:badSolidResolution', ...
        ['solid_phase_energy_resolution_J_cm2 must be a finite ', ...
         'nonnegative scalar.']);
end

function [SlEq, liquidFluxAdjustment] = projectPhaseEquilibrium( ...
        S, SlRaw, U, p, dt)
%PROJECTPHASEEQUILIBRIUM  Close S_l against S and U at the freezing point.
%
%   E_ice and E_liquid are the bulk energies at T_freeze for q_l=0 and 1.
%   Below E_ice the pack is all ice; above E_liquid it is all liquid; between
%   them latent heat fixes a mixed phase at T_freeze.  Projection preserves S
%   and U exactly.  liquidFluxAdjustment follows advanceBulkSnow's convention:
%   positive removes liquid, negative creates liquid by melting [m s-1].
%
%   The equilibrium relation is the inversion used by ClimaLand SnowModel's
%   bulk snow temperature, applied here as a COMPASS accepted-step closure.
%   Editor: Jizixing Zhu; Date: 2026-07-18.
    if S <= 0
        SlEq = 0;
        liquidFluxAdjustment = (SlRaw - SlEq) / dt;
        return;
    end
    rho = p.rho_liquid_kg_m3;
    dT = p.freezing_temperature_K - p.reference_temperature_K;
    capacityDepth = S + p.regularization_swe_m;
    Eice = rho * capacityDepth * p.specific_heat_ice_J_kg_K * dT ...
        - rho * S * p.latent_heat_fusion_J_kg;
    Eliquid = rho * capacityDepth * p.specific_heat_liquid_J_kg_K * dT;
    if Eliquid > Eice
        if U <= Eice
            qLiquid = 0;
        elseif U >= Eliquid
            qLiquid = 1;
        else
            qLiquid = (U - Eice) / (Eliquid - Eice);
        end
    else
        % For a trace pack, ClimaLand's fixed 0.1-m SWE regularization can
        % make the q-dependent sensible term exceed the pack latent heat,
        % reversing the two endpoint energies.  That artificial interval has
        % no stable mixed-phase branch.  Retain the cold all-ice branch until
        % its freezing-point upper energy is reached, then choose all liquid;
        % the existing terminal-drainage transaction removes a warm liquid
        % remnant without inventing a summertime snow tile.
        qLiquid = double(U >= Eice);
    end
    SlEq = min(max(qLiquid, 0), 1) * S;
    liquidFluxAdjustment = (SlRaw - SlEq) / dt;
end

function qmax = maximumLiquidFraction(d, p)
% ClimaLand maximum_liquid_mass_fraction.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    if d.temperature_K > p.freezing_temperature_K
        qmax = 0;
    else
        qmax = p.liquid_holding_capacity * p.rho_liquid_kg_m3 ...
            / max(d.density_kg_m3, eps);
    end
end

function f = clipWaterFlux(S, f0, dt)
% ClimaLand clip_water_flux.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    f = f0;
    if S - f0 * dt < 0, f = S / dt; end
end

function f = clipEnergyFlux(U, S, f0, water0, dt)
% ClimaLand clip_total_snow_energy_flux.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    f = f0;
    if U - f0 * dt > 0 || S - water0 * dt < 0
        f = U / dt;
    end
end

function f = clipLiquidFlux(Sl, S, f0, appliedWater, dt)
% ClimaLand clip_liquid_water_flux.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    predictedS = S - appliedWater * dt;
    predictedSl = Sl - f0 * dt;
    if predictedSl < 0
        f = Sl / dt;
    elseif predictedSl > predictedS
        f = (Sl - predictedS) / dt;
    else
        f = f0;
    end
end

function flux = phaseChangeFlux(U, S, ql, energyFlux, p, dt)
% ClimaLand phase_change_flux, including its 0.01-K reference offset.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    S = max(S, 0);
    cp = ql * p.specific_heat_liquid_J_kg_K ...
        + (1 - ql) * p.specific_heat_ice_J_kg_K;
    Ufreeze = p.rho_liquid_kg_m3 * (S + p.regularization_swe_m) ...
        * cp * (p.freezing_temperature_K - p.reference_temperature_K) ...
        - p.rho_liquid_kg_m3 * S * (1 - ql) ...
        * p.latent_heat_fusion_J_kg;
    predicted = U - energyFlux * dt;
    excess = predicted - Ufreeze;
    if excess > 0 || (excess < 0 && ql > 0)
        denom = (p.specific_heat_liquid_J_kg_K ...
            - p.specific_heat_ice_J_kg_K) ...
            * (p.freezing_temperature_K - p.reference_temperature_K) ...
            + p.latent_heat_fusion_J_kg;
        flux = -excess / dt / p.rho_liquid_kg_m3 / denom;
    else
        flux = 0;
    end
end

function value = scalarDriver(drivers, name, fallback)
% Read one finite scalar driver.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(drivers, name) && ~isempty(drivers.(name))
        value = double(drivers.(name));
    end
    assert(isscalar(value) && isfinite(value), ...
        'advanceBulkSnow:badDriver', '%s must be finite scalar.', name);
end
