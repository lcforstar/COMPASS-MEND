function registerAllDiagnostics(reg, cfg, profile)
%REGISTERALLDIAGNOSTICS  Register every variable in every group, gated by cfg.
%
% Order is conventional, not load-bearing. Each per-group register*.m guards
% itself with cfg.lFreezing / cfg.lNefc / cfg.lSoilairefc / cfg.lDualPor as
% appropriate, so this dispatcher can call them all unconditionally.
%
% Inputs:
%   reg  [inout]  diagnostics.Registry
%   cfg  [in]     legacy ModelSettings — config flags decide which groups get enabled
%
% Author: Zhu Jizixing, 2026-05-15
% Change: 2026-07-17 Jizixing Zhu -- register the accepted ClimaLand bulk
% snow state, phase-partition, flux, and closure outputs.
% Change: 2026-07-18 Jizixing Zhu -- register the no-snow zero-storage
% accepted surface-snow-soil energy ledger.
% Change: 2026-08-29 Codex -- register accepted native soil-organic C/N outputs.
% Change: 2026-08-20 Codex -- add a driver profile that retains hourly
% water/heat/N/C states and their conservation/numerical audit fields without
% allocating unrelated canopy, spectral, and remote-sensing time buffers.

    if nargin < 3 || isempty(profile)
        profile = 'full';
    end

    if strcmpi(strtrim(profile), 'driver')
        diagnostics.groups.registerSoilWater(reg, cfg);
        diagnostics.groups.registerSoilThermal(reg, cfg);
        diagnostics.groups.registerSurfaceEnergy(reg, cfg);
        diagnostics.groups.registerSoilAir(reg, cfg);
        diagnostics.groups.registerSoilNitrogen(reg, cfg);
        diagnostics.groups.registerSoilOrganic(reg, cfg);
        diagnostics.groups.registerNumerics(reg, cfg);
        diagnostics.groups.registerEnergyBalance(reg, cfg);
        diagnostics.groups.registerSolverInternal(reg, cfg);
        fprintf('[diag] driver profile registered %d variables in %d groups\n', ...
            numel(reg.list()), numel(reg.groupNames()));
        return;
    elseif ~strcmpi(strtrim(profile), 'full')
        error('diagnostics:registerAllDiagnostics:Profile', ...
            'Unknown diagnostic profile "%s"; expected full or driver.', profile);
    end

    diagnostics.groups.registerSoilWater(reg, cfg);
    diagnostics.groups.registerSoilThermal(reg, cfg);
    diagnostics.groups.registerSurfaceEnergy(reg, cfg);
    diagnostics.groups.registerSoilFreezing(reg, cfg);
    diagnostics.groups.registerSoilAir(reg, cfg);
    diagnostics.groups.registerSoilNitrogen(reg, cfg);
    diagnostics.groups.registerSoilOrganic(reg, cfg);
    diagnostics.groups.registerSoilOxygen(reg, cfg);
    diagnostics.groups.registerPlantCanopy(reg, cfg);
    diagnostics.groups.registerPlantNitrogen(reg, cfg);
    diagnostics.groups.registerPlantHydraulics(reg, cfg);
    diagnostics.groups.registerSurfaceWater(reg, cfg);
    diagnostics.groups.registerSnow(reg, cfg);
    diagnostics.groups.registerEvapotranspiration(reg, cfg);
    diagnostics.groups.registerRootUptake(reg, cfg);
    diagnostics.groups.registerDualPorosity(reg, cfg);
    diagnostics.groups.registerNumerics(reg, cfg);
    diagnostics.groups.registerRadiation(reg, cfg);
    diagnostics.groups.registerRadiationSpectral(reg, cfg);
    diagnostics.groups.registerEnergyBalance(reg, cfg);
    diagnostics.groups.registerCarbonFluxes(reg, cfg);
    diagnostics.groups.registerCanopyAerodynamics(reg, cfg);
    diagnostics.groups.registerCanopyProfiles(reg, cfg);
    diagnostics.groups.registerSolverInternal(reg, cfg);
    diagnostics.groups.registerDeepLearning(reg, cfg);
    if strcmpi(strtrim(getenv('STEMMUS_ONLINE_GS_DIAG')), '1')
        diagnostics.groups.registerOnlineGsValidation(reg, cfg);
    end

    fprintf('[diag] registered %d variables in %d groups\n', ...
        numel(reg.list()), numel(reg.groupNames()));
end
