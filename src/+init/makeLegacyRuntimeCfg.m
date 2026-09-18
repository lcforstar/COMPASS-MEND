function cfg = makeLegacyRuntimeCfg(built)
%MAKELEGACYRUNTIMECFG  Small legacy-shaped config view for diagnostics/allocators.
%
%   cfg = init.makeLegacyRuntimeCfg(built) assembles the minimal legacy-cfg struct
%   that transitional legacy consumers still expect (diagnostics.registerAllDiagnostics,
%   init.buildSoilModes, io.saveCompassPeriodicPartialOutput): switch flags + NIT +
%   checkpoint/diagnostic intervals + Veg.hydraulic_scheme. A compat shim during the
%   COMPASS.m -> new-chain migration, not a long-term container. Partial-output callers
%   augment the result with OutputPath + Time (env/time-derived, not in built).
%
%   Inputs:
%     built  struct  init.buildState output (control.switches/numerics, parameters.veg)
%   Returns:
%     cfg    struct  legacy-shaped config view
%
%   Editor: Jizixing Zhu; Date: 2026-07-01.
%   Change: 2026-07-01 Jizixing Zhu — extracted from init.prepareRuntime local to a package function so runOneStep periodic partial-output builds the same cfg (single source).
    cfg = struct();
    sw = built.control.switches;
    names = fieldnames(sw);
    for i = 1:numel(names)
        cfg.(names{i}) = sw.(names{i});
    end
    cfg.NIT = built.control.numerics.NIT;
    cfg.lSavePicardHistory = true;
    cfg.checkpointInterval = built.control.numerics.checkpointInterval;
    cfg.diagnosticInterval = built.control.numerics.diagnosticInterval;
    cfg.lPeriodicPartialOutput = 0;
    cfg.Veg = struct();
    cfg.Veg.hydraulic_scheme = built.parameters.veg.hydraulic.hydraulic_scheme;
end
