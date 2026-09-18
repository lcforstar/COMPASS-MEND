function soil_n = buildNitrogenSoilParams(built)
%BUILDNITROGENSOILPARAMS  Static soil-nitrogen solver parameters (D-soil).
%
%   soil_n = init.buildNitrogenSoilParams(built) repackages the static soil
%   nitrogen parameters the N solver needs, reproducing block 1 of the legacy
%   init.setInitialNitrogen (the static-transport / static-rate block). Caller
%   stores it at built.parameters.soil_n. Pure repackaging of already-built
%   per-layer soil parameters -- no physics, no dynamic state.
%
%   Inputs:
%     built  struct  built.parameters.soil (per-layer arrays incl. N rates +
%                    field_capacity) and built.domain.grid.nS
%
%   Returns soil_n (legacy NitrogenVariables block-1 names):
%     RHO_bulk [g/cm3], RHOKG (= KG.*bulk, adsorption buffer), DIF1111
%     [1 x NL x nS] (NH4=:,:,1 / NO3=:,:,2), DISP, K_mineralization, Vmax_nit,
%     Km_nit, Vmax_denit, Km_denit, K_nh4_dissolve, Kurea, pH, MSOC_volume,
%     Theta_f (field capacity).
%
%   NOT here (-> allocateWorkspace): dynamic N concentration state, the solver
%   work arrays (Retard/Dsh/Smin/Sden/...) and the tS diagnostic buffers.
%   Reference: init.setInitialNitrogen block 1 (src/+init/setInitialNitrogen.m:29-52).
% Editor: Jizixing Zhu; Date: 2026-06-27.

    s  = built.parameters.soil;
    nS = built.domain.grid.nS;
    NL = numel(s.disp);

    soil_n = struct();

    % static transport coefficients
    soil_n.RHO_bulk = s.bulk_density(:);
    soil_n.RHOKG    = s.kg(:) .* s.bulk_density(:);          % adsorption buffer [-]

    soil_n.DIF1111          = zeros(1, NL, nS);              % solute molecular diffusion layout
    soil_n.DIF1111(1, :, 1) = s.dif_nh4(:).';               % NH4 phase
    soil_n.DIF1111(1, :, 2) = s.dif_no3(:).';               % NO3 phase
    soil_n.DISP             = s.disp(:);

    % static biochemical rates
    soil_n.K_mineralization = s.k_mineralization(:);
    soil_n.Vmax_nit         = s.vmax_nit(:);
    soil_n.Km_nit           = s.km_nit(:);
    soil_n.Vmax_denit       = s.vmax_denit(:);
    soil_n.Km_denit         = s.km_denit(:);
    soil_n.K_nh4_dissolve   = s.k_nh4_dissolve(:);
    soil_n.Kurea            = s.kurea(:);
    soil_n.pH               = s.ph(:);
    soil_n.MSOC_volume      = s.msoc_volume(:);

    % field capacity (legacy VanGenuchten.Theta_f); from buildMoistureChar (B2)
    soil_n.Theta_f          = s.field_capacity(:);
end
