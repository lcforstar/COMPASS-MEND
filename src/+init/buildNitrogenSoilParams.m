function soil_n = buildNitrogenSoilParams(built,config)
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
%     config optional scalar raw site-config struct. The registered optional
%            [nitrogen].nh4_protected_per_clay_mg_kg is a finite nonnegative
%            real numeric scalar [mg N/kg dry soil per clay mass fraction].
%            Missing section/key defaults to zero; unknown keys are rejected.
%            nh4_protection_scope selects nitrification_only or the historical
%            nitrification_and_mend behavior when protection is enabled.
%
%   Returns soil_n (legacy NitrogenVariables block-1 names):
%     RHO_bulk [g/cm3], RHOKG (= KG.*bulk, adsorption buffer), DIF1111
%     [1 x NL x nS] (NH4=:,:,1 / NO3=:,:,2), DISP, K_mineralization, Vmax_nit,
%     Km_nit, Vmax_denit, Km_denit, K_nh4_dissolve, Kurea, pH, MSOC_volume,
%     Theta_f (field capacity).
%     nh4_threshold_enabled/per-clay coefficient are scalar; S_min_nh4_mg_kg
%     [mg N/kg dry soil] and M_min_nit_ug_cm3 [ug N/cm3 bulk soil] are NL-by-1,
%     bottom-up. These read-only thresholds do not add or repartition state.
%     Nonzero protection with explicit microbes requires source_hourly.
%
%   NOT here (-> allocateWorkspace): dynamic N concentration state, the solver
%   work arrays (Retard/Dsh/Smin/Sden/...) and the tS diagnostic buffers.
%   Reference: init.setInitialNitrogen block 1 (src/+init/setInitialNitrogen.m:29-52).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-10-23 Wenjie Liu -- carry the explicit biological pH-response
% switch independently of measured soil pH; old callers retain reference kinetics.
% Change: 2026-11-01 Wenjie Liu -- validate optional NH4 protection and derive
% bottom-up thresholds from mapped clay and the selected final bulk density.

    if nargin<2,config=struct();end
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
    soil_n.ph_response_enabled=true;
    if isfield(built,'control')&&isfield(built.control,'switches') ...
            &&isfield(built.control.switches,'lNitrogenPH')
        flag=built.control.switches.lNitrogenPH;
        assert(isscalar(flag)&&(islogical(flag)||isnumeric(flag)) ...
            &&isfinite(flag)&&ismember(flag,[0,1]), ...
            'buildNitrogenSoilParams:phResponse','lNitrogenPH must be zero or one.');
        soil_n.ph_response_enabled=logical(flag);
    end
    soil_n.MSOC_volume      = s.msoc_volume(:);

    % field capacity (legacy VanGenuchten.Theta_f); from buildMoistureChar (B2)
    soil_n.Theta_f          = s.field_capacity(:);

    % The default branch does not inspect clay or change legacy arithmetic.
    [perClay,scope]=readNh4Protection(config);
    soil_n.nh4_protected_per_clay_mg_kg=perClay;
    soil_n.nh4_threshold_enabled=perClay>0;
    soil_n.S_min_nh4_mg_kg=zeros(NL,1);
    soil_n.M_min_nit_ug_cm3=zeros(NL,1);
    if soil_n.nh4_threshold_enabled
        soil_n.nh4_protection_scope=scope;
        assert(isfield(s,'clay')&&isnumeric(s.clay)&&isreal(s.clay) ...
            &&numel(s.clay)==NL&&all(isfinite(s.clay(:))) ...
            &&all(s.clay(:)>=0&s.clay(:)<=1), ...
            'buildNitrogenSoilParams:nh4ProtectionClay', ...
            'NH4 protection requires NL mapped clay mass fractions in [0,1].');
        rho=soil_n.RHO_bulk;
        assert(numel(rho)==NL&&isreal(rho)&&all(isfinite(rho))&&all(rho>0), ...
            'buildNitrogenSoilParams:nh4ProtectionDensity', ...
            'NH4 protection requires NL finite positive final bulk densities.');
        if explicitMicrobesEnabled(built,config)
            assert(isSourceHourly(config), ...
                'buildNitrogenSoilParams:nh4ProtectionMode', ...
                'Nonzero NH4 protection with explicit soil microbes requires source_hourly.');
        end
        soil_n.S_min_nh4_mg_kg=perClay.*s.clay(:);
        soil_n.M_min_nit_ug_cm3=rho.*soil_n.S_min_nh4_mg_kg;
        assert(all(isfinite(soil_n.S_min_nh4_mg_kg)) ...
            &&all(isfinite(soil_n.M_min_nit_ug_cm3)), ...
            'buildNitrogenSoilParams:nh4ProtectionOverflow', ...
            'Derived NH4 protection profiles must remain finite.');
    end
end

function [value,scope]=readNh4Protection(config)
% Strict optional [nitrogen] contract; booleans are not numeric coefficients.
% Editor: Wenjie Liu; Date: 2026-11-01.
    assert(isstruct(config)&&isscalar(config), ...
        'buildNitrogenSoilParams:nitrogenConfig','Expected one raw site-config struct.');
    value=0;
    scope='nitrification_and_mend';
    if ~isfield(config,'nitrogen'),return;end
    cfg=config.nitrogen;
    assert(isstruct(cfg)&&isscalar(cfg), ...
        'buildNitrogenSoilParams:nitrogenConfig','[nitrogen] must be a scalar configuration section.');
    key='nh4_protected_per_clay_mg_kg';
    unknown=setdiff(fieldnames(cfg),{key,'nh4_protection_scope'});
    assert(isempty(unknown),'buildNitrogenSoilParams:nitrogenKeys', ...
        'Unknown [nitrogen] key(s): %s.',strjoin(unknown,', '));
    if isfield(cfg,'nh4_protection_scope')
        scope=cfg.nh4_protection_scope;
        assert(((ischar(scope)&&isrow(scope))||(isstring(scope)&&isscalar(scope))) ...
            &&any(strcmp(scope,{'nitrification_only','nitrification_and_mend'})), ...
            'buildNitrogenSoilParams:nh4ProtectionScope','Invalid nitrogen.nh4_protection_scope.');
        scope=char(scope);
    end
    if ~isfield(cfg,key),return;end
    value=cfg.(key);
    assert(isnumeric(value)&&isreal(value)&&isscalar(value) ...
        &&isfinite(value)&&value>=0, ...
        'buildNitrogenSoilParams:nh4ProtectionValue', ...
        'nitrogen.%s must be a finite nonnegative real numeric scalar.',key);
    value=double(value);
end

function enabled=explicitMicrobesEnabled(built,config)
% Prefer effective switches, which already cascade from the organic master.
% Editor: Wenjie Liu; Date: 2026-11-01.
    enabled=false;
    if isfield(built,'control')&&isfield(built.control,'switches') ...
            &&isfield(built.control.switches,'lExplicitSoilMicrobes')
        enabled=logical(built.control.switches.lExplicitSoilMicrobes);
    elseif isfield(config,'switches')
        sw=config.switches;
        enabled=isfield(sw,'soil_carbon_model')&&logical(sw.soil_carbon_model) ...
            &&isfield(sw,'explicit_soil_microbes')&&logical(sw.explicit_soil_microbes);
    end
end

function enabled=isSourceHourly(config)
% Editor: Wenjie Liu; Date: 2026-11-01.
    enabled=false;
    if isfield(config,'microbe')&&isstruct(config.microbe) ...
            &&isscalar(config.microbe)&&isfield(config.microbe,'reaction_method')
        method=config.microbe.reaction_method;
        enabled=((ischar(method)&&isrow(method))||(isstring(method)&&isscalar(method))) ...
            &&strcmp(method,'source_hourly');
    end
end
