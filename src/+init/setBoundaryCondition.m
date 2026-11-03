function BoundaryCondition = setBoundaryCondition(cfg, SoilProperties, GroundwaterSettings, MeteoData)
% setBoundaryCondition - Initialize boundary-condition settings.
%
% Inputs
%   cfg                 Model settings object.
%   SoilProperties      Soil hydraulic properties.
%   GroundwaterSettings Groundwater coupling settings.
%   MeteoData           Meteorological forcing; Pg is in Pa.

bc = cfg.BC;

% Water head boundary.
% NBCh/NBChB: 1 = Dirichlet, 2 = flux, 3 = atmospheric/free drainage.
if ~ismember(bc.NBCh_top, [1 2 3])
    error('setBoundaryCondition:InvalidBC', ...
        'Invalid top water boundary type NBCh_top=%d; expected 1, 2, or 3.', bc.NBCh_top);
end
if ~ismember(bc.NBCh_bottom, [1 2 3])
    error('setBoundaryCondition:InvalidBC', ...
        'Invalid bottom water boundary type NBCh_bottom=%d; expected 1, 2, or 3.', bc.NBCh_bottom);
end

BoundaryCondition.NBCh  = bc.NBCh_top;
BoundaryCondition.NBChB = bc.NBCh_bottom;
BoundaryCondition.BCh   = bc.h_top;
BoundaryCondition.BChB  = bc.h_bottom;

% Dynamic top water-head boundary flag used by the main loop.
if bc.NBCh_top ~= 1
    BoundaryCondition.NBChh = 2;
else
    BoundaryCondition.NBChh = 1;
end

% Temperature boundary.
% NBCT/NBCTB: 1 = Dirichlet, 2 = heat flux, 3 = atmospheric energy balance.
if cfg.lThmrlefc == 1
    if ~ismember(bc.NBCT_top, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid top temperature boundary type NBCT_top=%d; expected 1, 2, or 3.', bc.NBCT_top);
    end
    if ~ismember(bc.NBCT_bottom, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid bottom temperature boundary type NBCT_bottom=%d; expected 1, 2, or 3.', bc.NBCT_bottom);
    end

    BoundaryCondition.NBCT  = bc.NBCT_top;
    BoundaryCondition.NBCTB = bc.NBCT_bottom;
    BoundaryCondition.BCT   = bc.T_top;
    BoundaryCondition.BCTB  = bc.T_bottom;
end

% Air pressure boundary.
% NBCP: 1 = Dirichlet, 2 = specified flux/escape, 3 = atmospheric forcing.
% [2026-05-26] renamed NBCa_* -> NBCP_* (config field name).
if cfg.lSoilairefc == 1
    if ~ismember(bc.NBCP_top, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid top air-pressure boundary type NBCP_top=%d; expected 1, 2, or 3.', bc.NBCP_top);
    end
    if ~ismember(bc.NBCP_bottom, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid bottom air-pressure boundary type NBCP_bottom=%d; expected 1, 2, or 3.', bc.NBCP_bottom);
    end

    BoundaryCondition.NBCP  = bc.NBCP_top;
    BoundaryCondition.NBCPB = bc.NBCP_bottom;
    BoundaryCondition.BCP   = bc.P_top;
    BoundaryCondition.BCPB  = bc.P_bottom;
end

% Nitrogen boundary.
% [2026-05-26] renamed: NBCN_X -> NBC_X_top, NBCNB_X -> NBC_X_bottom,
%                       BCN_X -> BC_X_top,  BCNB_X -> BC_X_bottom (X in NH4/NO3/UREA).
% Internal BoundaryCondition.NBCN_NH4 / NBCNB_NH4 / BCN_NH4 / BCNB_NH4 names KEPT for
% downstream solver interface compatibility (NitrogenBoundary struct field names).
if cfg.lNefc == 1
    if ~ismember(bc.NBC_NH4_top, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid top NH4 boundary type NBC_NH4_top=%d; expected 1, 2, or 3.', bc.NBC_NH4_top);
    end
    if ~ismember(bc.NBC_NH4_bottom, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid bottom NH4 boundary type NBC_NH4_bottom=%d; expected 1, 2, or 3.', bc.NBC_NH4_bottom);
    end
    if ~ismember(bc.NBC_NO3_top, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid top NO3 boundary type NBC_NO3_top=%d; expected 1, 2, or 3.', bc.NBC_NO3_top);
    end
    if ~ismember(bc.NBC_NO3_bottom, [1 2 3])
        error('setBoundaryCondition:InvalidBC', ...
            'Invalid bottom NO3 boundary type NBC_NO3_bottom=%d; expected 1, 2, or 3.', bc.NBC_NO3_bottom);
    end

    BoundaryCondition.NBCN_NH4  = bc.NBC_NH4_top;       % keep internal name (downstream solver expects)
    BoundaryCondition.NBCNB_NH4 = bc.NBC_NH4_bottom;
    BoundaryCondition.BCN_NH4   = bc.BC_NH4_top;
    BoundaryCondition.BCNB_NH4  = bc.BC_NH4_bottom;
    BoundaryCondition.NBCN_NO3  = bc.NBC_NO3_top;
    BoundaryCondition.NBCNB_NO3 = bc.NBC_NO3_bottom;
    BoundaryCondition.BCN_NO3   = bc.BC_NO3_top;
    BoundaryCondition.BCNB_NO3  = bc.BC_NO3_bottom;
end

% Local atmospheric pressure used when a fixed air-pressure boundary is needed.
BoundaryCondition.BtmPg = 89821.33497 * 10; % Pa*10, STEMMUS internal pressure unit.

% Runtime ponding and runoff placeholders.
BoundaryCondition.DSTOR  = 0;
BoundaryCondition.DSTOR0 = 0;
BoundaryCondition.RS     = 0;
BoundaryCondition.DSTMAX = 0;
BoundaryCondition.Delt_t_outer = cfg.Time.Delt_t;
BoundaryCondition.IRPT1  = 0;
BoundaryCondition.IRPT2  = 0;

% Soil, groundwater, and meteorology derived boundary fields.
BoundaryCondition.fmax                = SoilProperties.fmax;
BoundaryCondition.Ks0                 = SoilProperties.Ks0;
BoundaryCondition.GroundwaterSettings = GroundwaterSettings;
BoundaryCondition.GroundwaterCoupling = GroundwaterSettings.GroundwaterCoupling;
BoundaryCondition.TopPg               = MeteoData.Pg(:) * 10; % Pa -> Pa*10.

% When the bottom air-pressure boundary is Dirichlet, calculate it from
% current surface pressure plus a small hydrostatic pressure increase.
BoundaryCondition.Tot_Depth = cfg.Grid.Tot_Depth;
BoundaryCondition.UseHydrostaticAirBottom = true;

% Fertilizer parameters.
BoundaryCondition.Fertilizer.t_fert_str = bc.fert_t_fert;
BoundaryCondition.Fertilizer.fert_kg_ha = str2double(strsplit(strtrim(bc.fert_amount), ';'));
BoundaryCondition.Fertilizer.DEPTHSTART = bc.fert_DEPTHSTART;
BoundaryCondition.Fertilizer.DEPTHEND   = bc.fert_DEPTHEND;

n = numel(BoundaryCondition.Fertilizer.fert_kg_ha);
types = strsplit(strtrim(bc.fert_type), ';');
if isscalar(types)
    types = repmat(types, 1, n);
end

parse_ratios = @(s) str2double(strsplit(strtrim(s), ';'));
if ~isempty(bc.fert_NH4_ratio)
    nh4_custom = parse_ratios(bc.fert_NH4_ratio);
else
    nh4_custom = zeros(1, n);
end
if ~isempty(bc.fert_NO3_ratio)
    no3_custom = parse_ratios(bc.fert_NO3_ratio);
else
    no3_custom = zeros(1, n);
end
if ~isempty(bc.fert_urea_ratio)
    urea_custom = parse_ratios(bc.fert_urea_ratio);
else
    urea_custom = ones(1, n);
end

NH4_r = zeros(1, n);
NO3_r = zeros(1, n);
urea_r = zeros(1, n);
for ii = 1:n
    switch lower(strtrim(types{ii}))
        case 'ammonium'
            r = [1.0, 0.0, 0.0];
        case 'nitrate'
            r = [0.0, 1.0, 0.0];
        case 'ammonium_nitrate'
            r = [0.5, 0.5, 0.0];
        case {'urea', 'lime_nitrogen'}
            r = [0.0, 0.0, 1.0];
        case 'custom'
            r = [nh4_custom(ii), no3_custom(ii), urea_custom(ii)];
        otherwise
            warning('setBoundaryCondition:UnknownFertilizerType', ...
                'Unknown fertilizer type at event %d: %s. Using urea.', ii, types{ii});
            r = [0.0, 0.0, 1.0];
    end

    NH4_r(ii) = r(1);
    NO3_r(ii) = r(2);
    urea_r(ii) = r(3);
end

BoundaryCondition.Fertilizer.NH4_ratio  = NH4_r;
BoundaryCondition.Fertilizer.NO3_ratio  = NO3_r;
BoundaryCondition.Fertilizer.urea_ratio = urea_r;

% Build fertilizer controller only when nitrogen is enabled.
if cfg.lNefc == 1 && ~isempty(bc.fert_t_fert)
    F = BoundaryCondition.Fertilizer;
    fert_params = struct( ...
        't_fert',     datetime(strsplit(strtrim(F.t_fert_str), ';'), ...
                               'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss'), ...
        'fert_kg_ha', F.fert_kg_ha, ...
        'NH4_ratio',  F.NH4_ratio, ...
        'NO3_ratio',  F.NO3_ratio, ...
        'urea_ratio', F.urea_ratio, ...
        'depth_cm',   F.DEPTHEND - F.DEPTHSTART);
    BoundaryCondition.Fertilizer.controller = soil.nitrogen.fertilizer.initFertilizerController(fert_params);
else
    BoundaryCondition.Fertilizer.controller = [];
end

end
