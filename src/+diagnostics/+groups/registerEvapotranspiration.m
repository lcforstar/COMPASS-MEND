function registerEvapotranspiration(reg, cfg) %#ok<INUSD>
%REGISTEREVAPOTRANSPIRATION  Canonical evaporation and transpiration outputs.
% Output: data/output_<tag>/evapotranspiration.nc
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-17 Jizixing Zhu -- distinguish rates [cm s-1] from
%         accepted-step amounts [mm] in names and metadata.
% Change: 2026-07-17 Jizixing Zhu -- distinguish accepted realized soil
%         evaporation from the upstream SCOPE atmospheric demand.
% Change: 2026-07-18 Jizixing Zhu -- separate positive ecological ET from
%         signed net vapor exchange and expose dew/deposition explicitly.
    g = 'evapotranspiration';
    registerRateAmount(reg, g, 'soil_evaporation', 'Soil evaporation');
    registerRateAmount(reg, g, 'soil_evaporation_demand', ...
        'Upstream SCOPE soil-evaporation demand');
    registerRateAmount(reg, g, 'soil_evaporation_deficit', ...
        'Atmospheric soil-evaporation demand not realized by the active set');
    registerRateAmount(reg, g, 'transpiration', 'Root water uptake exported to transpiration');
    registerRateAmount(reg, g, 'bare_surface_dew', ...
        'Bare-surface condensation routed as liquid supply');
    registerRateAmount(reg, g, 'canopy_dew', ...
        'Canopy condensation retained by the canopy water store');
    registerRateAmount(reg, g, 'snow_sublimation', ...
        'Positive snow vapor loss to the atmosphere');
    registerRateAmount(reg, g, 'snow_deposition', ...
        'Atmospheric vapor deposition into the snow store');
    registerRateAmount(reg, g, 'evapotranspiration_total', ...
        'Positive water loss: soil, wet canopy, transpiration, and snow');
    registerRateAmount(reg, g, 'net_surface_vapor_exchange', ...
        'Positive total vapor loss minus dew and snow deposition');
end

function registerRateAmount(reg, group, prefix, longName)
% Register one rate and its accepted-step integral.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    reg.register(struct('name',[prefix '_flux'], 'long_name',[longName ' flux'], ...
        'units','cm s-1', 'dims',{{'time'}}, 'group',group, 'transform','scalar'));
    reg.register(struct('name',[prefix '_amount'], 'long_name',[longName ' amount'], ...
        'units','mm', 'dims',{{'time'}}, 'group',group, 'transform','scalar'));
end
