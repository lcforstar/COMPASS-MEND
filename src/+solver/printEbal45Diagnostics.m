function printEbal45Diagnostics(KT, iter, canopy, meteo, rad, gap, thermal, fluxes, hydraulic, profiles)
%PRINTEBAL45DIAGNOSTICS Compact diagnostics for ebal_clm45 non-convergence.

trigger = false;
if isfield(iter, 'maxit') && isfield(iter, 'counter') && iter.counter >= iter.maxit
    trigger = true;
end
if isfield(iter, 'maxEBercu') && isfinite(iter.maxEBercu) && iter.maxEBercu > 1000
    trigger = true;
end
if isfield(thermal, 'Tsave') && isfinite(thermal.Tsave) && abs(thermal.Tsave) > 50
    trigger = true;
end
if isfield(thermal, 'Ts') && any(isfinite(thermal.Ts(:))) && ...
        (max(thermal.Ts(:), [], 'omitnan') > 50 || min(thermal.Ts(:), [], 'omitnan') < -80)
    trigger = true;
end
if isfield(profiles, 'Ta_layer') && any(isfinite(profiles.Ta_layer(:))) && ...
        (max(profiles.Ta_layer(:), [], 'omitnan') > 80 || min(profiles.Ta_layer(:), [], 'omitnan') < -80)
    trigger = true;
end
if ~trigger
    return;
end

fprintf('\n[EBAL45-DIAG] KT=%d\n', KT);
print_scalar('iter.counter', getfield_safe(iter, 'counter'));
print_scalar('iter.maxit', getfield_safe(iter, 'maxit'));
print_scalar('iter.maxEBercu', getfield_safe(iter, 'maxEBercu'));
print_scalar('iter.maxEBerch', getfield_safe(iter, 'maxEBerch'));
print_scalar('iter.maxEBers', getfield_safe(iter, 'maxEBers'));
print_scalar('canopy.LAI', getfield_safe(canopy, 'LAI'));
print_range('canopy.LAI_layer', getfield_safe(canopy, 'LAI_layer'));
print_scalar('meteo.Ta', getfield_safe(meteo, 'Ta'));
print_scalar('meteo.ea', getfield_safe(meteo, 'ea'));
print_scalar('meteo.u', getfield_safe(meteo, 'u'));
print_range('gap.Ps', getfield_safe(gap, 'Ps'));
print_range('rad.Rnuc', getfield_safe(rad, 'Rnuc'));
print_range('rad.Rnhc', getfield_safe(rad, 'Rnhc'));
print_range('thermal.Tcu', getfield_safe(thermal, 'Tcu'));
print_range('thermal.Tch', getfield_safe(thermal, 'Tch'));
print_range('thermal.Ts', getfield_safe(thermal, 'Ts'));
print_scalar('fluxes.Rnctot', getfield_safe(fluxes, 'Rnctot'));
print_scalar('fluxes.lEctot', getfield_safe(fluxes, 'lEctot'));
print_scalar('fluxes.Hctot', getfield_safe(fluxes, 'Hctot'));
print_scalar('fluxes.Tp_demand', getfield_safe(fluxes, 'Tp_demand'));
print_scalar('fluxes.Tp_t', getfield_safe(fluxes, 'Tp_t'));
print_scalar('fluxes.supply_ratio', getfield_safe(fluxes, 'supply_ratio'));
print_scalar('hydraulic.Tp_sun_demand', getfield_safe(hydraulic, 'Tp_sun_demand'));
print_scalar('hydraulic.Tp_sun', getfield_safe(hydraulic, 'Tp_sun'));
print_scalar('hydraulic.wt_sun', getfield_safe(hydraulic, 'wt_sun'));
print_scalar('hydraulic.psil_sun', getfield_safe(hydraulic, 'psil_sun'));
print_range('profiles.Tcu1d', getfield_safe(profiles, 'Tcu1d'));
print_range('profiles.lEcu1d', getfield_safe(profiles, 'lEcu1d'));
print_range('profiles.lEsun1d', getfield_safe(profiles, 'lEsun1d'));
print_range('profiles.Ta_layer', getfield_safe(profiles, 'Ta_layer'));
print_scalar('fluxes.Rnstot', getfield_safe(fluxes, 'Rnstot'));
print_scalar('fluxes.lEstot', getfield_safe(fluxes, 'lEstot'));
print_scalar('fluxes.Hstot', getfield_safe(fluxes, 'Hstot'));
print_scalar('fluxes.Gtot', getfield_safe(fluxes, 'Gtot'));
fprintf('[EBAL45-DIAG] end KT=%d\n\n', KT);
end

function value = getfield_safe(s, name)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = [];
end
end

function print_scalar(name, value)
if isempty(value)
    fprintf('  %-28s: <missing>\n', name);
elseif isnumeric(value) || islogical(value)
    fprintf('  %-28s: %.6g\n', name, double(value(1)));
else
    fprintf('  %-28s: <non-numeric>\n', name);
end
end

function print_range(name, value)
if isempty(value)
    fprintf('  %-28s: <missing>\n', name);
    return;
end
if ~(isnumeric(value) || islogical(value))
    fprintf('  %-28s: <non-numeric>\n', name);
    return;
end
v = double(value(:));
v = v(isfinite(v));
if isempty(v)
    fprintf('  %-28s: no finite values\n', name);
else
    fprintf('  %-28s: min %.6g  max %.6g  mean %.6g\n', name, min(v), max(v), mean(v));
end
end
