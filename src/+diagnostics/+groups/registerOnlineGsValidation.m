function registerOnlineGsValidation(reg, cfg) %#ok<INUSD>
%REGISTERONLINEGSVALIDATION  Temporary online Medlyn and thermal validation diagnostics.
% Output: data/output_<tag>/online_gs_validation.nc
%
% [DIAG 2026-06-19] This group is enabled only by STEMMUS_ONLINE_GS_DIAG=1.
% Delete after the annual Medlyn/SOX/thermal-acclimation validation is
% archived. These diagnostics are write-only and must not feed model physics.

    g = 'online_gs_validation';

    addScalar(reg, g, 'medlyn_gs_model_sun',    'Sunlit leaf Medlyn gs returned by biochemical', 'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_model_sun');
    addScalar(reg, g, 'medlyn_gs_model_sha',    'Shaded leaf Medlyn gs returned by biochemical', 'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_model_sha');
    addScalar(reg, g, 'medlyn_gs_from_rcw_sun', 'Sunlit leaf gs reconstructed from rcw',         'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_from_rcw_sun');
    addScalar(reg, g, 'medlyn_gs_from_rcw_sha', 'Shaded leaf gs reconstructed from rcw',         'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_from_rcw_sha');
    addScalar(reg, g, 'medlyn_gs_residual_sun', 'Sunlit Medlyn identity residual',               'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_residual_sun');
    addScalar(reg, g, 'medlyn_gs_residual_sha', 'Shaded Medlyn identity residual',               'mol m-2 s-1', 'fluxes.online_gs_diag.medlyn_gs_residual_sha');
    addScalar(reg, g, 'medlyn_a_sun',           'Sunlit net assimilation used by Medlyn',         'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_a_sun');
    addScalar(reg, g, 'medlyn_a_sha',           'Shaded net assimilation used by Medlyn',         'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_a_sha');
    addScalar(reg, g, 'medlyn_ag_sun',          'Sunlit gross assimilation',                      'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_ag_sun');
    addScalar(reg, g, 'medlyn_ag_sha',          'Shaded gross assimilation',                      'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_ag_sha');
    addScalar(reg, g, 'medlyn_cs_sun',          'Sunlit leaf-surface CO2',                        'ppm',          'fluxes.online_gs_diag.medlyn_cs_sun');
    addScalar(reg, g, 'medlyn_cs_sha',          'Shaded leaf-surface CO2',                        'ppm',          'fluxes.online_gs_diag.medlyn_cs_sha');
    addScalar(reg, g, 'medlyn_ci_sun',          'Sunlit intercellular CO2',                       'ppm',          'fluxes.online_gs_diag.medlyn_ci_sun');
    addScalar(reg, g, 'medlyn_ci_sha',          'Shaded intercellular CO2',                       'ppm',          'fluxes.online_gs_diag.medlyn_ci_sha');
    addScalar(reg, g, 'medlyn_vpd_sun',         'Sunlit leaf-surface VPD used by Medlyn',         'kPa',          'fluxes.online_gs_diag.medlyn_vpd_sun');
    addScalar(reg, g, 'medlyn_vpd_sha',         'Shaded leaf-surface VPD used by Medlyn',         'kPa',          'fluxes.online_gs_diag.medlyn_vpd_sha');
    addScalar(reg, g, 'medlyn_tleaf_sun',       'Sunlit leaf temperature used by biochemical',    'degC',         'fluxes.online_gs_diag.medlyn_tleaf_sun');
    addScalar(reg, g, 'medlyn_tleaf_sha',       'Shaded leaf temperature used by biochemical',    'degC',         'fluxes.online_gs_diag.medlyn_tleaf_sha');
    addScalar(reg, g, 'medlyn_q_sun',           'Sunlit absorbed PAR used by biochemical',        'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_q_sun');
    addScalar(reg, g, 'medlyn_q_sha',           'Shaded absorbed PAR used by biochemical',        'umol m-2 s-1', 'fluxes.online_gs_diag.medlyn_q_sha');
    addScalar(reg, g, 'medlyn_sfactor_sun',     'Sunlit hydraulic stress factor applied to capacity', '-',         'fluxes.online_gs_diag.medlyn_sfactor_sun');
    addScalar(reg, g, 'medlyn_sfactor_sha',     'Shaded hydraulic stress factor applied to capacity', '-',         'fluxes.online_gs_diag.medlyn_sfactor_sha');
    addScalar(reg, g, 'medlyn_g0',              'Medlyn intercept used online',                   'mol m-2 s-1',  'fluxes.online_gs_diag.medlyn_g0');
    addScalar(reg, g, 'medlyn_g1',              'Medlyn slope used online',                       'kPa^0.5',      'fluxes.online_gs_diag.medlyn_g1');

    addScalar(reg, g, 'thermal_tgrowth',         'Growth temperature used for thermal acclimation', 'degC',          'fluxes.online_gs_diag.thermal_tgrowth');
    addScalar(reg, g, 'thermal_tgrowth_clamped', 'Clamped growth temperature used for thermal acclimation', 'degC',  'fluxes.online_gs_diag.thermal_tgrowth_clamped');
    addScalar(reg, g, 'thermal_delta_s_vcmax',   'Vcmax deactivation entropy',                     'J mol-1 K-1',   'fluxes.online_gs_diag.thermal_delta_s_vcmax');
    addScalar(reg, g, 'thermal_delta_s_jmax',    'Jmax deactivation entropy',                      'J mol-1 K-1',   'fluxes.online_gs_diag.thermal_delta_s_jmax');
    addScalar(reg, g, 'thermal_topt_vcmax',      'Diagnostic Vcmax optimum implied by current response', 'degC',     'fluxes.online_gs_diag.thermal_topt_vcmax');
    addScalar(reg, g, 'thermal_topt_jmax',       'Diagnostic Jmax optimum implied by current response',  'degC',     'fluxes.online_gs_diag.thermal_topt_jmax');
    addScalar(reg, g, 'thermal_vcmax25_sun',     'Sunlit current 25C Vcmax capacity',               'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_vcmax25_sun');
    addScalar(reg, g, 'thermal_vcmax25_sha',     'Shaded current 25C Vcmax capacity',               'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_vcmax25_sha');
    addScalar(reg, g, 'thermal_jmax25_sun',      'Sunlit current 25C Jmax capacity',                'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_jmax25_sun');
    addScalar(reg, g, 'thermal_jmax25_sha',      'Shaded current 25C Jmax capacity',                'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_jmax25_sha');
    addScalar(reg, g, 'thermal_rd25_sun',        'Sunlit current 25C leaf dark respiration capacity', 'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_rd25_sun');
    addScalar(reg, g, 'thermal_rd25_sha',        'Shaded current 25C leaf dark respiration capacity', 'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_rd25_sha');
    addScalar(reg, g, 'thermal_vcmax_sun',       'Sunlit temperature- and stress-adjusted Vcmax',   'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_vcmax_sun');
    addScalar(reg, g, 'thermal_vcmax_sha',       'Shaded temperature- and stress-adjusted Vcmax',   'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_vcmax_sha');
    addScalar(reg, g, 'thermal_jmax_sun',        'Sunlit temperature- and stress-adjusted Jmax',    'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_jmax_sun');
    addScalar(reg, g, 'thermal_jmax_sha',        'Shaded temperature- and stress-adjusted Jmax',    'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_jmax_sha');
    addScalar(reg, g, 'thermal_rd_sun',          'Sunlit temperature-adjusted leaf dark respiration', 'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_rd_sun');
    addScalar(reg, g, 'thermal_rd_sha',          'Shaded temperature-adjusted leaf dark respiration', 'umol m-2 s-1', 'fluxes.online_gs_diag.thermal_rd_sha');
end

function addScalar(reg, group_name, name, long_name, units, source)
    reg.register(struct('name', name, 'long_name', long_name, ...
        'units', units, 'dims', {{'time'}}, 'group', group_name, ...
        'source', source, 'transform', 'scalar', ...
        'comment', '[DIAG 2026-06-19] Temporary online validation output; delete after annual validation is archived.'));
end
