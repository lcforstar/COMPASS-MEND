function directional = calc_brdf(options, directional, spectral, angles, rad, atmo, soil, leafopt, canopy, meteo, profiles, thermal)
%CALC_BRDF  Multi-angle optical/thermal BRDF products after the canopy energy balance.
%
%   directional = surface.radiation.calc_brdf(...) evaluates SCOPE-style observations at
%   the user angle list plus hotspot and principal-plane oversampling angles.
%
%   Returns:
%     directional [struct]  Angle list and per-angle products:
%       .psi / .tto     [deg]  Observer azimuth / zenith (extended list)
%       .brdf_          [-]    Spectral BRDF (rso) per angle [nwlS x n_angles]
%       .LoF_           [...]  Directional fluorescence radiance (if calc_fluor)
%       .Lot_ / .Lot    [...]  Directional thermal radiance (Planck or SB path)
%       .BrightnessT    [K]    Brightness temperature (Stefan-Boltzmann path)
%
%   Inputs:
%     options      [in] struct  SCOPE flags (calc_planck / calc_fluor / ...)
%     directional  [in] struct  Template with .tto, .psi, .noa from brdf_angles
%     spectral     [in] struct  SCOPE spectral bands
%     angles       [in] struct  Solar geometry (.tts [deg] required)
%     rad          [in] struct  Current-step optical radiation (RTMo seed)
%     atmo         [in] struct  Atmosphere table for RTMo
%     soil         [in] struct  Soil optical surface
%     leafopt      [in] struct  Leaf optical properties
%     canopy       [in] struct  Canopy structure
%     meteo        [in] struct  Meteorology for RTMo
%     profiles     [in] struct  Canopy profiles (fluorescence / xanthophyll)
%     thermal      [in] struct  ebal temperatures (.Tcu/.Tch/.Ts)
%
%   Method:
%     For each observer angle, re-run optical transfer (RTMo) then optional
%     thermal (RTMt_planck or RTMt_sb), fluorescence (RTMf), and xanthophyll
%     (RTMz). Hotspot oversampling is added only when a direct solar beam is
%     present. Ported from STEMMUS-SCOPE calc_brdf (same signature as COMPASS.m).
%
%   Reference:
%     van der Tol et al. (2009) Biogeosciences SCOPE; STEMMUS-SCOPE src/calc_brdf.m;
%     COMPASS.m directional call after ebal.
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-16 Jizixing Zhu — migrate STEMMUS calc_brdf into +Radiation
%   package with surface.radiation.* kernel calls and no global constants.

% Change: 2026-07-17 Jizixing Zhu — omit hotspot oversampling without a
%   direct solar beam and reject generated below-horizon view angles.

    constants = io.define_constants();
    tts = angles.tts;
    noa = directional.noa;

    direct_beam_active = isfinite(tts) ...
        && cos(tts * constants.deg2rad) > 1e-6 ...
        && isfinite(meteo.Rin) && meteo.Rin > 1e-9;
    if direct_beam_active
        psi_hoversampling = [0; 0; 0; 0; 0; 2; 358];
        tto_hoversampling = [tts; tts + 2; tts + 4; tts - 2; tts - 4; tts; tts];
        valid_hotspot_angle = isfinite(tto_hoversampling) ...
            & tto_hoversampling >= 0 & tto_hoversampling < 90;
        psi_hoversampling = psi_hoversampling(valid_hotspot_angle);
        tto_hoversampling = tto_hoversampling(valid_hotspot_angle);
    else
        psi_hoversampling = zeros(0, 1);
        tto_hoversampling = zeros(0, 1);
    end
    noah_o = size(tto_hoversampling, 1);

    psi_poversampling = [000 * ones(6, 1); 180 * ones(6, 1); ...
        090 * ones(6, 1); 270 * ones(6, 1)];
    tto_poversampling = [10:10:60, 10:10:60, 10:10:60, 10:10:60]';
    noap_o = size(tto_poversampling, 1);

    directional.psi = [directional.psi(:); psi_hoversampling; psi_poversampling];
    directional.tto = [directional.tto(:); tto_hoversampling; tto_poversampling];
    n_angles = noa + noah_o + noap_o;

    directional.brdf_ = zeros(length(spectral.wlS), n_angles);
    directional.Eoutte = zeros(1, n_angles);
    directional.BrightnessT = zeros(1, n_angles);
    directional.LoF_ = zeros(length(spectral.wlF), n_angles);
    directional.Lot_ = zeros(length(spectral.wlT), n_angles);
    directional.Lot = zeros(1, n_angles);

    directional_angles = angles;
    use_xanthophyll = optionOn(options, 'calc_xanthophyll_abs') ...
        || optionOn(options, 'calc_xanthophyllabs');

    for j = 1:n_angles
        directional_angles.tto = directional.tto(j);
        directional_angles.psi = directional.psi(j);
        [directional_rad, directional_gap] = surface.radiation.RTMo( ...
            spectral, atmo, soil, leafopt, canopy, directional_angles, ...
            meteo, rad, options);
        directional.brdf_(:, j) = directional_rad.rso;

        if optionOn(options, 'calc_planck')
            directional_rad = surface.radiation.RTMt_planck( ...
                spectral, directional_rad, soil, leafopt, canopy, ...
                directional_gap, directional_angles, ...
                thermal.Tcu, thermal.Tch, thermal.Ts(1), thermal.Ts(1), 1);
            directional.Lot_(:, j) = directional_rad.Lot_(spectral.IwlT);
        else
            directional_rad = surface.radiation.RTMt_sb( ...
                spectral, directional_rad, soil, leafopt, canopy, ...
                directional_gap, directional_angles, ...
                thermal.Tcu, thermal.Tch, thermal.Ts(1), thermal.Ts(1), 1);
            if isfield(directional_rad, 'Eoutte')
                directional.Lot(j) = directional_rad.Eoutte;
                directional.Eoutte(j) = directional_rad.Eoutte;
            end
            % Legacy STEMMUS formula used the pre-loop rad.Lot seed.
            if isfield(rad, 'Lot') && ~isempty(rad.Lot) && isfinite(rad.Lot)
                directional.BrightnessT(j) = ...
                    (pi * rad.Lot / constants.sigmaSB) ^ 0.25;
            end
        end

        if optionOn(options, 'calc_fluor')
            [directional_rad, ~] = surface.radiation.RTMf( ...
                spectral, directional_rad, soil, leafopt, canopy, ...
                directional_gap, directional_angles, profiles);
            if isfield(directional_rad, 'LoF_')
                directional.LoF_(:, j) = directional_rad.LoF_;
            end
        end
        if use_xanthophyll
            directional_rad = surface.radiation.RTMz( ...
                spectral, directional_rad, soil, leafopt, canopy, ...
                directional_gap, directional_angles, profiles); %#ok<NASGU>
        end
    end
end

function tf = optionOn(options, name)
% True when options.(name) is present and nonzero.
    tf = isstruct(options) && isfield(options, name) ...
        && ~isempty(options.(name)) && logical(options.(name));
end
