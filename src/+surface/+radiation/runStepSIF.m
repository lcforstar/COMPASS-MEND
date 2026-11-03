function model = runStepSIF(model)
%RUNSTEPSIF  Post-energy (or ebal-off) optical fluorescence and BRDF products.
%
%   Input:
%     model [inout] Model shell after RTMo (and ebal when calc_ebal=1).
%
%   Return:
%     model [out] Model with SIF / xanthophyll / Planck / BRDF auxiliaries
%           written to workspace.scope_step (diagnostic only).
%
%   Method:
%     calc_ebal=1 (default): RTMf, optional RTMz / RTMt_planck, optional
%       multi-angle BRDF (surface.radiation.calc_brdf), then Femtot — COMPASS.m case 1.
%     calc_ebal=0: fill canopy aPAR from RTMo gap statistics (legacy otherwise
%       branch), optional RTMf with unit eta profiles, then Femtot.
%
%   Reference:
%     COMPASS.m switch options.calc_ebal (case 1 + otherwise) and Femtot block;
%     van der Tol et al. (2009) SCOPE fluorescence / directional products.
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-16 Jizixing Zhu — migrate ebal-off aPAR branch and wire BRDF.
% Change: 2026-07-16 Jizixing Zhu — document that dispatcher and kernels now
%   contain no legacy soil bridge (superseded by explicit migrations).

    scope_static = model.built.p.scope_static;
    options = scope_static.options;
    scope_step = model.built.workspace.scope_step;
    spectral = scope_static.spectral;
    optical_parameters = scope_static.library.optipar;
    calc_ebal = optionOn(options, 'calc_ebal');
    % Default true when flag absent (historical SCOPE default).
    if ~isfield(options, 'calc_ebal')
        calc_ebal = true;
    end

    if calc_ebal
        model = runWithEbal(model, scope_step, scope_static, options, ...
            spectral, optical_parameters);
    else
        model = runWithoutEbal(model, scope_step, options, spectral, ...
            optical_parameters);
    end
end

function model = runWithEbal(model, scope_step, scope_static, options, ...
        spectral, optical_parameters)
% Full post-ebal radiation chain (COMPASS.m calc_ebal case 1).
    if optionOn(options, 'calc_fluor')
        if io.canopyLayerProfilesEnabled(options)
            [scope_step.rad, scope_step.profiles] = surface.radiation.RTMf( ...
                spectral, scope_step.rad, scope_step.Soil, ...
                scope_step.leafopt, scope_step.canopy, scope_step.gap, ...
                scope_step.angles, scope_step.profiles);
        else
            scope_step.rad = surface.radiation.RTMf( ...
                spectral, scope_step.rad, scope_step.Soil, ...
                scope_step.leafopt, scope_step.canopy, scope_step.gap, ...
                scope_step.angles, scope_step.profiles);
        end
    end
    if optionOn(options, 'calc_xanthophyll_abs') ...
            || optionOn(options, 'calc_xanthophyllabs')
        scope_step.rad = surface.radiation.RTMz( ...
            spectral, scope_step.rad, scope_step.Soil, scope_step.leafopt, ...
            scope_step.canopy, scope_step.gap, scope_step.angles, ...
            scope_step.profiles);
    end
    if optionOn(options, 'calc_planck')
        thermal = scope_step.thermal;
        scope_step.rad = surface.radiation.RTMt_planck( ...
            spectral, scope_step.rad, scope_step.Soil, scope_step.leafopt, ...
            scope_step.canopy, scope_step.gap, scope_step.angles, ...
            thermal.Tcu, thermal.Tch, thermal.Ts(2), thermal.Ts(1), 1);
    end
    if optionOn(options, 'calc_directional')
        scope_step.directional = runDirectional( ...
            scope_static, options, scope_step, spectral);
    end
    scope_step = appendFemtot(scope_step, options, spectral, optical_parameters);
    model.built.workspace.scope_step = scope_step;
end

function model = runWithoutEbal(model, scope_step, options, spectral, ...
        optical_parameters)
% ebal-off aPAR + optional fluorescence (COMPASS.m otherwise branch).
    canopy = scope_step.canopy;
    gap = scope_step.gap;
    rad = scope_step.rad;
    nl = canopy.nlayers;
    Fc = (1 - gap.Ps(1:end - 1))' / nl;
    lai = canopy.LAI;
    if ~isscalar(lai)
        lai = lai(1);
    end
    fluxes = readStruct(scope_step, 'fluxes');
    fluxes.aPAR = lai * (Fc * rad.Pnh + equations.meanleaf( ...
        canopy, rad.Pnu, 'angles_and_layers', gap.Ps));
    fluxes.aPAR_Cab = lai * (Fc * rad.Pnh_Cab + equations.meanleaf( ...
        canopy, rad.Pnu_Cab, 'angles_and_layers', gap.Ps));
    apar_wm2 = lai * (Fc * rad.Rnh_PAR + equations.meanleaf( ...
        canopy, rad.Rnu_PAR, 'angles_and_layers', gap.Ps));
    fluxes.aPAR_Wm2 = apar_wm2;
    fluxes.aPAR_Cab_eta = apar_wm2;
    scope_step.fluxes = fluxes;

    if optionOn(options, 'calc_fluor')
        profiles = readStruct(scope_step, 'profiles');
        profiles.etah = ones(nl, 1);
        profiles.etau = ones(13, 36, nl);
        scope_step.profiles = profiles;
        if io.canopyLayerProfilesEnabled(options)
            [scope_step.rad, scope_step.profiles] = surface.radiation.RTMf( ...
                spectral, scope_step.rad, scope_step.Soil, ...
                scope_step.leafopt, scope_step.canopy, scope_step.gap, ...
                scope_step.angles, scope_step.profiles);
        else
            scope_step.rad = surface.radiation.RTMf( ...
                spectral, scope_step.rad, scope_step.Soil, ...
                scope_step.leafopt, scope_step.canopy, scope_step.gap, ...
                scope_step.angles, scope_step.profiles);
        end
    end
    scope_step = appendFemtot(scope_step, options, spectral, optical_parameters);
    model.built.workspace.scope_step = scope_step;
end

function directional = runDirectional(scope_static, options, scope_step, spectral)
% Call package BRDF with the init-loaded angle template.
    if ~isfield(scope_static, 'directional') ...
            || isempty(scope_static.directional) ...
            || ~isfield(scope_static.directional, 'noa') ...
            || scope_static.directional.noa < 1
        error('runStepSIF:directionalTemplateMissing', ...
            ['calc_directional=true requires scope_static.directional ', ...
             '(tto/psi/noa from brdf_angles). Check init.buildScopeStatic.']);
    end
    template = scope_static.directional;
    % Fresh copy each step so oversampling does not grow the static template.
    directional0 = struct( ...
        'tto', template.tto(:), ...
        'psi', template.psi(:), ...
        'noa', template.noa);
    thermal = scope_step.thermal;
    if ~isfield(thermal, 'Tcu') || ~isfield(thermal, 'Tch') ...
            || ~isfield(thermal, 'Ts')
        error('runStepSIF:directionalNeedsThermal', ...
            'Directional BRDF after ebal needs thermal.Tcu/Tch/Ts.');
    end
    directional = surface.radiation.calc_brdf( ...
        options, directional0, spectral, scope_step.angles, scope_step.rad, ...
        scope_step.atmo, scope_step.Soil, scope_step.leafopt, ...
        scope_step.canopy, scope_step.meteo, scope_step.profiles, thermal);
end

function scope_step = appendFemtot(scope_step, options, spectral, optical_parameters)
% Total emitted fluorescence irradiance (legacy Femtot block after switch).
    if ~optionOn(options, 'calc_fluor')
        return;
    end
    fluorescence_yield = scope_step.leafbio.fqe;
    fluorescence_index = spectral.IwlF;
    fluxes = scope_step.fluxes;
    if optionOn(options, 'calc_PSI')
        scope_step.rad.Femtot = 1e3 * ( ...
            fluorescence_yield(2) ...
            * optical_parameters.phiII(fluorescence_index) ...
            * fluxes.aPAR_Cab_eta ...
            + fluorescence_yield(1) ...
            * optical_parameters.phiI(fluorescence_index) ...
            * fluxes.aPAR_Cab);
    else
        scope_step.rad.Femtot = 1e3 * fluorescence_yield ...
            * optical_parameters.phi(fluorescence_index) ...
            * fluxes.aPAR_Cab_eta;
    end
end

function tf = optionOn(options, name)
% True when options.(name) is present and nonzero.
    tf = isstruct(options) && isfield(options, name) ...
        && ~isempty(options.(name)) && logical(options.(name));
end

function value = readStruct(container, name)
% Read a struct field or return empty struct.
    value = struct();
    if isstruct(container) && isfield(container, name) ...
            && isstruct(container.(name))
        value = container.(name);
    end
end
