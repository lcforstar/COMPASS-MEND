function model = stageCanopyEnergy( ...
        model, kt, dt_outer_s, snow_surface_override, print_diagnostics)
%STAGECANOPYENERGY  Stage interception/snow, then run RTMo and energy balance.
%
%   Inputs:
%     model                 [inout] Model shell after stageCanopyInputs.
%     kt                    [-] One-based outer-step index.
%     dt_outer_s            [s] Actual outer-step duration.
%     snow_surface_override [struct|empty] Optional completed snow candidate
%                           for the single correction pass.
%     print_diagnostics     [logical] Print the ordinary-pass ebal diagnostic;
%                           false suppresses duplicate correction-pass output.
%
%   Return:
%     model [out] Model with interception, current snow-surface auxiliary,
%           optical radiation, and energy fields staged in workspace.
%
%   The explicit order is load-bearing: interception provides current ground
%   throughfall; a precipitation-only snow predictor updates the current-step
%   SCOPE surface; RTMo and ebal then consume the same snow view.  An optional
%   override is used by the single bounded melt/freeze correction pass.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- remove the one-hour snowfall lag by
%     staging interception and current-step snow before RTMo/ebal.
%   Change: 2026-07-17 Jizixing Zhu -- document predictor versus bounded-
%     correction inputs and remove the obsolete RTMo-before-interception text.

    if nargin < 4
        snow_surface_override = [];
    end
    if nargin < 5 || isempty(print_diagnostics)
        print_diagnostics = true;
    end
    model = vegetation.applyCanopyInterception(model, kt, dt_outer_s);
    built = model.built;
    scopeStep = built.workspace.scope_step;
    assert(isfield(scopeStep, 'Soil_bare'), ...
        'stageCanopyEnergy:missingBareSurface', ...
        'prepareScopeStep must provide Soil_bare before snow is applied.');
    if isempty(snow_surface_override)
        [snowSurface, phase, predictorFlux] = surface.snow.predictSurfaceState( ...
            built.Y.snow, built.parameters.snow, ...
            scopeStep.interception.throughfall_cms, ...
            scopeStep.meteo.Ta, scopeStep.meteo.RH, 0, dt_outer_s);
        correctionPass = false;
    else
        snowSurface = snow_surface_override;
        phase = surface.snow.partitionPrecipitation( ...
            scopeStep.interception.throughfall_cms, ...
            scopeStep.meteo.Ta, scopeStep.meteo.RH, ...
            built.parameters.snow);
        predictorFlux = struct();
        correctionPass = true;
    end
    scopeStep.Soil = surface.snow.applySurfaceProperties( ...
        scopeStep.Soil_bare, snowSurface, built.parameters.snow, ...
        scopeStep.scope_static.spectral);
    model.built.workspace.scope_step = scopeStep;
    model.built.workspace.snow = struct( ...
        'surface_state', snowSurface, ...
        'phase', phase, ...
        'predictor_flux', predictorFlux, ...
        'surface_correction_pass', correctionPass);
    model = surface.radiation.runStepRTMo(model);
    model = vegetation.energy.runStepEbal(model, kt, dt_outer_s);
    if print_diagnostics
        diagnostics.printStepEbalDiagnostics(model, kt);
    end
end
