function lai_state_next = accumulateLaiWaterStress(model, lai_state_next)
%ACCUMULATELAIWATERSTRESS  Update pending daily canopy water-stress for LAI.
%
%   Inputs:
%     model          [in] Model shell after canopy energy balance (has hydraulic
%                    wt_sun/wt_sha and preferably gap.Ps for sunlit fraction).
%     lai_state_next [inout] Pending heat-unit LAI state.
%
%   Return:
%     lai_state_next [out] Pending LAI state with wt_day_min updated.
%
%   Method:
%     Leaf-scale stress (wt_sun, wt_sha) is for photosynthesis only. Canopy LAI
%     is a canopy-scale state, so the hourly factor is area-weighted:
%       w_canopy = f_sun * wt_sun + (1 - f_sun) * wt_sha
%     with f_sun = sunlit leaf-area fraction from gap.Ps (LAI-layer weighted when
%     available). The running daily minimum of w_canopy is stored in wt_day_min
%     for the next daily heat-unit LAI advance. Non-model LAI modes leave state
%     unchanged.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu — use sun/shade area-weighted canopy stress
%     instead of min(wt_sun, wt_sha) for LAI aggregation.
    canopy_parameters = model.built.parameters.veg.canopy;
    if ~strcmpi(readField(canopy_parameters, 'LAI_mode', 'constant'), ...
            'model')
        return;
    end
    scope_step = model.built.workspace.scope_step;
    if ~isfield(scope_step, 'hydraulic') ...
            || ~isfield(scope_step.hydraulic, 'wt_sun') ...
            || ~isfield(scope_step.hydraulic, 'wt_sha')
        return;
    end

    wt_sun = double(scope_step.hydraulic.wt_sun);
    wt_sha = double(scope_step.hydraulic.wt_sha);
    if ~isscalar(wt_sun), wt_sun = mean(wt_sun(:), 'omitnan'); end
    if ~isscalar(wt_sha), wt_sha = mean(wt_sha(:), 'omitnan'); end
    if ~isfinite(wt_sun), wt_sun = 1; end
    if ~isfinite(wt_sha), wt_sha = 1; end
    wt_sun = min(1, max(0, wt_sun));
    wt_sha = min(1, max(0, wt_sha));

    f_sun = sunlitLeafAreaFraction(scope_step);
    w_canopy = f_sun * wt_sun + (1 - f_sun) * wt_sha;

    lai_state_next.wt_day_min = min( ...
        readField(lai_state_next, 'wt_day_min', 1.0), ...
        w_canopy);
end

function f_sun = sunlitLeafAreaFraction(scope_step)
% Return canopy sunlit leaf-area fraction in [0,1] for weighting wt_sun/wt_sha.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    f_sun = 0.5;
    if ~isfield(scope_step, 'gap') || ~isstruct(scope_step.gap) ...
            || ~isfield(scope_step.gap, 'Ps') || isempty(scope_step.gap.Ps)
        return;
    end
    Ps = double(scope_step.gap.Ps(:));
    Ps = Ps(isfinite(Ps));
    if isempty(Ps)
        return;
    end

    % Prefer LAI-layer weighted Ps (canopy layers only; exclude soil endpoint).
    if isfield(scope_step, 'canopy') && isstruct(scope_step.canopy) ...
            && isfield(scope_step.canopy, 'LAI_layer') ...
            && ~isempty(scope_step.canopy.LAI_layer)
        LAI_layer = double(scope_step.canopy.LAI_layer(:));
        nl = numel(LAI_layer);
        if nl >= 1 && numel(scope_step.gap.Ps) >= nl
            Ps_c = double(scope_step.gap.Ps(1:nl));
            Ps_c = min(1, max(0, Ps_c));
            LAI_layer = max(0, LAI_layer);
            wsum = sum(LAI_layer);
            if wsum > 1e-12
                f_sun = sum(LAI_layer .* Ps_c) / wsum;
                f_sun = min(1, max(0, f_sun));
                return;
            end
        end
    end

    % Fallback: mean of canopy Ps entries (drop last soil entry if present).
    if numel(Ps) >= 2
        Ps_c = Ps(1:end-1);
    else
        Ps_c = Ps;
    end
    Ps_c = min(1, max(0, Ps_c));
    f_sun = mean(Ps_c);
    f_sun = min(1, max(0, f_sun));
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
