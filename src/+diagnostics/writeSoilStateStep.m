function writeSoilStateStep(registry, built, kt)
%WRITESOILSTATESTEP  Export accepted water and frozen-soil state diagnostics.
%
%   COMPASS stores soil nodes from bottom to surface. Element row e joins
%   node e (column 1) to node e+1 (column 2). Registry transforms convert
%   these endpoint arrays to the public surface-to-bottom node order.
%
%   Effective saturation is evaluated with the van Genuchten definition:
%
%       Se = (theta_l - theta_r) / (theta_s - theta_r),  0 <= Se <= 1.
%
%   This function reads only the final refreshed cache associated with the
%   accepted h, U, and P_g state.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    if ~isfield(built.workspace, 'derived')
        return;
    end
    derived = built.workspace.derived;
    writeEndpoint(registry, 'soil_moisture', derived, 'Theta_UU', kt);
    writeEndpoint(registry, 'soil_unfrozen_water', derived, 'Theta_LL', kt);
    writeEndpoint(registry, 'soil_ice', derived, 'Theta_II', kt);
    writeEndpoint(registry, 'ice_ratio', derived, 'Ratio_ice', kt);
    writeEndpoint(registry, 'hydraulic_conductivity', derived, 'KL_h', kt);
    writeEndpoint(registry, 'frozen_hydraulic_conductivity', ...
        derived, 'KfL_h', kt);
    writeEndpoint(registry, 'soil_water_capacity', ...
        derived, 'DThetaLiquidDh', kt);

    if isfield(derived, 'Theta_LL') ...
            && isfield(built.parameters, 'soil')
        thetaResidual = built.parameters.soil.theta_r(:);
        thetaSaturated = built.parameters.soil.theta_s(:);
        denominator = max(thetaSaturated - thetaResidual, 1e-12);
        effectiveSaturation = (derived.Theta_LL - thetaResidual) ./ denominator;
        effectiveSaturation = max(0, min(1, effectiveSaturation));
        registry.write('effective_saturation_vg', effectiveSaturation, kt);
    end

    if isfield(built.Y, 'h_frez') && ~isempty(built.Y.h_frez)
        registry.write('frozen_matric_head', flip(built.Y.h_frez(:)), kt);
    elseif isfield(built.workspace.step, 'trial') ...
            && isfield(built.workspace.step.trial, 'h_frez')
        registry.write('frozen_matric_head', ...
            flip(built.workspace.step.trial.h_frez(:)), kt);
    end
end

function writeEndpoint(registry, outputName, source, sourceName, kt)
% Write one NL-by-2 bottom-up endpoint field when it is available.
    if isfield(source, sourceName) && ~isempty(source.(sourceName))
        registry.write(outputName, source.(sourceName), kt);
    end
end
