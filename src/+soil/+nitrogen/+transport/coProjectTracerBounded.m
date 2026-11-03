function [projected,diagnostic]=coProjectTracerBounded(rawTracer,projectionMap)
%COPROJECTTRACERBOUNDED Conservatively bound a source label by total N.
%
% RAWTRACER is a nodal passive-label solution [ug N cm-3 water].
% PROJECTIONMAP is the immutable map exported by the corresponding total-N
% solve.  This routine preserves the raw tracer column inventory under the
% total solve's exact trapezoidal weights while imposing
%       0 <= projected tracer <= projected total.
% If roundoff puts tracer inventory infinitesimally outside [0,total], the
% endpoint is used.  A material infeasibility rejects the candidate rather
% than creating/removing tagged nitrogen.
%
% The bounded redistribution uses the projected-total capacity.  Hence a
% zero tracer maps exactly to zero and a full-total tracer maps exactly to
% the total projected field, including when total-N positivity projection
% was active.
%
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-27 Wenjie Liu -- clamp the final roundoff correction back
%   to the closed bounds; report invariant magnitudes when a material
%   co-projection failure remains.
% Change: 2026-10-29 Wenjie Liu -- snap only exact zero/full inventories;
%   the feasibility tolerance must not erase a real tiny positive label or
%   promote a real unlabelled remainder to label on repeated accepted steps.

    rawTracer=rawTracer(:);
    required={'raw_total_ug_cm3','projected_total_ug_cm3', ...
        'storage_weights_cm','projected_total_mass_ug_cm2'};
    assert(isstruct(projectionMap)&&all(isfield(projectionMap,required)), ...
        'soil:nitrogen:tracerCoProjectionContract', ...
        'A complete frozen total projection map is required.');
    totalRaw=projectionMap.raw_total_ug_cm3(:);
    totalProjected=projectionMap.projected_total_ug_cm3(:);
    weights=projectionMap.storage_weights_cm(:);
    assert(numel(rawTracer)==numel(totalRaw)&&numel(rawTracer)==numel(weights) ...
            &&numel(totalProjected)==numel(rawTracer), ...
        'soil:nitrogen:tracerCoProjectionSize', ...
        'Tracer and total projection-map nodes are not aligned.');
    assert(all(isfinite(rawTracer))&&all(isfinite(totalRaw)) ...
            &&all(isfinite(totalProjected))&&all(isfinite(weights)) ...
            &&all(weights>=0)&&all(totalProjected>=0), ...
        'soil:nitrogen:tracerCoProjectionDomain', ...
        'Tracer co-projection inputs must be finite and physically valid.');

    rawMass=weights.'*rawTracer;
    totalMass=weights.'*totalProjected;
    scale=max(1,abs(weights).'*max(abs(rawTracer),abs(totalRaw)));
    tolerance=max(projectionMap.mass_tolerance_ug_cm2,512*eps(scale)*scale);
    assert(rawMass>=-tolerance&&rawMass<=totalMass+tolerance, ...
        'soil:nitrogen:tracerCoProjectionInfeasible', ...
        ['Raw tracer inventory lies materially outside [0,total]; a ' ...
         'bounded conservative co-projection is impossible.']);
    targetMass=min(max(rawMass,0),totalMass);

    if targetMass==0
        projected=zeros(size(rawTracer));
    elseif targetMass==totalMass
        projected=totalProjected;
        targetMass=totalMass;
    else
        projected=min(max(rawTracer,0),totalProjected);
        currentMass=weights.'*projected;
        if currentMass<targetMass
            capacity=max(totalProjected-projected,0);
            capacityMass=weights.'*capacity;
            assert(capacityMass>=targetMass-currentMass-tolerance, ...
                'soil:nitrogen:tracerCoProjectionCapacity', ...
                'Projected total has insufficient capacity for tracer mass.');
            if capacityMass>0
                projected=projected+capacity.*((targetMass-currentMass)/capacityMass);
            end
        elseif currentMass>targetMass
            if currentMass>0
                projected=projected.*(targetMass/currentMass);
            end
        end
        projected=min(max(projected,0),totalProjected);
        projected=correctRoundoffMass(projected,totalProjected,weights,targetMass,tolerance);
        % The mass correction is bounded analytically, but the final division
        % and addition can round a zero endpoint a few ulps outside [0,total].
        % Reapply the admissible interval; the resulting mass difference is
        % checked against the unchanged conservative tolerance below.
        projected=min(max(projected,0),totalProjected);
    end

    projectedMass=weights.'*projected;
    minimumProjected=min(projected);
    maximumExcess=max(projected-totalProjected);
    massResidual=projectedMass-targetMass;
    if any(projected<0)||any(projected>totalProjected) ...
            ||abs(massResidual)>tolerance
        error('soil:nitrogen:tracerCoProjectionFailure', ...
            ['Bounded tracer co-projection failed its mass/bounds invariant: ' ...
             'min=%.17g ug cm-3, max_excess=%.17g ug cm-3, ' ...
             'mass_residual=%.17g ug cm-2, tolerance=%.17g ug cm-2, ' ...
             'target_mass=%.17g ug cm-2, projected_mass=%.17g ug cm-2.'], ...
            minimumProjected,maximumExcess,massResidual,tolerance, ...
            targetMass,projectedMass);
    end
    diagnostic=struct( ...
        'method','bounded conservative co-projection using frozen total map', ...
        'raw_mass_ug_cm2',rawMass, ...
        'target_mass_ug_cm2',targetMass, ...
        'projected_mass_ug_cm2',projectedMass, ...
        'mass_residual_ug_cm2',projectedMass-targetMass, ...
        'minimum_raw_concentration',min(rawTracer), ...
        'maximum_raw_excess_over_total',max(rawTracer-totalRaw), ...
        'lower_bound_active',any(rawTracer<0), ...
        'upper_bound_active',any(rawTracer>totalProjected), ...
        'mass_tolerance_ug_cm2',tolerance);
end

function projected=correctRoundoffMass(projected,total,weights,target,tolerance)
% Put the final floating residual into one node with available weighted room.
    residual=target-weights.'*projected;
    if abs(residual)<=tolerance&&residual~=0
        if residual>0
            room=total-projected;
        else
            room=projected;
        end
        candidates=find(weights>0&room>0);
        for index=1:numel(candidates)
            node=candidates(index);
            availableMass=weights(node)*room(node);
            appliedMass=sign(residual)*min(abs(residual),availableMass);
            projected(node)=projected(node)+appliedMass/weights(node);
            residual=residual-appliedMass;
            if abs(residual)<=eps(max(1,abs(target)))
                break
            end
        end
    end
end
