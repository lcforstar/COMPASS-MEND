function partials = computePhaseStoragePartials(endpoint, state, switches, freezePartials, phys)
%COMPUTEPHASESTORAGEPARTIALS  Derive phase-volume partials for water storage.
%
%   partials = soil.water.computePhaseStoragePartials(endpoint, state,
%   switches, freezePartials, phys) converts local retention derivatives into
%   derivatives with respect to the solved signed pressure head and soil
%   temperature. The result is independent of the legacy DTheta_LLh/DTheta_UUh
%   naming and is suitable for unfrozen, non-isothermal, and freeze-thaw mass
%   storage assembly.
%
%   Inputs:
%     endpoint [in] struct of [NL x nD] phase values and local derivatives:
%                   .theta_l, .theta_u, .theta_i
%                   .dtheta_water_equiv_dbase_head [1/cm]
%                   .dtheta_liquid_dbase_head [1/cm]
%                   .dtheta_liquid_dfreeze_head [1/cm]
%                   .dtheta_ice_explicit_dT [1/K]
%     state    [in] struct of [NL x nD] raw variables:
%                   .h_raw [cm], .COR [-], .CORh [1/K]
%     switches [in] struct containing lFreezing and lSFCC
%     freezePartials [in] struct of [NL x nD] derivatives:
%                   .dh_frez_dh [cm/cm], .dh_frez_dT [cm/K]
%     phys     [in] physical constants with RHOL and RHOI [kg/m3]
%
%   Returns:
%     partials [out] struct of [NL x nD] derivatives:
%                   .dtheta_liquid_dh, .dtheta_liquid_dT
%                   .dtheta_water_equiv_dh, .dtheta_water_equiv_dT
%                   .dtheta_ice_dh, .dtheta_ice_dT
%                   .dtheta_gas_dh, .dtheta_gas_dT
%
%   Method:
%     The temperature-corrected base head is u=COR(T)*h, with
%     du/dh=COR and du/dT=-CORh*COR*h. Strict SFCC additionally uses the
%     implemented freezing head. In the unsaturated range, phase derivatives
%     close to total water. For h>=0, Ss*h is elastic storage outside the pore
%     phase partition, so liquid and ice derivatives close only to the capped
%     pore-water amount. The optional Dall'Amico phase relaxation remains
%     Picard-lagged: it changes phase values, while this Jacobian differentiates
%     the equilibrium phase closure used to build its target. Exact mixed-form
%     storage still uses the post-relaxation values.
%
%   References:
%     Zeng 2013, Coupled Dynamics in Soil, Eq. 4.17 and Eqs. 4.57-4.60.
%     Yu et al. 2018 JGR Atmospheres 123:7391, Appendix Eqs. A1-A5.
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-12 Jizixing Zhu -- separate positive-head elastic storage from
%           the liquid-ice pore-volume partition.

    if nargin < 5 || isempty(phys)
        phys = equations.physicalConstants();
    end

    hRaw = state.h_raw;
    COR = state.COR;
    Cw = state.CORh;
    duDh = COR;
    duDT = -Cw .* COR .* hRaw;

    cUbase = endpoint.dtheta_water_equiv_dbase_head;
    cLbase = endpoint.dtheta_liquid_dbase_head;
    cLfreeze = endpoint.dtheta_liquid_dfreeze_head;

    dUdh = cUbase .* duDh;
    dUdT = cUbase .* duDT;
    dPhaseDh = dUdh;
    dPhaseDT = dUdT;
    saturated = hRaw >= 0;
    dPhaseDh(saturated) = 0;
    dPhaseDT(saturated) = 0;
    freezing = isSwitchOn(switches, 'lFreezing');
    strictSFCC = freezing && isSwitchOn(switches, 'lSFCC');

    if ~freezing
        dLdh = dPhaseDh;
        dLdT = dPhaseDT;
        dIdh = zeros(size(dUdh));
        dIdT = zeros(size(dUdT));
    elseif strictSFCC
        dLdh = cLbase .* duDh + cLfreeze .* freezePartials.dh_frez_dh;
        dLdT = cLbase .* duDT + cLfreeze .* freezePartials.dh_frez_dT;
        rhoLOverI = phys.RHOL / phys.RHOI;
        dIdh = rhoLOverI .* (dPhaseDh - dLdh);
        dIdT = rhoLOverI .* (dPhaseDT - dLdT);
    else
        % The simplified branch prescribes ice explicitly as a smooth
        % temperature fraction and closes liquid water from total mass.
        dIdh = zeros(size(dUdh));
        dIdT = endpoint.dtheta_ice_explicit_dT;
        rhoIOverL = phys.RHOI / phys.RHOL;
        dLdh = dPhaseDh - rhoIOverL .* dIdh;
        dLdT = dPhaseDT - rhoIOverL .* dIdT;
    end

    dGdh = -dLdh - dIdh;
    dGdT = -dLdT - dIdT;
    requireFinite(dLdh, dLdT, dUdh, dUdT, dIdh, dIdT, dGdh, dGdT);

    partials = struct();
    partials.dtheta_liquid_dh = dLdh;
    partials.dtheta_liquid_dT = dLdT;
    partials.dtheta_water_equiv_dh = dUdh;
    partials.dtheta_water_equiv_dT = dUdT;
    partials.dtheta_ice_dh = dIdh;
    partials.dtheta_ice_dT = dIdT;
    partials.dtheta_gas_dh = dGdh;
    partials.dtheta_gas_dT = dGdT;
end

function requireFinite(varargin)
% Fail at the producer rather than hiding an invalid storage Jacobian.
    for i = 1:nargin
        assert(all(isfinite(varargin{i}(:))), ...
            'computePhaseStoragePartials:nonfinite', ...
            'All phase-storage partial derivatives must be finite.');
    end
end

function tf = isSwitchOn(switches, name)
% Read a logical legacy switch without inventing a fallback physical mode.
    tf = isfield(switches, name) && switches.(name) == 1;
end
