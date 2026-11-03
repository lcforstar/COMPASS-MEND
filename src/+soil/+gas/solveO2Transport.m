function [y_O2, C_O2_liq, f_aerobic, f_anaerobic, diagnostic] = ...
        solveO2Transport(y_O2_old, y_O2_reference, phaseOld, phaseNew, ...
        transport, R_O2_sink_reference, DeltZ, Delt_t, boundary)
%SOLVEO2TRANSPORT Advance one conservative soil-oxygen candidate.
%
%   This O2 adapter preserves the active production interface and physics.
%   The reusable conservative transport assembly is owned by
%   soil.gas.species.solveDryMoleFractionTransport; O2-specific diffusion, Henry partition,
%   and aerobic-response constants are supplied explicitly here.
%
%   Nodes are bottom-to-surface. Positive gas and liquid fluxes point upward.
%   y_O2 is dry-gas mole fraction [-], C_O2_liq is equilibrium
%   dissolved concentration [mol L^-1 water], and the aerobic/anaerobic
%   factors are [NL x 2]. DIAGNOSTIC is the exact discrete O2 ledger.
%
%   References: Millington and Quirk (1961), Massman (1998), Sander (2023).
%   Editor: Jizixing Zhu; Date: 2026-05-16.
%   Change: 2026-08-16 Jizixing Zhu -- use accepted/candidate phase storage,
%   shared gas dispersion and liquid flux, monotone conservative advection,
%   fixed-old-state Picard coupling, explicit boundary reactions, and an O2
%   conservation ledger; remove the positive diffusivity floor and state clamp.
%   Change: 2026-08-16 Jizixing Zhu -- use the Massman O2-in-air reference
%   diffusivity and its explicit local gas-pressure correction.
%   Change: 2026-08-20 Jizixing Zhu -- consume the typed applied species
%   boundary without changing the discrete transport equation.
%   Change: 2026-08-20 Jizixing Zhu -- apply zero free-gas bottom transport
%   with dissolved O2 free-drainage outflow as a monotone Robin row.
%   Change: 2026-08-21 Jizixing Zhu -- separate the species-independent
%   conservative transport kernel from the O2-specific parameter and
%   aerobic-response adapter without changing the active equation.
%   Change: 2026-08-21 Jizixing Zhu -- make dry-gas mole fraction the sole
%   prognostic O2 coordinate and use the converged shared carrier flux.

    properties = soil.gas.species.oxygenProperties();
    [y_O2, C_O2_liq, diagnostic] = ...
        soil.gas.species.solveDryMoleFractionTransport( ...
        y_O2_old, y_O2_reference, phaseOld, phaseNew, transport, ...
        R_O2_sink_reference, DeltZ, Delt_t, boundary, properties);

    f_aerobic_col = C_O2_liq ...
        ./ (C_O2_liq + properties.half_saturation_mol_L);
    f_aerobic = [f_aerobic_col, f_aerobic_col];
    f_anaerobic = 1 - f_aerobic;
end

% Historical concentration-coordinate implementation (superseded):
%SOLVEO2TRANSPORT Conservative positive gas/dissolved oxygen interval.
%   Bottom-up node C [mol/cm3 gas], endpoint liquid/gas fractions [NL,2],
%   requested layer sink [mol/cm3 bulk/s], velocity [cm/s] upward positive.
%   PHASE supplies theta_liquid_new, theta_liquid_old, theta_gas_old, and
%   temperature_old [NN,1]. Optional reserved_sink_mol_cm3_s gets priority.
%   Without PHASE only fixed, unfrozen phase tests are supported. Production
%   always supplies actual phase history. Top is Dirichlet; bottom defaults
%   to closed. Lumped Mnew*Cnew-Mold*Cold, existing molecular diffusion and
%   donor-upwind gas advection form a positive operator. An active set limits
%   demand and returns actual removal, unmet demand and boundary mass.
%   An infeasible reserved demand fails before commit: caller rebuilds C/N.
%   Source: retained COMPASS Henry/diffusion constants; repair contract at
%   docs/soil/biogeochemistry/OXYGEN_CONSERVATIVE_REPAIR_20260914.md.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-14 Codex -- conserve phase-resolved oxygen storage and
%   return applied consumption instead of clipping material negative states.
% Change: 2026-10-05 Wenjie Liu -- retain successful active-set solutions
%   exactly; if node-to-layer demand mapping produces a material negative,
%   preserve feasible reserved respiration and conservatively limit only
%   the remaining layer demand along its common nonnegative feasible ray.

% Change: 2026-10-23 Wenjie Liu -- restore the COMPASS-COMPASS_MATLAB
% oxygen transport/source boundary implementation; retain current carrier physics.
