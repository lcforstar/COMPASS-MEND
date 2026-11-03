function built = refreshDualPorositySink(built, dt_s)
%REFRESHDUALPOROSITYSINK Stage a dual-porosity candidate without changing Y.
%
%   built = soil.refreshDualPorositySink(built, dt_s) runs the HYDRUS Simunek (2003) dual-porosity
%   first-order exchange ONCE at the start of a trial substep and caches SinkIm_current for the
%   Richards RHS. Ports COMPASS.m:1646-1655. This is a substep-level PREPARE, NOT a Picard-internal
%   refresh: SinkIm_current is held CONSTANT across the Picard KIT loop so the mobile phase loses
%   exactly what the immobile phase gains (mass closure). A5 opens no Picard loop, is not placed in
%   a KIT loop, wires no runOneStep/attemptSolve, and calls no solve/refresh.
%
%   Inputs:
%     built  [inout] adapter container. Reads accepted mobile/immobile state
%                    and immutable parameters; writes only trial candidate/aux.
%     dt_s   [in]    [s] substep length. Passed through to the kernel (kernel applies its own
%                    <=0 default of 3600 s).
%
%   Returns:
%     built  with workspace.step.trial.soil_modes.dual_porosity candidate
%            and workspace.step.trial_aux.dual_porosity result. Accepted Y is
%            bitwise unchanged until solver.carrySubstepState accepts the tile.
%
%   Method: delegate to the value-semantic preferential-flow owner entrypoint.
%     The candidate and frozen exchange rate ride the normal solver snapshot;
%     rollback therefore discards both without restoring a prematurely edited Y.
%
%   Reference:
%     HYDRUS-1D WATFLOW.f DualPor line 1027-1075 (via soil.water.dualporosity.dualPorExchange).
%     Exchange coefficient comes from parameters.soil.modes.dual_porosity.
%   Oracle: tests/test_refreshDualPorositySink.m (legacy dualPorExchange called directly).
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-20 Jizixing Zhu -- replace the accepted-state write with a
%     soil-owned candidate and a trial-only Richards sink.
    sw = built.control.switches;
    if ~(isfield(sw, 'lDualPor') && sw.lDualPor == 1)
        return;   % dual-porosity off: no-op, Y.soil_modes untouched
    end

    accepted = built.Y.soil_modes.dual_porosity;
    p = built.parameters.soil;
    dpParameters = p.modes.dual_porosity;
    parameters = struct( ...
        'theta_s', p.theta_s, 'theta_r', p.theta_r, ...
        'theta_s_im', dpParameters.theta_s_im, ...
        'theta_r_im', dpParameters.theta_r_im, ...
        'OmegaExch', dpParameters.OmegaExch);
    [candidate, result] = ...
        soil.preferential_flow.stageDualPorosityCandidate( ...
            accepted, built.Y.Theta_L, parameters, dt_s);
    built.workspace.step.trial.soil_modes.dual_porosity = candidate;
    built.workspace.step.trial_aux.dual_porosity = result;
end
