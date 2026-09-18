function built = scatterNitrogenVariables(built, NV, kt)
%SCATTERNITROGENVARIABLES  Commit the mutated PlantCN slice into split containers.
%
%   built = init.scatterNitrogenVariables(built, NV, kt) writes the step-kt results of
%   the reused plant-C/N chain back into the new-chain containers. It is the commit half of the PR-D
%   bridge (init.assembleNitrogenVariables is the assemble half) and MUST be called ONLY
%   in the accepted-step branch of runOneStep: a rejected outer step must not advance
%   pools, accumulators, one-shot flags, or fire diagnostics. This mirrors the
%   phenNext/laiStateNext pending->commit discipline.
%
%   2-slot remap: the reused functions ran with a local kt2=2, so the last-write value of
%   every pool sits at slot 2 and every diagnostic row was written at slot 2; this maps
%   slot 2 back to the scalar Y.plant_cn pools and to the real-kt built.diagnostics rows.
%
%   Four commit categories (kept explicit + separate, per the PR-D spec):
%     1. pools       [g C or N / m2]  NV.<pool>(2) -> Y.plant_cn.{C,N,transfer}.* (last write)
%     2. accum       [scalar]         15 cum_* -> Y.plant_cn.accum.* (cum_U_N = plant-N inflow)
%     3. state       [scalar/struct]  6 one-shot flags / snapshots -> Y.plant_cn.state.*
%     4. diag rows   slot 2 -> kt     plantN_diag/nsc_diag/turnover_diag/frPHU_diag
%   Baselines (C/N_pool_init_const) are frozen constants and are NEVER written back.
%
%   Inputs:
%     built  struct   [inout] init container; .Y.plant_cn / .diagnostics updated
%     NV     struct   [in]    assembled and mutated PlantCN working frame
%     kt     [-]      [in]    real outer-step index (1-based); diag rows land at this row
%
%   Variable class (per compass-coding-standard §1):
%     NV -> reconstructed working frame (in);  Y.plant_cn pools/accum/state -> prognostic state (out);
%     diagnostics -> diagnostic (out, write-only).
%
%   Reference: COMPASS.m:1195-1421 (glue that produced these writes); docs/pr_d_plant_cn_spec.md
%     §2 (write map, accept-only). Companion: init.assembleNitrogenVariables.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-11 Jizixing Zhu -- remove Srt, SinkS, and root-uptake diagnostics; PlantCN commits only its own state.
    pc = built.Y.plant_cn;

    % same canonical remap as assembleNitrogenVariables (single source of truth for the pool set)
    cOrgans = {'labile','leaf','root_fine','livestem','deadstem','livecroot', ...
               'deadcroot','fruit','root_storage','stem_storage','xsmr'};   % 11 C pools
    nOrgans = cOrgans(1:end-1);                                             % 10 N pools
    xfer    = {'leaf','root','livestem'};                                   % 3 onset transfer pools

    % --- 1. pools: slot-2 (last write of the step) -> scalar Y.plant_cn pools ---
    for i = 1:numel(cOrgans)
        pc.C.(cOrgans{i}) = NV.(['C_' cOrgans{i}])(2);
    end
    for i = 1:numel(nOrgans)
        pc.N.(nOrgans{i}) = NV.(['N_' nOrgans{i}])(2);
    end
    for i = 1:numel(xfer)
        pc.transfer.C.(xfer{i}) = NV.(['C_' xfer{i} '_xfer'])(2);
        pc.transfer.N.(xfer{i}) = NV.(['N_' xfer{i} '_xfer'])(2);
    end

    % --- 2. accum: 15 cumulative flux accumulators (scalar) ---
    accF = fieldnames(pc.accum);
    for i = 1:numel(accF)
        pc.accum.(accF{i}) = NV.(accF{i});
    end

    % --- 3. state: 6 one-shot flags / year-start snapshots ---
    stF = fieldnames(pc.state);
    for i = 1:numel(stF)
        pc.state.(stF{i}) = NV.(stF{i});
    end

    built.Y.plant_cn = pc;   % baselines untouched (frozen; not part of pc mutation above)

    % --- 4. diagnostic rows: slot 2 -> real-kt built.diagnostics arrays ---
    built.diagnostics.plantN_diag(kt)   = NV.plantN_diag(2);
    built.diagnostics.nsc_diag(kt)      = NV.nsc_diag(2);
    built.diagnostics.turnover_diag(kt) = NV.turnover_diag(2);
    built.diagnostics.frPHU_diag(kt)    = NV.frPHU_diag(2);
end
