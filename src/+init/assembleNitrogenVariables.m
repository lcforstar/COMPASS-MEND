function NV = assembleNitrogenVariables(built, kt)
%ASSEMBLENITROGENVARIABLES  Rebuild the legacy NitrogenVariables plant-C/N slice from split containers.
%
%   NV = init.assembleNitrogenVariables(built, kt) gathers every NitrogenVariables
%   field the reused per-KT plant-C/N chain consumes (runPlantNDemand ->
%   applyRootUptake glue -> runPlantCNStep -> computeNEE, legacy COMPASS.m:1195-1421)
%   from the new-chain split containers, so those functions run VERBATIM (golden-safe).
%   This is the assemble half of the PR-D bridge; init.scatterNitrogenVariables is the
%   commit half. It reads only, mutates nothing.
%
%   2-slot remap (the load-bearing trick): legacy pools are tS x 1 time series read at
%   idx_prev=max(1,KT-1) and written at KT with the REAL KT. New-chain Y.plant_cn pools
%   are SCALARS. Each pool is assembled as a 2-element column [prev; prev]; the step-4
%   wrappers call the reused functions with a LOCAL kt2=2, so idx_prev=max(1,1)=1 reads
%   the accepted scalar (slot 1) and writes land at slot 2. Because runPlantCNStep uses
%   the SAME kt for pools AND diagnostic rows, the diag scratch is ALSO 2-slot (seeded
%   from the still-unwritten template row built.diagnostics.<name>(kt)); scatter maps the
%   slot-2 diag back to the real-kt diagnostics array. Diagnostics stay in built.diagnostics
%   (not migrated into Y.plant_cn) -- only the prognostic pools are scalar.
%
%   Inputs:
%     built  struct  init container. Reads:
%                      .Y.plant_cn.{C,N,transfer,baseline,accum,state}  (scalar pools + state)
%                      .diagnostics.{plantN_diag,nsc_diag,turnover_diag,frPHU_diag} (templates)
%     kt     [-]      real outer-step index (1-based); picks the unwritten diag template row
%
%   Returns NV struct (legacy NitrogenVariables plant subset):
%     C/N/transfer pools  [2 x 1, g C or N / m2]  2-slot columns [prev; prev]
%     15 cum_* accumulators [scalar]              plant C/N cumulative fluxes (cum_U_N = N inflow)
%     6 state fields                              turnover_state / onset_mobilized_flag /
%                                                 resp_night_hours / C_deadstem_soy /
%                                                 N_deadstem_soy / pruned_this_year
%     C_pool_init_const / N_pool_init_const [scalar]  frozen conservation baselines (read-only)
%     plantN_diag/nsc_diag/turnover_diag    [2 x 1 struct]  diag scratch (write target = slot 2)
%     frPHU_diag                            [2 x 1]         frPHU snapshot scratch
%
%   Variable class (per compass-coding-standard §1):
%     Y.plant_cn pools/accum/state -> prognostic state (in);  baseline -> parameter (in, read-only);
%     diag templates -> diagnostic scratch (out, write-only downstream);  NV -> reconstructed
%     working frame handed to the reused physics functions.
%
%   Reference: COMPASS.m:1195-1421 (glue that read/writes these fields); setInitialNitrogen.m
%     (legacy field origins); docs/pr_d_plant_cn_spec.md §0 (2-slot decision) / §1 (read map).
%     Companion: init.scatterNitrogenVariables (commit), init.assembleSoilStateN (SoilStateN).
%
%   LEGACY ADAPTER LAYER (Julia-port): this + scatterNitrogenVariables + assembleSoilStateN +
%   buildDemandMeteoView exist ONLY to feed the reused verbatim legacy physics functions their old
%   struct shapes. They are a TRANSITIONAL boundary adapter, not new state. When the plant-C/N
%   physics is ported to Julia (dY, Y, p, aux, t), these assemble/scatter bridges DELETE -- the
%   tendency reads Y/p directly. Do not build new physics on top of this shape.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-04 Jizixing Zhu -- doc: mark as the transitional legacy adapter layer (Julia-port).
    pc = built.Y.plant_cn;

    % --- canonical pool remap (drop the C_/N_ prefix; transfer drops _xfer, root_fine->root) ---
    cOrgans = {'labile','leaf','root_fine','livestem','deadstem','livecroot', ...
               'deadcroot','fruit','root_storage','stem_storage','xsmr'};   % 11 C pools
    nOrgans = cOrgans(1:end-1);                                             % 10 N pools (no xsmr = C debt only)
    xfer    = {'leaf','root','livestem'};                                   % 3 onset transfer pools (C + N)

    NV = struct();

    % --- 2-slot prognostic pools: [prev; prev], slot 1 = accepted scalar, slot 2 = write target ---
    for i = 1:numel(cOrgans)
        v = pc.C.(cOrgans{i});               % [g C/m2] scalar accepted value
        NV.(['C_' cOrgans{i}]) = [v; v];     % [2 x 1]
    end
    for i = 1:numel(nOrgans)
        v = pc.N.(nOrgans{i});               % [g N/m2]
        NV.(['N_' nOrgans{i}]) = [v; v];     % [2 x 1]
    end
    for i = 1:numel(xfer)
        vc = pc.transfer.C.(xfer{i});        % [g C/m2] onset transfer pool
        vn = pc.transfer.N.(xfer{i});        % [g N/m2]
        NV.(['C_' xfer{i} '_xfer']) = [vc; vc];
        NV.(['N_' xfer{i} '_xfer']) = [vn; vn];
    end

    % --- 15 cumulative flux accumulators (scalar, prognostic; cum_U_N = plant-N inflow) ---
    accF = fieldnames(pc.accum);
    for i = 1:numel(accF)
        NV.(accF{i}) = pc.accum.(accF{i});
    end

    % --- 6 persistent one-shot flags + year-start snapshots (scalar/struct, prognostic) ---
    stF = fieldnames(pc.state);
    for i = 1:numel(stF)
        NV.(stF{i}) = pc.state.(stF{i});
    end

    % --- frozen conservation baselines (read-only; never scattered back) ---
    NV.C_pool_init_const = pc.baseline.C_pool_init_const;   % [g C/m2]
    NV.N_pool_init_const = pc.baseline.N_pool_init_const;   % [g N/m2]

    % --- 2-slot diagnostic scratch: seed both slots from the unwritten template row at kt ---
    %   (runPlantCNStep writes row kt2=2; scatter maps slot 2 -> built.diagnostics(kt).)
    NV.plantN_diag   = repmat(built.diagnostics.plantN_diag(kt),   2, 1);   % [2 x 1 struct]
    NV.nsc_diag      = repmat(built.diagnostics.nsc_diag(kt),      2, 1);   % [2 x 1 struct]
    NV.turnover_diag = repmat(built.diagnostics.turnover_diag(kt), 2, 1);   % [2 x 1 struct]
    NV.frPHU_diag    = repmat(built.diagnostics.frPHU_diag(kt),    2, 1);   % [2 x 1]
end
