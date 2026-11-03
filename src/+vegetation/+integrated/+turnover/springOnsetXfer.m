function [p, diag] = springOnsetXfer(p, onset, fstor2tran, alloc, dt)
%SPRINGONSETXFER  Spring storage->xfer->display onset transfer (CLM5 perennial, faithful).
%
%   Faithful adaptation of the CLM5 FruitTree perennial onset (CNPhenologyMod.F90
%   :1946-1996): at budburst, a fraction fstor2tran of each organ's NSC storage is
%   mobilised into that organ's TRANSFER pool, which then drains into the displayed
%   organ pool evenly over the onset period (onset_counter = ndays_stor at budburst).
%   This bootstraps the canopy so an online leafc_sla LAI does not death-spiral from
%   a C_leaf=0 winter IC.
%
%   DIAGONAL DUAL-STORAGE (2026-06-07): CLM5 holds a SEPARATE storage pool per organ
%   (leafc_storage, frootc_storage, livecrootc_storage, livestemc_storage, ...), and
%   at onset each *_storage drains ONLY into the SAME organ's *_xfer -> display, never
%   into a combined reserve that is re-split (verified CNPhenologyMod.F90:923-941 +
%   CNCStateUpdate1Mod.F90:440-461: strictly diagonal plumbing). COMPASS collapses the
%   five CLM5 storage pools into TWO by spatial origin:
%     - C_root_storage  (belowground reserve, coarse+fine root) -> feeds ROOT only
%     - C_stem_storage  (aboveground woody reserve)             -> feeds LEAF + LIVESTEM
%   The leaf onset carbon is sourced from the aboveground (stem) reserve because the
%   spring leaf flush in deciduous trees is remobilised from proximal woody/branch NSC
%   via the phloem (Loescher et al. 1990 HortScience 25:274; Chapin et al. 1990 ARES
%   21:423); COMPASS has no separate leaf reserve pool. This keeps each reserve pool
%   self-contained: it fills from its own organ's allocation share (computeCarbonAllocation
%   Step 3c) and drains only to spatially-adjacent organs at onset, so the root reserve
%   can rebuild itself and is NEVER zero-ratcheted by a combined-reserve re-split (the
%   v59 bug: a combined total mobilised then re-split by the pre-mobilisation share fr_C,
%   which pinned C_root_storage at 0 once drained).
%
%   C and N are moved at the SAME fraction at every step (mobilise + drain), so each
%   transfer preserves the per-pool C:N and the plant C/N budget closes to machine
%   precision; the FlexibleCN clamp in the allocation step then pulls each organ toward
%   its C:N target band.
%
%   Returns:
%     p      [inout] struct, pools updated in place:
%              .C_root_storage/.N_root_storage   [g C|N/m²] belowground reserve (debited)
%              .C_stem_storage/.N_stem_storage   [g C|N/m²] aboveground reserve (debited)
%              .C_leaf/.N_leaf                   [g C|N/m²] displayed leaf (credited)
%              .C_root/.N_root                   [g C|N/m²] displayed fine root (credited)
%              .C_livestem/.N_livestem           [g C|N/m²] displayed live stem (credited)
%              .C_leaf_xfer/.N_leaf_xfer         [g C|N/m²] leaf transfer pool
%              .C_root_xfer/.N_root_xfer         [g C|N/m²] fine-root transfer pool
%              .C_livestem_xfer/.N_livestem_xfer [g C|N/m²] live-stem transfer pool
%     diag   struct: .C_onset_mobilized .C_onset_to_display  [g C/m²/step]
%
%   Inputs:
%     onset      struct .mobilize [bool] one-shot, true on the budburst step
%                       .counter  [s]    remaining onset duration (>0 -> drain active)
%     fstor2tran [-]    fraction of storage mobilised at budburst (CLM5 0.5)
%     alloc      struct .leaf .root .stem [-] budburst allocation split (sum = 1);
%                       used to apportion the ABOVEGROUND (stem) reserve between the
%                       leaf and livestem xfer pools (re-normalised over leaf+stem).
%     dt         [s]    outer timestep
%
%   Reference:
%     CLM5 FruitTree CNPhenologyMod.F90:1946-1996 (perennial onset: *_storage_to_xfer
%       = fstor2tran*organ_storage, per organ, C and N; onset_counter = ndays_stor*secspday),
%       :1890-1903 (xfer drain over onset_counter)
%       (D:/OneDrive/CTSM-release-clm5.0_FruitTree/src/biogeochem/)
%     CLM5 CNPhenologyMod.F90:923-941 (season-decid per-organ storage->xfer),
%       CNCStateUpdate1Mod.F90:440-461 (diagonal storage->xfer->display state update)
%       (D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/)
%     Lawrence DM et al. (2019) JAMES 11:4245 (CLM5 Tech Note, deciduous phenology)
%     Loescher WH et al. (1990) HortScience 25:274 (deciduous tree spring remob from woody NSC)
%
%   Author: Zhu Jizixing 2026-06-05; diagonal dual-storage rewrite 2026-06-07

    diag = struct('C_onset_mobilized', 0, 'C_onset_to_display', 0);

    % --- 1. Budburst: mobilise each reserve -> its OWN organ transfer pools (one-shot) ---
    % Diagonal: root reserve -> root xfer; stem reserve -> leaf + livestem xfer.
    % C and N at the SAME fraction (fstor2tran) -> mobilised C:N = each reserve's C:N.
    if onset.mobilize
        % Belowground (root) reserve -> fine-root xfer only.
        dC_root = fstor2tran * p.C_root_storage;       % [g C/m²]
        dN_root = fstor2tran * p.N_root_storage;       % [g N/m²]
        p.C_root_xfer   = p.C_root_xfer   + dC_root;
        p.N_root_xfer   = p.N_root_xfer   + dN_root;
        p.C_root_storage = p.C_root_storage - dC_root;
        p.N_root_storage = p.N_root_storage - dN_root;

        % Aboveground (stem) reserve -> leaf + livestem xfer, split by re-normalised alloc.
        dC_stem = fstor2tran * p.C_stem_storage;       % [g C/m²]
        dN_stem = fstor2tran * p.N_stem_storage;       % [g N/m²]
        ag_sum  = alloc.leaf + alloc.stem;
        if ag_sum > 0
            fr_leaf = alloc.leaf / ag_sum;             % [-] aboveground reserve -> leaf share
        else
            fr_leaf = 1;                               % degenerate: all to leaf bootstrap
        end
        p.C_leaf_xfer     = p.C_leaf_xfer     + fr_leaf       * dC_stem;
        p.C_livestem_xfer = p.C_livestem_xfer + (1 - fr_leaf) * dC_stem;
        p.N_leaf_xfer     = p.N_leaf_xfer     + fr_leaf       * dN_stem;
        p.N_livestem_xfer = p.N_livestem_xfer + (1 - fr_leaf) * dN_stem;
        p.C_stem_storage = p.C_stem_storage - dC_stem;
        p.N_stem_storage = p.N_stem_storage - dN_stem;

        diag.C_onset_mobilized = dC_root + dC_stem;
    end

    % --- 2. Onset drain: transfer -> display, even over the remaining counter ---
    % f = dt/counter empties each xfer pool exactly as counter -> 0 (CLM5 onset drain).
    if onset.counter > 0
        f = min(1, dt / onset.counter);           % [-] fraction of current xfer this step
        dCl = p.C_leaf_xfer * f;      p.C_leaf     = p.C_leaf     + dCl;  p.C_leaf_xfer     = p.C_leaf_xfer     - dCl;
        dNl = p.N_leaf_xfer * f;      p.N_leaf     = p.N_leaf     + dNl;  p.N_leaf_xfer     = p.N_leaf_xfer     - dNl;
        dCr = p.C_root_xfer * f;      p.C_root     = p.C_root     + dCr;  p.C_root_xfer     = p.C_root_xfer     - dCr;
        dNr = p.N_root_xfer * f;      p.N_root     = p.N_root     + dNr;  p.N_root_xfer     = p.N_root_xfer     - dNr;
        dCs = p.C_livestem_xfer * f;  p.C_livestem = p.C_livestem + dCs;  p.C_livestem_xfer = p.C_livestem_xfer - dCs;
        dNs = p.N_livestem_xfer * f;  p.N_livestem = p.N_livestem + dNs;  p.N_livestem_xfer = p.N_livestem_xfer - dNs;
        diag.C_onset_to_display = dCl + dCr + dCs;
    end
end
