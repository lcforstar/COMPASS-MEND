function C_labile_target_g = computeLabileSugarTarget(R_m_potential_kg, pftParams, dt)
%COMPUTELABILESUGARTARGET  Soluble-sugar (C_labile) working-pool floor for NSC allocation.
%
%   COMPASS treats C_labile as the soluble-sugar / immediately-mobile pool and
%   C_root_storage + C_stem_storage as the starch-like reserve. Without a floor,
%   the allocator sweeps all leftover C_labile into storage each step, leaving the
%   sugar pool at ~0 (unphysical: real trees keep soluble sugars above a critical
%   threshold while starch does the seasonal swinging). This function returns the
%   sugar floor the allocator must protect before sending excess carbon to storage.
%
%   Level 2 (respiration-buffer) closure, self-scaling across PFTs via the plant's
%   OWN potential maintenance respiration (no per-PFT absolute floor to calibrate):
%
%       C_labile_target = max(C_labile_abs_min, labile_resp_days * Rm_potential_day)
%
%   i.e. enough soluble sugar to fund 'labile_resp_days' of potential maintenance
%   respiration. Sized from the temperature/N-driven POTENTIAL R_m (before any
%   substrate down-regulation) so a low substrate signal cannot shrink its own
%   target. The Level-1 observation-constrained floor (measured soluble-sugar
%   concentration per organ) supersedes this when orchard sugar data are available.
%
%   Returns:
%     C_labile_target_g  [g C/m2]  soluble-sugar floor to protect this step
%
%   Inputs:
%     R_m_potential_kg   [kg C/(m2 s)]  potential maintenance respiration
%                                       (computeNPP R_m sum, BEFORE turn-down)
%     pftParams          struct, optional fields:
%                          .labile_resp_days  [d]      days of R_m held as sugar
%                                                      (default 1.0; calibration knob,
%                                                       not a measured constant)
%                          .C_labile_abs_min  [g C/m2] small nonzero floor (default 0.01)
%     dt                 [s]            time step (unused; kept for interface symmetry)
%
%   Method:
%     Rm_potential_day = R_m_potential_kg * 1000 * 86400        [g C/(m2 d)]
%     C_labile_target  = max(C_labile_abs_min, labile_resp_days * Rm_potential_day)
%
%   Physical basis (sugar buffered, starch swings):
%     Martinez-Vilalta et al. 2016 Ecol Monogr 86:495 Fig.8 + Discussion
%       (soluble-sugar depletion rare, starch depletion common; soluble sugars
%        kept above a critical threshold for immediate/osmotic functions)
%       (local PDF: data/references/fruit_tree/Martinez-Vilalta_et_al-2016-Ecological_Monographs.pdf)
%   Design note: docs/nsc_labile_sugar_floor_design.md (Level 1/2/3 hierarchy,
%     sensitivity set labile_resp_days = [0.25, 0.5, 1.0, 2.0]).
%   [SOURCE UNVERIFIED 2026-06-14] labile_resp_days is a parsimonious calibration
%     parameter, NOT a measured value; no published apple soluble-sugar residence
%     time located. Default 1.0 d pending the Level-1 observation-constrained floor.
%
%   Variable class (per compass-coding-standard §1):
%     R_m_potential_kg, pftParams, dt — forcing/parameters (in)
%     C_labile_target_g              — auxiliary (derived target, in)
%
%   Author: Zhu Jizixing 2026-06-14

    labile_resp_days = 1.0;     % [d]      default working-sugar residence (calibration knob)
    C_labile_abs_min = 0.01;    % [g C/m2] default small nonzero floor
    if isfield(pftParams, 'labile_resp_days') && ~isempty(pftParams.labile_resp_days)
        labile_resp_days = pftParams.labile_resp_days;
    end
    if isfield(pftParams, 'C_labile_abs_min') && ~isempty(pftParams.C_labile_abs_min)
        C_labile_abs_min = pftParams.C_labile_abs_min;
    end

    Rm_potential_day  = max(R_m_potential_kg, 0) * 1000 * 86400;        % [g C/(m2 d)]
    C_labile_target_g = max(C_labile_abs_min, labile_resp_days * Rm_potential_day);   % [g C/m2]
end
