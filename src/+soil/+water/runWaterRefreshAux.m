function built = runWaterRefreshAux(built, dt_s, L_f)
%RUNWATERREFRESHAUX  Stage-4 water-refresh entry point for the new chain.
%
%   built = soil.water.runWaterRefreshAux(built, dt_s, L_f) refreshes the soil-water
%   auxiliary cache (built.workspace.derived) by delegating to
%   soil.water.updateWaterContentAuxiliary -- the SAME refresh path initialization uses
%   (init.initializeSoilAux / init.refreshSoilAuxiliary). This is the stage-4 solver's
%   refresh boundary; it exists so the solve refreshes coefficients consistently with init.
%
%   Inputs:
%     built  [inout] canonical model container (built.Y, built.workspace, built.parameters,
%                    built.control)
%     dt_s   [in]    [s] outer update interval (0 during init / refresh-only use)
%     L_f    [in]    [J/kg] fusion latent heat (omit -> site freeze switch default)
%
%   Returns:
%     built  with refreshed built.workspace.derived (Theta_LL/UU/II, Se, KL_h, KfL_h,
%            KfL_T, DTheta_LLh/UUh, Ratio_ice, Gamma_hh, COR, CORh, ...)
%
%   RULE: new stage-4 code MUST NOT call the legacy UpdateSoilWaterContent -- that function
%   is a read-only oracle only (equivalence-checked by test_stageB_refresh_oracle). The
%   new refresh path is soil.water.updateWaterContentAuxiliary.
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-water-only-migration-oracle-design.md sec 3.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
    if nargin < 2, dt_s = 0; end
    if nargin < 3, L_f = []; end
    built = soil.water.updateWaterContentAuxiliary(built, dt_s, L_f);
end
