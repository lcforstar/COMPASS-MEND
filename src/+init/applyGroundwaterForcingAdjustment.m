function forcing = applyGroundwaterForcingAdjustment(forcing, gwDep, built)
%APPLYGROUNDWATERFORCINGADJUSTMENT  Dunnian-runoff forcing adapter for groundwater coupling.
%
%   forcing = init.applyGroundwaterForcingAdjustment(forcing, gwDep, built) applies
%   Dunnian saturation-excess runoff to the DAILY raw precipitation and, when that
%   changes precip, rebuilds the full source time/solar/hourly layers, aligns
%   their simulation-start suffix, then rebuilds derived F3/gas drivers.
%
%   Why a forcing adapter (not a beginRun raw edit or a per-step correction):
%   the new init path expands daily raw.Precip into the hourly driver series ONCE
%   at build time (init.buildForcing), so editing raw.Precip alone would never
%   reach the loop, and editing hourly.Precip (cm/s) would not match
%   updateDunnianRunoff's daily-mm interface. This adapter reproduces the legacy
%   COMPASS.m [9] pre-loop Dunnian call (COMPASS.m:325-326), which ran before the
%   in-loop daily->hourly downscale, by re-running F2/F3/gas on the corrected raw.
%
%   Inputs:
%     forcing  [inout] struct  built.forcing with complete raw source history;
%                              timeAxis/solar/hourly/derived/gas are rebuilt
%                              together when precipitation changes.
%     gwDep    [in]    [cm]    water-table depth from beginRun 3.3
%                              (groundwater.calculateGroundWaterDepth)
%     built    [in]    struct  control.time / domain.site / forcing.config for the
%                              F2/F3 rebuild (unchanged by this call)
%
%   Returns forcing with:
%     runoffDunnian  [mm]  saturation-excess runoff per raw step (diagnostic, only
%                          written here); raw.Precip and aligned drivers updated
%                          when gwDep <= 1 cm forces the all-or-nothing zeroing.
%
%   Reference: groundwater.updateDunnianRunoff (Mostafa et al. 2022 GMD preprint,
%   MODFLOW-STEMMUS BMI coupling) — full-saturation (gw_Dep <= 1 cm) all-or-nothing rule.
%
%   Editor: Jizixing Zhu; Date: 2026-07-01.
%   Change: 2026-10-28 Wenjie Liu -- rebuild complete source drivers and
%     reapply the shared time/solar/hourly alignment after raw rain changes.
    [runoffDunnian, correctedPrecip] = groundwater.updateDunnianRunoff(forcing.raw.Precip, gwDep);
    forcing.runoffDunnian = runoffDunnian;   % [mm] diagnostic, parallels raw.Precip

    % gw_Dep > 1 cm: not saturated -> precip unchanged -> hourly/derived/gas stay
    % valid; skip the rebuild so the common (unsaturated) case is a pure no-op.
    if isequaln(correctedPrecip, forcing.raw.Precip)
        return
    end

    % saturated: precip modified -> rebuild the dependent driver layers so the
    % zeroed precipitation reaches the loop (F2 hourly -> F3 derived -> gas).
    forcing.raw.Precip = correctedPrecip;
    forcing = init.buildForcingDrivers(forcing, built);
end
