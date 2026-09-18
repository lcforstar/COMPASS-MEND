function [carbon,organic,feedback]=applyNitrogenLimitation( ...
        carbon,organic,dinActualMgG,bulkDensityBottomUp,parameters,dt_s)
%APPLYNITROGENLIMITATION Apply MEND post-allocation microbial C:N correction.
%   CARBON is the top-down MEND candidate [mg C g-1 dry soil]; ORGANIC is
%   the bottom-up COMPASS organic-N candidate [g N cm-3 bulk soil].
%   DINACTUALMGG is the top-down COMPASS-applied microbial immobilization
%   [mg N g-1 dry soil]. MBD receives no new DIN. Carbon exceeding the
%   configured maximum microbial C:N is transferred from MBA/MBD to IC,
%   preserving carbon while exposing an overflow-CO2 result and limitation.
%   DIN uptake that would push MBA below its minimum C:N is returned to NH4.
%   Reference: MOD_MEND.F90 commit 92323c7 post-candidate C:N correction.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-13 Codex -- return post-uptake MBA nitrogen in excess
%   of the minimum C:N to NH4 through the accepted exchange ledger.
    bdTop=flip(bulkDensityBottomUp(:)); n=numel(bdTop);
    dinActualMgG=dinActualMgG(:);
    assert(numel(dinActualMgG)==n&&all(isfinite(dinActualMgG)) ...
        &&all(dinActualMgG>=0)&&isscalar(dt_s)&&isfinite(dt_s)&&dt_s>0, ...
        'soil:microbe:mendNitrogenFeedback', ...
        'Applied MEND DIN and step duration must be physical.');
    cnMin=parameters.CN_MB_min(:);
    if isscalar(cnMin), cnMin=repmat(cnMin,n,1); end
    cnMax=parameters.CN_MB_max(:);
    if isscalar(cnMax), cnMax=repmat(cnMax,n,1); end
    assert(numel(cnMin)==n&&numel(cnMax)==n ...
        &&all(isfinite(cnMin))&&all(isfinite(cnMax)) ...
        &&all(cnMin>0)&&all(cnMax>cnMin), ...
        'soil:microbe:mendNitrogenFeedback','MEND C:N bounds are invalid.');
    mbN=flip(organic.MB_N(:))./bdTop.*1000+dinActualMgG;
    bdN=flip(organic.BD_N(:))./bdTop.*1000;
    postUptakeRelease=max(mbN-carbon.MB(:)./cnMin,0);
    mbN=mbN-postUptakeRelease;
    overflowBD=max(carbon.BD-bdN.*cnMax,0);
    carbon.BD=carbon.BD-overflowBD;
    overflowMB=max(carbon.MB-mbN.*cnMax,0);
    carbon.MB=carbon.MB-overflowMB;
    overflow=overflowMB+overflowBD;
    carbon.IC=carbon.IC+overflow;
    carbon.POC=carbon.POM1+carbon.POM2;
    carbon.EP=carbon.ENZP1+carbon.ENZP2;
    carbon.Tot=carbon.POC+carbon.MOC+carbon.MB+carbon.BD+carbon.EP+carbon.EM;
    organic.MB_N=flip(mbN.*bdTop./1000);
    organic=soil.microbe.scaleReactionPhases(organic,'N');
    carbon=soil.microbe.scaleReactionPhases(carbon,'C');
    potential=max(carbon.MB+overflowMB,eps);
    feedback=struct('overflow_mba_mg_g',overflowMB, ...
        'overflow_mbd_mg_g',overflowBD,'overflow_co2_mg_g',overflow, ...
        'nh4_post_uptake_release_mg_g',postUptakeRelease, ...
        'growth_retained_fraction',carbon.MB./potential, ...
        'co2_c_rate_ug_cm3_s',flip(overflow.*bdTop.*1000./dt_s), ...
        'carbon_residual_mg_g',-overflowMB-overflowBD+overflow, ...
        'policy','MEND_MBD_then_MBA_CNmax_overflow');
    assert(max(abs(feedback.carbon_residual_mg_g))<=1e-12, ...
        'soil:microbe:mendNitrogenFeedbackClosure', ...
        'MEND nitrogen-limitation carbon correction failed closure.');
end
