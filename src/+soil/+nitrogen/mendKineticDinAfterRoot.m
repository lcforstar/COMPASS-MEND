function dinFree=mendKineticDinAfterRoot( ...
        remaining,reservation,nh4Node,no3Node,theta,rhokg,bd,S_min_nh4_mg_kg)
%MENDKINETICDINAFTERROOT Present dissolved DIN to the MEND kinetic formula.
%   The root-first shadow ledger REMAINING and RESERVATION use total NH4
%   (liquid plus equilibrium-sorbed) and liquid NO3 [mg N/g dry soil]. The
%   source MEND kinetic formula instead reads its free NH4 and NO3 pools.
%   Scale the post-root total NH4 by its accepted liquid/total fraction,
%   assuming the same instantaneous sorption equilibrium as COMPASS soil N.
%   Inventory allocation continues to use the unscaled REMAINING ledger.
%   Inputs are bottom-up accepted node states and layer parameters; output
%   is top-down NL-by-2 [mg N/g dry soil] of liquid NH4 and NO3.
%   Optional S_MIN_NH4_MG_KG is a finite nonnegative bottom-up NL-vector
%   [mg N/kg dry soil]. A nonzero threshold is subtracted from the total
%   post-root NH4 before applying the unchanged liquid/total fraction.
%   Empty or scalar/vector zero executes the original arithmetic. This
%   kinetic view never changes the physical post-root inventory ledger.
%   Reference: MOD_MEND.F90 commit 92323c7, Nim_NH4/Nim_NO3 use MNPOOL%NH4
%   and MNPOOL%NO3, whereas MNPOOL%NH4tot includes NH4ads.
%   Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-11-01 Wenjie Liu -- optionally protect total post-root NH4
% before presenting its dissolved fraction to MEND kinetic uptake.
    [totalNh4,~]=soil.organic_n.availableDinMgG( ...
        nh4Node,no3Node,theta,rhokg,bd);
    [liquidNh4,~]=soil.organic_n.availableDinMgG( ...
        nh4Node,no3Node,theta,zeros(size(rhokg)),bd);
    remainingNh4=remaining.nh4_mg_g(:);
    reservedNh4=reservation.nh4_mg_g(:);
    remainingNo3=remaining.no3_mg_g(:);
    assert(isequal(size(totalNh4),size(remainingNh4),size(reservedNh4), ...
        size(remainingNo3)) ...
        &&all(abs(totalNh4-remainingNh4-reservedNh4)<=1e-10) ...
        &&all(liquidNh4<=totalNh4+1e-12), ...
        'soil:nitrogen:mendKineticDinPhase', ...
        'The post-root DIN ledger does not match accepted NH4 phase storage.');
    liquidFraction=zeros(size(totalNh4));
    positive=totalNh4>0;
    liquidFraction(positive)=liquidNh4(positive)./totalNh4(positive);
    thresholdEnabled=false;
    if nargin>=8&&~isempty(S_min_nh4_mg_kg)
        assert(isnumeric(S_min_nh4_mg_kg)&&isreal(S_min_nh4_mg_kg) ...
            &&isvector(S_min_nh4_mg_kg) ...
            &&(numel(S_min_nh4_mg_kg)==numel(remainingNh4) ...
                ||(isscalar(S_min_nh4_mg_kg)&&S_min_nh4_mg_kg==0)) ...
            &&all(isfinite(S_min_nh4_mg_kg(:)))&&all(S_min_nh4_mg_kg(:)>=0), ...
            'soil:nitrogen:mendKineticProtection', ...
            'S_min_nh4_mg_kg must be a finite nonnegative NL-vector or zero scalar.');
        thresholdEnabled=any(S_min_nh4_mg_kg(:)>0);
    end
    if thresholdEnabled
        protectedMgG=1e-3*flip(S_min_nh4_mg_kg(:));
        accessibleNh4=max(remainingNh4-protectedMgG,0);
        dinFree=[accessibleNh4.*liquidFraction,remainingNo3];
    else
        dinFree=[remainingNh4.*liquidFraction,remainingNo3];
    end
    assert(all(isfinite(dinFree(:)))&&all(dinFree(:)>=0) ...
        &&all(dinFree(:,1)<=remainingNh4+1e-12), ...
        'soil:nitrogen:mendKineticDinPhase', ...
        'Free DIN substrate exceeded the post-root inventory.');
end
