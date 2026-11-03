function spec=organicPoolSpec()
%ORGANICPOOLSPEC Single CSV contract for cold and resolved organic priors.
% Fractions C/SOC and nfrac_* N/SOC are dimensionless mass ratios. Optional
% nc_* N/C ratios refer to the paired pool. Phase fractions refer to their
% parent C or N pool, not total SOC. AOM is outside SOC: element densities
% g element/cm3, surface stores g element/cm2 in the first (shallowest) row.
% Editor: Wenjie Liu; Date: 2026-10-17.
    spec.carbon={'POM1','POM2','MOC','QOC','MB','BD','ENZP1','ENZP2','EM','ENZN','IOM','DOC'};
    spec.nitrogen={'PON1','PON2','MON','QON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N','ENZN_N','DON'};
    spec.paired_carbon=spec.carbon([1:10,12]);
    spec.cold={'depth_cm','frac_POC','frac_MOC','frac_MB','frac_IOM'};
    spec.detailed=[{'depth_cm'},strcat('frac_',spec.carbon)];
    spec.phases={'mobile','attached','retained'};
end
