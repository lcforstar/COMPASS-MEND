function [nh4,no3]=availableDinMgG(nh4Node,no3Node,theta,rhokg,bd)
%AVAILABLEDINMGG Map accepted COMPASS DIN inventory to the microbial basis.
%   Inputs: NH4/NO3 [NN x 1, ug N/cm3 water], theta [NL x 2, cm3/cm3],
%   rhokg=rho*Kd [NL x 1, dimensionless], bd [NL x 1, g dry soil/cm3].
%   All inputs are bottom-up at the same accepted time t^n. Outputs are
%   top-down [NL x 1, mg N/g dry soil]. NH4 includes equilibrium sorption;
%   NO3 includes liquid only. Ice, granules and new fertilizer are excluded.
%   Source: soil.nitrogen.transport.calculateStorage, endpoint trapezoids.
%   Average(theta*C), not average(theta)*average(C), preserves that contract.
%   This adapter does not consume or modify any state. Invalid inputs error.
% Editor: Jizixing Zhu; Date: 2026-09-05.
    bd=bd(:); rhokg=rhokg(:); nh4Node=nh4Node(:); no3Node=no3Node(:);
    n=numel(bd);
    assert(isequal(size(theta),[n,2]) && numel(rhokg)==n ...
        &&numel(nh4Node)==n+1 &&numel(no3Node)==n+1, ...
        'organic_n:dinShape','DIN requires NL-by-2 endpoint water and NN concentrations.');
    assert(all(isfinite([theta(:);rhokg;bd;nh4Node;no3Node])) ...
        &&all(theta(:)>=0 & theta(:)<=1) &&all(rhokg>=0) ...
        &&all(bd>0) &&all(nh4Node>=0) &&all(no3Node>=0), ...
        'organic_n:dinValues','DIN input states and capacities must be physical.');
    nh4Endpoints=[nh4Node(1:end-1),nh4Node(2:end)];
    no3Endpoints=[no3Node(1:end-1),no3Node(2:end)];
    nh4=flip(mean((theta+rhokg).*nh4Endpoints,2)./bd)*1e-3;
    no3=flip(mean(theta.*no3Endpoints,2)./bd)*1e-3;
end
