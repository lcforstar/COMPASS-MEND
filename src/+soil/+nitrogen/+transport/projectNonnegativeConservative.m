function [projected, diagnostic] = projectNonnegativeConservative( ...
        concentration, thetaElement, storageCoefficient, dz)
%PROJECTNONNEGATIVECONSERVATIVE  Remove negative ADE values without changing mass.
%
%   The retained STEMMUS-N finite-element solve can produce round-off-scale
%   negative nodal concentrations during strong seasonal depletion.  A plain
%   MAX(C,0) creates nitrogen.  This routine projects the raw nodal solution
%   onto the nonnegative set while preserving its column inventory under the
%   same trapezoidal storage weights used by CALCULATESTORAGE.
%
%   Inputs:
%     concentration       [NN x 1] [ug cm-3 water]
%     thetaElement        [NL x 2] [-] effective liquid-water content
%     storageCoefficient  [NL x 1] [-] adsorbed storage term; zero for NO3
%     dz                  [NL x 1] [cm]
%
%   Returns:
%     projected  [NN x 1] [ug cm-3 water] nonnegative concentration
%     diagnostic struct with raw/projected column mass [ug cm-2]
%
%   The operation is a numerical admissibility projection after the unchanged
%   ADE linear solve; it does not alter reaction or transport equations.
%
%   Editor: Codex; Date: 2026-08-30.
% Change: 2026-09-14 Codex -- reject material negative nodal inventory
%   instead of redistributing it through the whole column.

    c=concentration(:);
    dz=dz(:);
    storageCoefficient=storageCoefficient(:);
    nl=numel(dz);
    assert(numel(c)==nl+1&&isequal(size(thetaElement),[nl,2]) ...
        &&numel(storageCoefficient)==nl, ...
        'soil:nitrogen:projectionSize', ...
        'Concentration and element storage geometry are not aligned.');
    assert(all(isfinite(c))&&all(isfinite(thetaElement),'all') ...
        &&all(isfinite(storageCoefficient))&&all(isfinite(dz)) ...
        &&all(thetaElement>=0,'all')&&all(storageCoefficient>=0) ...
        &&all(dz>0),'soil:nitrogen:projectionDomain', ...
        'Conservative projection inputs must be finite and physically valid.');

    weights=zeros(nl+1,1);
    for element=1:nl
        weights(element)=weights(element)+0.5*dz(element) ...
            *(thetaElement(element,1)+storageCoefficient(element));
        weights(element+1)=weights(element+1)+0.5*dz(element) ...
            *(thetaElement(element,2)+storageCoefficient(element));
    end
    rawMass=weights.'*c;
    positive=max(c,0);
    positiveMass=weights.'*positive;
    negativeMass=weights.'*max(-c,0);
    roundoffConcentration=128*eps(max(1,max(abs(c))));
    assert(min(c)>=-roundoffConcentration, ...
        'soil:nitrogen:materialNegativeConcentration', ...
        ['The ADE solution contains a material negative node; reject this ' ...
         'step instead of redistributing nitrogen across the column.']);
    scaleMass=max(1,abs(weights).'*abs(c));
    tolerance=128*eps(scaleMass)*scaleMass;
    assert(rawMass>=-tolerance,'soil:nitrogen:negativeColumnInventory', ...
        ['The ADE solution has a negative total column inventory; a ' ...
         'nonnegative conservative projection is impossible.']);

    targetMass=max(rawMass,0);
    if positiveMass>0
        projected=positive.*(targetMass/positiveMass);
    else
        projected=zeros(size(c));
    end
    projectedMass=weights.'*projected;
    diagnostic=struct('raw_mass_ug_cm2',rawMass, ...
        'projected_mass_ug_cm2',projectedMass, ...
        'mass_residual_ug_cm2',projectedMass-rawMass, ...
        'minimum_raw_concentration',min(c), ...
        'negative_mass_ug_cm2',negativeMass, ...
        'relative_redistribution',negativeMass/max(positiveMass,realmin), ...
        'roundoff_concentration_limit',roundoffConcentration, ...
        'active',any(c<0));
end
