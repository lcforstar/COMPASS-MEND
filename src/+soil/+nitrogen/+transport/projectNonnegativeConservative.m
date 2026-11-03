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
%   ADE linear solve; it does not alter reaction or transport equations. A
%   material-but-small oscillation may only be balanced against positive mass
%   in the two adjacent nodes. Larger/nonlocal corrections are rejected.
%
%   Editor: Codex; Date: 2026-08-30.
% Change: 2026-09-14 Codex -- reject material negative nodal inventory
%   instead of redistributing it through the whole column.
% Change: 2026-09-22 Wenjie Liu -- admit a tightly bounded, nearest-neighbour
%   conservative repair for local ADE undershoot; retain hard rejection for
%   larger or nonlocal defects.

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
    [minimumConcentration,minimumNode]=min(c);
    maximumRelativeRepair=1e-7;
    relativeRepair=negativeMass/max(positiveMass,realmin);
    if minimumConcentration < -roundoffConcentration ...
            && relativeRepair>maximumRelativeRepair
        error('soil:nitrogen:materialNegativeConcentration', ...
            ['The ADE solution contains a material negative node; reject this ', ...
             'step instead of redistributing nitrogen across the column. ', ...
             'minimum=%.17g ug cm-3 water at bottom-up node %d; ', ...
             'negative mass=%.17g ug cm-2; raw column mass=%.17g ug cm-2; ', ...
             'relative repair=%.17g exceeds %.17g; ', ...
             'roundoff limit=%.17g ug cm-3.'], ...
            minimumConcentration,minimumNode,negativeMass,rawMass, ...
            relativeRepair,maximumRelativeRepair,roundoffConcentration);
    end
    scaleMass=max(1,abs(weights).'*abs(c));
    tolerance=128*eps(scaleMass)*scaleMass;
    assert(rawMass>=-tolerance,'soil:nitrogen:negativeColumnInventory', ...
        ['The ADE solution has a negative total column inventory; a ' ...
         'nonnegative conservative projection is impossible.']);

    targetMass=max(rawMass,0);
    [projected,maxDonorDistance,localRepairMass]=localConservativeRepair( ...
        c,weights,targetMass,2,tolerance);
    projectedMass=weights.'*projected;
    diagnostic=struct('raw_mass_ug_cm2',rawMass, ...
        'projected_mass_ug_cm2',projectedMass, ...
        'mass_residual_ug_cm2',projectedMass-rawMass, ...
        'minimum_raw_concentration',min(c), ...
        'negative_mass_ug_cm2',negativeMass, ...
        'relative_redistribution',relativeRepair, ...
        'maximum_relative_repair',maximumRelativeRepair, ...
        'local_repair_mass_ug_cm2',localRepairMass, ...
        'maximum_donor_distance_nodes',maxDonorDistance, ...
        'projection_method','nearest-neighbour conservative active set', ...
        'roundoff_concentration_limit',roundoffConcentration, ...
        'active',any(c<0));
end

function [projected,maxDistance,repairedMass]=localConservativeRepair( ...
        raw,weights,targetMass,maxRadius,tolerance)
% Zero negative nodes and remove exactly the added mass from nearest donors.
    projected=raw(:);
    negativeNodes=find(projected<0);
    repairedMass=0;maxDistance=0;
    for k=1:numel(negativeNodes)
        node=negativeNodes(k);
        deficit=weights(node)*(-projected(node));
        projected(node)=0;
        repairedMass=repairedMass+deficit;
        remaining=deficit;
        for radius=1:maxRadius
            candidates=unique([node-radius,node+radius]);
            candidates=candidates(candidates>=1&candidates<=numel(projected));
            for j=1:numel(candidates)
                donor=candidates(j);
                available=weights(donor)*max(projected(donor),0);
                removed=min(remaining,available);
                if removed>0
                    projected(donor)=projected(donor)-removed/weights(donor);
                    remaining=remaining-removed;
                    maxDistance=max(maxDistance,radius);
                end
                if remaining<=tolerance,break;end
            end
            if remaining<=tolerance,break;end
        end
        if remaining>tolerance
            error('soil:nitrogen:nonlocalPositivityRepair', ...
                ['A negative ADE node needs %.17g ug cm-2 beyond its two-node ', ...
                 'neighbourhood; reject instead of applying nonlocal redistribution.'],remaining);
        end
    end
    residual=targetMass-weights.'*projected;
    if abs(residual)<=tolerance&&residual~=0
        donors=find(weights>0&projected>0);
        if ~isempty(donors)
            [~,at]=min(abs(donors-negativeNodes(1)));
            donor=donors(at);
            projected(donor)=projected(donor)+residual/weights(donor);
        end
    end
    assert(all(projected>=0),'soil:nitrogen:localProjectionBounds', ...
        'Local conservative ADE repair left a negative concentration.');
end
