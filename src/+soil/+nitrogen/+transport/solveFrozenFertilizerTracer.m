function [candidate,ledger]=solveFrozenFertilizerTracer(frozen,request)
%SOLVEFROZENFERTILIZERTRACER Advance a one-way dissolved fertilizer-N label.
%
% The tracer reuses one converged total-N ADE operator exported by
% SOLVETRANSPORTSTEP.  It cannot alter R/Rold, dispersion, liquid flux,
% boundary type, accepted gross reaction rates, or any total state.
% Atmospheric and bottom-boundary tracer concentrations are identically
% zero: only fertilizer already present in REQUEST.old or tagged internal
% release pathways can introduce label.
%
% Inputs
%   frozen  struct returned by BUILDFROZENADEOPERATOR for this interval.
%   request struct with
%     .old.nh4_ug_cm3, .old.no3_ug_cm3       [NN x 1]
%     .rates.nh4.sources / .sinks            named [NL x 2] gross rates
%     .rates.no3.sources / .sinks            named [NL x 2] gross rates
%   Missing named rates are zero.  Supported names match
%   frozen.gross_rates.  NO3 nitrification is derived as 0.98 times tagged
%   NH4 nitrification; a supplied copy must agree exactly within roundoff.
%   Every tagged gross rate must be between zero and its accepted total
%   gross rate [ug N cm-3 soil s-1].
%
% Returns
%   candidate.nh4_ug_cm3 / no3_ug_cm3        bounded projected label
%   candidate.raw_*_ug_cm3                   raw frozen-ADE label solution
%   ledger                                    per-path amounts [ug N cm-2],
%                                             boundary fluxes and closure
%
% The total positivity map is applied as a bounded conservative
% co-projection, guaranteeing 0 <= label <= total without changing either
% the total solution or tagged column inventory.  A materially infeasible
% request rejects the candidate.
%
% Editor: Wenjie Liu; Date: 2026-09-21.

    validateFrozen(frozen);
    [oldNH4,oldNO3,rates]=normalizeRequest(frozen,request);
    assertImplicitGrossContract(frozen);

    netNH4=sumRateGroup(rates.nh4.sources,frozen.NL) ...
        -sumRateGroup(rates.nh4.sinks,frozen.NL);
    netNO3=sumRateGroup(rates.no3.sources,frozen.NL) ...
        -sumRateGroup(rates.no3.sinks,frozen.NL);
    rhsNH4Pre=frozen.nh4.old_storage_lumped.*oldNH4./frozen.dt_s ...
        +assembleRateLoad(netNH4,frozen.dz_cm);
    rhsNO3Pre=frozen.no3.old_storage_lumped.*oldNO3./frozen.dt_s ...
        +assembleRateLoad(netNO3,frozen.dz_cm);
    rhsNH4=applyZeroTracerBoundary(rhsNH4Pre,frozen.boundary, ...
        frozen.ql_nodes_cm_s);
    rhsNO3=applyZeroTracerBoundary(rhsNO3Pre,frozen.boundary, ...
        frozen.ql_nodes_cm_s);

    [rawNH4,rawNO3]=soil.nitrogen.transport.solveSoluteLinearSystem( ...
        frozen.no3.matrix_diag_upper,frozen.nh4.matrix_diag_upper, ...
        oldNH4,oldNO3,frozen.NN,frozen.NL, ...
        frozen.no3.matrix_lower,frozen.nh4.matrix_lower,rhsNH4,rhsNO3);
    [projectedNH4,projectionNH4]=soil.nitrogen.transport.coProjectTracerBounded( ...
        rawNH4,frozen.nh4.projection_map);
    [projectedNO3,projectionNO3]=soil.nitrogen.transport.coProjectTracerBounded( ...
        rawNO3,frozen.no3.projection_map);

    % Boundary exchange belongs to the frozen ADE solve. The subsequent
    % mass-conservative admissibility projection is an internal redistribution
    % and must not manufacture a different external flux when co-projection
    % moves label capacity between interior and boundary nodes.
    boundaryFlux=realizedTracerBoundaryFlux(frozen,rhsNH4Pre,rhsNO3Pre, ...
        rawNH4,rawNO3);
    candidate=struct( ...
        'schema_version',1, ...
        'nh4_ug_cm3',projectedNH4, ...
        'no3_ug_cm3',projectedNO3, ...
        'raw_nh4_ug_cm3',rawNH4, ...
        'raw_no3_ug_cm3',rawNO3, ...
        'projection_nh4',projectionNH4, ...
        'projection_no3',projectionNO3);
    ledger=buildPathLedger(frozen,oldNH4,oldNO3,projectedNH4,projectedNO3, ...
        rates,boundaryFlux);
end

function validateFrozen(frozen)
    required={'schema_version','dt_s','NL','NN','dz_cm','ql_nodes_cm_s', ...
        'boundary','nh4','no3','gross_rates'};
    assert(isstruct(frozen)&&all(isfield(frozen,required)) ...
            &&frozen.schema_version==1&&frozen.NN==frozen.NL+1, ...
        'soil:nitrogen:frozenTracerContract', ...
        'A compatible frozen ADE operator is required.');
end

function assertImplicitGrossContract(frozen)
% Current production chemistry exposes accepted reactions as explicit gross rates.
    coefficient=[frozen.nh4.implicit_reaction_coefficient_s_1(:); ...
        frozen.no3.implicit_reaction_coefficient_s_1(:)];
    tolerance=256*eps(max(1,max(abs(coefficient))));
    assert(all(abs(coefficient)<=tolerance), ...
        'soil:nitrogen:frozenTracerImplicitReaction', ...
        ['Fertilizer tracing requires accepted explicit gross rates; a ' ...
         'nonzero implicit reaction coefficient has no named path ledger.']);
end

function [oldNH4,oldNO3,rates]=normalizeRequest(frozen,request)
    assert(isstruct(request)&&isfield(request,'old')&&isstruct(request.old), ...
        'soil:nitrogen:frozenTracerRequest', ...
        'Tracer request.old is required.');
    oldNH4=readOld(request.old,'nh4_ug_cm3',frozen.nh4.old_total_ug_cm3);
    oldNO3=readOld(request.old,'no3_ug_cm3',frozen.no3.old_total_ug_cm3);
    if isfield(request,'rates')&&~isempty(request.rates)
        requestedRates=request.rates;
    else
        requestedRates=struct();
    end
    rates=struct();
    rates.nh4.sources=normalizeRateGroup(requestedRates,'nh4','sources', ...
        frozen.gross_rates.nh4.sources,frozen.NL);
    rates.nh4.sinks=normalizeRateGroup(requestedRates,'nh4','sinks', ...
        frozen.gross_rates.nh4.sinks,frozen.NL);
    rates.no3.sources=normalizeRateGroup(requestedRates,'no3','sources', ...
        frozen.gross_rates.no3.sources,frozen.NL);
    rates.no3.sinks=normalizeRateGroup(requestedRates,'no3','sinks', ...
        frozen.gross_rates.no3.sinks,frozen.NL);

    derivedNitrification=0.98*rates.nh4.sinks.nitrification;
    supplied=rates.no3.sources.nitrification;
    scale=max(1,max(abs([derivedNitrification(:);supplied(:)])));
    tolerance=512*eps(scale)*scale;
    suppliedExplicitly=hasNestedField(requestedRates,'no3','sources','nitrification');
    if suppliedExplicitly
        assert(all(abs(supplied-derivedNitrification)<=tolerance,'all'), ...
            'soil:nitrogen:frozenTracerNitrificationMismatch', ...
            'Tagged NO3 nitrification must equal 0.98 times tagged NH4 nitrification.');
    end
    rates.no3.sources.nitrification=derivedNitrification;
    validateRateBound(derivedNitrification, ...
        frozen.gross_rates.no3.sources.nitrification,'no3.sources.nitrification');
end

function value=readOld(owner,name,total)
    assert(isfield(owner,name), ...
        'soil:nitrogen:frozenTracerOldState', ...
        'Both old NH4 and NO3 tracer concentrations are required.');
    value=owner.(name)(:);
    total=total(:);
    scale=max(1,max(abs(total)));
    tolerance=512*eps(scale)*scale;
    assert(numel(value)==numel(total)&&all(isfinite(value)) ...
            &&all(value>=-tolerance)&&all(value<=total+tolerance), ...
        'soil:nitrogen:frozenTracerOldBounds', ...
        'Old tracer concentration must lie within [0,total] at every node.');
    value=max(value,0);
    value=min(value,total);
end

function group=normalizeRateGroup(requested,species,direction,totalGroup,nl)
    names=fieldnames(totalGroup);
    supplied=struct();
    if isfield(requested,species)&&isstruct(requested.(species)) ...
            &&isfield(requested.(species),direction) ...
            &&isstruct(requested.(species).(direction))
        supplied=requested.(species).(direction);
    end
    unknown=setdiff(fieldnames(supplied),names);
    assert(isempty(unknown), ...
        'soil:nitrogen:frozenTracerUnknownPath', ...
        'Tracer request contains an unknown gross-rate pathway.');
    group=struct();
    for index=1:numel(names)
        name=names{index};
        if isfield(supplied,name)
            value=supplied.(name);
        else
            value=zeros(nl,2);
        end
        assert(isequal(size(value),[nl,2])&&all(isfinite(value),'all'), ...
            'soil:nitrogen:frozenTracerRateShape', ...
            'Every tracer gross rate must be a finite [NL x 2] array.');
        validateRateBound(value,totalGroup.(name), ...
            [species '.' direction '.' name]);
        group.(name)=max(value,0);
    end
end

function validateRateBound(value,total,name)
    scale=max(1,max(abs(total),[],'all'));
    tolerance=512*eps(scale)*scale;
    assert(all(value>=-tolerance,'all')&&all(value<=total+tolerance,'all'), ...
        'soil:nitrogen:frozenTracerRateBounds', ...
        'Tagged gross rate %s must lie within its accepted total gross rate.',name);
end

function present=hasNestedField(owner,a,b,c)
    present=isfield(owner,a)&&isstruct(owner.(a)) ...
        &&isfield(owner.(a),b)&&isstruct(owner.(a).(b)) ...
        &&isfield(owner.(a).(b),c);
end

function total=sumRateGroup(group,nl)
    total=zeros(nl,2);
    names=fieldnames(group);
    for index=1:numel(names)
        total=total+group.(names{index});
    end
end

function load=assembleRateLoad(rate,dz)
% Same linear finite-element source quadrature as ASSEMBLEMATRIXCOEFFICIENTS.
    nl=numel(dz);
    load=zeros(nl+1,1);
    for element=1:nl
        load(element)=load(element) ...
            +(2*rate(element,1)+rate(element,2))*dz(element)/6;
        load(element+1)=load(element+1) ...
            +(rate(element,1)+2*rate(element,2))*dz(element)/6;
    end
end

function rhs=applyZeroTracerBoundary(rhs,boundary,qlNodes)
% Reproduce APPLYBOUNDARYCONDITIONS with every external tracer concentration zero.
    if boundary.NBCNB==1
        rhs(1)=0;
    elseif boundary.NBCNB==2
        if ~(qlNodes(1)>0||(boundary.lThmrlefc==1&&boundary.QMB_now==0))
            rhs(1)=0;
        end
    elseif boundary.NBCNB==3
        rhs(1)=0;
    else
        error('soil:nitrogen:frozenTracerBoundaryType', ...
            'Unsupported bottom nitrogen boundary type %g.',boundary.NBCNB);
    end
    if boundary.NBCN==1
        rhs(end)=0;
    elseif boundary.NBCN~=2
        error('soil:nitrogen:frozenTracerBoundaryType', ...
            'Unsupported top nitrogen boundary type %g.',boundary.NBCN);
    end
end

function boundaryFlux=realizedTracerBoundaryFlux(frozen,rhsNH4Pre,rhsNO3Pre,nh4,no3)
    saved=zeros(2,3,5);
    saved=fillSavedRows(saved,4,rhsNH4Pre, ...
        frozen.nh4.pre_boundary_matrix_diag_upper);
    saved=fillSavedRows(saved,5,rhsNO3Pre, ...
        frozen.no3.pre_boundary_matrix_diag_upper);
    zero=0;
    [bottomNH4,topNH4,bottomNO3,topNO3]= ...
        soil.nitrogen.transport.calculateBoundaryFlux(saved,nh4,no3, ...
        frozen.NN,1,zero,zero,zero,zero,frozen.ql_nodes_cm_s, ...
        zero,zero,zero,zero,frozen.boundary.NBCN,frozen.boundary.NBCNB, ...
        1e-30,frozen.boundary.lThmrlefc,frozen.boundary.QMB_now);
    boundaryFlux=struct( ...
        'nh4',soil.nitrogen.boundaryFluxLedger(topNH4,bottomNH4), ...
        'no3',soil.nitrogen.boundaryFluxLedger(topNO3,bottomNO3));
end

function saved=fillSavedRows(saved,speciesIndex,rhs,matrix)
    nn=numel(rhs);
    saved(1,1,speciesIndex)=rhs(1);
    saved(1,2,speciesIndex)=matrix(1,1);
    saved(1,3,speciesIndex)=matrix(1,2);
    saved(2,1,speciesIndex)=rhs(nn);
    saved(2,2,speciesIndex)=matrix(nn-1,2);
    saved(2,3,speciesIndex)=matrix(nn,1);
end

function ledger=buildPathLedger(frozen,oldNH4,oldNO3,newNH4,newNO3,rates,boundaryFlux)
    ledger=struct();
    ledger.schema_version=1;
    ledger.units=struct('rate','ug N cm-3 soil s-1', ...
        'boundary_flux','ug N cm-2 s-1','amount','ug N cm-2');
    ledger.external_tracer_boundary_concentration_ug_cm3=0;
    ledger.nh4=oneSoluteLedger(frozen,oldNH4,newNH4,rates.nh4, ...
        boundaryFlux.nh4,frozen.nh4.old_storage_lumped, ...
        frozen.nh4.projection_map.storage_weights_cm);
    ledger.no3=oneSoluteLedger(frozen,oldNO3,newNO3,rates.no3, ...
        boundaryFlux.no3,frozen.no3.old_storage_lumped, ...
        frozen.no3.projection_map.storage_weights_cm);
    ledger.combined_din=struct( ...
        'storage_change_ug_cm2',ledger.nh4.storage_change_ug_cm2 ...
            +ledger.no3.storage_change_ug_cm2, ...
        'expected_change_ug_cm2',ledger.nh4.expected_change_ug_cm2 ...
            +ledger.no3.expected_change_ug_cm2, ...
        'residual_ug_cm2',ledger.nh4.residual_ug_cm2 ...
            +ledger.no3.residual_ug_cm2, ...
        'nitrification_n2o_loss_ug_cm2', ...
            0.02*integrateRate(rates.nh4.sinks.nitrification, ...
            frozen.dz_cm)*frozen.dt_s);
end

function item=oneSoluteLedger(frozen,oldValue,newValue,rates,boundaryFlux,oldWeights,newWeights)
    sourceAmounts=integrateRateGroup(rates.sources,frozen.dz_cm,frozen.dt_s);
    sinkAmounts=integrateRateGroup(rates.sinks,frozen.dz_cm,frozen.dt_s);
    sourceTotal=sumScalarFields(sourceAmounts);
    sinkTotal=sumScalarFields(sinkAmounts);
    boundaryAmount=scaleBoundary(boundaryFlux,frozen.dt_s);
    initialMass=oldWeights(:).'*oldValue(:);
    finalMass=newWeights(:).'*newValue(:);
    expected=sourceTotal-sinkTotal+boundaryAmount.net_into_column;
    item=struct( ...
        'initial_storage_ug_cm2',initialMass, ...
        'final_storage_ug_cm2',finalMass, ...
        'storage_change_ug_cm2',finalMass-initialMass, ...
        'sources_ug_cm2',sourceAmounts, ...
        'sinks_ug_cm2',sinkAmounts, ...
        'boundary_flux',boundaryFlux, ...
        'boundary_amount_ug_cm2',boundaryAmount, ...
        'expected_change_ug_cm2',expected, ...
        'residual_ug_cm2',(finalMass-initialMass)-expected);
end

function amounts=integrateRateGroup(group,dz,dt)
    amounts=struct();
    names=fieldnames(group);
    for index=1:numel(names)
        name=names{index};
        amounts.(name)=integrateRate(group.(name),dz)*dt;
    end
end

function value=integrateRate(rate,dz)
    value=dz(:).'*mean(rate,2);
end

function total=sumScalarFields(owner)
    total=0;
    names=fieldnames(owner);
    for index=1:numel(names)
        total=total+owner.(names{index});
    end
end

function amount=scaleBoundary(flux,dt)
    amount=struct( ...
        'raw_top_upward',dt*flux.raw_top_upward, ...
        'raw_bottom_upward',dt*flux.raw_bottom_upward, ...
        'top_input',dt*flux.top_input,'top_loss',dt*flux.top_loss, ...
        'bottom_input',dt*flux.bottom_input,'bottom_loss',dt*flux.bottom_loss, ...
        'net_into_column',dt*flux.net_into_column);
end
