function candidate=applyDormancy(accepted,parameters,dt_s)
%APPLYDORMANCY Apply MEND active/dormant switching and maintenance.
%
%   Substrate saturation omega=DOC/(K_d+DOC) drives active-to-dormant and
%   dormant-to-active fluxes. MEMC remains the sole owner of active growth,
%   active maintenance, active death, and enzyme production. This adapter owns
%   only state switching and dormant maintenance respiration. MEND has no
%   independent dormant-mortality flux, so none is introduced here.
%
%   Pools and integrated fluxes use [mg C g^-1 dry soil]. Rates are [d^-1].
% Editor: Jizixing Zhu; Date: 2026-08-28.

    required={'MB','BD','DOC','POC','IC'};
    for i=1:numel(required)
        assert(isfield(accepted,required{i})&&all(isfinite(accepted.(required{i}))) ...
            &&all(accepted.(required{i})>=0), ...
            'applyDormancy:state','MEND state %s is invalid.',required{i});
    end
    rateUnitSeconds=86400;
    if isfield(parameters,'rate_unit_seconds')
        assert(parameters.rate_unit_seconds==86400, ...
            'applyDormancy:badRateClock', ...
            'MEND/MEMC rates are day-based; rate_unit_seconds must equal 86400.');
        rateUnitSeconds=parameters.rate_unit_seconds;
    end
    state=accepted; elapsed=0; count=0;
    flux=struct('A_to_D',zeros(size(state.MB)), ...
        'D_to_A',zeros(size(state.MB)), ...
        'dormant_maintenance',zeros(size(state.MB)), ...
        'carbon_before',struct('MB',state.MB,'BD',state.BD));
    before=entityCarbon(state);
    while elapsed<dt_s-1e-12
        step=min(parameters.mend_substep_s,dt_s-elapsed);
        h=step/rateUnitSeconds;
        omega=state.DOC./(parameters.K_d+state.DOC);
        aToD=(1-omega).*parameters.m_r.*state.MB.*h;
        dToA=omega.*parameters.m_r.*state.BD.*h;
        maintenance=parameters.dormant_maintenance_ratio ...
            .*parameters.m_r.*state.BD.*h;
        activeOut=aToD; dormantOut=dToA+maintenance;
        activeScale=min(1,state.MB./max(activeOut,realmin));
        dormantScale=min(1,state.BD./max(dormantOut,realmin));
        aToD=aToD.*activeScale;
        dToA=dToA.*dormantScale;
        maintenance=maintenance.*dormantScale;
        % Change: 2026-09-11 Jizixing Zhu -- switching carries donor phase;
        % dormant respiration removes C from that same phase.
        if isfield(state,'MB_mobile')
            fractionA=zeros(size(state.MB));fractionD=zeros(size(state.BD));
            fractionR=zeros(size(state.BD));
            pa=state.MB>0;pd=state.BD>0;
            fractionA(pa)=aToD(pa)./state.MB(pa);
            fractionD(pd)=dToA(pd)./state.BD(pd);
            fractionR(pd)=maintenance(pd)./state.BD(pd);
            for phase={'mobile','attached','retained'}
                af=['MB_' phase{1}];df=['BD_' phase{1}];
                av=state.(af);dv=state.(df);
                state.(af)=av-av.*fractionA+dv.*fractionD;
                state.(df)=dv+av.*fractionA-dv.*(fractionD+fractionR);
            end
        end
        state.MB=state.MB-aToD+dToA;
        state.BD=state.BD+aToD-dToA-maintenance;
        state.IC=state.IC+maintenance;
        flux.A_to_D=flux.A_to_D+aToD;
        flux.D_to_A=flux.D_to_A+dToA;
        flux.dormant_maintenance=flux.dormant_maintenance+maintenance;
        elapsed=elapsed+step; count=count+1;
    end
    after=entityCarbon(state);
    residual=after-before;
    assert(max(abs(residual))<=1e-11.*max(1,max(before)), ...
        'applyDormancy:carbonClosure','MEND dormancy failed carbon closure.');
    candidate=struct('state',state,'integrated_flux_mg_g',flux, ...
        'substep_count',count,'elapsed_s',elapsed, ...
        'carbon_residual_mg_g',residual, ...
        'reaction_ownership','switching_dormant_maintenance');
end

function value=entityCarbon(state)
% Sum physical carbon plus cumulative IC; Tot is a source diagnostic only.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    value=state.POC+state.MOC+state.QOC+state.MB+state.BD ...
        +state.DOC+state.EP+state.EM+state.IC;
end
