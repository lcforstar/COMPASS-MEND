function state = boundLeafNitrogen(state, lnc, slatop, p)
%BOUNDLEAFNITROGEN Bound functional capacity to accepted green-leaf N.
% lnc [g N/m2 green leaf], slatop [m2 green leaf/g leaf C] and p are
% immutable inputs. Returned capacities [umol/m2/s] are projected jointly:
% N_lc+N_resp+Vcmax25/Fc25+Jmax25/Fj25 <= available functional leaf N.
% A lower accepted N inventory takes precedence over slow acclimation. This
% changes enzyme-investment diagnostics, never the plant's physical N pools.
% Structural N follows computeLuna/CLM5 SNC closure; this is a COMPASS
% conservation constraint on the rate-limited LUNA state, not a wheat fit.
% Editor: Wenjie Liu; Date: 2026-09-30.
    assert(isscalar(lnc)&&isfinite(lnc)&&lnc>=0&&isscalar(slatop) ...
        &&isfinite(slatop)&&slatop>0,'LUNA:invalidLeafNitrogenBudget', ...
        'Leaf area N and SLA must be finite and physically valid.');
    n=numel(state.Vcmax25_z);
    if 0.9*lnc>p.SNC/slatop,available=lnc-p.SNC/slatop;
    else,available=0.1*lnc;end
    % Retain the established LUNA high-N boundary in tiny-area/seed states.
    % The excess physical N stays in the plant pool and is not invented as
    % chlorophyll/enzyme investment by the auxiliary adapter.
    available=min([lnc,p.FNCa_max,max(0,available)]);
    if ~isfield(state,'PNresp_z')||numel(state.PNresp_z)~=n,state.PNresp_z=zeros(n,1);end
    state.FNCa_z=ones(n,1)*available;
    state.N_functional_used_z=zeros(n,1);
    state.capacity_n_scale_z=ones(n,1);
    for z=1:n
        props=[state.PNlc_z(z),state.PNresp_z(z)];
        assert(all(isfinite(props))&&all(props>=0), ...
            'LUNA:invalidNitrogenInvestment','Functional N proportions are invalid.');
        if sum(props)>1,props=props/sum(props);end
        light=props(1)*available;resp=props(2)*available;
        vc=state.Vcmax25_z(z);jm=state.Jmax25_z(z);
        assert(isfinite(vc)&&vc>=0&&isfinite(jm)&&jm>=0, ...
            'LUNA:invalidFunctionalCapacity','Photosynthetic capacity is invalid.');
        enzyme=vc/p.Fc25+jm/p.Fj25;
        budget=max(0,available-light-resp);
        factor=1;
        if enzyme>budget&&enzyme>0,factor=budget/enzyme;end
        state.Vcmax25_z(z)=vc*factor;
        state.Jmax25_z(z)=jm*factor;
        state.capacity_n_scale_z(z)=factor;
        state.N_functional_used_z(z)=light+resp+enzyme*factor;
        if available>0
            state.PNlc_z(z)=light/available;
            state.PNresp_z(z)=resp/available;
            state.PNcb_z(z)=state.Vcmax25_z(z)/p.Fc25/available;
            state.PNet_z(z)=state.Jmax25_z(z)/p.Fj25/available;
            state.PNstore_z(z)=max(0,1-state.N_functional_used_z(z)/available);
        else
            state.PNlc_z(z)=0;state.PNresp_z(z)=0;
            state.PNcb_z(z)=0;state.PNet_z(z)=0;state.PNstore_z(z)=0;
        end
    end
end
