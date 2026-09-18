function [snit,sden,sink,ledger] = limitReactionRates(snit,sden,sink,available,bd,dt,reserved)
%LIMITREACTIONRATES Apply common per-layer NH4/NO3 inventory rationing.
%   All layer arrays are bottom-up. Rates use [ug N cm^-3 soil s^-1];
%   available and reserved use [mg N g^-1 dry soil] for the same step start.
%   Root reservations are committed first; only the remaining inventory is
%   available to chemistry. Chemistry is scaled proportionally when it
%   exceeds that remainder, so root uptake is never diluted by reactions.
%   Surface volatilization is intentionally excluded (area-flux contract).
%   Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-12 Codex — reserve root DIN before chemistry scaling.
% Change: 2026-09-14 Codex -- make the layer-order and time-level contract
%   explicit for the fixed post-event inventory passed by solveNitrogen.
    if nargin<7||isempty(reserved), reserved=zeros(size(available)); end
    assert(isequal(size(snit),size(sden))&&ndims(sink)==3 ...
        &&size(sink,3)==2&&isequal(size(available),size(snit)) ...
        &&isequal(size(reserved),size(available)) ...
        &&numel(bd)==size(snit,1)&&isscalar(dt)&&dt>0, ...
        'soil:nitrogen:reactionRateShape','Reaction rationing shapes are invalid.');
    bd=bd(:); factors=ones(size(available));
    requested=zeros(size(available)); actual=requested; unmet=requested;
    for sp=1:2
        if sp==1, chem=snit; else, chem=sden; end
        root=reserved(:,sp);
        assert(all(root<=available(:,sp)+1e-12), ...
            'soil:nitrogen:rootReservationOverdraw', ...
            'Root reservation exceeds available DIN.');
        chemDemand=mean(chem,2).*dt./bd*1e-3;
        demand=chemDemand+root;
        requested(:,sp)=demand;
        limited=chemDemand>max(available(:,sp)-root,0)&chemDemand>0;
        factors(limited,sp)=max(available(limited,sp)-root(limited),0)./chemDemand(limited);
        if sp==1, snit=snit.*factors(:,sp); else, sden=sden.*factors(:,sp); end
        actual(:,sp)=chemDemand.*factors(:,sp)+root;
        unmet(:,sp)=demand-actual(:,sp);
    end
    ledger=struct('requested_mg_g',requested,'actual_mg_g',actual, ...
        'unmet_mg_g',unmet,'scale',factors,'remaining_mg_g',available-actual, ...
        'reserved_mg_g',reserved, ...
        'grid_order','bottom_up','inventory_time_level','fixed_step_start', ...
        'policy','root_first_then_same_stage_same_species_proportional');
    assert(all(ledger.remaining_mg_g(:)>=-1e-12), ...
        'soil:nitrogen:reactionOverdraw','Reaction rationing exceeded DIN inventory.');
end
