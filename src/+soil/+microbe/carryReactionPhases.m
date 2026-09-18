function state=carryReactionPhases(state,ledger,kind,bulkDensityBottomUp)
%CARRYREACTIONPHASES Carry MEND active/dormant switches with donor phase identity.
%   STATE is an accepted C (top-down mg/g) or N (bottom-up g/cm3) owner.
%   LEDGER contains accepted reaction substep before/after totals [mg/g] and
%   integrated active-to-dormant, dormant-to-active, and dormant-maintenance
%   amounts [mg/g]. Switching preserves mobile/attached/retained identity;
%   other within-group growth or loss retains the resulting phase shares.
%   No transport or attachment reaction is performed here.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-14 Codex -- invert density conversion before flipping N.
    suffix=''; factor=1;
    if strcmp(kind,'N')
        suffix='_N';
        bd=flip(bulkDensityBottomUp(:));
        factor=1000./bd;
    else
        assert(strcmp(kind,'C'),'soil:microbe:phaseKind', ...
            'Reaction phase kind must be C or N.');
    end
    names={'mobile','attached','retained'};
    first=['MB_' names{1} suffix];
    if ~isfield(state,first), return; end
    active=zeros(numel(state.(['MB' suffix])),3);
    dormant=active;
    for phase=1:3
        active(:,phase)=state.(['MB_' names{phase} suffix])(:);
        dormant(:,phase)=state.(['BD_' names{phase} suffix])(:);
    end
    if strcmp(kind,'N')
        active=flip(active).*factor;
        dormant=flip(dormant).*factor;
    end
    for j=1:numel(ledger)
        item=ledger{j};
        assert(all(abs(sum(active,2)-item.before.MB(:))<=1e-9.*max(1,item.before.MB(:))) ...
            &&all(abs(sum(dormant,2)-item.before.BD(:))<=1e-9.*max(1,item.before.BD(:))), ...
            'soil:microbe:phaseStart','Reaction phase totals do not match the accepted substep start.');
        a=item.a_to_d(:); d=item.d_to_a(:); loss=item.dormant_loss(:);
        assert(all([a;d;loss]>=-1e-12) ...
            &&all(a<=sum(active,2)+1e-9.*max(1,sum(active,2))), ...
            'soil:microbe:phaseSwitchOverdraw', ...
            'MEND active phase switch exceeds its donor.');
        fromActive=splitDonor(active,a);
        active=active-fromActive;
        dormant=dormant+fromActive;
        assert(all(d+loss<=sum(dormant,2)+1e-9.*max(1,sum(dormant,2))), ...
            'soil:microbe:phaseSwitchOverdraw', ...
            'MEND dormant switching or maintenance exceeds its donor.');
        fromDormant=splitDonor(dormant,d);
        dormant=dormant-fromDormant-splitDonor(dormant,loss);
        active=active+fromDormant;
        active=resizeGroup(active,item.after.MB(:),true);
        dormant=resizeGroup(dormant,item.after.BD(:),false);
    end
    if strcmp(kind,'N')
        active=flip(active./factor);
        dormant=flip(dormant./factor);
    end
    for phase=1:3
        state.(['MB_' names{phase} suffix])=active(:,phase);
        state.(['BD_' names{phase} suffix])=dormant(:,phase);
    end
end

function portion=splitDonor(group,amount)
    total=sum(group,2);
    share=zeros(size(total));
    positive=total>0;
    share(positive)=amount(positive)./total(positive);
    portion=group.*share;
end

function group=resizeGroup(group,target,allowNewGrowth)
    current=sum(group,2);
    positive=current>0;
    group(positive,:)=group(positive,:).*target(positive)./current(positive);
    empty=~positive&target>1e-12;
    assert(allowNewGrowth||~any(empty), ...
        'soil:microbe:phaseIdentity', ...
        'Dormant biomass appeared without a donor-phase switch.');
    group(empty,2)=target(empty);
    group(~positive&~empty,:)=0;
    assert(all(isfinite(group(:)))&&all(group(:)>=-1e-10), ...
        'soil:microbe:phaseState','Reaction phase state is invalid.');
end
