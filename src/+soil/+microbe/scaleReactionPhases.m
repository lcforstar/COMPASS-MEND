function state=scaleReactionPhases(state,kind)
%SCALEREACTIONPHASES Allocate within-group reaction changes, not switching.
% C fields use top-down mg/g; N fields use bottom-up g/cm3. Existing group
% phase shares are retained for growth/death/overflow. Empty groups seed
% attached phase as in the existing initialization closure. Call BEFORE
% switching, never to replace transfer of donor phase identity.
% Source: former bioparticleTransport reconciliation, moved to reaction owner.
% Editor: Jizixing Zhu; Date: 2026-09-11.
    suffix='';if strcmp(kind,'N'),suffix='_N';end
    for item={'MB','BD'}
        key=item{1}; aggregate=[key suffix];
        fields={[key '_mobile' suffix],[key '_attached' suffix],[key '_retained' suffix]};
        if ~isfield(state,fields{1}),continue;end
        total=state.(fields{1})+state.(fields{2})+state.(fields{3});
        positive=total>0;target=state.(aggregate);
        for j=1:3
            state.(fields{j})(positive)=state.(fields{j})(positive).*target(positive)./total(positive);
            state.(fields{j})(~positive)=0;
        end
        state.(fields{2})(~positive)=target(~positive);
    end
end
