function [state,diagnostic]=mixMicrobeCandidate(state,dzBottomUp,parameters,dt_s)
%MIXMICROBECANDIDATE Apply CoLM mixing to MEMC solid/microbial pools.
%
%   This is passive soil-profile redistribution, not a new MEMC reaction or an
%   active microbial motility equation. MEMC reaction equations remain frozen.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
    dz=flip(dzBottomUp(:))/100; source=zeros(size(dz));
    fields={'POC','MOC','QOC','MB','EP','EM'};
    changes=zeros(numel(fields),1);
    for i=1:numel(fields)
        name=fields{i}; before=sum(state.(name).*dz);
        state.(name)=soil.litter.colmVerticalTransport( ...
            state.(name),source,dz,dt_s,parameters.advection_m_s, ...
            parameters.diffusivity_m2_s);
        changes(i)=sum(state.(name).*dz)-before;
    end
    diagnostic=struct('fields',{fields},'storage_change',changes, ...
        'active_motility',false);
end
