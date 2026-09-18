function candidate = advanceAcceptedCandidate( ...
        accepted, theta, bulkDensity, dz_cm, q_cm_s, parameters, dt_s)
%ADVANCEACCEPTEDCANDIDATE Conservatively transport DOC and DON.
%
%   This DOM owner transports bottom-up element total DOC/DON [g/cm3 bulk
%   soil]. Linear-equilibrium states repartition by Kd; explicitly tagged
%   dynamic states retain their sorbed inventory and transport only free mass.
%   DON reactions belong to the paired microbial C/N ledger called from
%   soil.nitrogen.advanceAcceptedCandidate. This routine has no independent
%   DON-to-NH4 or DON-to-NO3 reaction and never mutates the accepted input.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
% Change: 2026-09-01 Codex -- remove the standalone DON ammonification
%   bypass so all soil DON processing is owned by the explicit microbial C-N
%   ledger; MEMC/MEND source equations are unchanged.
% Change: 2026-09-14 Jizixing Zhu -- support explicitly supplied dynamic
%   phase inventories; transport free DOM only and retain sorbed C/N exactly.
    dynamic=strcmp(parameters.sorption_model,'mend_dynamic');
    docPhase=soil.dom.phaseState(accepted,'doc',theta,bulkDensity,parameters);
    donPhase=soil.dom.phaseState(accepted,'don',theta,bulkDensity,parameters);
    if dynamic
        % The existing free-solute ADE needs positive liquid capacity.
        % Do not change dry/frozen transport physics through a water floor.
        assert(all(theta(:)>0),'soil:dom:dynamicDryTransport', ...
            'Dynamic DOM transport currently requires positive liquid water in every layer.');
        docInput=docPhase.free_g_cm3_bulk;donInput=donPhase.free_g_cm3_bulk;
        kdDoc=0;kdDon=0;
    else
        docInput=accepted.doc_total_g_cm3;donInput=accepted.don_total_g_cm3;
        kdDoc=parameters.kd_doc_cm3_g;kdDon=parameters.kd_don_cm3_g;
    end
    [doc,docTransport]=soil.dom.transportConservative( ...
        docInput,theta,bulkDensity, ...
        kdDoc,dz_cm,q_cm_s, ...
        parameters.diffusion_cm2_s,parameters.dispersion_cm,dt_s, ...
        parameters.transport_substep_s);
    [don,donTransport]=soil.dom.transportConservative( ...
        donInput,theta,bulkDensity, ...
        kdDon,dz_cm,q_cm_s, ...
        parameters.diffusion_cm2_s,parameters.dispersion_cm,dt_s, ...
        parameters.transport_substep_s);
    if dynamic
        doc=doc+accepted.doc_sorbed_g_cm3;
        don=don+accepted.don_sorbed_g_cm3;
        % Public storage diagnostics retain their total-inventory meaning.
        sorbedC=sum(accepted.doc_sorbed_g_cm3.*dz_cm(:));
        sorbedN=sum(accepted.don_sorbed_g_cm3.*dz_cm(:));
        docTransport.storage_before_g_cm2=docTransport.storage_before_g_cm2+sorbedC;
        docTransport.storage_after_g_cm2=docTransport.storage_after_g_cm2+sorbedC;
        donTransport.storage_before_g_cm2=donTransport.storage_before_g_cm2+sorbedN;
        donTransport.storage_after_g_cm2=donTransport.storage_after_g_cm2+sorbedN;
    end
    assert(all(don>=-1e-14)&&all(isfinite(don)), ...
        'soil:dom:badDonCandidate','DON candidate is invalid.');
    don=max(don,0);
    next=accepted;
    next.doc_total_g_cm3=doc;
    next.don_total_g_cm3=don;
    next.direct_don_to_no3=false;
    zeroRate=zeros(size(don));
    candidate=struct('state',next, ...
        'nh4_source_ug_cm3_s',zeroRate, ...
        'no3_source_ug_cm3_s',zeroRate, ...
        'mineralized_don_g_cm3',zeroRate, ...
        'transport',struct('doc',docTransport,'don',donTransport), ...
        'direct_don_to_no3',false, ...
        'reaction_ownership','explicit_microbial_cn_ledger');
end
