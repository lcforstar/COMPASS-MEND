function writeBioparticleExports(reg,built,kt)
%WRITEBIOPARTICLEEXPORTS Accepted mobile C/N bottom exports, never DOM losses.
% Step/cumulative amounts [g element cm-2], mean rates [g element cm-2 s-1].
% MB, BD and intracellular ENZN are disjoint cargo; total is derived only.
% Disabled transport/closed bottom gives physical zero. No state is mutated.
% Editor: Wenjie Liu; Date: 2026-09-29.
    enabled=built.parameters.soil_microbes.transport_enabled;
    amounts=zeros(3,2);cumulative=amounts;dt=1;
    groups={'MB','BD','ENZN'};elements={'carbon','nitrogen'};
    if enabled
        assert(isfield(built.workspace.soil_microbes,'bioparticle_transport'), ...
            'diagnostics:particleLedger','Enabled particle transport lacks accepted ledger.');
        work=built.workspace.soil_microbes.bioparticle_transport;dt=work.elapsed_s;
        assert(isfinite(dt)&&dt>0,'diagnostics:particleTime','Invalid accepted elapsed time.');
        memory=[];
        if isfield(built.Y.soil_microbes,'bottom_export_cumulative')
            memory=built.Y.soil_microbes.bottom_export_cumulative;
        end
        for i=1:3
            for j=1:2
                amounts(i,j)=work.diagnostics.(groups{i}).([elements{j} '_transport']).bottom_outflow_g_cm2;
                if ~isempty(memory),cumulative(i,j)=memory.(groups{i}).([elements{j} '_g_cm2']);end
            end
        end
    end
    amounts=[amounts;sum(amounts,1)];cumulative=[cumulative;sum(cumulative,1)];
    groups{4}='total';
    for i=1:4
        for j=1:2
            name=['bioparticle_' lower(groups{i}) '_' elements{j} '_bottom_outflow_'];
            reg.write([name 'step'],amounts(i,j),kt);
            reg.write([name 'flux'],amounts(i,j)/dt,kt);
            reg.write([name 'cumulative'],cumulative(i,j),kt);
        end
    end
end
