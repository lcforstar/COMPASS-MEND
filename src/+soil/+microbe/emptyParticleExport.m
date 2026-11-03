function memory=emptyParticleExport()
%EMPTYPARTICLEEXPORT Accepted-runtime-memory schema, g C/N per cm2 column.
% Each group is disjoint; ENZN is explicit intracellular cargo, not biomass
% counted a second time. Only candidate/commit updates cumulative amounts.
% Editor: Wenjie Liu; Date: 2026-09-29.
    memory=struct('schema_version',1);
    for g={'MB','BD','ENZN'}
        memory.(g{1})=struct('carbon_g_cm2',0,'nitrogen_g_cm2',0);
    end
end
