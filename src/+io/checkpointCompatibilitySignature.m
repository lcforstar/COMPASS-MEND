function signature=checkpointCompatibilitySignature(built)
%CHECKPOINTCOMPATIBILITYSIGNATURE Describe restart-critical BGC topology.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-08-28 Codex -- include native organic-C/N, dormancy, and
%     BIOPARTICLE topology in the restart-compatibility contract.
% Change: 2026-09-14 Jizixing Zhu -- bind source-hour restarts to their explicit
%   hourly policy/parameters in addition to the distinct pool schemas.
% Change: 2026-09-16 Jizixing Zhu -- bind source-hour equilibrium restarts
% to the DOC/DON Kd profiles that define their algebraic phase inventories.
    sw=built.control.switches;
    signature=struct('NL',built.domain.NL,'NN',built.domain.NN, ...
        'soil_organic',readFlag(sw,'lSoilOrganic'), ...
        'soil_dom',readFlag(sw,'lSoilDOM'), ...
        'explicit_soil_microbes',readFlag(sw,'lExplicitSoilMicrobes'), ...
        'litter_input',readFlag(sw,'lLitterInput'), ...
        'microbial_dormancy',readFlag(sw,'lMicrobialDormancy'), ...
        'microbial_transport',readFlag(sw,'lMicrobialTransport'), ...
        'organic_schema',stateSchema(built.Y,'soil_organic'), ...
        'dom_schema',stateSchema(built.Y,'soil_dom'), ...
        'microbe_schema',stateSchema(built.Y,'soil_microbes'));
    if isfield(built.parameters,'soil_microbes') ...
        &&isfield(built.parameters.soil_microbes,'hourly')
        signature.mend_hourly=built.parameters.soil_microbes.hourly;
        if strcmp(built.parameters.soil_dom.sorption_model,'linear_equilibrium')
            signature.dom_equilibrium=struct('sorption_model','linear_equilibrium', ...
                'kd_doc_cm3_g',built.parameters.soil_dom.kd_doc_cm3_g, ...
                'kd_don_cm3_g',built.parameters.soil_dom.kd_don_cm3_g);
        end
    end
end

function value=readFlag(s,name)
    value=false; if isfield(s,name), value=logical(s.(name)); end
end

function value=stateSchema(Y,name)
    value=0;
    if isfield(Y,name)&&isfield(Y.(name),'schema_version')
        value=Y.(name).schema_version;
    end
end
