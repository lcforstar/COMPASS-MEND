function writeOrganicBoundaryExchange(reg,built,kt)
%WRITEORGANICBOUNDARYEXCHANGE Accepted extra organic boundary routes only.
% Existing bottom-outflow writers remain authoritative. New top in/out and
% bottom-input values are positive g element/cm2, rates per elapsed second.
% Editor: Wenjie Liu; Date: 2026-09-29.
    routes={'top_inflow','top_outflow','bottom_inflow'};
    for carrier={'doc','don'}
        key=carrier{1};ledger=built.workspace.soil_dom.transport.(key);
        memory=[];
        if isfield(built.Y.soil_dom,'boundary_exchange_cumulative')
            memory=built.Y.soil_dom.boundary_exchange_cumulative.(key);
        end
        writeRoutes(key,ledger,memory,ledger.driver_elapsed_s);
    end
    enabled=built.parameters.soil_microbes.transport_enabled;
    for element={'carbon','nitrogen'}
        e=element{1};sumAmounts=zeros(1,3);sumMemory=sumAmounts;dt=1;
        for group={'MB','BD','ENZN'}
            g=group{1};key=['bioparticle_' lower(g) '_' e];ledger=[];memory=[];
            if enabled
                work=built.workspace.soil_microbes.bioparticle_transport;dt=work.elapsed_s;
                ledger=work.diagnostics.(g).([e '_transport']);
                if isfield(built.Y.soil_microbes,'boundary_exchange_cumulative')
                    memory=built.Y.soil_microbes.boundary_exchange_cumulative.([g '_' e]);
                end
            end
            [a,c]=writeRoutes(key,ledger,memory,dt);
            sumAmounts=sumAmounts+a;sumMemory=sumMemory+c;
        end
        for j=1:3
            writeValues(['bioparticle_total_' e '_' routes{j}],sumAmounts(j),sumMemory(j),dt);
        end
    end
    function [a,c]=writeRoutes(key,ledger,memory,dt)
        a=zeros(1,3);c=a;
        for j=1:3
            field=[routes{j} '_g_cm2'];
            if ~isempty(ledger),a(j)=ledger.(field);end
            if ~isempty(memory),c(j)=memory.(field);end
            writeValues([key '_' routes{j}],a(j),c(j),dt);
        end
    end
    function writeValues(key,amount,cumulative,dt)
        reg.write([key '_step'],amount,kt);
        reg.write([key '_flux'],amount/dt,kt);
        reg.write([key '_cumulative'],cumulative,kt);
    end
end
