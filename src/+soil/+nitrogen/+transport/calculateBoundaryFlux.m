function [cvBot,cvTop,cvBot1,cvTop1] = calculateBoundaryFlux(SAVE,Nn,Nno3,NN,KT, ...
    cvTop,cvBot,cvTop1,cvBot1,QL_nodes,cBot_NH4,cTop_NH4,cTop_NO3,cBot_NO3,NBCN,NBCNB,rMin,Thmrlefc,QMB)
%CALCULATEBOUNDARYFLUX  Compute realized NH4/NO3 top + bottom solute fluxes after solute solve.
%   All four returned fluxes use the transport coordinate convention:
%   positive is vertically upward. Thus top inflow is negative, bottom
%   inflow is positive, and net input to the column is cvBot-cvTop.
%
%   The candidate slot is reconstructed on every Picard call. It must not
%   retain the preceding accepted step or accumulate repeated iterations.
%
% Change: 2026-09-21 Wenjie Liu -- document the raw upward-positive contract
%   and make candidate reconstruction idempotent across Picard iterations.

    cvTop(KT)=0;
    cvBot(KT)=0;
    cvTop1(KT)=0;
    cvBot1(KT)=0;

%% NH4+__N
    if NBCN == 2
        if QL_nodes(NN)<0, cvTop(KT)=cvTop(KT) + QL_nodes(NN)*cTop_NH4(KT); end
    else
        cvTop(KT)=SAVE(2,1,4)-SAVE(2,2,4)*Nn(NN-1)-SAVE(2,3,4)*Nn(NN);
    end
    
    if NBCNB == 2
        if QL_nodes(1)>=0 , cvBot(KT)=cvBot(KT)+QL_nodes(1)*cBot_NH4(KT); end
        if QL_nodes(1)<0 ,  cvBot(KT)=cvBot(KT)+QL_nodes(1)*Nn(1); end
        if Thmrlefc == 1 && QMB(KT) == 0 ,  cvBot(KT)=0; end
    elseif   NBCNB == 3
        % cvBot(KT)=cvBot(KT)+QL_nodes(1)*Nn(1); 
        cvBot(KT)=-SAVE(1,1,4)+SAVE(1,2,4)*Nn(1)+SAVE(1,3,4)*Nn(2);
    else
        cvBot(KT)=-SAVE(1,1,4)+SAVE(1,2,4)*Nn(1)+SAVE(1,3,4)*Nn(2);
    end
    
    if(abs(cvTop(KT)) < rMin) ,cvTop(KT)=0; end
    if(abs(cvBot(KT)) < rMin) , cvBot(KT)=0; end

%% NO3-__N

    if NBCN == 2
        if QL_nodes(NN)<0, cvTop1(KT)=cvTop1(KT) + QL_nodes(NN)*cTop_NO3(KT); end
    else
        cvTop1(KT)=SAVE(2,1,5)-SAVE(2,2,5)*Nno3(NN-1)-SAVE(2,3,5)*Nno3(NN);
    end
    
    if NBCNB == 2
        if QL_nodes(1)>=0 , cvBot1(KT)=cvBot1(KT)+QL_nodes(1)*cBot_NO3(KT); end
        if QL_nodes(1)<0 ,  cvBot1(KT)=cvBot1(KT)+QL_nodes(1)*Nno3(1); end
        if Thmrlefc == 1 && QMB(KT) == 0 ,  cvBot1(KT)=0; end
    elseif   NBCNB == 3
        % cvBot1(KT)=cvBot1(KT)+QL_nodes(1)*Nno3(1); 
        cvBot1(KT)=-SAVE(1,1,5)+SAVE(1,2,5)*Nno3(1)+SAVE(1,3,5)*Nno3(2);
    else
        cvBot1(KT)=-SAVE(1,1,5)+SAVE(1,2,5)*Nno3(1)+SAVE(1,3,5)*Nno3(2);
    end
    
    if(abs(cvTop1(KT)) < rMin) , cvTop1(KT)=0; end
    if(abs(cvBot1(KT)) < rMin) , cvBot1(KT)=0; end
