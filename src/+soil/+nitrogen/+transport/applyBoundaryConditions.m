function [RHS,RHS1,C4,C4_a,C3,C3_a,cTop_NO3,cTop_NH4] = applyBoundaryConditions(RHS1,C3,C3_a,RHS,NBCN,NBCNB,KT,C4,NN,C4_a,QL_nodes,Thmrlefc,QMB,cBot_NH4,cTop_NH4,cBot_NO3,cTop_NO3,lBNF,BNF)
%APPLYBOUNDARYCONDITIONS  Apply NH4/NO3 transport boundary conditions to the linear system.
%   QL_nodes follows the global vertical coordinate: positive is upward.
%   Therefore downward top inflow has QL_nodes(NN)<0, whereas upward bottom
%   inflow has QL_nodes(1)>0. Public input/loss diagnostics are derived from
%   this raw convention after the solve; the ADE boundary equations below
%   retain their historical signs.
% Change: 2026-09-21 Wenjie Liu -- state the shared raw-flux sign contract.
%% NH4+__N
%%%%%%%%%% Apply the bottom boundary condition called for by NBCNB %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if NBCNB==1            %-----> Specify concentration at bottom to be ---BCNB;
        RHS(1)=cBot_NH4(KT);
        C3(1,1)=1;
        RHS(2)=RHS(2)-C3_a(1)*RHS(1);
        C3(1,2)=0;
        C3_a(1)=0;
    elseif NBCNB==2        %-----> Specify concentration flux at bottom to be ---BCNB (Positive upwards);
        if QL_nodes(1)>0 || (Thmrlefc==1 && QMB(KT)==0)
            RHS(1)=RHS(1)+QL_nodes(1)*cBot_NH4(KT);
        else
            C3(1,1)=-1;
            C3(1,2)=1;
            RHS(1)=0;
        end
    elseif NBCNB==3        %-----> NBChB=3,Gravity drainage at bottom--specify concentration flux= 0;
        C3(1,1)=-1;
        C3(1,2)=1;
        RHS(1)=0;
    end
    
%%%%%%%%%% Apply the surface boundary condition called for by NBCN  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if NBCN==1             %-----> Specified concentration at surface---equal to cTop_NH4;
        if lBNF == 0
        RHS(NN)=cTop_NH4(KT);
        else
         RHS(NN)=cTop_NH4(KT)+BNF;
        end
        C3(NN,1)=1;
        RHS(NN-1)=RHS(NN-1)-C3(NN-1,2)*RHS(NN);
        C3(NN-1,2)=0;
        C3_a(NN-1)=0;
    else
        if QL_nodes(NN)<0
            RHS(NN)=RHS(NN)-QL_nodes(NN)*cTop_NH4(KT);
        end
    end

   
    %% NO3-__N
%%%%%%%%%% Apply the bottom boundary condition called for by NBCNB %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if NBCNB==1            %-----> Specify concentration at bottom to be ---BCNB;
        RHS1(1)=cBot_NO3(KT);
        C4(1,1)=1;
        RHS1(2)=RHS1(2)-C4_a(1)*RHS1(1);
        C4(1,2)=0;
        C4_a(1)=0;
    elseif NBCNB==2        %-----> Specify concentration flux at bottom to be ---BCNB (Positive upwards);
        if QL_nodes(1)>0 || (Thmrlefc==1 && QMB(KT)==0)
            RHS1(1)=RHS1(1)+QL_nodes(1)*cBot_NO3(KT);
        else
            C4(1,1)=-1;
            C4(1,2)=1;
            RHS1(1)=0;
        end
    elseif NBCNB==3        %-----> NBChB=3,Gravity drainage at bottom--specify concentration flux= 0;       
        C4(1,1)=-1;
        C4(1,2)=1;
        RHS1(1)=0;
    end
    
%%%%%%%%%% Apply the surface boundary condition called for by NBCN  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if NBCN==1             %-----> Specified concentration at surface---equal to cTop_NH4;
        RHS1(NN)=cTop_NO3(KT);
        C4(NN,1)=1;
        RHS1(NN-1)=RHS1(NN-1)-C4(NN-1,2)*RHS1(NN);
        C4(NN-1,2)=0;
        C4_a(NN-1)=0;
    else
        if QL_nodes(NN) < 0  % If flux is downward
            RHS1(NN) = RHS1(NN) - QL_nodes(NN) * cTop_NO3(KT);%cTop_NO3(KT)/86400;%QL_nodes(NN) * cTop_NO3(KT);
        end
    end
end


% %  if(WLayer.and.hNew(NumNP).gt.0.) then ! mass balance in the surface layer
%         hT=hNew(NumNP)
%         do 15 jj=1,NS
%           if((hT+dt*(Prec-rSoil)).gt.0.)
%      !      cTop_NH4(jj)=(hT*cTop_NH4(jj)+dt*Prec*cT(jj))/(hT+dt*(Prec-rSoil))
% 15      continue
%       end if  解释代码物理含义来自于hydrus    考虑hnew有积水的时候再加，rsoil是蒸发
