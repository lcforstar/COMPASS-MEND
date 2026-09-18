function [Dsh, DIF2, DISP2, Dsh_h, Dsh_T, Dsh_A] = soluteDispersion( ...
        DIF1111, nS, NL, Theta_LL, Theta_s, QL, QL_h, QL_T, QL_A, DISP11, NN, DeltZ, ...
        R, Tr, TT, eT, Delt_t, Srt, Theta_L, lArtD, lUpW, RHO_b_Kd, PeCr, Retard)
%SOLUTEDISPERSION  Hydrodynamic dispersion coefficient (NH4 + NO3) per node.
%
% Computes the apparent dispersion D_sh = α_L × |q|/θ + D_eff for each
% solute species (s = 1 NH4, s = 2 NO3). The pore-water velocity |q|/θ is
% taken from the liquid water flux QL, with its three sub-components
% (QL_h, QL_T, QL_A) propagated for sensitivity output.
%
% Inputs (intent and units):
%   DIF1111  [in]    [m²/s, 1×NL×nS]      molecular diffusion D_0 per layer per species
%   nS       [in]    [-]                   number of solute species (= 2: NH4, NO3)
%   NL       [in]    [-]                   number of layers
%   Theta_LL [in]    [m³/m³, NL×2]        liquid content (top, bottom node)
%   Theta_s  [in]    [m³/m³, 1×NL]        saturated content (transposed in caller)
%   QL       [in]    [cm/s, NL×2]         liquid Darcy flux per layer node
%   QL_h     [in]    [cm/s, NL×2]         flux contribution from h-gradient
%   QL_T     [in]    [cm/s, NL×2]         flux from T-gradient
%   QL_A     [in]    [cm/s, NL×2]         flux from air-pressure gradient
%   DISP11   [in]    [cm, NL×1]            longitudinal dispersivity α_L per layer
%   NN       [in]    [-]                   number of nodes (= NL+1)
%   DeltZ    [in]    [cm, NL×1]           layer thickness
%   R        [in]    [J/(mol·K)]          gas constant = 8.314
%   Tr       [in]    [°C]                  reference temperature
%   TT       [in]    [°C, NN×1]           soil temperature at nodes
%   eT       [inout] [-, NN×1]            temperature modulator (overwritten)
%   Delt_t   [in]    [s]                   outer step length
%   Srt      [in]    [1/s, NL×2]          root water uptake sink (advective)
%   Theta_L  [in]    [m³/m³, NL×2]        previous-step liquid content
%   lArtD    [in]    [0/1]                 artificial dispersion switch
%   lUpW     [in]    [0/1]                 upwind scheme switch
%   RHO_b_Kd   [in]    [g/cm³, NL×1]        bulk density × Kd
%   PeCr     [in]    [-]                   critical Peclet number = 2
%   Retard   [in]    [-, NL×2]            retardation factor R_d
%
% Returns:
%   Dsh      [cm²/s, NL×2×nS]   total hydrodynamic dispersion per node + species
%   DIF2     [cm²/s, NL×2×nS]   effective molecular diffusion (Millington-Quirk)
%   DISP2    [cm²/s, NL×2×nS]   mechanical dispersion term
%   Dsh_h    [cm²/s, NL×2×nS]   sensitivity vs h-driven flux component
%   Dsh_T    [cm²/s, NL×2×nS]   sensitivity vs T-driven flux component
%   Dsh_A    [cm²/s, NL×2×nS]   sensitivity vs A-driven flux component
%
% Method:
%   Millington-Quirk diffusion correction: f_w = θ^(7/3) / θ_s²
%   D_eff(s) = f_w × D_0(s) × e_T
%   D_disp = α_L × |q|/θ
%   Dsh = D_eff + D_disp; artificial dispersion (lArtD=1) bumps Dsh when |Pe|>PeCr.
%
% Reference:
%   Millington & Quirk 1961 Trans Faraday Soc 57:1200 (tortuosity)
%   van Genuchten & Wagenet 1989 SSSAJ (ADE numerical form)
%   STEMMUS-N legacy CalculateSoluteDispersion.m

    EA=0;
    TD=zeros(NN);
    
    QL_nodes = equations.interpolateElementToNode(QL, DeltZ);
    QL_nodes_h = equations.interpolateElementToNode(QL_h, DeltZ);
    QL_nodes_T = equations.interpolateElementToNode(QL_T, DeltZ);
    QL_nodes_A = equations.interpolateElementToNode(QL_A, DeltZ);
    % QL_nodes = CalculateVelocities(QL, NL, NN, DeltZ, Theta_LL, Theta_L, Delt_t, Srt);
    
    fw = Theta_LL.^(7/3) ./ Theta_s'.^2;  % 向量化操作

    for MN=1:NN
        %Arrhenius FORMULA
        TD(MN)=(TT(MN)-Tr)/R/(TT(MN)+273.15)./(Tr+273.15);
        eT(MN)=exp(EA*TD(MN));
    end

    if eT(MN) > 1
        eT(MN) = 1;
    end    

    for js=1:nS
        for ML = 1:NL
            for ND = 1:2
                MN = ML + ND - 1;  % 节点编号
                DIF2(ML, ND, js) = DIF1111(1,ML,js)' * fw(ML,ND) * eT(MN);  % 利用广播
            end
        end
    end
   
    for js = 1:nS
        for ML = 1:NL
            J = ML;  % 当前单元的土壤类型索引
            for ND = 1:2
                MN = ML + ND - 1;  % 节点编号
                DISP2(ML, ND, js) = DISP11(J) * abs(QL_nodes(MN));  % 使用对应土壤类型的DISP1     
            end
        end
    end

    for js = 1:nS
        for ML = 1:NL
            J = ML;  % 当前单元的土壤类型索引
            for ND = 1:2
                MN = ML + ND - 1;  % 节点编号
                DISP2_h(ML, ND, js) = DISP11(J) * abs(QL_nodes_h(MN));  % 使用对应土壤类型的DISP1     
            end
        end
    end

    for js = 1:nS
        for ML = 1:NL
            J = ML;  % 当前单元的土壤类型索引
            for ND = 1:2
                MN = ML + ND - 1;  % 节点编号
                DISP2_T(ML, ND, js) = DISP11(J) * abs(QL_nodes_T(MN));  % 使用对应土壤类型的DISP1     
            end
        end
    end

    for js = 1:nS
        for ML = 1:NL
            J = ML;  % 当前单元的土壤类型索引
            for ND = 1:2
                MN = ML + ND - 1;  % 节点编号
                DISP2_A(ML, ND, js) = DISP11(J) * abs(QL_nodes_A(MN));  % 使用对应土壤类型的DISP1     
            end
        end
    end


    for js = 1:nS
        for ML = 1:NL
            for ND = 1:2
                Dsh(ML, ND, js) = Theta_LL(ML, ND) * DIF2(ML, ND, js) + DISP2(ML, ND, js);
            end
        end
    end

     for js = 1:nS
        for ML = 1:NL
            for ND = 1:2
                Dsh_h(ML, ND, js) = Theta_LL(ML, ND) * DIF2(ML, ND, js) + DISP2_h(ML, ND, js);
            end
        end
     end

      for js = 1:nS
        for ML = 1:NL
            for ND = 1:2
                Dsh_T(ML, ND, js) = Theta_LL(ML, ND) * DIF2(ML, ND, js) + DISP2_T(ML, ND, js);
            end
        end
      end

      for js = 1:nS
          for ML = 1:NL
              for ND = 1:2
                  Dsh_A(ML, ND, js) = Theta_LL(ML, ND) * DIF2(ML, ND, js) + DISP2_A(ML, ND, js);
              end
          end
      end


    if(~ lArtD && ~lUpW)
        DPom=zeros(NL,2,2);
        for ML = 1:NL
            for ND = 1:2
                DPom(ML, ND, 1) = max(Delt_t / (6  * ( Theta_LL(ML, ND) + RHO_b_Kd(J))),0);
                DPom(ML, ND, 2) = max(Delt_t / (6  * Theta_LL(ML, ND)),0);
            end
        end
        for js = 1:nS
            for ML = 1:NL
                for ND = 1:2
                    MN = ML + ND -1;
                    Dsh(ML, ND, js) = max(Dsh(ML, ND, js) - QL_nodes(MN) * QL_nodes(MN) * DPom(ML, ND, js) , Dsh(ML, ND, js)/2);
                end
            end
        end
        for js = 1:nS
            for ML = 1:NL
                for ND = 1:2
                    MN = ML + ND -1;
                    Dsh_h(ML, ND, js) = max(Dsh_h(ML, ND, js) - QL_nodes_h(MN) * QL_nodes_h(MN) * DPom(ML, ND, js) , Dsh_h(ML, ND, js)/2);
                end
            end
        end
        for js = 1:nS
            for ML = 1:NL
                for ND = 1:2
                    MN = ML + ND -1;
                    Dsh_T(ML, ND, js) = max(Dsh_T(ML, ND, js) - QL_nodes_T(MN) * QL_nodes_T(MN) * DPom(ML, ND, js) , Dsh_T(ML, ND, js)/2);
                end
            end
        end
        for js = 1:nS
            for ML = 1:NL
                for ND = 1:2
                    MN = ML + ND -1;
                    Dsh_A(ML, ND, js) = max(Dsh_A(ML, ND, js) - QL_nodes_A(MN) * QL_nodes_A(MN) * DPom(ML, ND, js) , Dsh_A(ML, ND, js)/2);
                end
            end
        end
    end


    if(lArtD)
        DD=zeros(NL,2,2);
        if(PeCr ~= 0 && abs(min(QL_nodes)) > 1.e-15)
            for ML = 1:NL
                for ND = 1:2
                    MN = ML + ND - 1;  % 节点编号（修复：原缺此行，误用上一循环残留的 MN=NN）
                    DD(ML, ND, 1) = QL_nodes(MN) * QL_nodes(MN) * Delt_t / Theta_LL(ML, ND) / Retard(ML, ND) / PeCr;
                    DD(ML, ND, 2) = QL_nodes(MN) * QL_nodes(MN) * Delt_t / Theta_LL(ML, ND) / 1 / PeCr;
                end
            end
        end
        for js = 1:nS
            for ML = 1:NL
                for ND = 1:2
                    if(DD(ML, ND, js) > Dsh(ML, ND, js))
                        Dsh(ML, ND, js) = DD(ML, ND, js);
                    end
                end
            end
        end
    end



















