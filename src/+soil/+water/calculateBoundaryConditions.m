function [WaterMatrix, RHS, AVAIL0, Throughfall, Infiltration, Infiltration_act, TopFlux_potential, runoffHort] = ...
    calculateBoundaryConditions(WaterMatrix, RHS, SoilVariables, BoundaryCondition, cfg, ...
    KT, hN, NBChh, h_target, Throughfall, EVAP, Delt_t)
% calculateBoundaryConditions Apply bottom and surface boundary conditions for Richards equation.
%
% Atmospheric production requests keep one active set per nonlinear solve.
% In the storage regime the endpoint head is liquid depth and the top row
% receives the HYDRUS surface-layer terms 1/dt and d0/dt. Capacity overflow
% is a fixed-head reaction evaluated by the scheduler. Direct legacy callers
% retain their old pre-clipped path for compatibility fixtures.
% Editor: Jizixing Zhu; Date: 2026-07-17.

NBCh  = BoundaryCondition.NBCh;
NBChB = BoundaryCondition.NBChB;
BCh   = BoundaryCondition.BCh;
BChB  = BoundaryCondition.BChB;
Ks0 = BoundaryCondition.Ks0;
GW  = BoundaryCondition.GroundwaterSettings;
NN  = cfg.Grid.NN;

C4   = WaterMatrix.C4;
C4_a = WaterMatrix.C4_a;
KL_eff = SoilVariables.KL_h;
if (isstruct(SoilVariables) && isfield(SoilVariables, 'KfL_h')) || ...
        (isobject(SoilVariables) && isprop(SoilVariables, 'KfL_h'))
    KL_eff = SoilVariables.KfL_h;
end

if ~GW.GroundwaterCoupling
    if NBChB == 1
        RHS(1)   = BChB;
        C4(1, 1) = 1;
        RHS(2)   = RHS(2) - C4_a(1) * RHS(1);
        C4(1, 2) = 0;
        C4_a(1)  = 0;
    elseif NBChB == 2
        RHS(1) = RHS(1) + BChB;
    elseif NBChB == 3
        RHS(1) = RHS(1) - KL_eff(1, 1);
    end
else
    indxBotm        = GW.indxBotmLayer;     % bottom-to-top (STEMMUS convention), dynamic
    indxBotmLayer_R = GW.indxBotmLayer_R;   % top-to-bottom (MODFLOW convention), dynamic
    topLevel        = GW.topLevel;
    soilLayerThickness = GW.soilLayerThickness;

    if NBChB == 1
        RHS(indxBotm)     = GW.headBotmLayer - topLevel + soilLayerThickness(indxBotmLayer_R);
        C4(indxBotm, 1)   = 1;
        RHS(indxBotm + 1) = RHS(indxBotm + 1) - C4(indxBotm, 2) * RHS(indxBotm);
        C4(indxBotm, 2)   = 0;
        C4_a(indxBotm)    = 0;
    elseif NBChB == 2
        RHS(indxBotm) = RHS(indxBotm) + BChB;
    elseif NBChB == 3
        RHS(indxBotm) = RHS(indxBotm) - KL_eff(indxBotm, 1);
    end
end

% Production atmospheric requests defer runoff to the capacity reaction.
% The Ks pre-clip remains only in the direct legacy compatibility path.
surfaceRegime = 'legacy';
if isfield(BoundaryCondition, 'SurfaceRegime')
    surfaceRegime = BoundaryCondition.SurfaceRegime;
end
if strcmp(surfaceRegime, 'legacy')
    runoffHort = max(0, Throughfall(KT) - Ks0);
    effectivePrecip = Throughfall(KT) - runoffHort;
else
    % Richards reaction, not Ks alone, determines infiltration capacity.
    runoffHort = 0;
    effectivePrecip = Throughfall(KT);
end
AVAIL0               = effectivePrecip + BoundaryCondition.DSTOR0 / Delt_t;
Infiltration(KT)     = effectivePrecip;
Infiltration_act(KT) = effectivePrecip;
TopFlux_potential(KT) = EVAP(KT) - AVAIL0;

if NBCh == 1
    RHS(NN)     = hN;
    C4(NN, 1)   = 1;
    RHS(NN-1)   = RHS(NN-1) - C4(NN-1, 2) * RHS(NN);
    C4(NN-1, 2) = 0;
    C4_a(NN-1)  = 0;
elseif NBCh == 2 || NBCh == 3
    if strcmp(surfaceRegime, 'storage')
        % QMT = E-P + (d1-d0)/dt with d1=h_surface.
        atmosphericFlux = EVAP(KT) - effectivePrecip;
        C4(NN, 1) = C4(NN, 1) + 1 / Delt_t;
        RHS(NN) = RHS(NN) - atmosphericFlux ...
            + BoundaryCondition.DSTOR0 / Delt_t;
    elseif NBChh == 1
        hN          = h_target;
        RHS(NN)     = hN;
        C4(NN, 1)   = 1;
        RHS(NN-1)   = RHS(NN-1) - C4(NN-1, 2) * RHS(NN);
        C4(NN-1, 2) = 0;
        C4_a(NN-1)  = 0;
    else
        if NBCh == 2
            RHS(NN) = RHS(NN) - BCh;
        else
            RHS(NN) = RHS(NN) - TopFlux_potential(KT);
        end
    end
end

RHS = RHS(:);
WaterMatrix.C4   = C4;
WaterMatrix.C4_a = C4_a;

% ---- 诊断日志（暴雨期边界矩阵组装），cfg.lDEBUG_RAIN_BC 控制 ----
if cfg.lDEBUG_RAIN_BC
    rain_thresh = 1e-6;  % cm/s ≈ 0.036 mm/h，过滤干天
    if Throughfall(KT) > rain_thresh
        fprintf(['[bc-debug KT=%d] TF=%.4e Ks0=%.4e runoffH=%.4e effPrec=%.4e ' ...
                 'EVAP=%.4e DSTOR0=%.4e AVAIL0=%.4e TopFlux_pot=%.4e ' ...
                 'NBCh=%d NBChh=%d hN=%.4e h_tgt=%.4e ' ...
                 'hh(NN)=%.4e RHS(NN-1)=%.4e RHS(NN)=%.4e ' ...
                 'C4(NN-1,2)=%.4e C4(NN,1)=%.4e C4_a(NN-1)=%.4e\n'], ...
                KT, Throughfall(KT), Ks0, runoffHort, effectivePrecip, ...
                EVAP(KT), BoundaryCondition.DSTOR0, AVAIL0, TopFlux_potential(KT), ...
                NBCh, NBChh, hN, h_target, ...
                SoilVariables.hh(NN), RHS(NN-1), RHS(NN), ...
                C4(NN-1, 2), C4(NN, 1), C4_a(NN-1));
    end
end
end
