function [SoilVariables] = applySoilHeteroWithInitialFreezing(L_f, SoilVariables, cfg)
% applySoilHeteroWithInitialFreezing - 应用初始土壤非均质性与冻结效应
%
% 输入:
%   L_f                  - 冻结潜热 (J/kg)
%   SoilVariables        - 土壤变量结构体
%   cfg                  - legacy ModelSettings 对象

    NN = cfg.Grid.NN;

    SoilVariables.ISFT    = 0;
    SoilVariables.Gamma_h  = ones(NN, 1);   % 滞后因子初值（循环内 updateGammaH 需要前值）
    SoilVariables.Gamma_hh = ones(NN, 1);

    for i = 1:NN
        SoilVariables.h_frez = init.updateHfreez(i, L_f, SoilVariables);
        SoilVariables.hh_frez(i) = SoilVariables.h_frez(i);
        SoilVariables.h(i)       = SoilVariables.h(i) - SoilVariables.h_frez(i);
        SoilVariables.hh(i)      = SoilVariables.h(i);
        SoilVariables.SAVEh(i)   = SoilVariables.h(i);
        SoilVariables.SAVEhh(i)  = SoilVariables.hh(i);

        [SoilVariables.Gamma_h, SoilVariables.Gamma_hh] = init.updateGammaH(i, SoilVariables);

        if cfg.lThmrlefc == 1
            SoilVariables.TT(i) = SoilVariables.T(i);
        end

        if cfg.lSoilairefc == 1
            SoilVariables.P_gg(i) = SoilVariables.P_g(i);
        end

        if i < NN
            SoilVariables.XWRE(i, 1) = 0;
            SoilVariables.XWRE(i, 2) = 0;
        end
    end

end
