function [SoilVariables, CHK, SAVE, QET, QEB, EnergyVariables] = solveThermalEquations( ...
    SoilVariables, TransportCoefficient, MatricCoeff, AirVars, cfg, ...
    BoundaryCondition, HBoundaryFlux, MeteoData, Precip, EVAP, ...
    Srt, L_f, T, h, P_g, L, Delt_t, KT, SAVE, QET, QEB, airHeatImpedance)
% solveThermalEquations - 鑳介噺鏂圭▼瀹屾暣姹傝В绠＄嚎锛堟浛浠?Enrgy_sub.m锛?%                         Full thermal solve pipeline (replaces Enrgy_sub.m)
%
% 鍑芥暟瀹氫箟 / Function Definition:
%   [SoilVariables, CHK, SAVE, QET, QEB] = solveThermalEquations( ...
%       SoilVariables, TransportCoefficient, MatricCoeff, AirVars, cfg, ...
%       BoundaryCondition, HBoundaryFlux, MeteoData, Precip, EVAP, ...
%       Srt, L_f, T, h, P_g, L, Delt_t, KT, SAVE, QET, QEB, airHeatImpedance)
%
% 杈撳叆 / Inputs:
%   SoilVariables    - 鍦熷￥鐘舵€佺粨鏋勪綋 / Soil state struct
%     .TT            - 褰撳墠 Newton 姝ユ俯搴?(NN脳1) [K]
%   TransportCoefficient - 浼犺緭绯绘暟缁撴瀯浣?/ Transport coefficient struct
%   MatricCoeff      - 姘村垎鐭╅樀绯绘暟 (soil.water.calculateMatricCoefficients 杈撳嚭)
%     .Khh/.KhT/.Kha/.Vvh/.VvT - 娑叉按瀵煎害/瀵规祦绯绘暟 (NL脳nD)
%   AirVars          - 姘旂浉姹傝В缁撴灉 (lSoilairefc=1 鏃舵湁鏁? 鍚﹀垯浼?[])
%   cfg              - 妯″瀷璁剧疆 / Model settings (cfg.Grid, cfg.lSoilairefc, ...)
%   BoundaryCondition - BC 缁撴瀯浣?/ Boundary condition struct
%     .NBCT          - 椤堕儴 BC 绫诲瀷 (1=Dirichlet Ts from ebal)
%     .BCT           - 椤堕儴 BC 鍊?/ Top BC value (Ts [K/掳C] from ebal)
%     .NBCTB         - 搴曢儴 BC 绫诲瀷 (1=Dirichlet Tbh)
%     .BCTB          - 搴曢儴 BC 鍊?/ Bottom BC value (Tbh [K/掳C])
%     .DSTOR0        - 鍒濆绉按娣卞害 [cm]
%   HBoundaryFlux    - 姘村垎杈圭晫閫氶噺 / Water boundary flux struct
%     .QMB           - 搴曢儴姘村垎閫氶噺 [cm/s]
%   MeteoData      - 姘旇薄寮鸿揩缁撴瀯浣?/ Meteorological forcing struct
%     .Ta_msr        - 姘旀俯鏃堕棿搴忓垪 [掳C] / Air temperature time series
%   Precip           - 闄嶆按鏃堕棿搴忓垪 [cm/s] / Precipitation time series
%   EVAP             - 钂稿彂鏃堕棿搴忓垪 [cm/s] / Evaporation time series
%   Srt              - 鏍圭郴鍚告按 (NL脳nD) [cm/s] / Root water uptake
%   L_f              - 鐔斿寲娼滅儹 [J/kg] / Latent heat of fusion (0 = no freeze)
%   T                - 涓婁竴鏃舵娓╁害 (NN脳1) / Previous outer time step temperature
%   h                - 涓婁竴鏃舵姘村ご (NN脳1) / Previous outer time step head
%   P_g              - 涓婁竴鏃舵姘斿帇 (NN脳1) / Previous outer time step pressure
%   L                - 姹藉寲娼滅儹 (NN脳1) [J/g] / Latent heat of vaporization
%   Delt_t           - 鏃舵闀?[s] / Time step size
%   KT               - 澶栧眰鏃舵绱㈠紩 / Outer time step index
%   SAVE             - 绯绘暟蹇収 (3脳3脳5) / Coefficient snapshot
%   QET              - 椤堕儴鐑€氶噺鏃堕棿搴忓垪 / Top heat flux time series
%   QEB              - 搴曢儴鐑€氶噺鏃堕棿搴忓垪 / Bottom heat flux time series
%   airHeatImpedance - 姘斿帇鐑樆鎶楃郴鏁帮紙榛樿 3.0锛? Air-pressure heat impedance scale
%
% 杈撳嚭 / Outputs:
%   SoilVariables    - 鏇存柊鍚庣粨鏋勪綋锛圱T 宸叉洿鏂帮級/ Updated struct (TT updated)
%   CHK              - Newton 鏀舵暃娈嬪樊 (NN脳1) / Newton convergence residuals
%   SAVE             - 鏇存柊鍚庡揩鐓?/ Updated coefficient snapshot
%   QET              - 鏇存柊鍚庨《閮ㄧ儹閫氶噺 / Updated top heat flux series
%   QEB              - 鏇存柊鍚庡簳閮ㄧ儹閫氶噺 / Updated bottom heat flux series
%
% 璋冪敤閾?/ Call chain:
%   calculateEnergyParameters 鈫?calculateMatricCoefficients 鈫?%   assembleCoefficientMatrices 鈫?calculateBoundaryConditions 鈫?%   solveTridiagonalMatrixEquations 鈫?calculateBoundaryFlux
%
% 寮曠敤 / References:
%   鍘?Enrgy_sub.m锛汼TEMMUS Technical Notes Section 5
%
% 淇敼鑰?/ Modified by:   鏈卞瘎瀛愭槦 (Zhu Jizixing)
% 淇敼鏃ユ湡 / Modified date: 2026-04-25

    % 1. 鑳介噺鏂圭▼绯绘暟 / Energy PDE coefficients
    % stor 淇锛?026-05-02锛夊悗锛孠_g 宸插甫 Lundin 鍐伴樆鎶?10^(-惟路R_ice)锛孠Ta/VTa
    % 閫氳繃 尾_g銆並aa 鑷姩缁ф壙鍐扮浉鏁堝簲銆傚啀鍋?exp(-惟路R_ice) 鍗冲弻閲嶅啺闃绘姉锛屼笌 Yu 2018
    % Disable the optional extra air-heat impedance by default.
    if nargin < 22 || isempty(airHeatImpedance)
        airHeatImpedance = 0;
    end

    SoilVariables.L = L(:);
    EnergyVariables = soil.thermal.calculateEnergyParameters( ...
        SoilVariables, TransportCoefficient, MatricCoeff, AirVars, cfg, Srt, L_f, airHeatImpedance);

    % 2. 涓夊瑙掔煩闃电粍瑁?/ Assemble tridiagonal matrices
    EnergyMatrices = soil.thermal.assembleThermalEquationMatrices(EnergyVariables, cfg);

    % 3. RHS 缁勮 / Assemble RHS
    [RHS, EnergyMatrices, SAVE] = soil.thermal.assembleCoefficientMatrices( ...
        EnergyMatrices, SoilVariables, cfg, Delt_t, T, h, P_g, SAVE);
    % 4. 杈圭晫鏉′欢 / Apply boundary conditions
    [RHS, EnergyMatrices] = soil.thermal.calculateBoundaryConditions( ...
        BoundaryCondition, EnergyMatrices, HBoundaryFlux, MeteoData, ...
        SoilVariables, Precip, EVAP, Delt_t, RHS, L, KT);

    % 5. Thomas 绠楁硶姹傝В / Solve via Thomas algorithm
    debugThermal = (isprop(cfg, 'lDEBUG_NAN_TRACE') && logical(cfg.lDEBUG_NAN_TRACE)) ...
        || (isprop(cfg, 'lDebug') && logical(cfg.lDebug));
    if debugThermal && KT <= 2
        fprintf('[Thermal KT=%d] diag=[%e,%e] rhs=[%e,%e] C5_a=[%e,%e] BCT=%e BCTB=%e\n', ...
            KT, min(EnergyMatrices.C5(:,1)), max(EnergyMatrices.C5(:,1)), ...
            min(RHS), max(RHS), min(EnergyMatrices.C5_a), max(EnergyMatrices.C5_a), ...
            BoundaryCondition.BCT, BoundaryCondition.BCTB);
        [minDiag, minDiagIdx] = min(EnergyMatrices.C5(:,1));
        if KT == 2 && minDiag < 0
            leftElem = max(minDiagIdx - 1, 1);
            rightElem = min(minDiagIdx, cfg.Grid.NL);
            fprintf(['[Thermal NEGDIAG] node=%d diag=%e C2=[%e,%e] C5sup=%e C5subPrev=%e ' ...
                'KTT_left=[%e,%e] KTT_right=[%e,%e] CTT_left=[%e,%e] CTT_right=[%e,%e] ' ...
                'CTTliq_left=[%e,%e] CTTliq_right=[%e,%e] DThetaT_left=[%e,%e] DThetaT_right=[%e,%e] ' ...
                'cunsat_left=[%e,%e] cunsat_right=[%e,%e] WW_left=[%e,%e] WW_right=[%e,%e] ' ...
                'KTT_range=[%e,%e] CTT_range=[%e,%e]\n'], ...
                minDiagIdx, minDiag, EnergyMatrices.C2(minDiagIdx,1), EnergyMatrices.C2(minDiagIdx,2), ...
                EnergyMatrices.C5(minDiagIdx,2), EnergyMatrices.C5_a(leftElem), ...
                EnergyVariables.KTT(leftElem,1), EnergyVariables.KTT(leftElem,2), ...
                EnergyVariables.KTT(rightElem,1), EnergyVariables.KTT(rightElem,2), ...
                EnergyVariables.CTT(leftElem,1), EnergyVariables.CTT(leftElem,2), ...
                EnergyVariables.CTT(rightElem,1), EnergyVariables.CTT(rightElem,2), ...
                EnergyVariables.CTT_liq(leftElem,1), EnergyVariables.CTT_liq(leftElem,2), ...
                EnergyVariables.CTT_liq(rightElem,1), EnergyVariables.CTT_liq(rightElem,2), ...
                EnergyVariables.DTheta_LLT(leftElem,1), EnergyVariables.DTheta_LLT(leftElem,2), ...
                EnergyVariables.DTheta_LLT(rightElem,1), EnergyVariables.DTheta_LLT(rightElem,2), ...
                EnergyVariables.c_unsat(leftElem,1), EnergyVariables.c_unsat(leftElem,2), ...
                EnergyVariables.c_unsat(rightElem,1), EnergyVariables.c_unsat(rightElem,2), ...
                EnergyVariables.WW(leftElem,1), EnergyVariables.WW(leftElem,2), ...
                EnergyVariables.WW(rightElem,1), EnergyVariables.WW(rightElem,2), ...
                min(EnergyVariables.KTT(:)), max(EnergyVariables.KTT(:)), ...
                min(EnergyVariables.CTT(:)), max(EnergyVariables.CTT(:)));
        end
    end
    [SoilVariables, CHK, ~, ~] = soil.thermal.solveTridiagonalMatrixEquations( ...
        EnergyMatrices, SoilVariables, RHS, cfg);
    TT = SoilVariables.TT(:);

    % 6. 闈炵墿鐞嗗€间繚鎶わ細NaN 鎴?|T|鈮?00掳C 鏃跺洖閫€鍒版棫鏃跺眰娓╁害 T
    %    Safeguard against NaN / unphysical TT (allow negative 掳C in soil)
    if any(isnan(TT)) || any(abs(TT(1:cfg.Grid.NL + 1)) >= 100)
        TT = T(:);
    end
    TT(TT <= -272) = -272;

    % 7. 杈圭晫鐑€氶噺锛堢敤 clamp 鍚庣殑 TT锛屾爣閲忊啋鍐欏叆鏃跺簭鏁扮粍锛?    %    Boundary heat fluxes using clamped TT (scalar 鈫?time series)
    [QET(KT), QEB(KT)] = soil.thermal.calculateBoundaryFlux(SAVE, TT, cfg);

    SoilVariables.TT = TT;
end
