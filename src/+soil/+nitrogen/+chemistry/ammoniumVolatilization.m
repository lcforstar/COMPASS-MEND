function [Svol] = ammoniumVolatilization(TT, N, NL, PH1, KT, Svol, Tm, u_2, NN, NoTime)
%AMMONIUMVOLATILIZATION  Surface NH3 volatilization rate via convective mass transfer.
%
% Computes NH3 efflux at the soil surface as the product of (a) the Henry-equilibrium
% NH3 partial pressure derived from soil NH4+ pool and pH, and (b) a convective
% mass-transfer coefficient h_m driven by air temperature and wind speed.
%
% Inputs:
%   TT     [in]  [°C, NN×1]    soil temperature at nodes
%   N      [in]  [mg/L, NN×1]  NH4 concentration at nodes (used at surface NN)
%   NL     [in]  [-]            number of soil layers
%   PH1    [in]  [-, NL×1]     soil pH (uses PH1(5) here — surface band)
%   KT     [in]  [-]            current outer step index
%   Svol   [inout] [μg/(cm²·s), NL×2]  volatilization rate matrix (overwritten at NL,2)
%   Tm     [in]  [°C, n_days×1 or hourly]  near-surface air temperature time series
%   u_2    [in]  [m/s, same indexing as Tm]  2-m wind speed time series
%   NN     [in]  [-]            number of nodes (NL+1)
%   NoTime [in]  [-, KT×1]     outer-step → forcing-day index map
%
% Returns:
%   Svol   [μg/(cm²·s), NL×2]  surface-node volatilization rate (NL, 2) updated
%
% Method:
%   Henry's law:   KH = (0.2138/(Tm+273.15)) × 10^(6.123 − 1825/(Tm+273.15))
%   Convective h_m = 0.000612 × u_2^0.8 × (Tm+273.15)^0.382 × L^(-0.2)
%   Flux:          Svol = h_m × [KH × NH4(NN) / (1 + 10^(−0.05 + 2788/(T+273.15) − pH(5))) − C_atm]
%   where L = 0.0025 m is the volatilization characteristic length, C_atm = 0 (clean air).
%
% Reference:
%   Sherlock & Goh 1985 Fert Res; alternative Freney 1985 form provided commented.
%   STEMMUS-N legacy CalculateAmmoniumVolatilization.m

    C_ammonium = 0;    % [μg/m³] atmospheric NH3 background, assume clean air
    L = 0.05;          % [m] volatilization characteristic length — 苹果园树冠下混合层 5cm
                       %     (2026-05-17: 改自 0.0025m 裸地值; Sherlock & Goh 1985 原文针对牧草地约 5cm)
                       %     L↑ → h_m ∝ L^-0.2 略↓ (×0.72), 但匹配实际果园物理高度
    
    %方法一：Freney(1985)
    % KH = (2.39*10^5/(Tm(KT)+273.15))*exp(-4151/(Tm(KT)+273.15)); % Henry’s law constant Kh
    % hm = 48.4*u_2(KT)^0.8*(Tm(KT)+273.15)^-1.4; % m/s   % hm is the convective mass transfer coefficient
    % Svol(NL,2)=hm*((KH*(N(NN))/(1+10^(0.09018+(2729.92/(TT(NN)+273.15))-PH1(5))))-C_ammonium); % μg/(cm2·s)
     
    %方法二
   
    if NoTime(KT) <= 0
        idx = 1;
    else
        idx = NoTime(KT) + 1;
    end
    if idx > numel(Tm)
        % Tm/u_2 can be daily forcing while NoTime is hourly. Convert to a
        % daily index for short diagnostic runs with only a few forcing days.
        idx = floor(max(NoTime(KT) - 1, 0) / 24) + 1;
    end
    idx = min(max(idx, 1), numel(Tm));

    KH = (0.2138/(Tm(idx)+273.15))*10^(6.123-1825/(Tm(idx)+273.15));
    hm = 0.000612*u_2(idx)^0.8*(Tm(idx)+273.15)^0.382*L^-0.2;
    % [2026-05-17] pH 取表层 PH1(1) (= 0-0.25cm), 之前用 PH1(5) (= 2cm) 跟表面挥发公式不一致
    %     Svol 是表面通量, 控制 NH4↔NH3 平衡的是表土 pH (Sherlock & Goh 1985 原文亦取顶层)
    Svol(NL,2)=hm*((KH*(N(NN))/(1+10^(-0.05+(2788/(TT(NN)+273.15))-PH1(1))))-C_ammonium);


  

