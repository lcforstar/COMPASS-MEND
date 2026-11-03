function [lE, H, ec, Cc, lambda, s, lE_dry, lE_wet] = heatFluxesClm45(ra_h, rs, Tc, ea, Ta, e_to_q, PSI, Ca, Ci, ra_w, rhoa, fwet)

    %{
     CLM4.5-style heat and water-vapor flux helper for ebal_clm45.

     Scientific role:
       This function keeps the original SCOPE heat-flux structure, but allows
       the aerodynamic resistance for sensible heat (ra_h) and the aerodynamic
       resistance for water vapor/CO2 coupling (ra_w) to differ. This matches
       the ebal_clm45 use of CLM4.5-style scalar roughness and separate rah/raw
       pathways while preserving the existing canopy-soil energy-balance logic.

     CLM4.5-style improvement:
       Jizixing Zhu (朱寄子星), 18 Apr 2026.

     Main updates:
       - Added optional ra_w input; if omitted, ra_w = ra_h for backward
         compatibility with the original heatfluxes call.
       - Added optional rhoa input so ebal_clm45 can use forcing-state air
         density, following the CLM-style forc_rho turbulent conductance.
       - Upgraded saturation vapor pressure and derivative from Magnus/s_fun
         style approximations to Calc_Esat_Flatau.

     Method document:
       docs/clm45_aerodynamic_resistance_method.md

     Editor: Jizixing Zhu; Date: 2026-07-20.
    %}

    % load physical constants
    phys = equations.physicalConstants();

    rhoa_std = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));
    if nargin < 11 || isempty(rhoa)
        rhoa = rhoa_std;
    else
        rhoa_bad = ~isfinite(rhoa) | rhoa <= 0;
        if any(rhoa_bad(:))
            rhoa(rhoa_bad) = rhoa_std;
        end
    end
    cp   = phys.Cp_air_SI;
    MH2O = phys.MW_wv;
    R    = phys.R;

    % author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
    % date:     7 Dec 2007
    % updated: 15 Apr 2009 CvdT     changed layout
    % updated: 14 Sep 2012 CvdT     added ec and Cc to output
    %
    % parent: ebal.m
    %
    % usage:
    % function [lE, H]  = heatfluxes(ra,rs,Tc,ea,Ta,e_to_q,PSI)
    %
    % this function calculates latent and sensible heat flux
    %
    % input:
    %   ra          aerodynamic resistance for heat         s m-1
    %   rs          stomatal resistance                     s m-1
    %   Tc          leaf temperature                        oC
    %   ea          vapour pressure above canopy            hPa
    %   Ta          air temperature above canopy            oC
    %   e_to_q      conv. from vapour pressure to abs hum   hPa-1
    %   PSI         leaf water potential                    J kg-1
    %   Ca          ambient CO2 concentration               umol m-3
    %   Ci          intercellular CO2 concentration         umol m-3
    %
    % output:
    %   lEc         latent heat flux of a leaf              W m-2
    %   Hc          sensible heat flux of a leaf            W m-2
    %   ec          vapour pressure at the leaf surface     hPa
    %   Cc          CO2 concentration at the leaf surface   umol m-3
    if nargin < 10 || isempty(ra_w)
        ra_w = ra_h;
    end
    % --- [v29 2026-05-09] CLM5 fwet 物理：干湿叶面分别算 lE ---
    % fwet ∈ [0,1]：湿叶面比例（从 COMPASS.m 主循环传入，h2ocan 水库决定）
    % 默认 fwet=0 兼容旧调用（如 soil flux、未升级的 ebal 等）
    if nargin < 12 || isempty(fwet)
        fwet = 0;
    end
    fwet = max(0, min(1, fwet));   % 保护：[0,1]

    lambda      = (2.501 - 0.002361 * Tc) * 1E6;  %      [J kg-1]  Evapor. heat (J kg-1)
    % --- [改进] 饱和 VP 升级: Magnus → Flatau；s 使用 Flatau 解析导数替代 s_fun 近似 ---
    %            Saturation VP upgraded: Magnus → Flatau (CLM5 high-precision, ice/water phase split)
    %            s from analytical des/dT, consistent with ei (accuracy ±0.02% vs ±0.3%)
    %   改进者 / Modified by:   朱寄子星 (Jizixing Zhu)
    %   改进日期 / Modified date: 2026-04-18
    [es_Pa, des_dT_Pa] = equations.Calc_Esat_Flatau(Tc, 101325);   % 只需 es/des
    es_hPa      = reshape(es_Pa, size(Tc)) * 0.01;                                   % Pa → hPa
    ei          = es_hPa .* exp(1E-3 * PSI * MH2O / R ./ (Tc + 273.15));
    s           = reshape(des_dT_Pa, size(Tc)) * 0.01;                               % 解析导数 / analytical derivative [hPa/K]
    qi          = ei .* e_to_q;
    qa          = ea .* e_to_q;

    % --- [v29] CLM5 风格干湿叶面拆分（CLM5 Tech Note Eq 5.97）---
    %   lE_dry = (1-fwet) × ρ × λ × (qi-qa) / (ra_w + rs)   ← 干叶气孔蒸腾
    %   lE_wet =    fwet  × ρ × λ × (qi-qa) /  ra_w         ← 湿叶 Penman（无 rs）
    %   lE     = lE_dry + lE_wet
    common_vap = rhoa .* lambda .* (qi - qa);                           % [W m-2 × s/m]
    lE_dry     = (1 - fwet) .* common_vap ./ (ra_w + rs);               % 干叶部分
    lE_wet     = fwet       .* common_vap ./ ra_w;                      % 湿叶部分（rs=0）
    lE         = lE_dry + lE_wet;                                        % [W m-2] 总潜热
    H          = (rhoa .* cp) ./ ra_h .* (Tc - Ta);                     % [W m-2] 感热（不分干湿）
    % ec / Cc 仍用气孔阻力路径计算（叶面 vapor/CO2 浓度，跟原物理一致）
    ec         = ea + (ei - ea) .* ra_w ./ (ra_w + rs);
    Cc         = Ca - (Ca - Ci) .* ra_w ./ (ra_w + rs);
end
