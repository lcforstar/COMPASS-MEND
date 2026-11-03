function r_ss = Calculate_Soil_Resistance(Theta_surf, Theta_r, Theta_s, Theta_f, rss_scale)
%CALCULATE_SOIL_RESISTANCE Field-capacity-normalized surface resistance.
% Inputs: surface/residual/saturated/field-capacity water fractions [m3/m3],
% scalar or conformable arrays; optional positive finite scalar multiplier.
% Theta_s is retained for the existing four-argument API. Output r_ss [s/m].
% Moisture weight is bounded to [0,1]; constitutive resistance to [50,4000].
% Requested Sellers-style coefficients [SOURCE UNVERIFIED 2026-10-17].
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- activate requested field-capacity law and scale.
% [改进 2026-04-18 朱寄子星] 从 +CalculateEvapotranspiration/ 迁移至 +equations/
%   按物理方程归类，与 zo_and_d.m 等阻力类方程并置。
%   调用改为 equations.Calculate_Soil_Resistance(...)
%
% CLM4.5-style improvement:
%   Jizixing Zhu (朱寄子星), 18 Apr 2026.
%
% 方法: S4 - Mahfouf and Noilhan (1991) / Sellers-style soil evaporation resistance
%   ebal_clm45 每个时间步用表层土壤含水量 Theta_surf 动态更新 rss，
%   用于替代原始固定 soil.rss。该公式使干土的水汽扩散阻力快速增大，
%   湿土的蒸发阻力降低，从而把土壤水分状态反馈到潜热通量。
%
% Method document:
%   docs/clm45_aerodynamic_resistance_method.md
    if nargin<5 || isempty(rss_scale), rss_scale=1.0; end
    validateattributes(rss_scale,{'numeric'},{'real','scalar','positive','finite'});
    w=min(1.0,max(0.0,(Theta_surf-Theta_r)./max(Theta_f-Theta_r,1e-6)));
    r_ss_raw=exp(8.206-4.255.*w);
    r_ss=min(4000,max(50,rss_scale.*r_ss_raw));
end
