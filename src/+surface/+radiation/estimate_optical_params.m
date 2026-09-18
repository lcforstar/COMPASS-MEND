function [omega, beta, beta0] = estimate_optical_params(wavelength_band)
%ESTIMATE_OPTICAL_PARAMS 估算光学参数
%   根据波段估算散射系数和向上散射参数
%
%   输入:
%     wavelength_band: 'visible' 或 'nir'（近红外）
%
%   输出:
%     omega: 散射系数
%     beta: 漫射光向上散射参数
%     beta0: 直射光向上散射参数

    switch wavelength_band
        case 'visible'  % 可见光波段 (< 0.7μm)
            omega = 0.2;   % 低散射（强吸收）
            beta = 0.6;   % 较多的向上散射
            beta0 = 0.7;  % 直射光更多向上散射
        case 'nir'      % 近红外波段 (≥ 0.7μm)
            omega = 0.8;   % 高散射（弱吸收）
            beta = 0.4;   % 较少的向上散射
            beta0 = 0.3;  % 直射光较少向上散射
        otherwise
            error('不支持的波段类型');
    end
end

