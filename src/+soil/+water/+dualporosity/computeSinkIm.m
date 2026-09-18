function SinkIm = computeSinkIm(SoilVariables, VanGenuchten, cfg)
% computeSinkIm - 计算 dual-porosity mobile→immobile 单位时间体积通量
%
% HYDRUS WATFLOW.f line 1037-1040 (iDualPor=1, water content driven):
%   Se    = (θ_m  - θ_r   ) / (θ_s    - θ_r   )
%   Se_im = (θ_im - θ_r_im) / (θ_s_im - θ_r_im)
%   SinkIm = ω · (Se - Se_im)        [单位 1/s, 体积分数/秒]
%
% 输出:
%   SinkIm (NL × nD)  ≥ 0 表示从 mobile 流向 immobile (变湿)
%                      < 0 表示从 immobile 流向 mobile (退水)
%
% 注:
%   该函数只算瞬时通量, 不更新状态.
%   Richards 矩阵把 -SinkIm 当作汇项加入 RHS (assembleCoefficientMatrices C9).
%   KT 末由 dualPorExchange 用最终通量更新 θ_im.

    NL = cfg.Grid.NL;
    nD = 2;
    SinkIm = zeros(NL, nD);

    use_dual = false;
    try
        use_dual = cfg.lDualPor == 1 && isfield(SoilVariables, 'Theta_im');
    catch
    end
    if ~use_dual, return; end

    omega = getCfgVal(cfg, 'OmegaExch', 0.01);
    omega_per_s = omega / 3600;

    Theta_m  = SoilVariables.Theta_LL;
    Theta_im = SoilVariables.Theta_im;
    theta_s  = VanGenuchten.Theta_s(:);
    theta_r  = VanGenuchten.Theta_r(:);

    if isfield(SoilVariables, 'theta_s_im') && ~isempty(SoilVariables.theta_s_im)
        theta_s_im = SoilVariables.theta_s_im;
        theta_r_im = SoilVariables.theta_r_im;
    else
        f_im = getCfgVal(cfg, 'fImmobile', 0.30);
        theta_s_im = f_im * ones(NL, 1);
        theta_r_im = zeros(NL, 1);
    end

    for i = 1:NL
        denom_m  = max(theta_s(i)    - theta_r(i),    1e-12);
        denom_im = max(theta_s_im(i) - theta_r_im(i), 1e-12);
        for j = 1:nD
            Se    = max(0, min(1, (Theta_m(i, j)  - theta_r(i))    / denom_m));
            Se_im = max(0, min(1, (Theta_im(i, j) - theta_r_im(i)) / denom_im));
            SinkIm(i, j) = omega_per_s * (Se - Se_im);  % 1/s
        end
    end
end

function v = getCfgVal(cfg, name, defaultValue)
    v = defaultValue;
    try
        if isprop(cfg, name) || isfield(cfg, name)
            v = cfg.(name);
        end
    catch
    end
end
