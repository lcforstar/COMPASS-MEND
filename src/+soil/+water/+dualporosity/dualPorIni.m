function [SoilVariables, VanGenuchten] = dualPorIni(SoilVariables, VanGenuchten, cfg)
% dualPorIni - 初始化 HYDRUS Šimůnek (2003) dual-porosity 状态（HYDRUS InitW 标准）
%
% =====================================================================
% HYDRUS-1D 双孔模型的语义约定（INPUT.f line 472-497, 508-532）
% ---------------------------------------------------------------------
% ParD(2, M) = θs_mobile     ← 用户在 SELECTOR.IN 里给的 θs，**直接是 mobile 域容量**
% ParD(8, M) = θs_immobile   ← 用户在 SELECTOR.IN 里**额外**给的 immobile 容量
% Theta(i)   = θ_obs         ← 用户在 PROFILE.DAT 里给的「实测总水量」（含 mobile + immobile）
%
% 关键点：HYDRUS 不砍 mobile θs。用户的 VG.θs 直接是 mobile 容量，
% immobile 容量是额外加在上面的「小空隙」。
%   TDR 实测理论上限 = θs_m + θs_im
%
% HYDRUS InitW (line 482-484) 拆分：
%   θ_m  = θ_obs × ParD(2) / (ParD(2) + ParD(8))
%        = θ_obs × θs_m / (θs_m + θs_im)
%   θ_im = θ_obs − θ_m       ← 简单减法，保证 θ_m + θ_im = θ_obs
%
% =====================================================================
% COMPASS 适配（v25, 2026-05-08，重写自 v16）
% ---------------------------------------------------------------------
% 之前 v16 的做法：θs_m_eff = VG.θs − ts_im（**砍 mobile 容量**），
% 导致冬季 θ_m 接近 mobile_eff 上限时 Se→1 → Picard 雅可比奇异 → 全年跑 KT=1
% Water 残差爆 1.5e+10。
%
% v25 修复：严格按 HYDRUS：
%   - VanGenuchten.Theta_s 不动（保持用户输入值）
%   - 新增 SoilVariables.theta_s_im (= ts_im) 给 dualPorExchange 用
%   - 用户初值 θ_obs 按 HYDRUS InitW 比例拆 mobile + immobile
%   - hh 不需要重算（VG 参数完全不变，θ_m 与 θ_obs 同 Se 比例）
%   - 移除 VG.Theta_s_total 备份字段（不再需要）
%
% 修改者: Claude / 朱寄子星

    if ~hasFlag(cfg, 'lDualPor', 0)
        return;
    end

    NL = cfg.Grid.NL;
    nD = 2;

    % --- 1. immobile 域容量参数 ----------------------------------------
    %   ts_im_val = ParD(8) HYDRUS 语义
    %   tr_im_val = ParD(7) HYDRUS 语义（默认 0）
    ts_im_val = getCfgVal(cfg, 'ThetaSImmobile', 0.05);
    tr_im_val = getCfgVal(cfg, 'ThetaRImmobile', 0.0);
    theta_s_im = ts_im_val * ones(NL, 1);
    theta_r_im = tr_im_val * ones(NL, 1);

    % --- 2. VG 参数完全不动 -------------------------------------------
    % HYDRUS 把用户输入的 VG.Theta_s 直接当 mobile 域容量，无需修改。
    theta_s_mobile = VanGenuchten.Theta_s(:);   % 保持原值（用户输入）
    theta_r_mobile = VanGenuchten.Theta_r(:);

    % 物理合理性检查：mobile θs > 0
    if any(theta_s_mobile <= 0)
        error('dualPorIni:InvalidThetaS', ...
            'VanGenuchten.Theta_s 必须 > 0，发现 %d 层异常', sum(theta_s_mobile <= 0));
    end

    % --- 3. 按 HYDRUS InitW 拆 θ_obs → θ_m + θ_im ---------------------
    % HYDRUS INPUT.f line 482-484:
    %   ThMobile = θ_obs × ParD(2) / (ParD(2) + ParD(8))
    %   ThNewIm  = θ_obs − ThMobile
    Theta_obs = SoilVariables.Theta_LL;       % NL × nD, 用户初值（实测 TDR）
    Theta_m   = zeros(NL, nD);
    Theta_im  = zeros(NL, nD);
    for i = 1:NL
        ts_m  = theta_s_mobile(i);
        ts_im = theta_s_im(i);
        denom = max(ts_m + ts_im, 1e-12);
        for j = 1:nD
            tobs = Theta_obs(i, j);
            % 钉到 [0, ts_m + ts_im] 上限（理论 TDR 最大）
            tobs = max(0, min(ts_m + ts_im, tobs));
            Theta_m(i, j)  = tobs * ts_m / denom;
            Theta_im(i, j) = tobs - Theta_m(i, j);   % 守恒：θ_m + θ_im = θ_obs
        end
    end

    % --- 4. 写入 SoilVariables -----------------------------------------
    % VG 参数完全不动，只更新 Theta_LL 与 Theta_im
    SoilVariables.Theta_LL = Theta_m;
    SoilVariables.Theta_UU = Theta_m;
    SoilVariables.Theta_L  = Theta_m;
    SoilVariables.Theta_U  = Theta_m;

    SoilVariables.Theta_im       = Theta_im;
    SoilVariables.Theta_im_old   = Theta_im;
    SoilVariables.Theta_m_old    = Theta_m;
    SoilVariables.theta_r_im     = theta_r_im;
    SoilVariables.theta_s_im     = theta_s_im;
    SoilVariables.SinkIm_current = zeros(NL, nD);

    % hh 不需要重算：VG 参数不变，θ_m 与 θ_obs 的 Se 比值相同。
    % HYDRUS InitW 在 line 487-493 也只用 Kappa 选择曲线，不重算 hh。
end

function tf = hasFlag(cfg, name, defaultValue)
    tf = defaultValue;
    try
        if isprop(cfg, name) || isfield(cfg, name)
            tf = cfg.(name);
        end
    catch
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
