function SoilVariables = dualPorExchange(SoilVariables, VanGenuchten, cfg, Delt_t)
% dualPorExchange - HYDRUS Šimůnek (2003) dual-porosity 一阶交换 (iDualPor=1)
%
% 严格对照 HYDRUS-1D WATFLOW.f line 1027-1075 (subroutine DualPor):
%   - 用 ThOld（步初 mobile θ_m）和 ThOldIm（步初 immobile θ_im）算 SinkIm
%   - SinkIm 在整个 Picard 迭代里保持常数（确保 Richards RHS 扣的水量
%     与 immobile 加的水量完全一致 → 守恒）
%   - 立即更新 ThNewIm = ThOldIm + SinkIm·Δt（截断到 [θr_im, θs_im]）
%
% 调用时机：每个 substep 开始时调用一次（Picard 循环之前）。
% 之后 solveSoilMoistureBalance 直接读 SoilVariables.SinkIm_current 进 RHS。
%
% 修改者: Claude / 朱寄子星

    if ~hasFlag(cfg, 'lDualPor', 0)
        return;
    end

    NL = cfg.Grid.NL;
    nD = 2;
    omega = getCfgVal(cfg, 'OmegaExch', 0.01);   % 1/h
    omega_per_s = omega / 3600;                    % 1/s

    if nargin < 4 || isempty(Delt_t) || Delt_t <= 0
        Delt_t = 3600;
    end

    % 关键修复：用步初值 Theta_L（不是 Theta_LL，后者是 Picard 当前迭代值）
    % HYDRUS WATFLOW.f line 1039: Se = (ThOld - thr) / (ths - thr)
    if isfield(SoilVariables, 'Theta_L') && ~isempty(SoilVariables.Theta_L) ...
            && size(SoilVariables.Theta_L,1) == NL
        Theta_m = SoilVariables.Theta_L;
    else
        Theta_m = SoilVariables.Theta_LL;
    end
    Theta_im = SoilVariables.Theta_im;        % ThOldIm（步初 immobile）
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

    Theta_im_new = Theta_im;
    SinkIm_field = zeros(NL, nD);

    for i = 1:NL
        ts = theta_s(i);
        tr = theta_r(i);
        ts_im = theta_s_im(i);
        tr_im = theta_r_im(i);
        denom_m  = max(ts - tr, 1e-12);
        denom_im = max(ts_im - tr_im, 1e-12);

        for j = 1:nD
            tm  = Theta_m(i, j);
            tim = Theta_im(i, j);

            Se    = max(0, min(1, (tm  - tr   ) / denom_m));
            Se_im = max(0, min(1, (tim - tr_im) / denom_im));

            % HYDRUS line 1040: SinkIm = ω·(Se - Se_im)
            SinkIm = omega_per_s * (Se - Se_im);

            % HYDRUS line 1041-1042: 限幅 DeltaTh
            DeltaTh = (Se - Se_im) / (denom_m + denom_im) * denom_m * denom_im;

            % HYDRUS line 1057-1066: TrMax 限幅
            if SinkIm > 0  % 流向 immobile
                TrMaxIm = (ts_im - tim) / Delt_t;
                TrMaxIm = min(DeltaTh / Delt_t, TrMaxIm);
                if SinkIm > TrMaxIm, SinkIm = TrMaxIm; end
            elseif SinkIm < 0  % 流向 mobile
                TrMaxIm = -(ts - tm) / Delt_t;
                TrMaxIm = max(DeltaTh / Delt_t, TrMaxIm);
                if SinkIm < TrMaxIm, SinkIm = TrMaxIm; end
            end

            % HYDRUS line 1068: 更新 immobile 含水量
            tim_new = tim + SinkIm * Delt_t;
            tim_new = max(tr_im, min(ts_im, tim_new));

            Theta_im_new(i, j) = tim_new;
            SinkIm_field(i, j) = SinkIm;
        end
    end

    SoilVariables.Theta_im_old   = Theta_im;
    SoilVariables.Theta_m_old    = Theta_m;
    SoilVariables.Theta_im       = Theta_im_new;
    % 关键：缓存 SinkIm 给 Richards RHS 使用，确保两边用同一个值（守恒）
    SoilVariables.SinkIm_current = SinkIm_field;
    % NOT 修改 Theta_LL/Theta_UU: Richards 已通过 SinkIm 汇项自动扣除
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
