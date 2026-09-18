function [VanGenuchten, SoilProperties] = setVanGenuchtenParameters(SoilProperties)
%setVanGenuchtenParameters  装配 VG 参数 + 校验 POR ≥ VG.Theta_s + 干空气余量
%
% 参考：Van Genuchten MTh, Leij FJ, Yates SR (1991) The RETC code for
%       quantifying the hydraulic functions of unsaturated soils. EPA/600/2-91/065.
%
% 输入:
%   SoilProperties - 含 Theta_s/Theta_r/Field_mc/Alpha/n/m/Porosity 字段
%
% 输出:
%   VanGenuchten   - VG 参数结构体
%   SoilProperties - 若 Porosity 与 Theta_s 不自洽则被强制抬升后回写
%
% 物理约束：Porosity ≥ VG.Theta_s + 干空气余量
%   - POR  : 几何孔隙度 = 1 - bd/2.65（含干空气）
%   - VG.θs: 饱和含水量上限（液态水容量）
%   - 差额 ≥ trapped_air_margin 才能保证 Theta_V = POR - θL - θII > 0
%   - 双孔模式 (lDualPor=1) 下 dualPorIni 后续会改写 VG.Theta_s 为 mobile-eff，
%     此处仍按用户输入的原始 VG.Theta_s 校验是合理的（dualPorIni 后约束更宽松）。
%
% 用户输入优先：
%   - 若 POR(来自 Excel bd) >= VG.θs + margin → 不动，使用用户值
%   - 否则强制抬升 POR = VG.θs + margin，并 warning 通知用户检查 Excel 容重列
%
% 修改者：朱寄子星 / Claude
% 修改日期：2026-05-08

    VanGenuchten.Theta_s = SoilProperties.Theta_s;
    VanGenuchten.Theta_r = SoilProperties.Theta_r;
    VanGenuchten.Theta_f = SoilProperties.Field_mc;
    VanGenuchten.Alpha   = SoilProperties.Alpha;
    VanGenuchten.n       = SoilProperties.n;
    VanGenuchten.m       = SoilProperties.m;

    % ----- POR vs VG.Theta_s 自洽性强约束 -----
    if isfield(SoilProperties, 'Porosity') && ~isempty(SoilProperties.Porosity)
        trapped_air_margin = 0.03;   % 3% 干空气余量
        POR_required = VanGenuchten.Theta_s(:) + trapped_air_margin;
        POR_user     = SoilProperties.Porosity(:);
        n_bad = sum(POR_user < POR_required);
        if n_bad > 0
            diff_vec = POR_required - POR_user;
            warning('init:PorosityVGConflict', ...
                ['用户输入的 Porosity (来自 Excel 容重) 在 %d / %d 层 < VG.Theta_s + %.2f\n', ...
                 '  最大缺口: %.4f cm³/cm³\n', ...
                 '  原始 POR 范围: [%.4f, %.4f]\n', ...
                 '  自动抬升 POR = max(POR_user, VG.Theta_s + %.2f)，\n', ...
                 '  请检查 InitialConditions.xlsx 容重列与 soilproperty.xlsx theta_s 是否一致。'], ...
                n_bad, length(POR_user), trapped_air_margin, max(diff_vec), ...
                min(POR_user), max(POR_user), trapped_air_margin);
            SoilProperties.Porosity = max(POR_user, POR_required);
        end
    end

end
