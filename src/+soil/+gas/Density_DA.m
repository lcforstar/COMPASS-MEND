function SoilVariables = Density_DA(SoilVariables, cfg, Delt_t)
% Density_DA - 干空气密度及其偏导数 / Dry air density and partial derivatives
%
% 输入 / Inputs:
%   SoilVariables - 包含 T, h, hh, TT, P_g, P_gg, RHOV, DRHOVh, DRHOVT
%   cfg           - ModelSettings，含 Grid.NL/NN/DeltZ
%   Delt_t        - 当前时间步长 (s)
%
% 输出 / Output:
%   SoilVariables.Xaa, .XaT, .Xah, .DRHODAt, .DRHODAz, .RHODA

NL    = cfg.Grid.NL;
NN    = cfg.Grid.NN;
DeltZ = cfg.Grid.DeltZ;
phys  = equations.physicalConstants();
RDA   = phys.RDA;
Rv    = phys.Rv;

T      = SoilVariables.T;
h      = SoilVariables.h;
hh     = SoilVariables.hh;
TT     = SoilVariables.TT;
P_g    = SoilVariables.P_g;
P_gg   = SoilVariables.P_gg;
RHOV   = SoilVariables.RHOV;
DRHOVh = SoilVariables.DRHOVh;
DRHOVT = SoilVariables.DRHOVT;

if cfg.lDEBUG_NAN_TRACE
    fprintf('    [Density_DA] Input:\n');
    fprintf('      P_g NaN=%d, P_gg NaN=%d, TT NaN=%d\n', sum(isnan(P_g(:))), sum(isnan(P_gg(:))), sum(isnan(TT(:))));
    fprintf('      RHOV NaN=%d, DRHOVh NaN=%d, DRHOVT NaN=%d\n', sum(isnan(RHOV(:))), sum(isnan(DRHOVh(:))), sum(isnan(DRHOVT(:))));
end

Xaa     = zeros(1, NN);
XaT     = zeros(1, NN);
Xah     = zeros(1, NN);
DRHODAt = zeros(1, NN);
RHODA   = zeros(1, NN);
DRHODAz = zeros(1, NL);

for MN = 1:NN
    T_K = TT(MN) + 273.15;
    Xaa(MN) = 1 / (RDA * T_K);
    XaT(MN) = -(P_gg(MN) / (RDA * T_K^2) + Rv * DRHOVT(MN) / RDA);
    Xah(MN) = -Rv * DRHOVh(MN) / RDA;

    if cfg.lDEBUG_NAN_TRACE
        if isnan(Xaa(MN)) || isnan(XaT(MN)) || isnan(Xah(MN))
            fprintf('      [NaN@%d] Xaa/XaT/Xah. TT=%.2f, P_gg=%.0f, DRHOVT=%.2e\n', MN, TT(MN), P_gg(MN), DRHOVT(MN));
        end
    end

    DRHODAt(MN) = Xaa(MN)*(P_gg(MN)-P_g(MN))/Delt_t + XaT(MN)*(TT(MN)-T(MN))/Delt_t + Xah(MN)*(hh(MN)-h(MN))/Delt_t;
    RHODA(MN)   = P_gg(MN) / (RDA * T_K) - RHOV(MN) * Rv / RDA;

    if cfg.lDEBUG_NAN_TRACE
        if isnan(RHODA(MN))
            fprintf('      [NaN@%d] RHODA. P_gg/RDA=%.2e, RHOV=%.2e\n', MN, P_gg(MN)/RDA, RHOV(MN));
        end
    end
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('    [Density_DA] Output: Xaa NaN=%d, Xah NaN=%d, RHODA NaN=%d\n', sum(isnan(Xaa(:))), sum(isnan(Xah(:))), sum(isnan(RHODA(:))));
end

for ML = 1:NL
    XaaBAR = (Xaa(ML+1) + Xaa(ML)) / 2;
    XahBAR = (Xah(ML+1) + Xah(ML)) / 2;
    XaTBAR = (XaT(ML+1) + XaT(ML)) / 2;
    DRHODAz(ML) = XaaBAR*(P_gg(ML+1)-P_gg(ML))/DeltZ(ML) + XaTBAR*(TT(ML+1)-TT(ML))/DeltZ(ML) + XahBAR*(hh(ML+1)-hh(ML))/DeltZ(ML);

    if cfg.lDEBUG_NAN_TRACE
        if isnan(DRHODAz(ML))
            fprintf('      [NaN@%d] DRHODAz. Xaa/Xah/XaT=%.2e/%.2e/%.2e\n', ML, XaaBAR, XahBAR, XaTBAR);
        end
    end
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('    [Density_DA] DRHODAz NaN=%d\n', sum(isnan(DRHODAz(:))));
end

SoilVariables.Xaa     = Xaa;
SoilVariables.XaT     = XaT;
SoilVariables.Xah     = Xah;
SoilVariables.DRHODAt = DRHODAt;
SoilVariables.DRHODAz = DRHODAz;
SoilVariables.RHODA   = RHODA;
end
