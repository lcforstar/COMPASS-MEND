function [ETCON, EHCAP, TETCON, EfTCON, ZETA] = calculateSoilThermalProperites(SoilVariables, ThermalConductivity, VanGenuchten, cfg, L)
%CALCULATESOILTHERMALPROPERITES  Compute multiphase soil thermal properties.
%
%   Inputs provide liquid, gas, and actual ice volume fractions [-], node
%   temperature [degC], gas densities [g cm^-3], soil composition, and the
%   selected thermal-conductivity method. L is vaporization heat [J g^-1].
%
%   ETCON, TETCON, and EfTCON are alternative effective conductivities;
%   EHCAP is volumetric sensible heat capacity [J cm^-3 K^-1], and ZETA is
%   the dimensionless pore-temperature-gradient ratio.
%
%   Actual ice volume is used in phase-weighted conductivity and heat
%   capacity. The density-weighted liquid equivalent belongs only to the
%   conserved water-storage definition and is not a thermal volume fraction.
%
%   Reference: Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.43-4.44,
%   p. 75 (PDF p. 89); project three-phase extension uses theta_i directly.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- use actual ice and gas phase volumes
%           in thermal conductivity and sensible heat capacity.

NL        = cfg.Grid.NL;
nD        = cfg.Grid.nD;
NN        = cfg.Grid.NN;
ThermCond = cfg.lThermCond;

XWILT    = SoilVariables.XWILT;
XK       = SoilVariables.XK;
TT       = SoilVariables.TT;
POR      = SoilVariables.POR;
Theta_LL = SoilVariables.Theta_LL;
Theta_V  = SoilVariables.Theta_V;
Theta_II = SoilVariables.Theta_II;
if isfield(SoilVariables, 'XSOC')
    XSOC = SoilVariables.XSOC;
else
    XSOC = zeros(NL, 1);
end
DRHOVT   = SoilVariables.DRHOVT;
RHOV     = SoilVariables.RHOV;
RHODA    = SoilVariables.RHODA;

Theta_s = VanGenuchten.Theta_s;

HCAP  = ThermalConductivity.HCAP;
SF    = ThermalConductivity.SF;
TCA   = ThermalConductivity.TCA;
GA1   = ThermalConductivity.GA1;
GA2   = ThermalConductivity.GA2;
GB1   = ThermalConductivity.GB1;
GB2   = ThermalConductivity.GB2;
HCD   = ThermalConductivity.HCD;
ZETA0 = ThermalConductivity.ZETA0;
CON0  = ThermalConductivity.CON0;
PS1   = ThermalConductivity.PS1;
PS2   = ThermalConductivity.PS2;

phys = equations.physicalConstants();
c_a  = phys.c_a;
c_V  = equations.getPhysicalConstant('c_V');

ETCON  = zeros(NL, nD);
EHCAP  = zeros(NL, nD);
TETCON = zeros(NL, nD);
EfTCON = zeros(NL, nD);
ZETA   = zeros(NL, nD);
D_A    = zeros(1, NN);

TCON    = zeros(1, 7);
TCON(1) = 1.37e-3 * 4.182;

for i = 1:NL
    for j = 1:nD
        MN  = i + j - 1;
        XXX = Theta_LL(i, j);
        % Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.43-4.44, p. 75.
        % Every phase weight is its actual volumetric fraction.
        XII = Theta_II(i, j);
        if Theta_LL(i, j) < XK(i)
            XXX = XK(i);
        end

        % --- shape factor SF(2) ---
        if XXX < XWILT(i)
            SF(2) = GA2 + GB2(i) * XXX;
        else
            SF(2) = GA1 + GB1(i) * XXX;
        end

        D_A(MN) = 0.229 * (1 + TT(MN) / 273)^1.75;
        TCON(2) = TCA + D_A(MN) * L(MN) * DRHOVT(MN);
        TCON(6) = 5.2e-3 * 4.182;  % ice
        SF(6)   = 0.125;

        TARG(2) = TCON(2) / TCON(1) - 1;
        TARG(6) = TCON(6) / TCON(1) - 1;
        GRAT(2) = 0.667 / (1 + TARG(2) * SF(2)) + 0.333 / (1 + TARG(2) * (1 - 2 * SF(2)));
        GRAT(6) = 0.667 / (1 + TARG(6) * SF(6)) + 0.333 / (1 + TARG(6) * (1 - 2 * SF(6)));

        % --- de Vries ETCON / ZETA ---
        ETCON(i,j) = (PS1(i) + XXX*TCON(1) + (POR(i)-XXX-XII)*GRAT(2)*TCON(2) + XII*GRAT(6)*TCON(6)) / ...
                     (PS2(i) + XXX + (POR(i)-XXX-XII)*GRAT(2) + XII*GRAT(6));
        ZETA(i,j)  = GRAT(2) / (GRAT(2)*(POR(i)-XXX-XII) + XXX + PS2(i) + XII*GRAT(6));

        % --- EHCAP ---
        % Zeng (2013), Coupled Dynamics in Soil, Eq. 4.43, p. 75.
        EHCAP(i,j) = HCD(i) + HCAP(1) * Theta_LL(i, j);
        EHCAP(i,j) = EHCAP(i,j) ...
            + (RHODA(MN) * c_a + RHOV(MN) * c_V) * Theta_V(i, j) ...
            + HCAP(6) * XII;
        if Theta_LL(i, j) ~= XXX
            ZETA(i,j)  = ZETA0(i) + (ZETA(i,j)  - ZETA0(i)) * Theta_LL(i,j) / XXX;
            ETCON(i,j) = CON0(i)  + (ETCON(i,j) - CON0(i))  * Theta_LL(i,j) / XXX;
        end
        if TT(MN) < 0
            ZETA(i,j) = 1;
        end

        % Tian (2016) simplified de Vries method.
        if isfield(ThermalConductivity, 'TPS1')
            TPS1  = ThermalConductivity.TPS1;
            TPS2  = ThermalConductivity.TPS2;
            TCON0 = ThermalConductivity.TCON0;
            TSF(2) = 0.333 * (1 - (POR(i)-XXX-XII) / POR(i));
            TSF(6) = 0.333 * (1 - XII / POR(i));
            TSF(7) = 0.5;
            TCON(7) = 0.25e-2;
            TARG(7) = TCON(7) / TCON(1) - 1;
            TTGRAT(2) = 0.667/(1+TARG(2)*TSF(2)) + 0.333/(1+TARG(2)*(1-2*TSF(2)));
            TTGRAT(6) = 0.667/(1+TARG(6)*TSF(6)) + 0.333/(1+TARG(6)*(1-2*TSF(6)));
            TTGRAT(7) = 0.667/(1+TARG(7)*TSF(7)) + 0.333/(1+TARG(7)*(1-2*TSF(7)));
            if XSOC(i) == 0
                TETCON(i,j) = (TPS1(i)+XXX*TCON(1)+(POR(i)-XXX-XII)*TTGRAT(2)*TCON(2)+XII*TTGRAT(6)*TCON(6)) / ...
                              (TPS2(i)+XXX+(POR(i)-XXX-XII)*TTGRAT(2)+XII*TTGRAT(6));
            else
                TETCON(i,j) = (TPS1(i)+XXX*TCON(1)+(POR(i)-XXX-XII)*TTGRAT(2)*TCON(2)+XII*TTGRAT(6)*TCON(6)+XSOC(i)*TTGRAT(7)*TCON(7)) / ...
                              (TPS2(i)+XXX+(POR(i)-XXX-XII)*TTGRAT(2)+XII*TTGRAT(6)+XSOC(i)*TTGRAT(7));
            end
            if Theta_LL(i,j) ~= XXX
                TETCON(i,j) = TCON0(i) + (TETCON(i,j) - TCON0(i)) * Theta_LL(i,j) / XXX;
            end
        end

        % Farouki heat capacity for lThermCond == 4.
        if ThermCond == 4 && isfield(ThermalConductivity, 'FEHCAP')
            FEHCAP = ThermalConductivity.FEHCAP;
            EHCAP(i,j) = (FEHCAP(i)*(1-POR(i)) + 1000*4217.7*Theta_LL(i,j) + 917*2117.27*XII) * 1e-6;
        end

        % Peters-Lidard/Johansen conductivity for methods 2 and 4.
        if (ThermCond == 2 || ThermCond == 4) && isfield(ThermalConductivity, 'TCON_dry')
            TCON_dry = ThermalConductivity.TCON_dry;
            TCON_s_i = ThermalConductivity.TCON_s(i);
            TCON_i_i = 2.2;
            TCON_L   = 0.57;
            if Theta_II(i,j) <= 0
                if Theta_LL(i,j) / POR(i) > 0.1
                    K_e = log10(Theta_LL(i,j) / POR(i)) + 1;
                elseif Theta_LL(i,j) / POR(i) > 0.05
                    K_e = 0.7 * log10(Theta_LL(i,j) / POR(i)) + 1;
                else
                    K_e = 0;
                end
                TCON_sat = TCON_s_i^(1-Theta_s(i)) * TCON_L^Theta_s(i);
                EfTCON(i,j) = K_e * (TCON_sat - TCON_dry(i)) + TCON_dry(i);
            else
                K_e = (Theta_LL(i,j) + Theta_II(i,j)) / POR(i);
                TCON_sat = TCON_s_i^(1-Theta_s(i)) * TCON_i_i^(Theta_s(i)-Theta_LL(i,j)) * TCON_L^Theta_LL(i,j);
                EfTCON(i,j) = K_e * (TCON_sat - TCON_dry(i)) + TCON_dry(i);
            end
        end
    end
end
end
