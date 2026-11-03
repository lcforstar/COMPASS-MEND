function TransportCoefficient = calculateVaporVariables(SoilVariables, SoilProperties, VanGenuchten, TransportCoefficient, cfg)
%CALCULATEVAPORVARIABLES  Compute pore-gas diffusivity and thermal enhancement.
%
%   Inputs contain liquid and gas volume fractions [-], porosity [-], the
%   liquid-continuity threshold XK [-], temperature [degC], and the de Vries
%   pore-temperature-gradient ratio ZETA [-]. The returned D_V and D_A are
%   [cm2 s^-1], and Eta is dimensionless.
%
%   The Philip-de Vries liquid-island cross section is theta_g+theta_l below
%   the liquid-continuity threshold. Above the threshold, its liquid part is
%   reduced by f(theta_g)=theta_g/(porosity-theta_k). Actual liquid volume is
%   used; water-equivalent storage is not a vapor-transfer cross section.
%
%   Reference: Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.29-4.40,
%   pp. 71-74 (PDF pp. 85-88), especially Eqs. 4.35-4.36.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- evaluate the liquid-island area from
%           actual liquid and gas volumes and honor DVT_Switch safely.

NL           = cfg.Grid.NL;
nD           = cfg.Grid.nD;
NN           = cfg.Grid.NN;
ThmrlCondCap = cfg.lThmrlCondCap;

Theta_LL = SoilVariables.Theta_LL;   % NL×nD
Theta_g  = SoilVariables.Theta_g;    % NL×nD
POR      = SoilVariables.POR;        % NL×1
XK       = SoilVariables.XK;         % NL×1
TT       = SoilVariables.TT;         % 1×NN

fc      = SoilProperties.fc;         % Scalar texture parameter
Theta_s = VanGenuchten.Theta_s;      % NL×1
ZETA    = TransportCoefficient.ZETA; % NL×nD

D_V = zeros(NL, nD);
Eta = zeros(NL, nD);
D_A = zeros(1, NN);

if cfg.lDEBUG_NAN_TRACE
    fprintf('[calculateVaporVariables] Input: Theta_LL NaN=%d, Theta_g NaN=%d, ZETA NaN=%d\n', ...
        sum(isnan(Theta_LL(:))), sum(isnan(Theta_g(:))), sum(isnan(ZETA(:))));
end

for ML = 1:NL
    for ND = 1:nD
        MN = ML + ND - 1;

        if ThmrlCondCap
            thetaLiquid = Theta_LL(ML, ND);
            thetaGas = Theta_g(ML, ND);
            % Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.35-4.36,
            % p. 73: cross section = theta_g + f(theta_g)*theta_l.
            if thetaLiquid < XK(ML)
                EnhnLiqIsland = thetaGas + thetaLiquid;
            else
                EnhnLiqIsland = thetaGas ...
                    * (1 + thetaLiquid / max(POR(ML) - XK(ML), eps));
            end

            f0 = thetaGas^(7/3) / Theta_s(ML)^2;
            D_A(MN) = 0.229 * (1 + TT(MN) / 273)^1.75;

            if SoilVariables.DVT_Switch == 1 && thetaGas > 1e-2
                Eta(ML, ND) = ZETA(ML, ND) * EnhnLiqIsland ...
                    / (f0 * thetaGas);
            else
                Eta(ML, ND) = 0;
            end
        else
            f0 = Theta_g(ML, ND)^0.67;
            D_A(MN) = 1e4 * 2.12e-5 * (1 + TT(MN) / 273.15)^2;
            Eta(ML, ND) = 8 + 3 * Theta_LL(ML, ND) / Theta_s(ML) - ...
                7 * exp(-((1 + 2.6 / sqrt(fc)) * Theta_LL(ML, ND) / Theta_s(ML))^3);
        end

        D_V(ML, ND) = f0 * Theta_g(ML, ND) * D_A(MN);

        if cfg.lDEBUG_NAN_TRACE
            if isnan(D_A(MN)) || isnan(D_V(ML,ND))
                fprintf('  [CRITICAL@ML=%d,ND=%d] D_A=%.2e, D_V=%.2e, Theta_g=%.4f\n', ...
                    ML, ND, D_A(MN), D_V(ML,ND), Theta_g(ML,ND));
            end
        end
    end
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[calculateVaporVariables] Output: D_V NaN=%d, D_A NaN=%d, Eta NaN=%d\n', ...
        sum(isnan(D_V(:))), sum(isnan(D_A(:))), sum(isnan(Eta(:))));
end

TransportCoefficient.D_V = D_V;
TransportCoefficient.Eta = Eta;
TransportCoefficient.D_A = D_A;
end
