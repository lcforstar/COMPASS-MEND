function TransportCoefficient = calculateThermalLiquidConductivity(SoilVariables, cfg, TransportCoefficient)
%CALCULATETHERMALLIQUIDCONDUCTIVITY Compute temperature-driven liquid conductivity.
%
%   Inputs:
%     SoilVariables         [in] struct containing hh [cm], KL_h [cm/s],
%                           and optional freeze-impeded KfL_h [cm/s]
%     cfg                   [in] model settings containing Grid.NL and
%                           lKLT_Switch
%     TransportCoefficient [inout] transport auxiliary structure
%
%   Output:
%     TransportCoefficient.KL_T [cm2 s-1 K-1, NL x 2] thermal liquid
%                               conductivity at layer endpoints
%
%   Method:
%     COMPASS uses one empirical logarithmic temperature slope Cw for both
%     the retention correction COR=exp[-Cw*(T-Tr)] and thermal liquid flow.
%     Thus KL_T=KL_eff*h*d(ln COR)/dT=-KL_eff*Cw*h. This replaces the legacy
%     sand-specific GWT/Gamma0 slope without changing the head-gradient term.
%
%   Reference:
%     Zeng and Su, STEMMUS Governing Equations and Constitutive Equations,
%     Section 3.2.1, p.19: hTemCorr=h*exp[-Cpsi*(T-Tr)], Cpsi=0.0068 K^-1.
%
%   Editor: Jizixing Zhu; Date: 2026-07-13.

NL         = cfg.Grid.NL;
hh         = SoilVariables.hh;
KL_h       = SoilVariables.KL_h;
Cw         = equations.getPhysicalConstant('Cw');

KL_T = zeros(NL, 2);
for ML = 1:NL
    for ND = 1:2
        MN = ML + ND - 1;
        if cfg.lKLT_Switch == 1 && hh(MN) < 0
            KL_eff = KL_h(ML,ND);
            if (isstruct(SoilVariables) && isfield(SoilVariables, 'KfL_h')) || ...
                    (isobject(SoilVariables) && isprop(SoilVariables, 'KfL_h'))
                KL_eff = SoilVariables.KfL_h(ML,ND);
            end
            KL_T(ML,ND) = -KL_eff * Cw * hh(MN);
        end
    end
end

TransportCoefficient.KL_T = KL_T;
end
