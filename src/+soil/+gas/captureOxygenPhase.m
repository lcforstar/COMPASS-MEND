function phase = captureOxygenPhase(built, useTrial)
%CAPTUREOXYGENPHASE Freeze the original species thermodynamic phase contract.
% Bottom-up endpoint theta [-], node T [C], pressure [Pa], porosity [-].
% This value-only record owns no physical state and changes no carrier field.
% Source: COMPASS-COMPASS_MATLAB advanceOxygenSubstep phase adapters.
% Editor: Wenjie Liu; Date: 2026-10-23.
    phys=equations.physicalConstants();
    if nargin<2, useTrial=false; end
    if useTrial
        Y=built.workspace.step.trial;
        thetaL=built.workspace.derived.Theta_LL;
        thetaG=built.workspace.derived.Theta_g;
    else
        Y=built.Y;
        thetaL=Y.Theta_L;
        thetaG=repmat(built.parameters.soil.porosity(:),1,2)-thetaL-Y.Theta_I;
    end
    T=Y.T(:);
    phase=struct('theta_l',thetaL,'theta_g',thetaG, ...
        'temperature_C',T,'porosity',built.parameters.soil.porosity(:), ...
        'gas_pressure_Pa',Y.p_gas(:)/10, ...
        'vapor_pressure_Pa',built.workspace.derived.RHOV(:).*phys.Rv.*(T+phys.C_to_K)/10);
end
