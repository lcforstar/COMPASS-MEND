function theta_i = iceContent(T_C, h, fieldCapacity, switches)
%ICECONTENT  Compute soil ice volume fraction for the auxiliary refresh.
%
%   theta_i = soil.water.iceContent(T_C, h, fieldCapacity, switches) computes
%   the simplified sinusoidal freezing partition used when freezing is enabled
%   and strict SFCC is disabled.
%
%   Inputs:
%     T_C            [in] [degC] node temperature
%     h              [in] [cm]   matric head
%     fieldCapacity  [in] [-]    upper bound for ice content
%     switches       [in] struct built.control.switches
%
%   Output:
%     theta_i        [out] [-]   ice volume fraction
%
%   Reference:
%     Lundin 1990 Cold Reg. Sci. Technol. 18:47-59, sinusoidal transition
%     structure inherited from STEMMUS UpdateSoilWaterContent.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if ~isFreezingOn(switches) || h <= -1e7
        theta_i = 0;
        return;
    end
    if isfield(switches, 'lSFCC') && switches.lSFCC == 1
        theta_i = 0;
        return;
    end

    T_K = T_C + 273.15;
    Tf1 = 274.15;  % [K] fully liquid above this temperature
    Tf2 = 270.15;  % [K] fully frozen below this temperature
    if T_K > Tf1
        theta_i = 0;
    elseif T_K >= Tf2
        frac = 0.5 * (1 - sin(pi() * (T_K - 0.5 * Tf1 - 0.5 * Tf2) / (Tf1 - Tf2)));
        theta_i = max(0, min(1, frac)) * fieldCapacity;
    else
        theta_i = fieldCapacity;
    end
end

function tf = isFreezingOn(switches)
% Return true when freezing physics is enabled.
    tf = isfield(switches, 'lFreezing') && switches.lFreezing == 1;
end
