function theta_l = liquidContent(theta_u, theta_i, soilPar, state, switches, phys)
%LIQUIDCONTENT  Compute liquid water content under configured freeze physics.
%
%   theta_l = soil.water.liquidContent(theta_u, theta_i, soilPar, state,
%   switches, phys) partitions total water into liquid water for one
%   layer-endpoint auxiliary state.
%
%   Inputs:
%     theta_u   [in] [-] total water content from the retention curve
%     theta_i   [in] [-] ice volume fraction before mass-balance correction
%     soilPar   [in] struct scalar soil parameters: theta_s/theta_r/alpha/n/m
%     state     [in] struct scalar state: h, h_frez, T_C, gamma_hh
%     switches  [in] struct built.control.switches
%     phys      [in] struct physical constants
%
%   Output:
%     theta_l   [out] [-] liquid water content
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    strictSFCC = isfield(switches, 'lSFCC') && switches.lSFCC == 1 ...
        && isfield(switches, 'lThmrlefc') && switches.lThmrlefc == 1;
    if strictSFCC
        if state.h >= -1 && state.h_frez >= 0
            theta_l = soilPar.theta_s;
        else
            theta_m = state.gamma_hh * soilPar.theta_s;
            theta_l = state.gamma_hh * soilPar.theta_r + ...
                (theta_m - state.gamma_hh * soilPar.theta_r) / ...
                (1 + abs(soilPar.alpha * (state.h + state.h_frez))^soilPar.n)^soilPar.m;
        end
    else
        T_K = state.T_C + phys.T_freeze;
        if T_K > phys.T_freeze + 1 || state.h >= -1e-6
            theta_l = soilPar.theta_s;
        else
            theta_l = theta_u - theta_i * phys.RHOI / phys.RHOL;
        end
    end
end
