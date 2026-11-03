function [es_Pa, des_dT, qs, dqs_dT] = Calc_Esat_Flatau(T_degC, P_atm_Pa)
%CALC_ESAT_FLATAU  Saturation vapour pressure and humidity derivatives.
%
%   [es_Pa, des_dT, qs, dqs_dT] = equations.Calc_Esat_Flatau(T_degC,
%   P_atm_Pa) evaluates the Flatau polynomial used by CLM-style atmospheric
%   forcing calculations. Water coefficients are used for T >= 0 degC and ice
%   coefficients for T < 0 degC.
%
%   Inputs:
%     T_degC    [degC] scalar or vector air temperature
%     P_atm_Pa  [Pa] scalar or vector atmospheric pressure
%
%   Returns:
%     es_Pa     [Pa] saturation vapour pressure
%     des_dT    [Pa/K] derivative of saturation vapour pressure
%     qs        [kg/kg] saturation specific humidity
%     dqs_dT    [kg/kg/K] derivative of saturation specific humidity
%
%   Reference:
%     Flatau et al. 1992 J. Appl. Meteorol. 31:1507-1513; coefficients used by
%     CLM-style saturation vapour pressure routines.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    T_degC = T_degC(:);
    P_atm_Pa = P_atm_Pa(:);

    nT = length(T_degC);
    nP = length(P_atm_Pa);
    n_pts = max(nT, nP);

    if nT == nP
        % Matching scalar/vector lengths need no expansion.
    elseif nT == 1
        T_degC = repmat(T_degC, n_pts, 1);
    elseif nP == 1
        P_atm_Pa = repmat(P_atm_Pa, n_pts, 1);
    else
        error('Calc_Esat_Flatau:dimensionMismatch', ...
            'T_degC and P_atm_Pa must have matching lengths, or one input must be scalar.');
    end

    % Flatau polynomial coefficients, in hPa, for water and ice surfaces.
    a_wat = [6.11213476; 4.44007856e-1; 1.43064234e-2; 2.64461437e-4; ...
             3.05903558e-6; 1.96237241e-8; 8.92344772e-11; -3.73208410e-13; 2.09339997e-16];
    b_wat = [4.44017302e-1; 2.86064092e-2; 7.94683137e-4; 1.21211669e-5; ...
             1.03354611e-7; 4.04125005e-10; -7.88037859e-13; -1.14596802e-14; 3.81294516e-17];

    a_ice = [6.11123516; 5.03109514e-1; 1.88369801e-2; 4.20547422e-4; ...
             6.14396778e-6; 6.02780717e-8; 3.87940929e-10; 1.49436277e-12; 2.62655803e-15];
    b_ice = [5.03277922e-1; 3.77289173e-2; 1.26801703e-3; 2.49468427e-5; ...
             3.13703411e-7; 2.57180651e-9; 1.33268878e-11; 3.94116744e-14; 4.98070196e-17];

    mask_ice = T_degC < 0;

    a_sel = repmat(a_wat.', n_pts, 1);
    b_sel = repmat(b_wat.', n_pts, 1);
    if any(mask_ice)
        n_ice = nnz(mask_ice);
        a_sel(mask_ice, :) = repmat(a_ice.', n_ice, 1);
        b_sel(mask_ice, :) = repmat(b_ice.', n_ice, 1);
    end

    poly_es  = zeros(n_pts, 1);
    poly_des = zeros(n_pts, 1);
    for k = 1:9
        term_T   = T_degC.^(k - 1);
        poly_es  = poly_es  + a_sel(:, k) .* term_T;
        poly_des = poly_des + b_sel(:, k) .* term_T;
    end

    es_Pa  = 100 .* poly_es;
    des_dT = 100 .* poly_des;

    if nargout >= 3
        term_denom = P_atm_Pa - 0.378 .* es_Pa;
        qs = (0.622 .* es_Pa) ./ term_denom;
        dqs_dT = (0.622 .* P_atm_Pa) ./ (term_denom.^2) .* des_dT;
    end
end