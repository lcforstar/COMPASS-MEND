function [c_nh4, c_no3] = calcNitrogenConcentration(depths_input, c_nh4_input, c_no3_input, ...
                                                    depths_target, bulk_density, ModelSettings) %#ok<INUSD>
    %CALCNITROGENCONCENTRATION  Interpolate nitrogen concentrations and convert units.
    %
    % Syntax:
    %   [c_nh4, c_no3] = calcNitrogenConcentration(depths_input, c_nh4_input, ...
    %       c_no3_input, depths_target, bulk_density, ModelSettings)
    %
    % Inputs:
    %   depths_input   - measurement depths [cm]
    %   c_nh4_input    - measured ammonium concentration [mg/kg, mass basis]
    %   c_no3_input    - measured nitrate concentration  [mg/kg, mass basis]
    %   depths_target  - target depths [cm]
    %   bulk_density   - bulk density at each depth [g/cm3]
    %   ModelSettings  - model settings (currently unused; kept for signature)
    %
    % Returns:
    %   c_nh4          - ammonium concentration [ug/cm3, volume basis]
    %   c_no3          - nitrate concentration  [ug/cm3, volume basis]
    %
    % Method:
    %   1. linear-interpolate the mass concentration [mg/kg]
    %   2. unit conversion: mg/kg * rho_b [g/cm3] = ug/cm3
    %      (1 mg/kg * 1 g/cm3 = 1 ug/cm3)
    %   3. non-negative protection: negatives replaced by neighbouring positives or 0
    %
    % Units:
    %   input : mg/kg (mass) = mg N per kg dry soil = ppm (mass)
    %   output: ug/cm3 (volume) = ug N per cm3 dry soil
    %   conversion: c[ug/cm3] = c[mg/kg] * rho_b[g/cm3]
    %
    % Editor: Jizixing Zhu; Date: 2026-06-27.
    % Change: 2026-06-27 Jizixing Zhu -- comments/messages translated to English.

    % force column vectors
    depths_input = depths_input(:);
    c_nh4_input = c_nh4_input(:);
    c_no3_input = c_no3_input(:);
    depths_target = depths_target(:);
    bulk_density = bulk_density(:);

    fprintf('  nitrogen concentration init...\n');

    % =========================================================================
    % Step 1: linear-interpolate the mass concentration [mg/kg]
    % =========================================================================
    c_nh4_mass = init.interpLinearWithExtrapolation(depths_input, c_nh4_input, depths_target);
    c_no3_mass = init.interpLinearWithExtrapolation(depths_input, c_no3_input, depths_target);

    % =========================================================================
    % Step 2: unit conversion (mg/kg -> ug/cm3): c[ug/cm3] = c[mg/kg] * rho_b[g/cm3]
    % =========================================================================
    c_nh4 = c_nh4_mass .* bulk_density;
    c_no3 = c_no3_mass .* bulk_density;

    % =========================================================================
    % Step 3: non-negative protection
    % =========================================================================
    c_nh4 = apply_non_negative_protection(c_nh4);
    c_no3 = apply_non_negative_protection(c_no3);

    % =========================================================================
    % Step 4: summary
    % =========================================================================
    fprintf('  nitrogen concentration init done\n');
    fprintf('      NH4+ range: [%.3f, %.3f] ug/cm3\n', min(c_nh4(c_nh4 >= 0)), max(c_nh4(c_nh4 >= 0)));
    fprintf('      NO3-  range: [%.3f, %.3f] ug/cm3\n', min(c_no3(c_no3 >= 0)), max(c_no3(c_no3 >= 0)));

end

% =========================================================================
% Helper: non-negative protection
% =========================================================================
function c_clean = apply_non_negative_protection(c_input)
    %APPLY_NON_NEGATIVE_PROTECTION  Replace negative values.
    %
    % Strategy: search forward and backward for the nearest non-negative value.
    %   - if found on both sides, use their average
    %   - if found on one side, use that value
    %   - if none found, use 0
    %
    % Editor: Jizixing Zhu; Date: 2026-06-27.

    c_clean = c_input;

    for i = 1:length(c_clean)
        if c_clean(i) < 0
            % forward search
            forward_val = [];
            for j = i + 1:length(c_clean)
                if c_clean(j) >= 0
                    forward_val = c_clean(j);
                    break;
                end
            end

            % backward search
            backward_val = [];
            for j = i - 1:-1:1
                if c_clean(j) >= 0
                    backward_val = c_clean(j);
                    break;
                end
            end

            % pick the available value(s)
            if ~isempty(forward_val) && ~isempty(backward_val)
                c_clean(i) = (forward_val + backward_val) / 2;
            elseif ~isempty(forward_val)
                c_clean(i) = forward_val;
            elseif ~isempty(backward_val)
                c_clean(i) = backward_val;
            else
                c_clean(i) = 0;  % no positive value found -> 0
            end
        end
    end
end
