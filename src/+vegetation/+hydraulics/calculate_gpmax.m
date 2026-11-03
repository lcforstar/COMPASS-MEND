function gpmax = calculate_gpmax(H, Ksmax, CA, SCA, rho_w)
% CALCULATE_GPMAX Calculate specific maximum hydraulic conductance
%
% Syntax:
%   gpmax = calculate_gpmax(DBH, H, Ksmax, b, c, rho_w)
%
% Inputs:
%   H      - Tree height (m)
%   Ksmax  - Stem area-specific maximum hydraulic conductance kg m − 1 s − 1 MPa− 1
%   
%   rho_w  - Water density (kg/m³), default = 1000 kg/m³
%
% Output:
%   gpmax  - Specific maximum hydraulic conductance
%
% References:
%   Based on equations from Manzoni et al. (2013)
%   gpmax = (Ksmax * SCA) / (H * CA * rho_w)

    % Set default water density if not provided
    if nargin < 5
        rho_w = 1000; % kg/m³
    end
    
    % Calculate specific maximum hydraulic conductance using equation (8)
    gpmax = (Ksmax * SCA) / (H * CA * rho_w);
    
  
end

