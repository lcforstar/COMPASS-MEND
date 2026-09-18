function gamma = computeGamma(h, constants)
%COMPUTEGAMMA  Compute liquid-water availability factor for vapor enhancement.
%
%   gamma = soil.water.computeGamma(h, constants) evaluates the STEMMUS
%   dry-transition factor from matric head. It is an auxiliary field used by
%   coupled water-vapor terms and can be recomputed from the current head
%   state.
%
%   Inputs:
%     h          [in]  numeric [cm] node matric head
%     constants  [in]  struct with hd/hm [cm] dry-end thresholds
%
%   Output:
%     gamma  [out] numeric [-] liquid-water availability factor in [0,1]
%
%   Method:
%     Uses the legacy STEMMUS thresholds:
%       |h| >= |hd| -> 0, |h| < |hm| -> 1, logarithmic interpolation in between.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if nargin < 2 || isempty(constants)
        constants = soil.water.constants();
    end

    gamma = ones(size(h));
    dry = abs(h) >= abs(constants.hd);
    transition = ~dry & abs(h) >= abs(constants.hm);
    gamma(dry) = 0;
    gamma(transition) = log(abs(constants.hd) ./ abs(h(transition))) ./ ...
        log(abs(constants.hd) / abs(constants.hm));
end
