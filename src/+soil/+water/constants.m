function constants = constants()
%CONSTANTS  Return fixed constants for the STEMMUS soil-water formulation.
%
%   constants = soil.water.constants() returns empirical dry-end thresholds
%   used by the soil-water auxiliary equations. These are not universal
%   physical constants: they are process-model constants for the STEMMUS
%   dry-transition hydraulic behavior.
%
%   Outputs:
%     constants.hd          [cm]     dry-end matric-head cutoff
%     constants.hm          [cm]     dry-transition matric-head threshold
%
%   Notes:
%     hd = -1e7 cm is intentionally retained from the stabilized legacy path.
%     A previous local experiment with hd = -3e4 cm caused a near-singular
%     first-step matrix through the dry-transition gamma factor, so this value
%     should not be changed without a dedicated water/thermal stability test.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    constants.hd = -1e7;
    constants.hm = -9899;
end
