function SoilVariables = Density_V(SoilVariables, cfg)
%DENSITY_V  Compute phase-consistent pore-vapor density and state derivatives.
%
%   SoilVariables [inout] provides trial temperature TT [degC] and matric
%   head hh [cm]. Trial freezing head hh_frez [cm] (or h_frez as a legacy
%   fallback) selects the ice-equilibrium branch. The function returns RHOV
%   [g cm^-3], DRHOVh [g cm^-3 cm^-1], and DRHOVT [g cm^-3 K^-1].
%
%   Unfrozen nodes use the clipped Kelvin relation over liquid water. Frozen
%   nodes use the IAPWS ice-Ih sublimation pressure directly; consequently,
%   equilibrium vapor density is locally independent of matric head on that
%   branch. Temperatures outside [-100, 150] degC are rejected rather than
%   silently projected onto a constant-density plateau.
%
%   References:
%     Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.23-4.24, pp. 70-71.
%     IAPWS R14-08(2011), Eq. 6 and Table 3, pp. 5 and 7; local copy:
%     data/references/soil_gas_transport/IAPWS_2011_MeltSub.pdf.
%     Piecewise derivation: docs/compass_soil_water_heat_physics.md,
%     Eqs. 3.8-3.10 and 4.1.5p.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- differentiate the actual clipped
%           Kelvin and saturation-density functions.
%   Change: 2026-07-13 Jizixing Zhu -- replace the -20 degC vapor-density
%           plateau with an IAPWS ice-equilibrium branch and analytic slope.
%   Change: 2026-07-13 Jizixing Zhu -- vectorize ice saturation density and
%           its temperature derivative for simultaneous frozen nodes.

NN  = cfg.Grid.NN;
TT  = SoilVariables.TT(:).';
hh  = SoilVariables.hh(:).';
g   = equations.getPhysicalConstant('g');
Rv  = equations.getPhysicalConstant('Rv');
if isfield(SoilVariables, 'hh_frez') && ~isempty(SoilVariables.hh_frez)
    hFrez = SoilVariables.hh_frez(:).';
elseif isfield(SoilVariables, 'h_frez') && ~isempty(SoilVariables.h_frez)
    hFrez = SoilVariables.h_frez(:).';
else
    hFrez = zeros(size(hh));
end

assert(numel(TT) == NN && numel(hh) == NN && numel(hFrez) == NN, ...
    'Density_V:invalidNodeCount', ...
    'TT, hh, and freezing head must each contain cfg.Grid.NN values.');
assert(all(isfinite(TT(:))) && all(TT(:) >= -100) && all(TT(:) <= 150), ...
    'Density_V:temperatureOutsideEarthRange', ...
    'Soil temperature must remain within [-100, 150] degC.');

T_K = TT + 273.15;
frozen = hFrez < -1e-12;
unfrozen = ~frozen;

% Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.23-4.24, pp. 70-71.
% Evaluate the liquid-water Kelvin branch in one vector operation.
HRraw = exp(hh .* g ./ (Rv .* T_K));
HR = min(1, max(0.041, HRraw));
chiKelvin = double(HRraw > 0.041 & HRraw < 1);
rhoWater = 1e-6 .* exp(31.3716 - 6014.79 ./ T_K ...
    - 7.92495e-3 .* T_K) ./ T_K;
dRhoWaterDT = rhoWater .* (6014.79 ./ T_K.^2 ...
    - 7.92495e-3 - 1 ./ T_K);

RHOV = rhoWater .* HR;
RHOVSat = rhoWater;
DRHOVSatT = dRhoWaterDT;
DRHOVh = rhoWater .* chiKelvin .* HR .* g ./ (Rv .* T_K);
DRHOVT = HR .* dRhoWaterDT ...
    - rhoWater .* chiKelvin .* HR .* hh .* g ./ (Rv .* T_K.^2);

if any(frozen)
    % IAPWS R14-08(2011), Eq. 6, p. 5.
    % Ice, liquid, and vapor share one chemical potential while ice is
    % active, so rho_v=rho_vi(T) and its local h derivative is zero.
    [rhoIce, dRhoIceDT] = iceSaturationDensity(T_K(frozen), Rv);
    RHOV(frozen) = rhoIce;
    RHOVSat(frozen) = rhoIce;
    DRHOVSatT(frozen) = dRhoIceDT;
    DRHOVh(frozen) = 0;
    DRHOVT(frozen) = dRhoIceDT;
end

SoilVariables.RHOV   = RHOV;
SoilVariables.RHOV_s = RHOVSat;
SoilVariables.DRHOV_sT = DRHOVSatT;
SoilVariables.DRHOVh = DRHOVh;
SoilVariables.DRHOVT = DRHOVT;
end

function [rhoIce, dRhoIceDT] = iceSaturationDensity(T_K, RvCgs)
% Return IAPWS ice-Ih saturation vapor density and its analytic derivative.
    tripleT = 273.16;                              % [K]
    tripleP = 611.657;                             % [Pa]
    theta = T_K / tripleT;                         % [-]
    a1 = -21.2144006; b1 = 0.00333333333;          % [-]
    a2 = 27.3203819;  b2 = 1.20666667;             % [-]
    a3 = -6.10598130; b3 = 1.70333333;             % [-]

    logPressureRatio = a1 .* theta.^(b1 - 1) ...
        + a2 .* theta.^(b2 - 1) + a3 .* theta.^(b3 - 1);
    dLogPressureDT = (a1 .* (b1 - 1) .* theta.^(b1 - 2) ...
        + a2 .* (b2 - 1) .* theta.^(b2 - 2) ...
        + a3 .* (b3 - 1) .* theta.^(b3 - 2)) ./ tripleT;
    pressureIce = tripleP * exp(logPressureRatio); % [Pa]
    RvSi = RvCgs / 1e4;                            % [J kg^-1 K^-1]
    rhoIce = 1e-3 .* pressureIce ./ (RvSi .* T_K); % [g cm^-3]
    dRhoIceDT = rhoIce .* (dLogPressureDT - 1 ./ T_K); % [g cm^-3 K^-1]
end
