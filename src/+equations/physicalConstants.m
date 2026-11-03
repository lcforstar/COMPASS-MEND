function phys = physicalConstants()
%PHYSICALCONSTANTS  Fixed physical constants shared by COMPASS model modules.
%
%   phys = physicalConstants() returns scalar physical constants in a plain
%   value struct. The function is side-effect free: MATLAB callers get a small
%   immutable-by-convention struct, and the Julia port can map the same values
%   to module-level const bindings or an immutable PhysicalConstants struct.
%
%   Returns:
%     phys  struct  global physical constants, gas constants, thermal constants,
%                   radiation constants, model constants and unit conversions.
%
%   Notes:
%     Field names mirror legacy PhysicalConstants where possible, so modules can
%     migrate from the legacy handle singleton without renaming physical symbols.
%     Barometric fallback constants are included explicitly because forcing
%     pressure estimation previously hard-coded them in buildForcing.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-13 Jizixing Zhu -- centralize the matric-head
%           temperature correction slope used by retention and liquid flow.
%   Change: 2026-07-20 Jizixing Zhu -- cache the immutable value struct so
%           hot process kernels no longer depend on the legacy core singleton.

    persistent cached
    if ~isempty(cached)
        phys = cached;
        return
    end

    phys = struct();

    % --- basic physics ---
    phys.g       = 980.616;       % [cm s-2] gravitational acceleration
    phys.g_SI    = 9.80616;       % [m s-2] gravitational acceleration
    phys.Re      = 6.37122e6;     % [m] Earth radius
    phys.k_VK    = 0.4;           % [-] von Karman constant
    phys.N_Avog  = 6.02214e26;    % [kmol-1] Avogadro constant
    phys.k_Boltz = 1.38065e-23;   % [J K-1] Boltzmann constant

    % --- molecular weights and gas constants ---
    phys.MW_da = 28.966;          % [kg kmol-1] dry-air molecular weight
    phys.MW_wv = 18.016;          % [kg kmol-1] water-vapour molecular weight

    phys.R_gas_kmol = phys.N_Avog * phys.k_Boltz;    % [J kmol-1 K-1]
    phys.R          = phys.R_gas_kmol / 1000;        % [J mol-1 K-1]
    phys.R_da_SI    = phys.R_gas_kmol / phys.MW_da;  % [J kg-1 K-1]
    phys.RDA        = phys.R_da_SI * 1e4;            % [cm2 s-2 degC-1]
    phys.R_wv_SI    = phys.R_gas_kmol / phys.MW_wv;  % [J kg-1 K-1]
    phys.Rv         = phys.R_wv_SI * 1e4;            % [cm2 s-2 degC-1]
    phys.epsilon    = phys.MW_wv / phys.MW_da;       % [-] MW_wv / MW_da

    % Compatibility aliases for explicit-unit users.
    phys.MW_da_kgkmol = phys.MW_da;       % [kg kmol-1]
    phys.MW_wv_kgkmol = phys.MW_wv;       % [kg kmol-1]

    % --- heat capacity and latent heat ---
    phys.Tr       = 20;          % [degC] reference temperature
    phys.T_freeze = 273.15;      % [K] freezing-point temperature

    phys.Cp_liq_SI = 4188.0;     % [J kg-1 K-1] liquid-water heat capacity
    phys.c_V       = 1.870;      % [J g-1 degC-1] water-vapour heat capacity
    phys.Cp_air_SI = 1004.64;    % [J kg-1 K-1] dry-air heat capacity
    phys.Cp_ice_SI = 2117.27;    % [J kg-1 K-1] ice heat capacity
    phys.c_L       = phys.Cp_liq_SI * 0.001;  % [J g-1 degC-1]
    phys.c_a       = phys.Cp_air_SI * 0.001;  % [J g-1 degC-1]

    phys.L0         = 2501.0;    % [J g-1] vaporization latent heat
    phys.lambda_vap = 2.501e6;   % [J kg-1] vaporization latent heat
    phys.lambda_fus = 3.337e5;   % [J kg-1] fusion latent heat
    phys.lambda_sub = 2.8347e6;  % [J kg-1] sublimation latent heat

    % --- densities, conductivity and material properties ---
    phys.RHOL       = 1;         % [g cm-3] liquid-water density
    phys.RHOI       = 0.92;      % [g cm-3] ice density
    phys.rho_liq_SI = 1000.0;    % [kg m-3] liquid-water density
    phys.rho_ice_SI = 917.0;     % [kg m-3] ice density

    phys.Lambda1 = 0.228 / 100;     % [W m-1 degC-1] soil thermal coefficient
    phys.Lambda2 = -2.406 / 100;    % [W m-1 degC-1] soil thermal coefficient
    phys.Lambda3 = 4.909 / 100;     % [W m-1 degC-1] soil thermal coefficient
    phys.lambda_liq_cond = 0.57;    % [W m-1 K-1] liquid-water thermal conductivity
    phys.lambda_ice_cond = 2.29;    % [W m-1 K-1] ice thermal conductivity
    phys.lambda_air_cond = 0.023;   % [W m-1 K-1] air thermal conductivity

    phys.MU_W0  = 2.4152e-4;        % [g cm-1 s-1] liquid-water viscosity coefficient
    phys.MU1    = 4742.8;           % [J mol-1] viscosity temperature coefficient
    phys.b      = 4e-6;             % [cm] liquid-water viscosity parameter
    phys.MU_a   = 1.846e-4;         % [g cm-1 s-1] air viscosity coefficient
    phys.CKTN   = 50 + 2.575 * 20;  % [-] viscosity temperature correction term
    phys.Gamma0 = 71.89;            % [g s-2] water surface tension at 25 degC
    phys.Cw     = 0.0068;            % [K-1] matric-head temperature correction slope
    phys.Hc     = 0.02;             % [-] Henry-law constant used by legacy gas code

    % --- radiation and meteorology ---
    phys.Gsc     = 1360;          % [W m-2] solar constant used by legacy code
    phys.S_solar = 1361.0;        % [W m-2] WMO-style solar constant
    phys.Sigma_E = 4.90e-9;       % [MJ m-2 degC-4] Stefan-Boltzmann legacy unit
    phys.sigma   = 5.67e-8;       % [W m-2 K-4] Stefan-Boltzmann constant
    phys.Pstd    = 101325.0;      % [Pa] standard atmospheric pressure
    phys.Gamma_d = 0.0098;        % [K m-1] dry adiabatic lapse rate
    phys.O2_mol_fraction_dry_air = 0.2095;   % [-] dry-air O2 mole fraction
    phys.Oa_mmolmol_default      = phys.O2_mol_fraction_dry_air * 1000; % [mmol/mol]

    % --- model-specific constants and unit conversions ---
    phys.rBot      = 0;            % [LT-1] bottom flux resistance
    phys.W0        = 1.001e3;      % [-] Milly wetting-heat coefficient
    phys.GWT       = 7;            % [-] temperature dependence factor
    phys.Gamma_w   = 981;          % [g cm-2 s-2] unit weight of liquid water
    phys.PeCr      = 2;            % [-] critical Peclet number
    phys.ExcesInt  = 0;            % [cm] previous-step excess interception
    phys.aInterc   = 0.025;        % [cm d-1] interception model coefficient

    phys.Unit_C    = 1;            % [-] water-head unit conversion
    phys.UnitC     = 100;          % [cm m-1] length conversion
    phys.hPa_to_Pa = 100.0;        % [Pa hPa-1]
    phys.mm_to_m   = 0.001;        % [m mm-1]
    phys.cm_to_m   = 0.01;         % [m cm-1]
    phys.g_to_kg   = 0.001;        % [kg g-1]
    phys.C_to_K    = 273.15;       % [K] Celsius-to-Kelvin offset
    phys.J_to_MJ   = 1.0e-6;       % [MJ J-1]

    % --- pressure fallback constants ---
    phys.g_barometric   = 9.80665;    % [m s-2] gravity in legacy barometric formula
    phys.MW_da_kg_mol   = 0.0289644;  % [kg mol-1] dry-air molecular weight
    phys.R_universal_SI = 8.31432;    % [J mol-1 K-1] universal gas constant

    cached = phys;
end
