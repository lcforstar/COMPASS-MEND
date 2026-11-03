function [RWU, hydraulic] = assembleRWU(hydraulic, soil_state, plantHydro, domain, DeltZ)
%ASSEMBLERWU CLM5-style root water uptake sink assembly (with HR + gravity).
%
% Physics (CLM5 SoilWaterPlantSinkMod.F90:317-325 HydStress):
%   flux_j = gsr_j * (psi_soil_j - psi_root - soil_gravity_j)         [m/s]
%   gravity_j = rho_w g * depth_j
%   RWU_j  = flux_j / DeltZ_node_j * 100                              [1/s]
%   positive = soil->root, negative = root->soil (HR)
%   mass closure: sum(flux) = sumgsr * (psis - psir) = SPAC actual_total
%
% Inputs:
%   hydraulic    SPAC output: psir [MPa], gsr [NN x 1 m/s/MPa], Tp_t [cm/s]
%   soil_state   target soil slice with .h [cmH2O, NN x 1]
%   plantHydro   plant hydraulic bundle (rootContact diagnostic optional)
%   domain       .nl, .nn, .dz_cm [cm]
%   DeltZ        [cm, NL x 1] layer thickness (bottom-up: DeltZ(1)=deepest)
%
% Outputs:
%   RWU          [1/s, NN x 1] Richards sink units
%   hydraulic    adds rwu_diag and RWU_pos/neg node diagnostics
%
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-16 Jizixing Zhu — read soil_state.h and domain.nl/nn/dz_cm (drop hh/Grid bridge).

    if ~isstruct(domain) || ~isfield(domain, 'nl') || ~isfield(domain, 'nn') || ~isfield(domain, 'dz_cm')
        error('assembleRWU:badDomain', 'domain must provide nl/nn/dz_cm.');
    end
    NL_rwu = domain.nl;
    NN_rwu = domain.nn;
    dz_cm = domain.dz_cm;
    if ~isstruct(soil_state) || ~isfield(soil_state, 'h')
        error('assembleRWU:badSoil', 'soil_state.h [cmH2O] is required.');
    end
    hhMPA   = soil_state.h(:) * 9.80665e-5;       % cmH2O -> MPa
    psir_sc = hydraulic.psir;                          % MPa
    gsr_rwu = hydraulic.gsr;                           % NNx1 [m/s/MPa]

    % Node depths [m] below surface; bottom-up grid: node 1 deepest, NN surface
    DeltZ_cum = flip(cumsum(flip(dz_cm(:))));  % NL x 1 [cm], layer top depth
    depth_node_m = zeros(NN_rwu, 1);
    depth_node_m(NN_rwu) = 0;                          % surface node depth = 0
    for MN = 1:NN_rwu-1
        depth_node_m(MN) = DeltZ_cum(MN) / 100;        % cm -> m
    end
    soil_gravity_MPa = 1000 * 9.80665 * depth_node_m / 1e6;  % NNx1 [MPa]

    % ---- 驱动力（CLM5 风格，含重力）----
    grad_node = hhMPA - psir_sc - soil_gravity_MPa;    % NNx1 [MPa]

    % ---- 节点支配长度 DeltZ_node [cm]: 边界节点占一半相邻层，内部节点各占一半 ----
    DeltZ_node = zeros(NN_rwu, 1);
    DeltZ_node(1) = DeltZ(1) / 2;
    for MN = 2:NL_rwu
        DeltZ_node(MN) = (DeltZ(MN-1) + DeltZ(MN)) / 2;
    end
    DeltZ_node(NN_rwu) = DeltZ(NL_rwu) / 2;

    target_uptake = max(0, hydraulic.Tp_t) / 100;      % cm/s -> m/s
    zero_transpiration = target_uptake <= 1e-12;
    disable_root_return = strcmpi(strtrim(getenv('STEMMUS_DISABLE_ROOT_RETURN')), '1');

    % ---- 节点交换通量 + 诊断聚合 ----
    flux_node_full = gsr_rwu .* grad_node;             % NNx1 [m/s], 含正负
    if disable_root_return && ~zero_transpiration
        % Diagnostic-only A/B switch. Remove root->soil return while keeping
        % total net uptake equal to the SPAC transpiration demand.
        flux_pos_only = max(flux_node_full, 0);
        pos_total = sum(flux_pos_only);
        if pos_total > 0
            flux_node_full = flux_pos_only * (target_uptake / pos_total);
        else
            flux_node_full = zeros(NN_rwu, 1);
        end
    end
    F_pos = sum(max(flux_node_full, 0));               % [m/s] 土→根总流入
    F_neg = sum(min(flux_node_full, 0));               % [m/s] 根→土回流 ≤0
    F_net = F_pos + F_neg;                             % [m/s] 净 = SPAC actual_total

    % ---- 诊断 ----
    rc_enabled = false;
    if isfield(plantHydro, 'rootContact') && isstruct(plantHydro.rootContact) ...
            && isfield(plantHydro.rootContact, 'enabled')
        rc_enabled = logical(plantHydro.rootContact.enabled);
    end
    hydraulic.rwu_diag = struct( ...
        'root_contact_enabled', rc_enabled, ...
        'min_hh_mpa', min(hhMPA), ...
        'target_uptake_m_s', target_uptake, ...
        'max_soil_gravity_mpa', max(soil_gravity_MPa), ...
        'F_pos_m_s', F_pos, ...
        'F_neg_m_s', F_neg, ...
        'F_net_m_s', F_net, ...
        'qflx_phs_neg_m_s', F_neg, ...
        'capacity_integral_m_s', F_pos, ...
        'root_return_disabled', disable_root_return, ...
        'zero_transpiration', zero_transpiration);

    % ---- 组装 RWU 单位转换 ----
    if zero_transpiration
        RWU_capacity = zeros(NN_rwu, 1);
    else
        % CLM5 单步隐式 sink: 正负原样进 Richards rmx (Srt → C9 装配支持负值)
        % Σ flux = SPAC actual_total 数学恒等 (见 Plantrootconductivity.m 同样的 gravity)
        RWU_capacity = flux_node_full ./ DeltZ_node;   % NNx1 [m/s/cm], 含正负
    end

    hydraulic.rwu_diag.active_integral_m_s = sum(RWU_capacity .* DeltZ_node);
    hydraulic.RWU_pos_node = max(RWU_capacity, 0);     % 干层吸水强度 [m/s/cm]
    hydraulic.RWU_neg_node = min(RWU_capacity, 0);     % 湿层 HR 回灌强度 [m/s/cm]

    RWU = RWU_capacity * 100;   % [m/s/cm] -> [1/s] for Richards solver Srt
end
