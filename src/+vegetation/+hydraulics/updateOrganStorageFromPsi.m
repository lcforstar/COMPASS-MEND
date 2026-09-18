function [W_new, psi_v2_MPa, diag] = updateOrganStorageFromPsi( ...
    W_old, hydraulic, SoilVariables, plantHydro, LAI_now, Delt_t) %#ok<INUSL>
%UPDATEORGANSTORAGEFROMPSI  Diagnose organ water storage from SPAC ψ.
%
%   Active when hydraulic_scheme is storage_from_psi or storage_dynamic.
%   Leaf/root water potentials are solved OUTSIDE this function by the steady
%   SPAC path in ebal (calSPAC: soil-root conductance, shared supply, path
%   drops and vulnerability). This routine only maps those algebraic ψ values
%   to organ water storage W for diagnostics (VOD / RWC). It does not solve
%   SPAC again and does not feed W back into the current ebal Picard.
%
%   Method (constant-capacitance linear storage, not full PV curve):
%     RWC = W / W_cap = 1 + a * ψ_m ,   ψ_m <= 0  [m water column]
%     W_eq = W_cap * (1 + a * ψ_m) ,    clamped to [0, W_cap]
%   Equivalently capacitance C = dW/dψ_m = a * W_cap  [m per m-head].
%   Parameters:
%     W_cap (= plantHydro.W_capacity)  [m]  saturated organ storage depth
%     a     (= plantHydro.dyn_a_retention) [m^-1]  dRWC/dψ; empty ψ = -1/a
%   This is the same family as constant-C organ storage in land models
%   (e.g. m = m_max + C*ψ*scale in ORCHIDEE-CAN-NHA; Yao et al. 2021 GMD
%   discussion, §2.1.4.1). It is NOT a full Tyree-Hammel pressure-volume
%   analysis (no turgor-loss breakpoint; not -1/ψ vs RWC below TLP).
%
%   Stem ψ is not a SPAC state here: approximate mid-path value
%     ψ_stem ≈ 0.5*(ψ_root + ψ_leaf) - 0.5*ρgH_eff
%   with H_eff = dz_rs + dz_sl when available.
%
%   Optional stem lag (dyn_tau_stem_s = τ > 0), diagnostic only:
%     dW_stem/dt = (W_eq_stem - W_stem)/τ
%     W_stem^{n+1} = W_old + (W_eq - W_old)*(1 - exp(-Δt/τ))
%   Root and leaf stay at W_eq. τ = 0 => W = W_eq (instantaneous).
%   Lagged stem W does not alter ebal ψ (avoids noon ψ-W mismatch).
%
%   Returns:
%     W_new      [m, 3x1]   [W_root; W_stem; W_leaf]
%     psi_v2_MPa [MPa, 3x1] Echo of input/derived ψ used for the map
%     diag       [struct]   W, ψ, and RWC = W/W_cap diagnostics
%
%   Inputs:
%     W_old        [m, 3x1]  Previous W (used only if stem τ > 0)
%     hydraulic    [struct]  Post-ebal SPAC fields .psil, .psir [MPa]
%     SoilVariables [struct] Unused (call-site compatibility)
%     plantHydro   [struct]  .dyn_a_retention [3x1 m^-1], .W_capacity [3x1 m],
%                            optional .dz_rs/.dz_sl [m], .dyn_tau_stem_s [s]
%     LAI_now      [m2/m2]   Unused (reserved for biomass-based W_cap later)
%     Delt_t       [s]       Outer-step length for stem lag
%
%   Reference:
%     Path ψ: vegetation.hydraulics.calSPAC (steady SPAC, not time-integrated C).
%     Linear C storage: land-model constant capacitance; Yao et al., GMD
%       discussions 2021, ORCHIDEE-CAN-NHA r7236, §2.1.4.1 (m_max + C*ψ*scale).
%     Full leaf PV curves (Tyree & Hammel 1972; PROMETHEUS protocols) differ
%       and are not implemented here.
%     CLM5 PHS: algebraic plant ψ, no tissue capacitance (contrast).
%
%   Former name: Plant_hydraulics_v2_wrapper (deprecated alias still available).
%   Editor: Jizixing Zhu; Date: 2026-05-28.
%   Change: 2026-07-17 Jizixing Zhu — rename from Plant_hydraulics_v2_wrapper.
%   Change: 2026-07-17 Jizixing Zhu — rewrite docs: constant-C diagnostic map,
%     SPAC ψ ownership, parameter meanings, not full PV; NHA-style citation.

    m_to_MPa = 9.81e-3;       % [MPa m^-1] 1 m water column ≈ 9.81e-3 MPa

    % --- SPAC ψ from ebal (path-connected root-leaf), not re-solved here -----
    psi_leaf_MPa = hydraulic.psil;   % [MPa] leaf (sun/sha already aggregated)
    psi_root_MPa = hydraulic.psir;   % [MPa] root
    % Stem ψ: calSPAC has no explicit stem state; mid-path + half gravity drop.
    gravity_MPa = 1000 * 9.80665 * 1.0 / 1e6;   % [MPa] ρ_w g H for H=1 m
    if isfield(plantHydro, 'dz_rs') && isfield(plantHydro, 'dz_sl')
        H_eff = plantHydro.dz_rs + plantHydro.dz_sl;   % [m] root-stem + stem-leaf
        gravity_MPa = 1000 * 9.80665 * H_eff / 1e6;
    end
    psi_stem_MPa = 0.5 * (psi_root_MPa + psi_leaf_MPa) - 0.5 * gravity_MPa;

    % Retention slope a is in m^-1, so convert ψ to meters of water head.
    psi_root_m = psi_root_MPa / m_to_MPa;   % [m]
    psi_stem_m = psi_stem_MPa / m_to_MPa;
    psi_leaf_m = psi_leaf_MPa / m_to_MPa;

    % --- Constant-C inverse map: W_eq = W_cap * (1 + a * ψ_m) ----------------
    % W_eq is the storage in instantaneous equilibrium with this step's SPAC ψ.
    % Not a full PV curve; no TLP piecewise branch.
    a = plantHydro.dyn_a_retention(:);    % [m^-1]  dRWC/dψ_m
    Wcap = plantHydro.W_capacity(:);      % [m]    saturated organ storage

    W_eq = zeros(3, 1);
    W_eq(1) = Wcap(1) * (1 + a(1) * psi_root_m);
    W_eq(2) = Wcap(2) * (1 + a(2) * psi_stem_m);
    W_eq(3) = Wcap(3) * (1 + a(3) * psi_leaf_m);
    W_eq = max(0, min(W_eq, Wcap));       % [m] physical clamp

    % --- Optional stem first-order lag (VOD diurnal; does not change ebal ψ) -
    % dW_stem/dt = (W_eq_stem - W_stem)/τ
    % W_stem^{n+1} = W_old + (W_eq - W_old)*(1 - exp(-Δt/τ))
    tau_stem = 0;
    if isfield(plantHydro, 'dyn_tau_stem_s') && ~isempty(plantHydro.dyn_tau_stem_s)
        tau_stem = plantHydro.dyn_tau_stem_s;   % [s]
    end

    W_new = W_eq;   % root + leaf: instantaneous equilibrium with ψ
    if tau_stem > 0 && ~isempty(W_old) && numel(W_old) >= 2 && isfinite(W_old(2))
        relax = 1 - exp(-Delt_t / tau_stem);
        W_new(2) = W_old(2) + (W_eq(2) - W_old(2)) * relax;
        W_new(2) = max(0, min(W_new(2), Wcap(2)));
    end

    % Diagnostic ψ copy (same numbers as SPAC inputs / stem mid-path approx).
    psi_v2_MPa = [psi_root_MPa; psi_stem_MPa; psi_leaf_MPa];

    diag.W_root              = W_new(1);
    diag.W_stem              = W_new(2);
    diag.W_leaf              = W_new(3);
    diag.psi_root_v2         = psi_root_MPa;
    diag.psi_stem_v2         = psi_stem_MPa;
    diag.psi_leaf_v2         = psi_leaf_MPa;
    diag.RWC_root            = W_new(1) / max(Wcap(1), 1e-9);   % [-] W/W_cap
    diag.RWC_stem            = W_new(2) / max(Wcap(2), 1e-9);
    diag.RWC_leaf            = W_new(3) / max(Wcap(3), 1e-9);
    diag.psi_leaf_v1         = hydraulic.psil;   % same as leaf SPAC ψ
    diag.delta_psi_leaf_v1v2 = 0;                % map does not re-solve leaf ψ
end
