function out = computeCanopyWaterEvap(interception, fluxes, Tss, TT, Ps_max_diurnal, dt)
%COMPUTECANOPYWATEREVAP  Post-ebal wet-canopy evaporation deduction + water balance + soil EVAP.
%
%   out = vegetation.computeCanopyWaterEvap(interception, fluxes, Tss, TT, Ps_max_diurnal, dt)
%   is the PRE-substep surface-water step (legacy COMPASS.m:1425-1550): it evaporates the
%   canopy interception reservoir with the ebal wet-leaf latent flux, closes the canopy
%   water balance (residual), and computes the soil evaporation rate EVAP from the ebal soil
%   latent flux. PURE: no state mutation -- the caller stages the returned reservoir contents
%   into Y.canopy_water (state, accept-only) + built.diagnostics.canopy_water
%   (history row); the caller then advances the ClimaLand bulk snow transaction.
%
%   Two schemes (interception.scheme):
%     'single'        one reservoir; max(0,lE_wet) evaporated (no dew, v29 behaviour).
%     'four_quadrant' sunlit reservoir (no dew) + shaded reservoir (lE<0 -> dew, capped at
%                     max_storage_sha); lE_wet split by fwet x reservoir weight (v30).
%
%   Inputs:
%     interception    [in]  struct  the pre-evap canopy interception (vegetation.computeCanopyInterception):
%                                    .scheme .h2ocan/.h2ocan_sun/.h2ocan_sha [mm] .max_storage_sha [mm]
%                                    .fwet_sun/.fwet_sha [-] .CI_potential [cm/s]
%                                    .drip_mm/.drip_sun_to_sha/.drip_to_ground [mm]
%     fluxes          [in]  struct  ebal fluxes (SCOPE-guarded): .lE_wet_canopy [W/m2], .lEstot [W/m2]
%     Tss             [in]  [degC]  guarded soil-surface temperature (SCOPE-GUARD)
%     TT              [in]  [degC, NN x 1]  soil node temperature (EVAP fires only if any(TT>5))
%     Ps_max_diurnal  [in]  [-]     diurnal sunlit fraction (reservoir split weight)
%     dt              [in]  [s]     outer step length (legacy Delt_t_outer)
%
%   Returns out:
%     .state_h2ocan/.state_h2ocan_sun/.state_h2ocan_sha [mm]  post-evap reservoir contents
%     .evap_wet_mm [mm]        net wet-canopy evaporation this step (may be < 0 = shaded dew)
%     .lE_CI_W [W/m2]          latent flux equivalent of evap_wet_mm
%     .canopy_water_residual [mm]  canopy water-balance residual (|.| should be < 1e-6)
%     .CanopyInterception [cm/s]   evap_wet_mm as a rate (diagnostic)
%     .EVAP [cm/s]             nonnegative soil evaporation demand
%     .DEW [cm/s]              nonnegative bare-surface liquid supply
%   Change: 2026-07-17 Jizixing Zhu -- replace the stale scalar-SWE snow
%     reference with the accepted ClimaLand snow transaction.
%
%   Reference:
%     COMPASS.m:1425-1550 (legacy wet-canopy + soil EVAP block).
%     Canopy lambda_v: fixed 2.45e6 J/kg (~20 degC), same constant as the
%       legacy interception-evaporation conversion (not T-dependent here).
%     Soil lambda_s: (2.501-0.002361*Tss)*1e6 J/kg, Henderson-Sellers 1984
%       form as in heatfluxes.m / heatfluxes_clm45.m (lambda vs Tc).
%     Intentional split: canopy path keeps the fixed SCOPE-style constant;
%       soil EVAP uses Tss-dependent lambda. Unifying canopy to lambda(T) is
%       optional (~1-2% mass effect) and not required for leaf-psi work.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-05 Jizixing Zhu -- doc: caller stages state into Y.canopy_water + history into built.diagnostics.canopy_water (was runtime.canopy_water). Body unchanged (increment 3).
% Change: 2026-07-17 Jizixing Zhu -- use the area-weighted bare-ground
%           latent component when the snow-aware SCOPE adapter provides it.
% Change: 2026-07-17 Jizixing Zhu -- document fixed canopy lambda_v vs T-dependent soil lambda_s.
% Change: 2026-07-18 Jizixing Zhu -- keep signed vapor physics while routing
%           evaporation as demand and bare condensation as liquid supply.
    % Canopy wet-leaf: FIXED latent heat [J/kg]. Not (2.501-0.002361*T)*1e6.
    % ~20 degC representative value; matches legacy COMPASS interception-evap.
    % heatfluxes uses T-dependent lambda for ebal leaf LE; soil EVAP below uses
    % Tss-dependent lambda_s. Deliberate: keep mass conversion simple/stable.
    lambda_v = 2.45e6;
    scheme   = interception.scheme;
    CI_mm    = interception.CI_potential * dt * 10;   % [mm] potential interception this step

    lE_wet = 0;
    hasLEwet = isfield(fluxes, 'lE_wet_canopy') && isfinite(fluxes.lE_wet_canopy);
    if hasLEwet, lE_wet = fluxes.lE_wet_canopy; end

    if strcmpi(scheme, 'single')
        % --- v29 single reservoir: max(0,lE_wet) evaporated (no dew) ---
        h2o_before = interception.h2ocan;    % [mm] post-interception, pre-evap
        if hasLEwet
            evap_wet_mm  = min(max(0, lE_wet) / lambda_v * dt, h2o_before);   % [mm], capped at reservoir
            state_h2ocan = max(0, h2o_before - evap_wet_mm);
            lE_CI_W      = evap_wet_mm * lambda_v / dt;
        else
            evap_wet_mm  = 0;
            state_h2ocan = h2o_before;
            lE_CI_W      = 0;
        end
        state_h2ocan_sun = 0;
        state_h2ocan_sha = 0;
        evap_wet_sun_mm  = 0;
        evap_wet_sha_mm  = 0;
        dew_sha_mm       = 0;

        % conservation: reconstruct the pre-interception reservoir and close the balance
        h2o_at_start   = h2o_before - (CI_mm - interception.drip_mm);
        delta_storage  = state_h2ocan - h2o_at_start;
        residual_mm    = CI_mm - (interception.drip_mm + evap_wet_mm) - delta_storage;
    else
        % --- v30 four-quadrant: split lE_wet by fwet x reservoir weight; sha allows dew ---
        wt_sun = interception.fwet_sun * Ps_max_diurnal;
        wt_sha = interception.fwet_sha * (1 - Ps_max_diurnal);
        wt_tot = wt_sun + wt_sha;
        if wt_tot > 1e-12
            lE_wet_sun = lE_wet * wt_sun / wt_tot;
            lE_wet_sha = lE_wet * wt_sha / wt_tot;
        else
            lE_wet_sun = 0;
            lE_wet_sha = 0;
        end

        % sunlit: no dew (daytime PAR keeps leaf above dew point)
        h2o_sun_before   = interception.h2ocan_sun;
        evap_wet_sun_mm  = min(max(0, lE_wet_sun) / lambda_v * dt, h2o_sun_before);
        state_h2ocan_sun = max(0, h2o_sun_before - evap_wet_sun_mm);

        % shaded: lE<0 -> dew adds to reservoir, capped at max_storage_sha
        h2o_sha_before      = interception.h2ocan_sha;
        evap_sha_signed_mm  = lE_wet_sha / lambda_v * dt;   % [mm], negative = dew
        if evap_sha_signed_mm >= 0
            evap_wet_sha_mm  = min(evap_sha_signed_mm, h2o_sha_before);
            state_h2ocan_sha = max(0, h2o_sha_before - evap_wet_sha_mm);
            dew_sha_mm       = 0;
        else
            dew_sha_mm       = -evap_sha_signed_mm;
            state_h2ocan_sha = h2o_sha_before + dew_sha_mm;
            if state_h2ocan_sha > interception.max_storage_sha
                dew_sha_mm       = interception.max_storage_sha - h2o_sha_before;   % actual dew accepted
                state_h2ocan_sha = interception.max_storage_sha;
            end
            evap_wet_sha_mm  = -dew_sha_mm;   % negative = condensation
        end

        state_h2ocan = state_h2ocan_sun + state_h2ocan_sha;
        evap_wet_mm  = evap_wet_sun_mm + evap_wet_sha_mm;   % includes dew (negative)
        lE_CI_W      = evap_wet_mm * lambda_v / dt;

        % dual-reservoir conservation (legacy COMPASS.m:1509-1521)
        h2o_start_sun = h2o_sun_before - (CI_mm * Ps_max_diurnal - interception.drip_sun_to_sha);
        h2o_start_sha = h2o_sha_before - (CI_mm * (1 - Ps_max_diurnal) + ...
            interception.drip_sun_to_sha - interception.drip_to_ground);
        delta_storage = (state_h2ocan_sun - h2o_start_sun) + (state_h2ocan_sha - h2o_start_sha);
        water_out     = interception.drip_to_ground + evap_wet_sun_mm + max(0, evap_wet_sha_mm);
        residual_mm   = CI_mm + max(0, dew_sha_mm) - water_out - delta_storage;
    end

    % --- soil evaporation from the ebal soil latent flux (guarded) ---
    % Soil path IS temperature-dependent (unlike canopy lambda_v above).
    % Same linear form as heatfluxes.m: lambda = (2.501 - 0.002361*T)*1e6 [J/kg].
    lambda_s = (2.501 - 0.002361 * Tss) * 1e6;
    bareVaporMass = NaN;
    if isfield(fluxes, 'bare_ground_upward_water_vapor_flux_kg_m2_s')
        bareVaporMass = ...
            fluxes.bare_ground_upward_water_vapor_flux_kg_m2_s;
    elseif isfield(fluxes, 'qv_bare_net_kg_m2_s')
        bareVaporMass = fluxes.qv_bare_net_kg_m2_s;
    elseif isfield(fluxes, 'evaporation_bare_kg_m2_s')
        bareVaporMass = fluxes.evaporation_bare_kg_m2_s;
    end
    bareDewMass = 0;
    if isfinite(bareVaporMass)
        bareDewMass = max(-bareVaporMass, 0);
    elseif isfield(fluxes, 'bare_dew_kg_m2_s') ...
            && isfinite(fluxes.bare_dew_kg_m2_s)
        bareDewMass = max(fluxes.bare_dew_kg_m2_s, 0);
    end
    lEstot = NaN;
    if isfield(fluxes, ...
            'bare_ground_upward_latent_heat_flux_contribution_W_m2')
        lEstot = fluxes. ...
            bare_ground_upward_latent_heat_flux_contribution_W_m2;
    elseif isfield(fluxes, 'lEstot_bare')
        lEstot = fluxes.lEstot_bare;
    elseif isfield(fluxes, 'lEstot')
        lEstot = fluxes.lEstot;
    end
    if isfinite(bareVaporMass) && any(TT > 5)
        EVAP = max(bareVaporMass, 0) * 0.1; % kg m-2 s-1 -> cm s-1
    elseif isfinite(lambda_s) && lambda_s > 0 && isfinite(lEstot) && lEstot < 800 && lEstot > -300 && any(TT > 5)
        EVAP = max(lEstot, 0) / lambda_s * 0.1;   % [cm/s] (kg/m2/s -> mm/s -> x0.1 cm/s)
    else
        EVAP = 0;
    end

    out = struct( ...
        'state_h2ocan',          state_h2ocan, ...
        'state_h2ocan_sun',      state_h2ocan_sun, ...
        'state_h2ocan_sha',      state_h2ocan_sha, ...
        'evap_wet_mm',           evap_wet_mm, ...
        'lE_CI_W',               lE_CI_W, ...
        'canopy_water_residual', residual_mm, ...
        'CanopyInterception',    evap_wet_mm / dt / 10, ...   % [cm/s]
        'EVAP',                  EVAP, ...
        'DEW',                   bareDewMass * 0.1);
end
