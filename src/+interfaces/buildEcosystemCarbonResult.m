function carbon = buildEcosystemCarbonResult(built, plantCandidate, nitrogenCandidate, dt_s)
%BUILDECOSYSTEMCARBONRESULT Actual accepted-interval atmospheric C exchange.
% Inputs: built.Y accepted start state; PlantCN two-slot endpoint; soil-N
% candidate containing final AOM and MEND owners; actual interval [s].
% Outputs: GPP/Ra/Rh/NEE [kg C m^-2 s^-1], Resp [umol CO2 m^-2 s^-1].
% Cumulative differences include microbial N-feedback overflow and crop
% cleanup once. Litter/root/manure transfers are internal, not respiration.
% No empirical soil-respiration formula, state mutation, or pool inversion.
% Editor: Wenjie Liu; Date: 2026-09-30.
    assert(isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'buildEcosystemCarbonResult:duration','An actual positive interval is required.');
    old = built.Y.plant_cn;
    gpp_g_m2 = increment(old.accum.cum_GPP, plantCandidate.cum_GPP, 'GPP');
    ra_g_m2 = increment(old.accum.cum_R_a, plantCandidate.cum_R_a, 'plant respiration');
    aom_g_cm2 = 0; mend_g_cm2 = 0;
    sw = built.control.switches;
    if enabled(sw,'lSoilOrganic')
        assert(isstruct(nitrogenCandidate) && isfield(nitrogenCandidate,'organic') ...
            && ~isempty(nitrogenCandidate.organic), ...
            'buildEcosystemCarbonResult:missingAOM','Enabled AOM requires its final candidate.');
        aom_g_cm2 = increment(built.Y.soil_organic.cumulative.co2_C_g_cm2, ...
            nitrogenCandidate.organic.cumulative.co2_C_g_cm2,'AOM respiration');
    end
    if enabled(sw,'lExplicitSoilMicrobes')
        assert(isstruct(nitrogenCandidate) && isfield(nitrogenCandidate,'microbes') ...
            && ~isempty(nitrogenCandidate.microbes), ...
            'buildEcosystemCarbonResult:missingMEND','Enabled MEND requires its final candidate.');
        previous = built.Y.soil_microbes; next = nitrogenCandidate.microbes.state;
        assert(strcmp(previous.grid_order,'top_down') && strcmp(next.grid_order,'top_down'), ...
            'buildEcosystemCarbonResult:grid','MEND IC must use the declared top-down grid.');
        amount = next.IC(:) - previous.IC(:);
        assert(numel(amount)==built.domain.NL && all(isfinite(amount)) && all(amount>=0), ...
            'buildEcosystemCarbonResult:IC','Cumulative MEND respiratory carbon must not decrease.');
        mend_g_cm2 = sum(flip(amount).*built.parameters.soil.bulk_density(:) ...
            .*built.domain.dz(:))/1000; % mg/g drysoil -> g/cm2 ground
    end
    gpp = gpp_g_m2/(1000*dt_s); ra = ra_g_m2/(1000*dt_s);
    aom = aom_g_cm2*10/dt_s; mend = mend_g_cm2*10/dt_s;
    reco = ra+aom+mend;
    carbon = struct('NEE',reco-gpp,'Resp',reco*1e9/12, ...
        'GPP',gpp,'Ra_plant',ra,'Rh_aom',aom,'Rh_mend',mend, ...
        'Rh',aom+mend,'Reco',reco,'source','accepted_owner_cumulative_v1', ...
        'dt_s',dt_s);
end

function value=increment(previous,next,key)
    assert(isscalar(previous)&&isscalar(next)&&isfinite(previous)&&isfinite(next) ...
        &&next>=previous,'buildEcosystemCarbonResult:cumulative', ...
        'Accepted cumulative %s must be finite and monotone.',key);
    value=next-previous;
end
function value=enabled(sw,key)
    value=isfield(sw,key)&&sw.(key)==1;
end
