function summary = finalizeFluxSummary(summary)
%FINALIZEFLUXSUMMARY  Convert accepted weighted sums to mean boundary rates.
%
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-19 Jizixing Zhu -- finalize the accepted time-mean soil-surface
%           temperature used by the affine top-energy response.

    for name = {'water','heat','dry_air'}
        field = name{1};
        process = summary.(field);
        if process.enabled && process.weight_sum > 0
            process.top_flux_into_domain = ...
                process.top_flux_into_domain / process.weight_sum;
            process.bottom_flux_into_domain = ...
                process.bottom_flux_into_domain / process.weight_sum;
            if strcmp(field, 'heat') && process.top_response_weight_sum > 0
                process.top_response_evaluation_temperature_C = ...
                    process. ...
                    top_response_evaluation_temperature_weighted_sum_C ...
                    / process.top_response_weight_sum;
                process.top_response_evaluated_flux_into_domain = ...
                    process.top_response_evaluated_flux_weighted_sum ...
                    / process.top_response_weight_sum;
            end
        end
        summary.(field) = process;
    end
end
