function driver = buildWaterBoundaryForcing(rawTable)
%BUILDWATERBOUNDARYFORCING  Convert prescribed water-BC table to plain driver data.
%
%   driver = init.buildWaterBoundaryForcing(rawTable) copies the optional
%   water_boundary.csv table into a scalar plain-data struct. The public table
%   contract is time_s plus top_value and bottom_value. The selected endpoint
%   mode assigns value units and interpolation; no interpolation or trial
%   evaluation occurs here.
%
%   Inputs:
%     rawTable [in] table from raw.water_boundary, with time_s [s], top_value,
%                   and bottom_value. Value units are [cm] for head modes and
%                   [cm/s] for flux modes;
%                   [] means no prescribed water-boundary series.
%
%   Returns:
%     driver scalar struct with one column vector per table variable plus
%            n_rows [-]. Empty input returns time_s=[] and n_rows=0.
%
%   Architecture reference:
%     ClimaLand.jl src/shared_utilities/drivers.jl:44-83,175-187.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-11 Jizixing Zhu -- document the single flat [boundary] series-column input fields.
%   Change: 2026-07-11 Jizixing Zhu -- document the fixed two-value CSV contract derived from variable boundary modes.
%   Change: 2026-07-11 Jizixing Zhu -- make the same CSV columns the only value source for constant and variable prescribed modes.

    if isempty(rawTable)
        driver = struct('time_s', zeros(0, 1), 'n_rows', 0);
        return;
    end
    if ~istable(rawTable)
        error('buildWaterBoundaryForcing:badInput', ...
            'raw.water_boundary must be a table or empty.');
    end

    driver = struct();
    names = rawTable.Properties.VariableNames;
    for i = 1:numel(names)
        name = names{i};
        values = rawTable.(name);
        if ~isnumeric(values) || ~isvector(values)
            error('buildWaterBoundaryForcing:badColumn', ...
                'water_boundary column "%s" must be a numeric vector.', name);
        end
        driver.(name) = double(values(:));
    end
    driver.n_rows = height(rawTable);
end
