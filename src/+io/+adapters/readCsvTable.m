function tbl = readCsvTable(filepath)
%READCSVTABLE  Read a COMPASS site CSV into a table (stage-1 adapter, pure read).
%
%   tbl = io.adapters.readCsvTable(filepath) reads a comma-separated table,
%   skipping '#' comment lines, and returns a MATLAB table with one column per
%   header field. Used for soil.csv (one row per horizon),
%   soil_initial_profile.csv (one row per depth), plant_initial_pools.csv
%   (pool/value/unit long form), and root.csv.
%   No interpolation / grid mapping / unit conversion (those live in buildState).
%
%   Inputs:
%     filepath  [char]  path to a .csv file (first non-'#' line = header names)
%
%   Returns:
%     tbl       table   variables named by the header row; numeric columns
%                       become double, text columns become string/cellstr.
%
%   Convention: lines starting with '#' (after optional leading spaces) are
%   comments and are skipped, so column metadata can live as '# ...' headers.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-16 Jizixing Zhu -- document the descriptive site input filenames.

    if exist(filepath, 'file') ~= 2
        error('readCsvTable:notFound', 'CSV file not found: %s', filepath);
    end
    opts = detectImportOptions(filepath, 'CommentStyle', '#');
    opts.VariableNamingRule = 'preserve';   % keep header names verbatim (depth_cm, theta_s, ...)
    tbl = readtable(filepath, opts);
end
