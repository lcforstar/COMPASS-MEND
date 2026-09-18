function tbl = readTableFile(filepath)
%READTABLEFILE  Read a site-package table from CSV or Excel.
%
%   tbl = io.adapters.readTableFile(filepath) reads a simple table file into a
%   MATLAB table without grid construction, interpolation, or unit conversion.
%   It is used for optional site inputs that may be supplied as CSV, XLSX, or
%   XLS files. CSV files use the COMPASS comment-aware CSV adapter; Excel files
%   read sheet 1 with preserved variable names.
%
%   Inputs:
%     filepath  [char]  path to a .csv, .xlsx, or .xls table file
%
%   Returns:
%     tbl       table   raw table at native file granularity.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if exist(filepath, 'file') ~= 2
        error('readTableFile:notFound', 'Table file not found: %s', filepath);
    end

    [~, ~, ext] = fileparts(filepath);
    switch lower(ext)
        case '.csv'
            tbl = io.adapters.readCsvTable(filepath);
        case {'.xlsx', '.xls'}
            tbl = readtable(filepath, 'Sheet', 1, 'NumHeaderLines', 1, ...
                'VariableNamingRule', 'preserve');
        otherwise
            error('readTableFile:unsupportedExtension', ...
                'Unsupported table extension "%s" for %s.', ext, filepath);
    end
end
