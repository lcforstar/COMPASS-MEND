function [data, metadata] = read_excel_data(filename, sheetname, remove_nan_rows)
    %READ_EXCEL_DATA 智能读取 Excel 文件，自动检测和转换数据类型
    %
    % 语法：
    %   [data, metadata] = read_excel_data(filename, sheetname)
    %   [data, metadata] = read_excel_data(filename, sheetname, remove_nan_rows)
    %
    % 输入参数：
    %   filename        - Excel 文件路径 (char/string)
    %   sheetname       - 工作表名称，默认为 'sheet1' (char/string)
    %   remove_nan_rows - 是否删除包含 NaN 的行，默认为 true (logical)
    %
    % 输出参数：
    %   data      - 读取的数据，自动转换类型
    %   metadata  - 元数据结构，包含：
    %       .original_size       - 原始数据大小 [rows, cols]
    %       .final_size          - 清理后数据大小 [rows, cols]
    %       .removed_rows        - 移除的空行索引
    %       .removed_cols        - 移除的空列索引
    %       .removed_nan_rows    - 移除的含 NaN 行索引
    %       .type_conversion     - 类型转换信息
    %
    % 功能说明：
    %   1. 使用 readtable 作为主要方法读取 Excel
    %   2. 自动检测和移除完全空的行和列
    %   3. 自动转换数据类型：
    %      - 纯数值字符串 → double
    %      - 文本 → cell array
    %      - 原生数值 → 保留原样
    %   4. 删除包含 NaN 值的行（可选）
    %   5. 处理 NaN、inf 和空值
    %
    % 示例：
    %   [soil_data, metadata] = read_excel_data('soilproperty.xlsx', 'sheet1');
    %   [init_data, metadata] = read_excel_data('InitialConditions.xlsx', 'sheet1');

    % 处理输入参数
    if nargin < 2
        sheetname = 'sheet1';
    end
    if nargin < 3
        remove_nan_rows = true;  % 默认删除包含 NaN 的行
    end

    % 验证文件存在
    if ~isfile(filename)
        error('文件不存在: %s', filename);
    end

    % 初始化元数据
    metadata = struct();
    metadata.removed_rows = [];
    metadata.removed_cols = [];
    metadata.removed_nan_rows = [];
    metadata.type_conversion = {};

    % 尝试用 readtable 读取（首选方法）
    try
        raw_table = readtable(filename, 'Sheet', sheetname, 'VariableNamingRule', 'preserve', 'UseExcel', false);
        data_cell = table2cell(raw_table);
        col_names = raw_table.Properties.VariableNames;
    catch
        % 备选方案1：使用 xlsread（用于处理混合数据）
        try
            [num, txt, raw] = xlsread(filename, sheetname);
            if isempty(raw)
                error('无法读取数据');
            end
            data_cell = raw;
            col_names = [];
        catch
            % 备选方案2：使用 readcell（最后的手段）
            try
                data_cell = readcell(filename, 'Sheet', sheetname, 'UseExcel', false);
                col_names = [];
            catch ME
                error('无法读取文件 %s: %s', filename, ME.message);
            end
        end
    end

    % 记录原始大小
    [orig_rows, orig_cols] = size(data_cell);
    metadata.original_size = [orig_rows, orig_cols];

    % 第1步：检测并移除完全空的行
    empty_row_mask = false(orig_rows, 1);
    for i = 1:orig_rows
        row_empty = true;
        for j = 1:orig_cols
            cell_val = data_cell{i, j};
            if ~isempty(cell_val) && ~(isnan_value(cell_val))
                row_empty = false;
                break;
            end
        end
        empty_row_mask(i) = row_empty;
    end
    metadata.removed_rows = find(empty_row_mask);
    data_cell(empty_row_mask, :) = [];

    % 第2步：检测并移除完全空的列
    [curr_rows, curr_cols] = size(data_cell);
    empty_col_mask = false(1, curr_cols);
    for j = 1:curr_cols
        col_empty = true;
        for i = 1:curr_rows
            cell_val = data_cell{i, j};
            if ~isempty(cell_val) && ~(isnan_value(cell_val))
                col_empty = false;
                break;
            end
        end
        empty_col_mask(j) = col_empty;
    end
    metadata.removed_cols = find(empty_col_mask);
    data_cell(:, empty_col_mask) = [];

    % 更新列名
    if ~isempty(col_names)
        col_names(empty_col_mask) = [];
    end

    % 第3步：自动类型转换
    [final_rows, final_cols] = size(data_cell);
    data = cell(final_rows, final_cols);

    for j = 1:final_cols
        for i = 1:final_rows
            cell_val = data_cell{i, j};

            % 处理空值
            if isempty(cell_val)
                data{i, j} = NaN;
                continue;
            end

            % 尝试转换为数值
            if isnumeric(cell_val)
                % 已是数值，直接保留
                data{i, j} = cell_val;
            elseif ischar(cell_val) || isstring(cell_val)
                % 字符或字符串，尝试转换为数值
                str_val = char(cell_val);

                % 尝试解析为数值
                if is_numeric_string(str_val)
                    try
                        num_val = str2double(str_val);
                        if ~isnan(num_val)
                            data{i, j} = num_val;
                            % 记录类型转换
                            if j == 1 && i <= 3
                                metadata.type_conversion{end+1} = ...
                                    sprintf('[%d,%d] ''%s'' → double: %g', i, j, str_val, num_val);
                            end
                            continue;
                        end
                    catch
                        % 转换失败，保留为字符
                    end
                end

                % 无法转换为数值，保留为字符
                data{i, j} = str_val;
            else
                % 其他类型，直接保留
                data{i, j} = cell_val;
            end
        end
    end

    % 第4步：删除包含 NaN 的行（如果启用）
    if remove_nan_rows
        nan_row_mask = false(final_rows, 1);
        for i = 1:final_rows
            has_nan = false;
            for j = 1:final_cols
                cell_val = data{i, j};
                if isnumeric(cell_val) && isnan(cell_val)
                    has_nan = true;
                    break;
                end
            end
            nan_row_mask(i) = has_nan;
        end
        metadata.removed_nan_rows = find(nan_row_mask);
        data(nan_row_mask, :) = [];
        final_rows = size(data, 1);
    end

    % 第5步：转换为表格，附加列名
    if ~isempty(col_names)
        data = cell2table(data, 'VariableNames', col_names);
    else
        % 创建默认列名
        col_labels = arrayfun(@(x) sprintf('Col%d', x), 1:final_cols, 'UniformOutput', false);
        data = cell2table(data, 'VariableNames', col_labels);
    end

    % 更新元数据
    metadata.final_size = [final_rows, final_cols];
    metadata.removed_count = length(metadata.removed_rows) + length(metadata.removed_cols) + length(metadata.removed_nan_rows);

    % 打印清理摘要
    if metadata.removed_count > 0
        fprintf('数据清理摘要:\n');
        fprintf('  原始大小: %d × %d\n', orig_rows, orig_cols);
        fprintf('  最终大小: %d × %d\n', final_rows, final_cols);
        if ~isempty(metadata.removed_rows)
            fprintf('  移除的空行: %d\n', length(metadata.removed_rows));
        end
        if ~isempty(metadata.removed_cols)
            fprintf('  移除的空列: %d\n', length(metadata.removed_cols));
        end
        if ~isempty(metadata.removed_nan_rows)
            fprintf('  移除的 NaN 行: %d\n', length(metadata.removed_nan_rows));
        end
        if ~isempty(metadata.type_conversion)
            fprintf('  类型转换示例: %s\n', metadata.type_conversion{1});
        end
    end
end

% =========================================================================
% 辅助函数：检查是否为 NaN 值
% =========================================================================
function tf = isnan_value(val)
    %ISNAN_VALUE 检查值是否为 NaN（处理各种数据类型）

    if isnumeric(val)
        tf = isnan(val);
    elseif ischar(val) || isstring(val)
        str = lower(char(val));
        tf = strcmp(str, 'nan') || strcmp(str, 'n/a') || strcmp(str, '#n/a');
    else
        tf = false;
    end
end

% =========================================================================
% 辅助函数：检查是否为纯数值字符串
% =========================================================================
function tf = is_numeric_string(str)
    %IS_NUMERIC_STRING 检查字符串是否表示数值

    if ~ischar(str) && ~isstring(str)
        tf = false;
        return;
    end

    % 移除前后空格
    str = strtrim(char(str));

    % 空字符串不是数值
    if isempty(str)
        tf = false;
        return;
    end

    % 检查特殊值
    if strcmp(lower(str), 'inf') || strcmp(lower(str), '-inf') || ...
       strcmp(lower(str), 'nan') || strcmp(lower(str), '+inf')
        tf = true;
        return;
    end

    % 尝试转换
    num = str2double(str);
    tf = ~isnan(num);
end
