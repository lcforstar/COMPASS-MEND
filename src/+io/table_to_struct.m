function opts = table_to_struct(tbl, value_col, fieldname_col, convert_zero_to_empty)
    % TABLE_TO_STRUCT 将表格转换为结构体
    %
    % 语法：
    %   opts = table_to_struct(tbl)              % 默认：第1列为值，第2列为字段名
    %   opts = table_to_struct(tbl, 1, 2)      % 指定列索引
    %   opts = table_to_struct(tbl, 1, 2, true) % 并将值为0的转换为[]
    %
    % 输入参数：
    %   tbl                    - 输入表格 (table)
    %   value_col              - 值所在列 (默认: 1)
    %   fieldname_col          - 字段名所在列 (默认: 2)
    %   convert_zero_to_empty  - 是否将值为0的转换为[] (默认: false)
    %
    % 输出参数：
    %   opts         - 输出结构体 (struct)
    %
    % 示例：
    %   opts = table_to_struct(options_table)  % 第1列为值，第2列为字段名
    %   opts = table_to_struct(options_table, 2, 1)  % 第2列为值，第1列为字段名
    %   opts = table_to_struct(options_table, 1, 2, true)  % 0转换为[]

    % 默认参数
    if nargin < 2
        value_col = 1;
    end
    if nargin < 3
        fieldname_col = 2;
    end
    if nargin < 4
        convert_zero_to_empty = false;
    end

    % 初始化结构体
    opts = struct();

    % 获取表格的行数
    n_rows = height(tbl);

    % 遍历每一行，构建结构体
    for i = 1:n_rows
        % 获取字段名和值
        fieldname = tbl{i, fieldname_col};
        value = tbl{i, value_col};

   
        % 确保字段名是字符串
        if iscell(fieldname)
            fieldname = fieldname{1};
        end
        if isstring(fieldname)
            fieldname = char(fieldname);
        end

        % 清理字段名：移除空格，替换特殊字符
        fieldname = strtrim(fieldname);
        fieldname = regexprep(fieldname, '\s+', '_');  % 空格转下划线
        fieldname = regexprep(fieldname, '[^\w]', '');  % 删除特殊字符

        % 如果启用了转换功能，将值为0的转换为[]
        if convert_zero_to_empty
            % 处理不同类型的0值
            if iscell(value) && ~isempty(value)
                val_unwrapped = value{1};
            else
                val_unwrapped = value;
            end

            % 检查是否为0（数字或字符串）
            is_zero = false;
            if isnumeric(val_unwrapped) && val_unwrapped == 0
                is_zero = true;
            elseif ischar(val_unwrapped) && strcmp(val_unwrapped, '0')
                is_zero = true;
            elseif isstring(val_unwrapped) && strcmp(string(val_unwrapped), "0")
                is_zero = true;
            end

            if is_zero
                value = [];
            end
        end

        % 赋值到结构体
        opts.(fieldname) = value;
    end

    fprintf('[Table to Struct] 转换完成: %d 行，字段数: %d\n', n_rows, length(fieldnames(opts)));
end
