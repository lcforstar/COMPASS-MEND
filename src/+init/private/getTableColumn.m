function v = getTableColumn(tbl, name, dflt)
% Read a named table column as a column vector, or return dflt if absent.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ismember(name, tbl.Properties.VariableNames)
        v = tbl.(name)(:);
    else
        v = dflt(:);
    end
end
