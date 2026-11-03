function a=componentTotals(a,checkExisting)
%COMPONENTTOTALS Derive legacy totals from prognostic litter components.
% Components: metabolic, cellulose, lignin, CWD. Surface [g/cm2], element
% [g/cm3], bottom-up. Totals are caches, never additional C/N stores.
% Editor: Jizixing Zhu; Date: 2026-09-08.
    names={'metabolic','cellulose','lignin','cwd'};
    fields={'surface_C','surface_N','element_C','element_N'};
    for j=1:numel(fields)
        key=fields{j}; value=zeros(size(a.(key)));
        for i=1:numel(names)
            x=a.components.(names{i}).(key);
            assert(isequal(size(x),size(value))&&all(isfinite(x(:))) ...
                &&all(x(:)>=0),'litter:badComponent','Invalid component state.');
            value=value+x;
        end
        if checkExisting
            assert(all(abs(a.(key)(:)-value(:))<=1e-12*max(1,max(value(:)))), ...
                'litter:componentMismatch','Litter totals disagree with components.');
        end
        a.(key)=value;
    end
end
