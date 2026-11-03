function flags=numericalOptions(control)
%NUMERICALOPTIONS Actual upwind/artificial-dispersion bits for the DIN ADE.
% Parsed site numerics take precedence over historical switch-only inputs.
% Omitting both sources retains the old zero default. No state is mutated.
% Editor: Wenjie Liu; Date: 2026-11-03.
    flags=struct('lUpW',0,'lArtD',0);
    for key={'lUpW','lArtD'}
        name=key{1};
        if isfield(control,'numerics') && isfield(control.numerics,name)
            value=control.numerics.(name);
        elseif isfield(control,'switches') && isfield(control.switches,name)
            value=control.switches.(name);
        else
            continue;
        end
        assert((isnumeric(value)||islogical(value))&&isscalar(value)&&isreal(value) ...
            &&(value==0||value==1),'soil:nitrogen:transportOption','%s must be 0 or 1.',name);
        flags.(name)=double(value);
    end
end
