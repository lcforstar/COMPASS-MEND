function registerNumerics(reg, cfg) %#ok<INUSD>
%REGISTERNUMERICS  Solver diagnostics (group: numerics).
% Output: data/output_<tag>/numerics.nc
    g = 'numerics';

    reg.register(struct('name','step_failed', 'long_name','Picard convergence failure flag', ...
        'units','1', 'dims',{{'time'}}, 'group',g, 'transform','scalar', ...
        'comment','0 = OK; 1 = Picard did not converge within max iters or sub-step shrunk'));
end
