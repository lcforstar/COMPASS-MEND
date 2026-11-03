function out = solveCanopySoilEnergy(in)
%SOLVECANOPYSOILENERGY  Run one canopy-soil energy-balance solve from target input.
%
%   out = vegetation.energy.solveCanopySoilEnergy(in) selects vegetation.energy.ebal or
%   vegetation.energy.ebalProfile from control.canopy_air_profile_eb and returns the
%   kernel output. No model mutation. Dispatch only.
%
%   Returns:
%     out [struct]  iteration, fluxes, radiation, gap, thermal, profiles,
%                   soil_surface, root_water_uptake, frac, hydraulic
%
%   Inputs:
%     in [in] struct  From vegetation.energy.assembleEbalInput (EBAL_IO_CONTRACT §5.2)
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-16 Jizixing Zhu — Phase B: thin dispatch to ebal(in).
%   Change: 2026-07-16 Jizixing Zhu — rename kernels to ebal / ebalProfile.

    assert(isstruct(in) && isfield(in, 'soil') && isfield(in, 'scope_step'), ...
        'solveCanopySoilEnergy:badInput', ...
        'in must be a target-name struct from assembleEbalInput.');

    if logical(in.control.canopy_air_profile_eb)
        out = vegetation.energy.ebalProfile(in);
    else
        out = vegetation.energy.ebal(in);
    end
end
