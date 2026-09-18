function latent_heat_J_kg = freezeLatentHeat(switches)
%FREEZELATENTHEAT  Resolve the per-step water-fusion latent heat.
%
%   Input:
%     switches [in] Model physics switches including lFreezing.
%
%   Return:
%     latent_heat_J_kg [J/kg] Fusion latent heat, zero when freezing is off.
%
%   The environment override is retained for controlled diagnostics and
%   reproduces the previous outer-step staging behavior.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    constants = equations.physicalConstants();
    if isfield(switches, 'lFreezing') && switches.lFreezing == 1
        latent_heat_J_kg = constants.lambda_fus;
    else
        latent_heat_J_kg = 0;
    end
    override = str2double(io.getEnv('LFUSION_OVERRIDE'));
    if isfinite(override)
        latent_heat_J_kg = override;
    end
end

