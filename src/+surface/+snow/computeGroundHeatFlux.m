function flux_W_m2 = computeGroundHeatFlux(built, snowState, parameters)
%COMPUTEGROUNDHEATFLUX  Conductive heat from bulk snow to the top soil node.
%
%   Positive flux is downward into soil. This is the sign-reversed form of
%   ClimaLand integrated/soil_snow_model.jl::update_soil_snow_ground_heat_flux!,
%   including its 0.1-m minimum exchange depth for the single-layer snowpack.
%   COMPASS has an explicit z=0 soil-surface node; this scalar is the reference
%   evaluation of the affine interface response at that node.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-07-17 Jizixing Zhu -- map the COMPASS top element to the
%           ClimaLand top-cell-centre geometry instead of treating the
%           surface Dirichlet node as a cell centre.
%   Change: 2026-07-19 Jizixing Zhu -- delegate to the z=0 surface-node
%           affine response used by the coupled U-primary boundary row.

    response = surface.snow.computeGroundHeatResponse(built, snowState, parameters);
    flux_W_m2 = response.reference_flux_into_soil_W_m2;
end
