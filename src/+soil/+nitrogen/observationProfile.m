function [nh4,no3,audit]=observationProfile(Y,parameters,dz)
%OBSERVATIONPROFILE Accepted dissolved+adsorbed DIN on dry-soil basis.
% Y is one accepted column; nh4/no3 liquid NN-by-1 [ug/cm3 water]. The
% authoritative N effective-liquid snapshot is NL-by-2, bottom-up. Parameters
% RHO_bulk [g/cm3] and RHOKG [-] are NL-by-1; dz [cm] is positive NL-by-1.
% Outputs are bottom-up NN-by-1 [ug/g = mg/kg dry soil]. Each incident half
% element contributes to nodal capacity and dry mass, so summing profile
% times dry mass recovers the FE dissolved+adsorbed inventory exactly.
% Ice-entrapped N, unreacted granules and urea are excluded from this observable;
% they remain in the physical conservation ledger. No fallback to diagnostic
% total-N slots or initialization-only Y.theta; missing snapshots are errors.
% Does not change state, chemistry, transport or boundary conditions.
% Editor: Wenjie Liu; Date: 2026-09-29.
    dz=dz(:); nl=numel(dz); nn=nl+1;
    assert(all(isfinite(dz))&&all(dz>0)&&nl>=1, ...
        'soil:nitrogen:observationGeometry','Positive finite layer thickness required.');
    rho=parameters.RHO_bulk(:); buffer=parameters.RHOKG(:);
    assert(numel(rho)==nl&&numel(buffer)==nl ...
        &&all(isfinite(rho))&&all(rho>0)&&all(isfinite(buffer))&&all(buffer>=0), ...
        'soil:nitrogen:observationParameters','Invalid layer density or sorption buffer.');
    assert(isfield(Y,'soil_nitrogen')&&isfield(Y.soil_nitrogen,'theta_liquid_prev'), ...
        'soil:nitrogen:observationSnapshot','Accepted N liquid snapshot is required.');
    theta=Y.soil_nitrogen.theta_liquid_prev;
    assert(isequal(size(theta),[nl,2])&&all(isfinite(theta),'all') ...
        &&all(theta>=0,'all')&&numel(Y.nh4)==nn&&numel(Y.no3)==nn ...
        &&all(isfinite(Y.nh4(:)))&&all(Y.nh4(:)>=0) ...
        &&all(isfinite(Y.no3(:)))&&all(Y.no3(:)>=0), ...
        'soil:nitrogen:observationState','Invalid accepted N concentrations or water snapshot.');
    audit.node_volume_cm=lump(ones(nl,2),dz);
    audit.dry_mass_g_cm2=lump(repmat(rho,1,2),dz);
    audit.nh4_capacity_cm=lump(theta+buffer,dz);
    audit.no3_capacity_cm=lump(theta,dz);
    nh4=Y.nh4(:).*audit.nh4_capacity_cm./audit.dry_mass_g_cm2;
    no3=Y.no3(:).*audit.no3_capacity_cm./audit.dry_mass_g_cm2;
    audit.nh4_inventory_ug_cm2=sum(nh4.*audit.dry_mass_g_cm2);
    audit.no3_inventory_ug_cm2=sum(no3.*audit.dry_mass_g_cm2);
    audit.definition='dissolved_plus_equilibrium_adsorbed_excludes_ice_granules_urea';
    audit.schema_version=1;
end

function value=lump(endpoint,dz)
    contribution=endpoint.*(dz/2);
    value=[contribution(:,1);0]+[0;contribution(:,2)];
end
