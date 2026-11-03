function [value,detail]=docFreeLayerMean(built,topCm,bottomCm)
%DOCFREELAYERMEAN Accepted dissolved DOC, excluding QOC, in mg C/kg dry soil.
% Uses the production phase reader and accepted liquid-water cache. A layer
% interval [topCm,bottomCm] is integrated by exact overlap thickness times
% dry bulk density, not by averaging concentrations or layer centers.
% Read-only: D+Q remains the sole DOM state owner. MOC+Q is MAOC only as a
% diagnostic sum, never a second state or an extra SOC contribution.
% Editor: Wenjie Liu; Date: 2026-10-12.
    assert(isscalar(topCm)&&isscalar(bottomCm)&&isfinite(topCm) ...
        &&isfinite(bottomCm)&&topCm>=0&&bottomCm>topCm, ...
        'diagnostics:docFreeSupport','Require a finite positive depth interval.');
    z=built.domain.z_node(:);rho=built.parameters.soil.bulk_density(:);
    assert(all(diff(z)<0)&&topCm>=z(end)&&bottomCm<=z(1), ...
        'diagnostics:docFreeSupport','Depth interval must lie within the bottom-up column.');
    overlap=max(0,min(z(1:end-1),bottomCm)-max(z(2:end),topCm));
    theta=mean(built.workspace.derived.Theta_LL,2);
    phase=soil.dom.phaseState(built.Y.soil_dom,'doc',theta,rho,built.parameters.soil_dom);
    dryMass=sum(overlap.*rho);
    assert(dryMass>0,'diagnostics:docFreeSupport','Depth interval has zero dry-soil mass.');
    value=sum(overlap.*phase.free_g_cm3_bulk)/dryMass*1e6;
    detail=struct('doc_contract_version','doc_free_d_v2','units','mg C/kg dry soil', ...
        'depth_support_cm',[topCm bottomCm],'aggregation','overlap_dry_soil_mass_weighted', ...
        'doc_free_mg_kg',value,'qoc_mg_kg', ...
        sum(overlap.*phase.sorbed_g_cm3_bulk)/dryMass*1e6, ...
        'doc_plus_qoc_mg_kg',sum(overlap.*phase.total_g_cm3_bulk)/dryMass*1e6);
    if isfield(built.Y,'soil_microbes')
        moc=flip(built.Y.soil_microbes.MOC(:)).*rho/1000;
        detail.moc_mg_kg=sum(overlap.*moc)/dryMass*1e6;
        detail.maoc_mg_kg=detail.moc_mg_kg+detail.qoc_mg_kg;
    end
end
