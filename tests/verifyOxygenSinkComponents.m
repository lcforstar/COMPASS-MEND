function report=verifyOxygenSinkComponents()
% One native fixed-O2 tile exercises the independent nitrification ledger
% even when adding respiration consumes all representable digits of its rate.
% A zero-O2 tile must still reject the positive respiratory request.
% Editor: Wenjie Liu; Date: 2026-11-03.
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
    built.domain=struct('NL',1,'NN',2,'dz',10);
    built.parameters.soil.theta_s=.5;
    built.parameters.soil_n=struct('ph_response_enabled',false,'Theta_f',.25, ...
        'RHOKG',1,'Vmax_nit',1,'Km_nit',1,'pH',[7;7]);
    built.Y.soil_gas.y_O2=[.2;.2];
    phase=struct('theta_g',[.2 .2],'theta_l',[.3 .3],'porosity',.5, ...
        'temperature_C',[20;20],'gas_pressure_Pa',[101325;101325],'vapor_pressure_Pa',[0;0]);
    tr=struct('free_dry_air_mass_flux_g_cm2_s',0,'gas_dispersion_cm2_s',0,'liquid_flux_cm_s',0);
    endpoint=struct('kind','fixed_dry_mole_fraction','composition_source','prescribed', ...
        'requested_external_dry_mole_fraction',.2,'requested_flux_into_domain_mol_cm2_s',0);
    applied=struct('applied_kind','fixed_dry_mole_fraction','dry_mole_fraction',.2);
    tile=struct('t_start_s',0,'dt_s',3600,'phase_old',phase,'phase_new',phase, ...
        'transport',tr,'boundary_request',struct('t_end_s',3600,'top',endpoint,'bottom',endpoint), ...
        'boundary_applied',struct('species','O2','top',applied,'bottom',applied), ...
        'nh4',[1;1],'matric_head_cm',[-50;-50],'dry_air_storage_after',[.001 .001], ...
        'QL',0,'QL_h',0,'QL_T',0,'QL_a',0);
    hour=struct('schema_version',1,'tile_count',1,'tiles',{{tile}},'elapsed_s',3600);
    cap=[1e-23 1e-23]; % ug N cm^-3 s^-1 at both endpoints
    respiration=4e-13*12e6; % ug C cm^-3 s^-1
    [total,nit]=soil.gas.calcO2Sink(cap,respiration,1);
    assert(nit>1e-30 && total-respiration*1e-6/12==0);
    gas=soil.gas.replayOxygenInterval(built,hour,respiration,cap);
    assert(gas.feasible && gas.accepted_tiles==1 && isequal(gas.Snit,cap));
    assert(abs(gas.ledger.nitrification_actual_by_layer_mol_cm3-nit*3600)<=eps(nit*3600));
    assert(abs(gas.ledger.balance_residual_mol_cm2)<1e-14);
    empty=built; empty.Y.soil_gas.y_O2=[0;0];
    before=empty;
    hour.tiles{1}.boundary_applied.top.dry_mole_fraction=0;
    hour.tiles{1}.boundary_applied.bottom.dry_mole_fraction=0;
    hour.tiles{1}.boundary_request.top.requested_external_dry_mole_fraction=0;
    hour.tiles{1}.boundary_request.bottom.requested_external_dry_mole_fraction=0;
    rejected=soil.gas.replayOxygenInterval(empty,hour,respiration,cap);
    assert(~rejected.feasible && ~rejected.has_data && rejected.elapsed_s==0 ...
        && strcmp(rejected.failure_reason,'positive_respiration_without_oxygen') && isequaln(empty,before));
    report=struct('status','PASS','nitrification_sink_mol_cm3_s',nit, ...
        'old_reconstructed_sink_mol_cm3_s',total-respiration*1e-6/12,'rejection',rejected.failure_reason);
end
