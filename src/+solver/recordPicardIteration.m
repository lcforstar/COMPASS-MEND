function diag = recordPicardIteration(diag, kit, waterError, energyError, soilVariables, bcSwitch, stage)
%RECORDPICARDITERATION Store one Picard iteration diagnostic sample.

diag.water_error(kit) = waterError;
diag.energy_error(kit) = energyError;
diag.hh_surface(kit) = soilVariables.hh(end);
diag.TT_surface(kit) = soilVariables.TT(end);
diag.bc_switch(kit) = logical(bcSwitch);
diag.stage{kit} = stage;
end
