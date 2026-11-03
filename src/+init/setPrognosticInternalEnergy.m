function built = setPrognosticInternalEnergy(built, energyStorage)
%SETPROGNOSTICINTERNALENERGY  Write Y.U, trial.U, previous.U from endpoint U.
%
%   built = init.setPrognosticInternalEnergy(built, energyStorage) computes
%   nodal internal energy from the endpoint energyStorage (NL×2) and writes
%   it to the three state mirrors: Y.U, workspace.step.trial.U, and
%   workspace.step.previous.U.  Regional [NN × nCell] shapes are handled by
%   applying the computation independently to each column cell.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- validate endpoint/cell dimensions and
%           reject ambiguous regional replication.

    NL = built.domain.NL;
    dz = built.domain.dz;
    nCell = size(built.Y.T, 2);  % 1 for single-column, >1 for regional

    if size(energyStorage, 1) ~= NL || size(energyStorage, 2) ~= 2 ...
            || any(~isfinite(energyStorage(:)))
        error('setPrognosticInternalEnergy:invalidEndpointEnergy', ...
            'energyStorage must be finite with leading shape [NL x 2].');
    end
    if nCell > 1 && (ndims(energyStorage) ~= 3 ...
            || size(energyStorage, 3) ~= nCell)
        error('setPrognosticInternalEnergy:ambiguousRegionalEnergy', ...
            'Regional energyStorage must have shape [NL x 2 x nCell].');
    end

    U_all = zeros(built.domain.NN, nCell);
    for c = 1:nCell
        ep = energyStorage;
        if nCell > 1 && ndims(energyStorage) == 3
            ep = energyStorage(:, :, c);
        end
        U_nodal = soil.thermal.assembleNodalInternalEnergy( ...
            ep, zeros(NL,2), zeros(NL,2), zeros(NL,2), dz);
        U_all(:, c) = U_nodal;
    end

    built.Y.U = U_all;
    if isfield(built.workspace.step, 'trial') ...
            && isfield(built.workspace.step.trial, 'U')
        built.workspace.step.trial.U = U_all;
    end
    if isfield(built.workspace.step, 'previous') ...
            && isfield(built.workspace.step.previous, 'U')
        built.workspace.step.previous.U = U_all;
    end
end
