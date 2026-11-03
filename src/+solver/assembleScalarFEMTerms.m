function terms = assembleScalarFEMTerms(storage, conductance, advection, gradientSource, distributedSource, dz)
%ASSEMBLESCALARFEMTERMS  Assemble one scalar two-node finite-element operator.
%
%   terms = solver.assembleScalarFEMTerms(storage, conductance, advection,
%       gradientSource, distributedSource, dz) applies the common endpoint
%   quadrature used by the COMPASS water, heat, and dry-air equations.
%
%   Inputs:
%     storage          [coefficient, NL x 2] endpoint storage coefficient
%     conductance      [coefficient, NL x 2] endpoint diffusion/conduction coefficient
%     advection        [coefficient, NL x 2] endpoint first-order transport coefficient
%     gradientSource   [source, NL x 2] endpoint gradient/gravity source
%     distributedSource[source, NL x 2] endpoint volumetric source
%     dz               [cm, NL x 1] element thickness, bottom-up
%
%   Returns:
%     terms.mass       [assembled coefficient, NN x 2] lumped storage matrix
%     terms.operator   [assembled coefficient, NN x 2] diagonal/upper operator
%     terms.lower      [assembled coefficient, NL x 1] lower operator
%     terms.gradient   [assembled source, NN x 1] signed endpoint-gradient source
%     terms.load       [assembled source, NN x 1] consistent distributed load
%
%   Method:
%     This is the shared numerical skeleton only. Physical coefficient
%     formulas, equation RHS combinations, boundary policies, and state
%     updates remain in their water/heat/dry-air modules. The 1/2, 1/3,
%     and 1/6 weights reproduce the existing COMPASS scalar FEM assembly.
%
%   Reference:
%     soil.water.assembleWaterEquationMatrices
%     soil.thermal.assembleThermalEquationMatrices
%     soil.gas.assembleDryAirEquationMatrices
%   Editor: Jizixing Zhu; Date: 2026-07-13.

    dz = dz(:);
    NL = numel(dz);
    NN = NL + 1;
    storage = requireEndpoints(storage, NL, 'storage');
    conductance = requireEndpoints(conductance, NL, 'conductance');
    advection = requireEndpoints(advection, NL, 'advection');
    gradientSource = requireEndpoints(gradientSource, NL, 'gradientSource');
    distributedSource = requireEndpoints(distributedSource, NL, 'distributedSource');

    mass = zeros(NN, 2);
    operator = zeros(NN, 2);
    lower = zeros(NL, 1);
    gradient = zeros(NN, 1);
    load = zeros(NN, 1);

    for element = 1:NL
        width = dz(element);
        mass(element, 1) = mass(element, 1) + storage(element, 1) * width / 2;
        mass(element + 1, 1) = mass(element + 1, 1) + storage(element, 2) * width / 2;

        diffusion = (conductance(element, 1) + conductance(element, 2)) / (2 * width);
        advLeft = advection(element, 1) / 3 + advection(element, 2) / 6;
        advRight = advection(element, 1) / 6 + advection(element, 2) / 3;
        operator(element, 1) = operator(element, 1) + diffusion - advLeft;
        operator(element, 2) = operator(element, 2) - diffusion - advRight;
        operator(element + 1, 1) = operator(element + 1, 1) + diffusion + advRight;
        lower(element) = -diffusion + advLeft;

        gradientMean = (gradientSource(element, 1) + gradientSource(element, 2)) / 2;
        gradient(element) = gradient(element) - gradientMean;
        gradient(element + 1) = gradient(element + 1) + gradientMean;

        loadLeft = (2 * distributedSource(element, 1) + distributedSource(element, 2)) * width / 6;
        loadRight = (distributedSource(element, 1) + 2 * distributedSource(element, 2)) * width / 6;
        load(element) = load(element) + loadLeft;
        load(element + 1) = load(element + 1) + loadRight;
    end

    terms = struct('mass', mass, 'operator', operator, 'lower', lower, ...
        'gradient', gradient, 'load', load);
end

function value = requireEndpoints(value, NL, name)
% Normalize empty optional terms to zero and reject ambiguous shapes.
    if isempty(value)
        value = zeros(NL, 2);
    end
    assert(isequal(size(value), [NL, 2]) && all(isfinite(value(:))), ...
        'assembleScalarFEMTerms:badShape', ...
        '%s must be a finite NL-by-2 endpoint array.', name);
end
