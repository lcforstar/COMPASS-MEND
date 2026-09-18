function [zom, d] = zo_and_d(soil, canopy)
%ZO_AND_D Compute momentum roughness length and zero-plane displacement.
%
%   [zom, d] = equations.zo_and_d(soil, canopy)
%
%   Inputs
%   ------
%   soil.z0m_ground
%       Optional bare-ground momentum roughness length [m]. The default is
%       0.01 m, matching the historical COMPASS CLM5 and optimal branches.
%   canopy.hc
%       Canopy height [m], scalar or N-by-1.
%   canopy.LAI
%       Leaf area index [m2 m-2], scalar or N-by-1.
%   canopy.Method
%       One of clm5, verhoef, optimal, standard, or table. If omitted, the
%       historical COMPASS default is optimal.
%   canopy.SAI
%       Stem area index [m2 m-2]. Required by clm5.
%   canopy.z0mr, canopy.displar
%       CLM5 PFT-specific ratios z0m/h and d/h [-]. Required by clm5; these
%       values must come from the selected PFT rather than a universal crop
%       fallback.
%
%   Outputs
%   -------
%   zom
%       Momentum roughness length [m], N-by-1.
%   d
%       Zero-plane displacement height [m], N-by-1.
%
%   Methods
%   -------
%   clm5
%       Exact CLM5 sparse-to-dense canopy geometry from CLM5 Technical Note
%       equations 5.125-5.127 and CanopyFluxesMod.F90. PAI = LAI + SAI is
%       capped at 2 for the aerodynamic vegetation fraction. At PAI = 0,
%       d = 0 and zom equals the bare-ground roughness.
%   standard
%       Exact FAO-56 height ratios: d = (2/3)h and zom = 0.123h.
%   optimal
%       Original COMPASS composite parameterization by Jizixing Zhu. Its
%       sparse and dense displacement branches meet at LAI = 4. Exact zero
%       LAI is evaluated without a positive numerical surrogate.
%   verhoef
%       Verhoef et al. (1997), retained for SCOPE and V2 production parity.
%   table
%       Historical categorical lookup retained for compatibility.
%
%   A zero canopy height always represents bare ground: d = 0 and zom is
%   soil.z0m_ground (default 0.01 m), independent of the selected method.
%
%   References
%   ----------
%   Lawrence et al. (2018), CLM5 Technical Note, equations 5.125-5.127.
%   Zeng and Wang (2007), J. Hydrometeorology, doi:10.1175/JHM607.1.
%   Allen et al. (1998), FAO Irrigation and Drainage Paper 56.
%   Verhoef et al. (1997), HESS, doi:10.5194/hess-1-81-1997.
%
% Editor: Jizixing Zhu; Date: 2026-04-18.
% Change: 2026-08-07 Jizixing Zhu -- use exact FAO-56 coefficients, replace
% the mislabeled CLM5 equations with the documented CLM5 PAI formulation,
% evaluate exact zero LAI in the COMPASS optimal method, and restore the
% declared bare-ground roughness when canopy height is zero.

    if isfield(canopy, 'Method') && ~isempty(canopy.Method)
        method = canopy.Method;
    elseif isfield(soil, 'Method') && ~isempty(soil.Method)
        method = soil.Method;
    else
        method = 'optimal';
    end

    requireField(canopy, 'hc');
    requireField(canopy, 'LAI');
    h_c = canopy.hc(:);
    lai = canopy.LAI(:);
    n_pts = max(numel(h_c), numel(lai));
    h_c = expandNumericInput(h_c, n_pts, 'canopy.hc');
    lai = expandNumericInput(lai, n_pts, 'canopy.LAI');

    if any(h_c < 0)
        error('zo_and_d:NegativeCanopyHeight', ...
            'canopy.hc must be nonnegative.');
    end
    if any(lai < 0)
        error('zo_and_d:NegativeLAI', ...
            'canopy.LAI must be nonnegative.');
    end

    z0m_ground = 0.01;
    if isfield(soil, 'z0m_ground') && ~isempty(soil.z0m_ground)
        if ~(isnumeric(soil.z0m_ground) && isscalar(soil.z0m_ground) ...
                && isfinite(soil.z0m_ground) && soil.z0m_ground > 0)
            error('zo_and_d:InvalidGroundRoughness', ...
                'soil.z0m_ground must be a finite positive scalar in metres.');
        end
        z0m_ground = double(soil.z0m_ground);
    end

    d = zeros(n_pts, 1);
    zom = zeros(n_pts, 1);

    switch lower(strtrim(method))
        case 'clm5'
            requireField(canopy, 'SAI');
            requireField(canopy, 'z0mr');
            requireField(canopy, 'displar');

            sai = expandNumericInput(canopy.SAI, n_pts, 'canopy.SAI');
            z0mr = expandNumericInput(canopy.z0mr, n_pts, 'canopy.z0mr');
            displar = expandNumericInput(canopy.displar, n_pts, 'canopy.displar');
            if any(sai < 0)
                error('zo_and_d:NegativeSAI', ...
                    'canopy.SAI must be nonnegative.');
            end
            if any(z0mr <= 0)
                error('zo_and_d:InvalidCLM5RoughnessRatio', ...
                    'canopy.z0mr must be positive.');
            end
            if any(displar < 0 | displar >= 1)
                error('zo_and_d:InvalidCLM5DisplacementRatio', ...
                    'canopy.displar must satisfy 0 <= displar < 1.');
            end

            pai_critical = 2.0;
            pai = lai + sai;
            effective_vegetation_fraction = ...
                (1 - exp(-min(pai, pai_critical))) ./ ...
                (1 - exp(-pai_critical));
            dense_canopy_z0m = z0mr .* h_c;
            active = h_c > 0 & pai > 0;
            zom(:) = z0m_ground;
            zom(active) = exp( ...
                effective_vegetation_fraction(active) .* ...
                    log(dense_canopy_z0m(active)) + ...
                (1 - effective_vegetation_fraction(active)) .* ...
                    log(z0m_ground));
            d(active) = effective_vegetation_fraction(active) .* ...
                displar(active) .* h_c(active);

        case 'standard'
            d = (2 / 3) .* h_c;
            zom = 0.123 .* h_c;

        case 'table'
            if isfield(canopy, 'SurfaceType')
                surface_type = canopy.SurfaceType;
            elseif isfield(soil, 'SurfaceType')
                surface_type = soil.SurfaceType;
            else
                surface_type = 'grass';
            end
            z0_value = getTableRoughness(surface_type);
            zom(:) = z0_value;
            d = 0.64 .* h_c;
            if z0_value < 0.01
                d(:) = 0;
            end

        case 'optimal'
            z_0c = zeros(n_pts, 1);
            low_canopy = h_c <= 1;
            medium_canopy = h_c > 1 & h_c < 10;
            high_canopy = h_c >= 10;
            z_0c(low_canopy) = 0.13 .* h_c(low_canopy);
            z_0c(medium_canopy) = ...
                0.139 .* h_c(medium_canopy) ...
                - 0.009 .* h_c(medium_canopy).^2;
            z_0c(high_canopy) = 0.05 .* h_c(high_canopy);

            drag_coefficient = repmat(1.4 / 1000, n_pts, 1);
            vegetated_height = h_c > 0.05;
            exponent_argument = 0.909 - 3.03 .* ...
                (z_0c(vegetated_height) ./ h_c(vegetated_height));
            drag_coefficient(vegetated_height) = ...
                (-1 + exp(exponent_argument)).^4 ./ 4;

            dense = lai >= 4;
            sparse = lai < 4;
            d(dense) = h_c(dense) - z_0c(dense) ./ 0.3;
            logarithm_argument = ...
                1 + (drag_coefficient(sparse) .* lai(sparse)).^0.25;
            d(sparse) = 1.1 .* h_c(sparse) .* log(logarithm_argument);
            d = min(d, 0.99 .* h_c);

            vegetation_roughness = z0m_ground + 0.3 .* h_c .* ...
                (drag_coefficient .* lai).^0.5;
            zom = min(0.3 .* (h_c - d), vegetation_roughness);

        case 'verhoef'
            constants = io.define_constants();
            requireField(canopy, 'CR');
            requireField(canopy, 'CD1');
            requireField(canopy, 'Psicor');
            requireField(soil, 'CSSOIL');

            kappa = constants.kappa;
            cr = expandNumericInput(canopy.CR, n_pts, 'canopy.CR');
            cd1 = expandNumericInput(canopy.CD1, n_pts, 'canopy.CD1');
            psicor = expandNumericInput(canopy.Psicor, n_pts, 'canopy.Psicor');
            cssoil = expandNumericInput(soil.CSSOIL, n_pts, 'soil.CSSOIL');

            square_root_drag = sqrt(cd1 .* lai ./ 2);
            inverse_wind_ratio = max( ...
                3.3, (cssoil + cr .* lai ./ 2).^(-0.5));
            active = lai > 1e-7 & h_c > 1e-7;
            d(active) = h_c(active) .* ...
                (1 - (1 - exp(-square_root_drag(active))) ./ ...
                    square_root_drag(active));
            zom = (h_c - d) .* exp(-kappa .* inverse_wind_ratio + psicor);

        otherwise
            error('zo_and_d:UnknownMethod', ...
                ['Unknown method "%s". Expected clm5, verhoef, optimal, ' ...
                 'standard, or table.'], method);
    end

    zom = max(zom, 1e-4);
    d = max(d, 0);

    bare_ground = h_c == 0;
    zom(bare_ground) = z0m_ground;
    d(bare_ground) = 0;
end

function requireField(value, field_name)
    if ~isfield(value, field_name) || isempty(value.(field_name))
        error('zo_and_d:MissingInput', ...
            'Required input %s is missing or empty.', field_name);
    end
end

function values = expandNumericInput(values, n_pts, name)
    if ~(isnumeric(values) && isreal(values) && ~isempty(values))
        error('zo_and_d:InvalidInput', ...
            '%s must be a nonempty real numeric scalar or vector.', name);
    end
    values = values(:);
    if isscalar(values)
        values = repmat(values, n_pts, 1);
    elseif numel(values) ~= n_pts
        error('zo_and_d:InputSizeMismatch', ...
            '%s must be scalar or contain %d values.', name, n_pts);
    end
    if any(~isfinite(values))
        error('zo_and_d:NonfiniteInput', ...
            '%s must contain only finite values.', name);
    end
end

function z0 = getTableRoughness(surface_type)
    switch lower(strtrim(surface_type))
        case {'grass', 'grass_mowed'}
            z0 = 0.001;
        case {'bare_soil', 'soil'}
            z0 = 0.004;
        case 'thick_grass'
            z0 = 0.09;
        case {'forest', 'forest_level'}
            z0 = 0.95;
        case {'coniferous', 'coniferous_forest'}
            z0 = 1.10;
        case 'alfalfa'
            z0 = 0.03;
        case 'potato'
            z0 = 0.04;
        case 'cotton'
            z0 = 0.13;
        case {'orchard', 'citrus', 'apple'}
            z0 = 0.35;
        otherwise
            z0 = 0.01;
    end
end
