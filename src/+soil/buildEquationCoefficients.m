function equations = buildEquationCoefficients(built, opts)
%BUILDEQUATIONCOEFFICIENTS  Build coupled water, dry-air, and heat coefficients.
%
%   equations = soil.buildEquationCoefficients(built, opts) reads the current
%   trial state and the already refreshed soil auxiliary fields and returns
%   coefficient banks only. It does not assemble finite-element matrices,
%   construct a right-hand side, apply boundary conditions, enter a linear
%   solve, or mutate built.
%
%   Inputs:
%     built [in] model container with trial state [h: cm, T: degC,
%                p_gas: internal pressure units] and refreshed auxiliary data
%     opts  [in] coefficient options:
%                .DSTMAX [cm] surface-water storage used by dry-air transport
%                .Srt [1/s, NL x nD] root-water sink
%                .L_f [J/kg] fusion latent heat; zero disables latent effects
%                .airHeatImpedance [-] retained dry-air heat impedance factor
%
%   Returns:
%     equations.water   [struct] conservative water-equation coefficients
%     equations.dry_air [struct] retained dry-air coefficients, or empty
%                                fields when lSoilairefc is disabled
%     equations.heat    [struct] retained heat coefficients, or empty fields
%                                when lThmrlefc is disabled
%     equations.flux_family [struct] one shared liquid/vapor/dry-air/energy
%                                element-flux family used by all equation banks
%
%   Method:
%     The water bank uses explicit phase-storage partials. Compatibility names
%     are derived from those physical partials only for the retained dry-air
%     and heat coefficient kernels. Matrix/RHS/boundary assembly belongs to
%     the later equation-system assembly layer. The lSoilairefc and
%     lThmrlefc master switches gate their complete coefficient banks.
%
%   References:
%     soil.water.calculateWaterEquationCoefficients
%     soil.gas.calculateDryAirParameters
%     soil.thermal.calculateEnergyParameters
%   Editor: Jizixing Zhu; Date: 2026-07-12.
%   Change: 2026-07-12 Jizixing Zhu -- gate dry-air and heat coefficient
%           kernels with lSoilairefc and lThmrlefc, respectively.
%   Change: 2026-07-12 Jizixing Zhu -- use the explicit water-equation
%           coefficient function name.
%   Change: 2026-07-13 Jizixing Zhu -- marshal gas pressure mobility through
%           the Kg field without changing coefficient values.
%   Change: 2026-07-14 Jizixing Zhu -- marshal dW_dh and dW_dTheta into the
%           heat coefficient frame for the complete wetting-heat derivative.
%   Change: 2026-08-14 Jizixing Zhu -- create
%           one complete element flux family and pass it to dry-air and heat
%           coefficient consumers without changing the water FEM rule.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = normalizeOptions(built, opts);
    frame = assembleCoefficientFrame(built);
    SV = frame.SoilVariables;
    TC = frame.TransportCoefficient;
    cfg = frame.cfg;

    [water, SV, TC] = soil.water.calculateWaterEquationCoefficients(SV, TC, cfg);
    fluxFamily = soil.coupling.evaluateCoupledElementFlux(SV, TC, cfg);
    dryAir = struct();
    if cfg.lSoilairefc
        dryAir = soil.gas.calculateDryAirParameters( ...
            SV, TC, cfg, opts.DSTMAX, fluxFamily);
    end
    heat = struct();
    if cfg.lThmrlefc
        heat = soil.thermal.calculateEnergyParameters( ...
            SV, TC, water, dryAir, cfg, opts.Srt, opts.L_f, ...
            opts.airHeatImpedance, fluxFamily);
    end

    equations = struct('water', water, 'dry_air', dryAir, 'heat', heat, ...
        'flux_family', fluxFamily);
end

function frame = assembleCoefficientFrame(built)
% Assemble one read-only compatibility frame for the retained coefficient kernels.
    tr = built.workspace.step.trial;
    Y = built.Y;
    d = built.workspace.derived;
    fl = built.workspace.flux;
    S = built.parameters.soil;
    sw = built.control.switches;

    SV = struct();
    SV.hh = tr.h(:); SV.h = Y.h(:);
    SV.TT = tr.T(:); SV.T = Y.T(:);
    SV.P_gg = tr.p_gas(:);
    SV.h_frez = tr.h_frez(:); SV.hh_frez = tr.h_frez(:);
    SV.COR = d.COR(:); SV.CORh = d.CORh(:);
    SV.Theta_LL = d.Theta_LL; SV.Theta_UU = d.Theta_UU; SV.Theta_U = d.Theta_U;
    SV.Theta_V = d.Theta_V; SV.Theta_g = d.Theta_g; SV.Theta_II = d.Theta_II;
    SV.Theta_L = Y.Theta_L; SV.Ratio_ice = d.Ratio_ice;
    SV.DThetaLiquidDh = d.DThetaLiquidDh; SV.DThetaLiquidDT = d.DThetaLiquidDT;
    SV.DThetaWaterEquivDh = d.DThetaWaterEquivDh; SV.DThetaWaterEquivDT = d.DThetaWaterEquivDT;
    SV.DThetaIceDh = d.DThetaIceDh; SV.DThetaIceDT = d.DThetaIceDT;
    SV.DThetaGasDh = d.DThetaGasDh; SV.DThetaGasDT = d.DThetaGasDT;
    SV.DTheta_LLh = d.DThetaLiquidDh; SV.DTheta_UUh = d.DThetaWaterEquivDh;
    SV.SAVEDTheta_LLh = SV.DTheta_LLh; SV.SAVEDTheta_UUh = SV.DTheta_UUh;
    SV.KfL_h = d.KfL_h;
    SV.RHOV = d.RHOV(:); SV.DRHOVh = d.DRHOVh(:); SV.DRHOVT = d.DRHOVT(:);
    SV.RHODA = d.RHODA(:); SV.Xaa = d.Xaa(:); SV.XaT = d.XaT(:); SV.Xah = d.Xah(:);
    SV.DRHODAz = d.DRHODAz(:); SV.L = d.L(:);
    SV.POR = S.porosity(:);
    if readSwitch(sw, 'lKLa_Switch')
        SV.KLa_Switch = 1;
    else
        SV.KLa_Switch = [];
    end

    TC = struct('D_V', d.D_V, 'Eta', d.Eta, 'KL_T', d.KL_T, ...
        'D_Ta', d.D_Ta, 'Kg', d.Kg, 'D_Vg', d.D_Vg, ...
        'V_A', fl.V_A, 'Lambda_eff', d.Lambda_eff, ...
        'c_unsat', d.c_unsat, 'WW', d.WW, ...
        'dW_dh', d.dW_dh, 'dW_dTheta', d.dW_dTheta, ...
        'DTheta_LLT', d.DThetaLiquidDT);

    cfg = struct();
    activeBottomNode = 1;
    if isfield(built.workspace, 'boundary') ...
            && isfield(built.workspace.boundary, 'groundwater_driver')
        driver = built.workspace.boundary.groundwater_driver;
        if isfield(driver, 'GroundwaterCoupling') && driver.GroundwaterCoupling == 1
            assert(isfield(driver, 'indexBotmLayer') ...
                    && isscalar(driver.indexBotmLayer) ...
                    && driver.indexBotmLayer >= 1 ...
                    && driver.indexBotmLayer <= built.domain.NN, ...
                'soil:buildEquationCoefficients:badGroundwaterNode', ...
                'Enabled groundwater driver must provide a valid active node.');
            activeBottomNode = driver.indexBotmLayer;
        end
    end
    cfg.Grid = struct('NL', built.domain.NL, 'NN', built.domain.NN, ...
        'nD', built.domain.grid.nD, 'DeltZ', built.domain.dz, ...
        'Tot_Depth', built.domain.grid.Tot_Depth, ...
        'indxBotmLayer', activeBottomNode);
    for name = {'lhThmrl','lThmrlefc','lKLa_Switch','lDVa_Switch', ...
            'lSoilairefc','lSWCC','lFreezing','lSFCC'}
        cfg.(name{1}) = readSwitch(sw, name{1});
    end
    frame = struct('SoilVariables', SV, 'TransportCoefficient', TC, 'cfg', cfg);
end

function opts = normalizeOptions(built, opts)
% Fill coefficient-only defaults without deriving new physical values.
    defaults = struct('DSTMAX', 0, 'Srt', built.workspace.flux.Srt, ...
        'L_f', built.workspace.step.L_f, 'airHeatImpedance', 0);
    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end
end

function value = readSwitch(switches, name)
% Read one legacy integer switch as a scalar logical value.
    value = isfield(switches, name) && switches.(name) == 1;
end
