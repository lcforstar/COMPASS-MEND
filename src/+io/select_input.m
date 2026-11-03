function [soil, leafbio, canopy, meteo, angles, SpaceTimeInfo] = select_input(ScopeParameters, Theta_LL, digitsVector, canopy, options, SoilProperties)
    % 从ScopeParameters中选择输入数据
    % ScopeParametersNames 索引（与 loadParameters.m 保持一致）：
    %  1-Cab  2-Cca  3-Cdm  4-Cw  5-Cs  6-N  7-rho_thermal  8-tau_thermal
    %  9-Vcmo  10-m  11-Type  12-kV  13-Rdparam  14-Tparam  15-fqe  16-spectrum
    %  17-rss  18-rs_thermal  19-cs  20-rhos  21-lambdas  22-LAI  23-hc
    %  24-zo  25-d  26-LIDFa  27-LIDFb  28-leafwidth  29-z  30-Rin  31-Ta  32-Ts
    %  33-Rli  34-p  35-ea  36-u  37-Ca  38-Oa  39-rb  40-Cd  41-CR  42-CD1
    %  43-Psicor  44-CSSOIL  45-rbs  46-rwc  47-startDOY  48-endDOY
    %  49-LAT  50-LON  51-timezn  52-tts  53-tto  54-psi  55-SMC
    %  56-Tyear  57-beta  58-kNPQs  59-qLs  60-stressfactor  61-Cant
    %  62-BSMBrightness  63-BSMlat  64-BSMlon  65-BallBerry0  66-SAI
    theta_s0 = max(SoilProperties.Theta_s_type);

    soil.spectrum      = ScopeParameters.spectrum(digitsVector(16));
    soil.rss           = ScopeParameters.rss(digitsVector(17));
    soil.rs_thermal    = ScopeParameters.rs_thermal(digitsVector(18));
    soil.cs            = ScopeParameters.cs(digitsVector(19));
    soil.rhos          = ScopeParameters.rhos(digitsVector(20));
    soil.CSSOIL        = ScopeParameters.CSSOIL(digitsVector(44));
    soil.lambdas       = ScopeParameters.lambdas(digitsVector(21));
    soil.rbs           = ScopeParameters.rbs(digitsVector(45));
    soil.SMC           = Theta_LL(end - 1, 1);
    soil.BSMBrightness = ScopeParameters.BSMBrightness(digitsVector(62));
    soil.BSMlat        = ScopeParameters.BSMlat(digitsVector(63));
    soil.BSMlon        = ScopeParameters.BSMlon(digitsVector(64));

    leafbio.Cab     = ScopeParameters.Cab(digitsVector(1));
    leafbio.Cca     = ScopeParameters.Cca(digitsVector(2));
    if options.Cca_function_of_Cab
        leafbio.Cca = 0.25 * ScopeParameters.Cab(digitsVector(1));
    end
    leafbio.Cdm     = ScopeParameters.Cdm(digitsVector(3));
    leafbio.Cw      = ScopeParameters.Cw(digitsVector(4));
    leafbio.Cs      = ScopeParameters.Cs(digitsVector(5));
    leafbio.Cant    = ScopeParameters.Cant(digitsVector(61));
    leafbio.N       = ScopeParameters.N(digitsVector(6));
    leafbio.Vcmo    = ScopeParameters.Vcmo(digitsVector(9));
    leafbio.m       = ScopeParameters.m(digitsVector(10));
    leafbio.BallBerry0 = ScopeParameters.BallBerry0(digitsVector(65)); % JAK 2016-10. Accidentally left out of v1.70
    leafbio.Type    = ScopeParameters.Type(digitsVector(11));
    leafbio.Tparam  = ScopeParameters.Tparam(digitsVector(14), :); % this is correct (: instead of 14)
    fqe             = ScopeParameters.fqe(digitsVector(15));
    leafbio.Rdparam = ScopeParameters.Rdparam(digitsVector(13));

    leafbio.rho_thermal = ScopeParameters.rho_thermal(digitsVector(7));
    leafbio.tau_thermal = ScopeParameters.tau_thermal(digitsVector(8));

    leafbio.Tyear         = ScopeParameters.Tyear(digitsVector(56));
    leafbio.beta          = ScopeParameters.beta(digitsVector(57));
    leafbio.kNPQs         = ScopeParameters.kNPQs(digitsVector(58));
    leafbio.qLs           = ScopeParameters.qLs(digitsVector(59));
    leafbio.stressfactor  = ScopeParameters.stressfactor(digitsVector(60));

    canopy.LAI  = ScopeParameters.LAI(digitsVector(22));
    canopy.hc  = ScopeParameters.hc(digitsVector(23));
    if isfield(ScopeParameters, 'SAI')
        if numel(digitsVector) >= 66
            canopy.SAI = ScopeParameters.SAI(min(digitsVector(66), numel(ScopeParameters.SAI)));
        else
            canopy.SAI = ScopeParameters.SAI(min(digitsVector(22), numel(ScopeParameters.SAI)));
        end
    else
        canopy.SAI = 0.05;
    end
    canopy.LIDFa = ScopeParameters.LIDFa(digitsVector(26));
    canopy.LIDFb  = ScopeParameters.LIDFb(digitsVector(26)); % this is correct (26 instead of 27)
    canopy.leafwidth  = ScopeParameters.leafwidth(digitsVector(28));
    canopy.rb   = ScopeParameters.rb(digitsVector(39));
    canopy.Cd  = ScopeParameters.Cd(digitsVector(40));
    canopy.CR = ScopeParameters.CR(digitsVector(41));
    canopy.CD1  = ScopeParameters.CD1(digitsVector(42));
    canopy.Psicor  = ScopeParameters.Psicor(digitsVector(43));
    canopy.rwc  = ScopeParameters.rwc(digitsVector(46));
    canopy.kV = ScopeParameters.kV(digitsVector(12));
    canopy.zo  = ScopeParameters.zo(digitsVector(24));
    canopy.d = ScopeParameters.d(digitsVector(25));

    meteo.z  = ScopeParameters.z(digitsVector(29));
    meteo.Rin   = ScopeParameters.Rin(digitsVector(30));
    meteo.Ta = ScopeParameters.Ta(digitsVector(31));
    meteo.Ts = ScopeParameters.Ts(digitsVector(32));  % 实测地表温度 (°C)
    meteo.Rli  = ScopeParameters.Rli(digitsVector(33));
    meteo.p  = ScopeParameters.p(digitsVector(34));
    meteo.ea  = ScopeParameters.ea(digitsVector(35));
    meteo.u   = ScopeParameters.u(digitsVector(36));
    meteo.Ca = ScopeParameters.Ca(digitsVector(37));
    meteo.Oa  = ScopeParameters.Oa(digitsVector(38));

    SpaceTimeInfo.startDOY = ScopeParameters.startDOY(digitsVector(47));
    SpaceTimeInfo.endDOY = ScopeParameters.endDOY(digitsVector(48));
    SpaceTimeInfo.LAT = ScopeParameters.LAT(digitsVector(49));
    SpaceTimeInfo.LON = ScopeParameters.LON(digitsVector(50));
    SpaceTimeInfo.timezn = ScopeParameters.timezn(digitsVector(51));

    angles.tts = ScopeParameters.tts(digitsVector(52));
    angles.tto = ScopeParameters.tto(digitsVector(53));
    angles.psi = ScopeParameters.psi(digitsVector(54));


    %% derived input
    if options.soil_heat_method == 1
        soil.GAM =  equations.Soil_Inertia1(soil.SMC, theta_s0);
    else
        soil.GAM  = equations.Soil_Inertia0(soil.cs, soil.rhos, soil.lambdas);
    end
    if options.calc_rss_rbs
        SoilProperties.SMC = soil.SMC;  % 传入当前表层含水量
        [soil.rss, soil.rbs] = equations.calc_rssrbs(SoilProperties, canopy.LAI, soil.rbs);
    end

    if leafbio.Type
        leafbio.Type = 'C4';
    else
        leafbio.Type = 'C3';
    end
    canopy.hot  = canopy.leafwidth / canopy.hc;
    if options.calc_zo
        [canopy.zo, canopy.d]  = equations.zo_and_d(soil, canopy);
    end

    if options.calc_PSI == 1
        leafbio.fqe(1) = fqe / 5;
        leafbio.fqe(2) = fqe;
    else
        leafbio.fqe = fqe;
    end
