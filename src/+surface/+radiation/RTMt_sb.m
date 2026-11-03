function [rad] = RTMt_sb(spectral, rad, soil, leafopt, canopy, gap, angles, Tcu, Tch, Tsu, Tsh, obsdir)
    % RTMo（光学）：基础散射/吸收模型。从这开始——双流近似或离散纵标法（取决于具体实现）。
    % RTMf（荧光）：在长波段添加再发射项。SIF反演必需，但计算成本低，因为叶绿素荧光仅占吸收PAR的~1%。
    % 热辐射模型——关键权衡：
    % 
    % RTMt_sb：Stefan-Boltzmann近似（T⁴标度）。快速但假设灰体行为——对通量计算可接受
    % RTMt_planck：全光谱积分。物理上更精确但：
    % 
    % 慢~10-100倍（取决于光谱分辨率）
    % 需要光谱发射率ε(λ,T)——实际很难获取足够精度的数据
    % 多数应用没有验证数据支撑这种复杂度
    % 
    % 
    % 
    % 实际问题：除非做高光谱热红外反演或耦合大气辐射传输，Stefan-Boltzmann版本够用。叶温和发射率参数化的不确定性通常远大于灰体近似误差。
    % RTMz批判：叶黄素循环（玉米黄质←→紫黄质转化）是光保护机制，但：
    % 
    % 时间尺度不匹配：V→Z转化需数分钟；标准RT假设瞬时稳态
    % 与NPQ耦合：没有独立测量无法分离动态玉米黄质效应和其他非光化学淬灭机制
    % 除非有PRI（光化学反射指数）数据或受控实验，否则增加参数但不增加约束能力
    % function 'RTMt_sb' calculates total outgoing radiation in hemispherical
    % direction and total absorbed radiation per leaf and soil component.
    % Radiation is integrated over the whole thermal spectrum with
    % Stefan-Boltzman's equation. This function is a simplified version of
    % 'RTMt_planck', and is less time consuming since it does not do the
    % calculation for each wavelength separately.
    %
    % Authors: Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
    % date:     5  Nov 2007
    % update:   13 Nov 2007
    %           16 Nov 2007 CvdT    improved calculation of net radiation
    %           27 Mar 2008 JT      added directional calculation of radiation
    %           24 Apr 2008 JT      Introduced dx as thickness of layer (see parameters)
    %           31 Oct 2008 JT      introduced optional directional calculation
    %           31 Oct 2008 JT      changed initialisation of F1 and F2 -> zeros
    %           07 Nov 2008 CvdT    changed layout
    %           16 Mar 2009 CvdT    removed Tbright calculation
    %              Feb 2013 WV      introduces structures for version 1.40
    %           19 Apr 2026 (Jizixing Zhu / 朱寄子星)
    %                               non-uniform layer thickness: replaced scalar dx/iLAI
    %                               with canopy.layer_weight-driven dx_layer/LAI_layer;
    %                               updated two-stream iteration and directional terms
    %
    % Table of contents of the function
    %   0       preparations
    %       0.0     globals
    %       0.1     initialisations
    %       0.2     parameters
    %       0.3     geometric factors of Observer
    %       0.4     geometric factors associated with extinction and scattering
    %       0.5     geometric factors to be used later with rho and tau
    %       0.6     fo for all leaf angle/azumith classes
    %   1       calculation of upward and downward fluxes
    %   2       total net fluxes
    %   Appendix A. Stefan-Boltzmann
    %
    % usage:
    % [rad] = RTMt_sb(options,spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Tsu,Tsh)
    %
    % Most input and output are structures. These structures are further
    % specified in a readme file. The temperatures Tcu, Tch, Tsu and Tsh are
    % variables.
    %
    % Input:
    %   options       calculation options
    %   spectral    information about wavelengths and resolutions
    %   rad         a large number of radiative fluxes: spectrally distributed
    %               and integrated, and canopy radiative transfer coefficients
    %   soil        soil properties
    %   leafopt     leaf optical properties
    %   canopy      canopy properties (such as LAI and height)
    %   gap         probabilities of direct light penetration and viewing
    %   angles      viewing and observation angles
    %   Tcu         Temperature of sunlit leaves    (oC), [13x36xnl]
    %   Tch         Temperature of shaded leaves    (oC), [13x36xnl]
    %   Tsu         Temperature of sunlit soil      (oC), [1]
    %   Tsh         Temperature of shaded soil      (oC), [1]
    %
    % Output:
    %   rad         a large number of radiative fluxes: spectrally distributed
    %               and integrated, and canopy radiative transfer coefficients.
    %               Here, thermal fluxes are added

    % Editor: Jizixing Zhu; Date: 2026-04-19.
    % Change: 2026-07-16 Jizixing Zhu — use Omega*LAI for thermal radiative transfer and return leaf fluxes per physical leaf area.

    %% 0.1 parameters

    IT          = find(spectral.wlS == 10000);   % Take 10 microns as representative wavelength for the thermal

    % load Constants
    Constants = io.define_constants();

    deg2rad     = Constants.deg2rad;
    nl          = canopy.nlayers;
    lidf        = canopy.lidf;
    litab       = canopy.litab;
    lazitab     = canopy.lazitab;
    nlazi       = length(lazitab);
    tto         = angles.tto;
    psi         = angles.psi;
    Ps          = gap.Ps;
    K           = gap.K;

    rho         = leafopt.refl(IT);       % [nwl]               Leaf/needle reflection
    tau         = leafopt.tran(IT);       % [nwl]               Leaf/needle transmission
    rs          = soil.refl(IT);          % [nwl]               Soil reflectance
    epsc        = 1 - rho - tau;              % [nwl]               Emissivity vegetation
    epss        = 1 - rs;                   % [nwl]               Emissivity soil
    crit        = max(1E-2);              % [1]                 Desired minimum accuracy
    LAI_physical = canopy.LAI;
    if isfield(canopy, 'clumping_index') && ~isempty(canopy.clumping_index)
        clumpingIndex = canopy.clumping_index;
    else
        clumpingIndex = 1.0;
    end
    if ~isscalar(clumpingIndex) || ~isfinite(clumpingIndex) || ...
            clumpingIndex <= 0 || clumpingIndex > 1
        error('RTMt_sb:badClumpingIndex', ...
            'canopy.clumping_index must be a finite scalar in (0,1].');
    end
    LAI_optical = clumpingIndex * LAI_physical;
    % [改进 2026-04-19 朱寄子星] 非均匀层厚：有 canopy.layer_weight 时按实际层厚加权，
    % 否则退回均匀 1/nl。dx_layer [nl×1] 归一化，LAI_layer = LAI * dx_layer。
    % Non-uniform layer thickness: use canopy.layer_weight if present, else 1/nl fallback.
    if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
        dx_layer = canopy.layer_weight(:);
    else
        dx_layer = ones(nl, 1) / nl;
    end
    LAI_layer   = LAI_optical * dx_layer;

    %% 0.2 initialiations
    Rnhc        = zeros(nl, 1);  % [nl]
    Rnuc        = zeros(size(Tcu));   % [13,36,nl]

    %% 0.3 geometric factors of observer
    if obsdir
        cos_tto = cos(tto * deg2rad);   % [1]                 cos observation angle
        sin_tto = sin(tto * deg2rad);   % [1]                 sin observation angle
    end

    %% 0.4 geometric factors associated with extinction and scattering
    cos_ttl     = cos(litab * deg2rad);
    if obsdir
        sin_ttl = sin(litab * deg2rad);
        cos_ttlo = cos((lazitab - psi) * deg2rad);
    end
    bfli        = cos_ttl.^2;
    bf          = bfli' * lidf;

    %% 0.5 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
    ddb         = 0.5 * (1 + bf);         %                     f1^2 + f2^2
    ddf         = 0.5 * (1 - bf);         %                     2*f1*f2
    if obsdir
        dob     = 0.5 * (K + bf);         %                     fo*f1
        dof     = 0.5 * (K - bf);         %                     fo*f1
    end

    %% 0.6 fo for all leaf angle/azumith classes
    if obsdir
        Co      = cos_ttl * cos_tto;    % [nli]               pag 305
        So      = sin_ttl * sin_tto;    % [nli]               pag 305
        cos_deltao  = Co * ones(1, nlazi) + So * cos_ttlo'; % [nli, nlazi] projection of leaves in  in direction of sun (pag 125-126)
        fo      = cos_deltao / abs(cos_tto); % [nli, nlazi]    leaf area projection factors in direction of observation
    end

    %% 1. calculation of upward and downward fluxes pag 305
    sigb        = ddb * rho + ddf * tau;  % [nwlt]              Diffuse backscatter scattering coefficient
    sigf        = ddf * rho + ddb * tau;  % [nwlt]              Diffuse forward     scattering coefficient
    if obsdir
        vb      = dob * rho + dof * tau;  % [nwlt]              Directional backscatter scattering coefficient for diffuse  incidence
        vf      = dof * rho + dob * tau;  % [nwlt]              Directional forward     scattering coefficient for diffuse  incidence
    end
    a           = 1 - sigf;             % [nwlt]              Attenuation
    m           = sqrt(a * a - sigb * sigb); % [nwlt]
    rinf        = (a - m) / sigb;         % [nwlt]              Reflection coefficient for infinite thick canopy
    rinf2       = rinf * rinf;          % [nwlt]

    fHs         = (1 - rinf2) * (1 - rs) / (1 - rinf * rs);
    % [改进 2026-04-19 朱寄子星] fHc 改为 [nl×1] 向量：每层用 LAI_layer(j) 替代标量 iLAI；
    % 两流迭代步长同步改为 LAI_layer(j)，保持逐层能量守恒。
    % fHc is now [nl×1]; two-stream step uses per-layer LAI_layer(j).
    fHc_layer   = LAI_layer * m * (1 - rinf);
    fbottom     = (rs - rinf) / (1 - rs * rinf);

    % 1.1 radiance by components
    C2K = Constants.C2K;
    sigmaSB = Constants.sigmaSB;
    Hcsu3       = Stefan_Boltzmann(Tcu, C2K, sigmaSB); %                   Radiance by sunlit leaves
    Hcsh        = Stefan_Boltzmann(Tch, C2K, sigmaSB); %                   Radiance by shaded leaves
    Hssu        = Stefan_Boltzmann(Tsu, C2K, sigmaSB); %                   Radiance by sunlit soil
    Hssh        = Stefan_Boltzmann(Tsh, C2K, sigmaSB); %                   Radiance by shaded soil

    % 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)
    v1 = repmat(1 / size(Hcsu3, 2), 1, size(Hcsu3, 2)); % vector for computing the mean
    Hcsu2 = reshape(Hcsu3, size(Hcsu3, 1), []);   % create a block matrix from the 3D array
    Hcsu = (v1 * reshape(Hcsu2' * lidf, size(Hcsu3, 2), []))'; % compute column means for each level

    Hc          = Hcsu .* Ps(1:nl) + Hcsh .* (1 - Ps(1:nl));      % hemispherical emittance by leaf layers
    Hs          = Hssu .* Ps(nl + 1) + Hssh .* (1 - Ps(nl + 1));      % hemispherical emittance by soil surface

    % 1.3 Diffuse radiation
    cont        = 1;                                        % continue iteration (1:yes, 0:no)
    counter     = 0;                                        % number of iterations
    F1          = zeros(nl + 1, 1);
    F2          = zeros(nl + 1, 1);
    F1top       = 0;
    while cont
        F1topn  = -rinf * F2(1);
        F1(1)   = F1topn;
        for j   = 1:nl
            F1(j + 1) = F1(j) * (1 - m * LAI_layer(j)) + fHc_layer(j) * Hc(j);
        end
        F2(nl + 1) = fbottom * F1(nl + 1) + fHs * Hs;
        for j   = nl:-1:1
            F2(j)   = F2(j + 1) * (1 - m * LAI_layer(j)) + fHc_layer(j) * Hc(j);
        end
        cont    = abs(F1topn - F1top) > crit;
        F1top   = F1topn;
        counter = counter + 1;
    end

    Emin        = (F1 + rinf * F2) / (1 - rinf2);
    Eplu        = (F2 + rinf * F1) / (1 - rinf2);

    % 1.4 Directional radiation
    if obsdir
        % [改进 2026-04-19 朱寄子星] 方向辐射贡献改用 LAI_layer 逐层加权替代标量 iLAI。
        % Directional thermal terms weighted by LAI_layer instead of scalar iLAI.
        piLo1       = epsc * K * (LAI_layer .* Hcsh)' * (gap.Po(1:nl) - gap.Pso(1:nl));                % directional   emitted     radation by shaded leaves
        % JAK 2015-01: replaced earlier loop by this: all-at-once with more efficient mean
        absfo_rep = repmat(abs(fo), 1, nl);
        piLo2 = epsc * (LAI_layer .* (v1 * reshape((Hcsu2 .* absfo_rep)' * lidf, size(Hcsu3, 2), []))') .* gap.Pso(1:nl); % compute column means for each level

        piLo3       = (vb * Emin(1:nl) + vf * Eplu(1:nl))' * (LAI_layer .* gap.Po(1:nl));      % directional   scattered   radiation by vegetation for diffuse incidence
        piLo4       = epss * Hssh * (gap.Po(nl + 1) - gap.Pso(nl + 1));                        % directional   emitted     radiation by shaded soil
        piLo5       = epss * Hssu * gap.Pso(nl + 1);                                   % directional   emitted     radiation by sunlit soil
        piLo6       = rs * Emin(nl + 1) * gap.Po(nl + 1);                                % directional   scattered   radiation by soil       for diffuse incidence       [1]

        piLot       = piLo1 + sum(piLo2) + piLo3 + piLo4 + piLo5 + piLo6;
    else
        piLot       = NaN;
    end
    Lot             = piLot / pi;

    %% 2. total net fluxes
    % net radiation per component, in W m-2 (leaf or soil surface)
    for j = 1:nl
        % Convert the RT result per effective leaf area to an average per
        % physical leaf area before ebal integrates over physical LAI.
        Rnuc(:, :, j) = clumpingIndex * ...
            (Emin(j) + Eplu(j + 1) - 2 * Hcsu3(:, :, j)) * epsc;    % sunlit leaf
        Rnhc(j) = clumpingIndex * ...
            (Emin(j) + Eplu(j + 1) - 2 * Hcsh(j)) * epsc;           % shaded leaf
    end
    Rnus            = (Emin(nl + 1) - Hssu) * epss;                       % sunlit soil
    Rnhs            = (Emin(nl + 1) - Hssh) * epss;                       % shaded soil

    %% 3. Write the output to the rad structure

    rad.Emint   = Emin;
    rad.Eplut   = Eplu;
    rad.Eoutte  = Eplu(1) - Emin(1); % 1)
    rad.Lot     = Lot;
    rad.Rnuct   = Rnuc;
    rad.Rnhct   = Rnhc;
    rad.Rnust   = Rnus;
    rad.Rnhst   = Rnhs;
    return

    % 1) CvdT, 11 December 2015.
    % We subtract Emin(1), because ALL incident (thermal) radiation from Modtran
    % has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!

    %% Appendix A. Stefan-Boltzmann
function H      =   Stefan_Boltzmann(T_C, C2K, sigmaSB)

    H       = sigmaSB * (T_C + C2K).^4;
    return
