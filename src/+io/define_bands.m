function [spectral] = define_bands()

    % Define spectral regions for SCOPE v_1.40
    % All spectral regions are defined here as row vectors
    % WV Jan. 2013
    % Editor: Jizixing Zhu; Date: 2026-07-03.
    % Change: 2026-07-03 Jizixing Zhu — add spectral.IwlP/IwlT (optical/thermal band index sets into wlS), restored from standard SCOPE define_bands; consumed by RTMo (RTMo.m:424-425) and the static leafopt_s/rs_static assembly.
    % Change: 2026-07-04 Jizixing Zhu -- add spectral.IwlF (fluorescence emission band indices into wlS, ((640:850)-399)), restored from legacy run_stemmus_initialize.m:370; consumed by RTMf + Femtot in surface.radiation.runStepSIF (post-ebal migration PR-A).

    % 3 spectral regions for SCOPE

    reg1 =   400:1:2400;
    reg2 =  2500:100:15000;
    reg3 = 16000:1000:50000;

    spectral.wlS  = [reg1 reg2 reg3];

    % Other spectral (sub)regions

    spectral.wlP   = reg1;                            % PROSPECT data range
    spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
    spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
    spectral.wlO   = reg1;                            % optical part
    spectral.wlT   = [reg2 reg3];                     % thermal part
    spectral.wlZ   = 500:1:600;                       % xanthophyll region
    wlS            = spectral.wlS;
    spectral.wlPAR = wlS(wlS >= 400 & wlS <= 700);  % PAR range

    % Optical/thermal band index sets into wlS (wlS = [wlP wlT], no overlap).
    % RTMo splits hemispherical-out into optical (Eouto) vs thermal (Eoutt) ranges
    % via these; the static leafopt_s/rs_static assembly fills optical vs thermal bands.
    spectral.IwlP  = 1 : numel(spectral.wlP);                      % PROSPECT optical band indices
    spectral.IwlT  = numel(spectral.wlP) + 1 : numel(spectral.wlS); % thermal band indices
    spectral.IwlF  = ((640:850) - 399)';                          % fluorescence emission band indices into wlS (SIF; RTMf/Femtot); legacy run_stemmus_initialize.m:370

    % Data used by aggreg routine to read in MODTRAN data

    spectral.SCOPEspec.nreg = 3;
    spectral.SCOPEspec.start = [400  2500  16000];
    spectral.SCOPEspec.end   = [2400 15000  50000];
    spectral.SCOPEspec.res   = [1   100   1000];

end
