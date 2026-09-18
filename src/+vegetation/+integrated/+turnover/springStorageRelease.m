function [N_storage_new, N_labile_new, dN_release] = springStorageRelease( ...
    N_storage, N_labile, frPHU, spring_release_endphu, spring_release_days, dt)
%SPRINGSTORAGERELEASE  Daily N_storage release to N_labile during spring window.
%
%   Spring window: 0 < frPHU < spring_release_endphu (typically 0.20).
%   Release rate: 1/spring_release_days (constant, e.g. 60 days → 1.67%/day).
%
%   Returns:
%     N_storage_new
%     N_labile_new
%     dN_release   [g N/m² ground]
%
%   Inputs:
%     N_storage, N_labile     [g N/m² ground]
%     frPHU                   [-]   spring iff > 0 && < spring_release_endphu
%     spring_release_endphu   [-]   Apple 0.20 (frPHU window upper bound)
%     spring_release_days     [-]   Apple 60 days (decay constant for N→labile)
%     dt                      [s]
%
%   Method:
%     d_N_release = N_storage / spring_release_days × dt/86400
%     N_storage -= d_release
%     N_labile  += d_release
%
%   Physical: Spring burst (Apple buds break around DOY 104; storage that was
%   retransloc'd last autumn now feeds new leaf + root growth).
%
%   Reference:
%     Tagliavini M et al (2003) — Apple stem N as winter storage; spring remobilization
%     Aerts 1996 — perennial N cycling
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.4 S2.4)

    if frPHU <= 0 || frPHU >= spring_release_endphu || N_storage <= 0
        N_storage_new = N_storage;
        N_labile_new  = N_labile;
        dN_release    = 0;
        return;
    end

    rate_per_day  = 1.0 / max(spring_release_days, 1);   % [1/day]
    dN_release    = N_storage * rate_per_day * (dt / 86400);
    dN_release    = min(dN_release, N_storage);          % cannot release more than available
    N_storage_new = N_storage - dN_release;
    N_labile_new  = N_labile  + dN_release;
end
