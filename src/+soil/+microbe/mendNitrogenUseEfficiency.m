function phi=mendNitrogenUseEfficiency(activeCarbon,activeNitrogen,cnMin,cnMax,w)
%MENDNITROGENUSEEFFICIENCY Intrinsic NUE from active-microbe C:N.
%   Inputs are top-down layer inventories [mg element g-1 dry soil], C:N
%   bounds [mg C mg-1 N], and dimensionless exponent W. PHI is NL-by-1 in
%   [0,1]. Zero C and N gives the C-limited endpoint; N-free positive C
%   gives the N-limited endpoint. The same MEND closure drives DON-N
%   mineralization and the COMPASS-coupled DIN kinetic request.
%   Source: MOD_MEND.F90 fNLimit_MB, commit 92323c7; Wang et al. (2020),
%   Eq. 38.
%   Editor: Jizixing Zhu; Date: 2026-09-13.
% Change: 2026-10-27 Wenjie Liu -- allow a scalar or layer-resolved exponent.
    activeCarbon=activeCarbon(:); activeNitrogen=activeNitrogen(:);
    n=numel(activeCarbon);
    cnMin=expand(cnMin,n); cnMax=expand(cnMax,n); w=expand(w,n);
    assert(numel(activeNitrogen)==n ...
        &&all(isfinite([activeCarbon;activeNitrogen])) ...
        &&all(activeCarbon>=-1e-12)&&all(activeNitrogen>=-1e-12) ...
        &&all(cnMin>0)&&all(cnMax>cnMin) ...
        &&all(w>=1), ...
        'soil:microbe:mendNueParameters', ...
        'MEND NUE requires nonnegative MBA C/N and valid C:N bounds.');
    carbon=max(activeCarbon,0); nitrogen=max(activeNitrogen,0);
    cn=inf(n,1); hasN=nitrogen>0;
    cn(hasN)=carbon(hasN)./nitrogen(hasN);
    cn(carbon==0&~hasN)=cnMin(carbon==0&~hasN);
    phi=((min(cnMax,max(cnMin,cn))-cnMin)./(cnMax-cnMin)).^w;
end

function value=expand(value,n)
    value=value(:);
    if isscalar(value), value=repmat(value,n,1); end
    assert(numel(value)==n&&all(isfinite(value)), ...
        'soil:microbe:mendNueParameters', ...
        'MEND NUE C:N bounds must be scalar or NL-by-1.');
end
