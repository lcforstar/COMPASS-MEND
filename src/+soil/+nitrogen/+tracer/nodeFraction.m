function fraction=nodeFraction(tag,total,label)
%NODEFRACTION Bounded passive source share of a nonnegative physical donor.
% Inputs are aligned native nodal inventories or concentrations in the same
% units. A zero donor must carry exactly zero source label. Returns a fraction.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-11-02 Wenjie Liu -- share the original attribution guard with
% the joint passive DIN solve without changing its numerical behavior.
    tag=tag(:);total=total(:);tol=2e-11.*max(1,abs(total));
    assert(all(isfinite([tag;total]))&&all(total>=0) ...
        &&all(tag>=-tol)&&all(tag<=total+tol), ...
        'soil:nitrogen:tracer:nodeFraction', ...
        'Invalid fertilizer fraction in %s.',label);
    fraction=zeros(size(total));positive=total>0;
    fraction(positive)=tag(positive)./total(positive);
    assert(all(tag(~positive)==0), ...
        'soil:nitrogen:tracer:zeroNodeDonor', ...
        'Zero total %s donor contains fertilizer label.',label);
    fraction=min(1,max(0,fraction));
end
