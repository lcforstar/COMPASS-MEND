function a = fruitTreeAllocation(gdd_budburst, harvested, p)
%FRUITTREEALLOCATION  5-phase GDD carbon-allocation fractions for a perennial fruit tree.
%
%   Faithful port of the CLM5-FruitTree 5-phase perennial-woody allocation
%   (Dombrowski et al. 2022 GMD 15:5167), keyed to thermal time since bud-break
%   g = GDD-since-budburst. Returns per-organ NPP fractions (leaf/root/stem/repr,
%   sum = 1). Replaces the earlier 2-phase-linear simplification, which merged the
%   huiripe (grnrp) split and kept a_root elevated through fruit-fill (-> root too
%   high, fruit too low vs FruitTree-sim).
%
%   Returns:
%     a   struct, fields each [-] (fraction of NPP), a.leaf+a.root+a.stem+a.repr = 1
%       .leaf .root .stem .repr
%
%   Inputs:
%     gdd_budburst [degC.d]  GDD since bud-break (NaN/<=0 before bud-break -> dormant)
%     harvested    [0/1]     1 after harvest -> storage-refill branch
%     p            struct  CLM5-FruitTree citrus (nc idx 35) params:
%                  .arooti 0.7 .arootf 0.2 .arootf2 0.08 .astemf 0.22 .fleafi 0.85
%                  .aleafstor 0.3 .ft_bfact -0.5 (J-Norman shape)
%                  .grnfill 400 (huigrain) .grnrp 1100 (huiripe) .ft_GDDlfmat 1350 (huilfmat)
%                  .hybgdd 2880 (gddmaturity)   [all GDD offsets from bud-break]
%
%   Method (CLM5-FruitTree NutrientCompetitionFlexibleCNMod.F90:1645-1754 +
%   CNPhenologyMod.F90:1950-1954). Phase boundaries (g): grnfill(400) | grnrp(1100) |
%   lfmat(1350) | hybgdd(2880):
%     P1 vegetative  [0,grnfill)   : aroot arooti->arootf, leaf=(1-aroot)*fleafi*Cur, stem=rem, repr=0
%     P2 fruit-fill  [grnfill,grnrp): aroot arooti2->arootf2 (denom grnrp), leaf prefactor = SAVED
%                                     aleafi (compounded taper, :1695), stem->astemf, repr=remainder
%     P3 ripening    [grnrp,lfmat) : aroot~arootf2, leaf prefactor=fleafi again (:1720), stem=astemf, repr=rem
%     P4 leaf-mature [lfmat,hybgdd): aroot=arootf2, leaf=1e-5, stem=astemf, repr=remainder (rises to ~0.7)
%     P5 maturity/harvest          : display=0 (all organs); C -> C_labile surplus -> Step 3c reserve refill
%                                     (CLM5 P5 storage split mapped onto COMPASS's surplus-fed dual reserve)
%   Cur(g) = J-Norman taper (= vegetation.carbon.leafGDDtaper, exp form, bfact). Root re-anchors
%   arooti2 (=aroot at grnfill) and arooti3 (=aroot at grnrp) are DETERMINISTIC, computed inline.
%
%   STEM NOTE: CLM5's stem decline is a daily-Markovian power law astem*=（1-min(.,1)^allconss)
%   (declfact=4, allconss=1.5) that compounds to the astemf floor by ~grnrp. The per-step Markovian
%   does not translate to COMPASS's hourly + dual C/N caller without persistent state, and its
%   envelope reaches astemf by ~grnrp; it is reproduced here by a linear astemi->astemf decline over
%   P2 then held at the astemf floor (P3/P4) -- same floor, same arrival window; the narrow stem band
%   (0.22-0.29) makes the fruit=remainder path-insensitive. [SOURCE: F90:1685-1745]
%
%   Author: Zhu Jizixing 2026-06-13 (stage-13a 5-phase faithful port; replaces 2-phase linear)

    grnfill = p.grnfill;        % [degC.d] huigrain offset (fruit-fill start, 400)
    grnrp   = p.grnrp;          % [degC.d] huiripe offset  (ripening start, 1100)  -- the missing split
    lfmat   = p.ft_GDDlfmat;    % [degC.d] huilfmat offset  (leaf taper -> 0, 1350)
    hybgdd  = p.hybgdd;         % [degC.d] gddmaturity offset (physiological maturity, 2880)
    bfact   = p.ft_bfact;       % [-] J-Norman shape factor (-0.5)
    arooti=p.arooti; arootf=p.arootf; arootf2=p.arootf2; astemf=p.astemf; fleafi=p.fleafi;
    cl01 = @(x) min(max(x,0),1);

    % --- Pre-budburst (dormant): no growth (dormancy gate upstream handles it) ---
    if ~isfinite(gdd_budburst) || gdd_budburst <= 0
        a = struct('leaf', 0, 'root', 0, 'stem', 0, 'repr', 0);
        return;
    end
    g = gdd_budburst;

    % --- P5: maturity / post-harvest -> reserve refill via ZERO display allocation ---
    % CLM5-FruitTree routes the P5 share (aleafstor + (1-aleafstor)/2 split) into the leaf/
    % stem/root STORAGE sub-pools (NutrientCompetitionFlexibleCNMod.F90:1742-1747), i.e. it is
    % reserve loading, not vegetative growth. COMPASS has a surplus-fed root/stem dual reserve refilled
    % target-driven from C_labile SURPLUS (computeCarbonAllocation.m Step 3c, 2026-06-10): storage
    % tops up only from the carbon left after display demand. So the architecture-faithful post-
    % harvest behaviour is ZERO display allocation -> all fixed C stays in C_labile surplus ->
    % Step 3c refills C_root_storage/C_stem_storage toward their wood-based targets. The former
    % display split (leaf=aleafstor, stem=root=(1-aleafstor)/2) instead grew the vegetative display
    % pools post-harvest, inflating live tissue + maintenance respiration AND (via display_total at
    % Step 3c:324) starving the very surplus the reserve refill depends on. Mirrors the dormant
    % branch above, which already returns all-zeros for the same "C -> reserve, not display" intent.
    if harvested || g >= hybgdd
        a.leaf = 0; a.root = 0; a.stem = 0; a.repr = 0;   % [-] display fractions; reserve refill handled by Step 3c surplus
        return;
    end

    % --- Deterministic phase-boundary anchors (CLM5 saves these at the transition; here analytic) ---
    arooti2 = cl01(arooti - (arooti - arootf));                              % aroot at end P1 (g=grnfill) = arootf
    aleafi  = max(0.05, (1 - arooti2) * fleafi * ...
                  vegetation.carbon.leafGDDtaper(grnfill, lfmat, bfact));    % leaf at end P1 (P2 prefactor)
    astemi  = 1 - arooti2 - aleafi;                                          % stem at end P1
    arooti3 = cl01(arooti2 - (arooti2 - arootf2));                           % aroot at end P2 (g=grnrp) = arootf2

    Cur = vegetation.carbon.leafGDDtaper(g, lfmat, bfact);                   % [-] J-Norman 1->0

    if g < grnfill
        % --- P1 vegetative / canopy build ---
        a.root = cl01(arooti - (arooti - arootf) * min(1, g/grnfill));
        a.leaf = max(0.05, (1 - a.root) * fleafi * Cur);
        a.stem = 1 - a.root - a.leaf;
        a.repr = 0;
    elseif g < grnrp
        % --- P2 fruit-fill begins: root -> arootf2 (denom grnrp); leaf prefactor = saved aleafi ---
        a.root = cl01(arooti2 - (arooti2 - arootf2) * min(1, g/grnrp));
        a.leaf = max(1.0e-5, (1 - a.root) * aleafi * Cur);
        a.stem = max(astemf, astemi - (astemi - astemf) * min(1, (g - grnfill)/(grnrp - grnfill)));
        a.repr = max(0, 1 - a.root - a.stem - a.leaf);
    elseif g < lfmat
        % --- P3 ripening: root ~arootf2 (denom lfmat); leaf prefactor = fleafi again ---
        a.root = cl01(arooti3 - (arooti3 - arootf2) * min(1, g/lfmat));
        a.leaf = max(1.0e-5, (1 - a.root) * fleafi * Cur);
        a.stem = astemf;
        a.repr = max(0, 1 - a.root - a.stem - a.leaf);
    else
        % --- P4 leaf maturity / leaf-fall: root flat arootf2, leaf ~0, fruit = remainder (rises) ---
        a.root = arootf2;
        a.leaf = 1.0e-5;
        a.stem = astemf;
        a.repr = max(0, 1 - a.root - a.stem - a.leaf);
    end

    % Renormalize defensively so the four fractions sum to exactly 1
    tot = a.leaf + a.root + a.stem + a.repr;
    if tot > 0
        a.leaf = a.leaf / tot;  a.root = a.root / tot;
        a.stem = a.stem / tot;  a.repr = a.repr / tot;
    end
end
