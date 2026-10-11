function [D_N, diag] = computeOrganDemand(NPP, current_DOY, pftParams, phen_alloc)
%COMPUTEORGANDEMAND  Per-organ N demand from LPJ allocation + target C:N ratios.
%
%   Returns:
%     D_N    struct fields [g N/(m²·s)]: .leaf, .root, .stem, .fruit, .total
%            Per-organ N demand rate. Aggregated by computeOrganDemand from
%            potential growth C × CN_target reciprocal.
%     diag   struct: .G_pot_leaf/root/stem/fruit (potential growth C, kg C/(m²·s))
%                    .fruit_active (0/1 phenology gate)
%
%   Inputs:
%     NPP           [kg C/(m²·s)]  Net primary production from computeNPP
%     current_DOY   [-]           Day of year (1-366) for fruit phenology gate
%     pftParams     struct        Required fields (PFT 表):
%                                   .has_fruit       [0/1]  fruit pool switch
%                                   .f_alloc_leaf, .f_alloc_root, .f_alloc_stem,
%                                   .f_alloc_fruit   [-]    sums to 1
%                                   .CN_L            [-]    leaf C:N (target)
%                                   .CN_target_root, .CN_target_stem,
%                                   .CN_target_fruit [-]
%
%   Method (LPJ Sitch 2003 + Wells & Eissenstat 2002 apple alloc):
%
%     1. Potential growth C per organ (Sitch 2003 Table 1):
%        G_pot,i = f_alloc,i × NPP        [kg C/(m²·s)]
%
%        Fruit gate (M2-B.1b simplified, M2-B.5 phenology升级):
%        if has_fruit == 1 AND fruit_active(DOY) == 1:
%           G_pot,fruit = f_alloc_fruit × NPP
%        else:
%           G_pot,fruit = 0, and fraction redistributed to leaf+root
%
%        Fruit active window for Apple: DOY [150, 290] (June 1 — Oct 17)
%        Covers fruit set + expansion + maturation (Mao 2018 long-term apple).
%
%     2. Per-organ N demand (Sitch 2003 Eq. 12; LPJ allocation by C:N):
%        D_N,i = G_pot,i × 1000 / CN_target,i       [g N/(m²·s)]
%                       ↑↑↑↑
%                      kg→g
%
%        D_N,total = sum over leaf/root/stem/fruit
%
%   Notes:
%     - When NPP = 0 (winter dormancy / freeze), all D_N = 0.
%     - When fruit gate closed, allocated f_alloc_fruit fraction is redistributed
%       to leaf and root proportionally (keeps C balance, avoids "lost C").
%     - This is M2-B.1b simplified; M2-B.5 will use GDD-based phenology with
%       species-specific fruit timing.
%
%   Reference:
%     Sitch S et al. (2003) GCB 9:161-185 Eq. 12 (C → N demand by C:N).
%     Wells CE, Eissenstat DM (2002) Marked differences in survivorship among
%       apple roots. Ecology 83(11):3083-3091.
%     Mao L et al. (2018) Long-term apple orchard N budget. (CN_target_fruit reference)
%
%   Phase C M2-B.1b (2026-05-19): new module. ClimaLand alignment — pure function,
%   no state mutation, called per KT in main loop AFTER computeNPP.
% Change: 2026-06-28 Jizixing Zhu — allometry branch treats empty [] froot_leaf/stem_leaf
%   as 0 (mirror computeCarbonAllocation; avoids 1+f1+[] -> [] -> 1/[] dimension error).

% [13a] phen_alloc optional (4th arg): struct with .gdd_budburst [°C·d], .phen [-]
% (dormancy), .harvested [0/1]. Absent -> dormant defaults.
if nargin < 4 || isempty(phen_alloc), phen_alloc = struct(); end
if ~isfield(phen_alloc, 'gdd_budburst'), phen_alloc.gdd_budburst = NaN; end
if ~isfield(phen_alloc, 'phen'),         phen_alloc.phen         = 1;   end
if ~isfield(phen_alloc, 'harvested'),    phen_alloc.harvested    = 0;   end
gdd_budburst = phen_alloc.gdd_budburst;   % [°C·d] GDD since bud-break
phen         = phen_alloc.phen;           % [-] phenology activity (0=dormant .. 1=active)
harvested    = phen_alloc.harvested;      % [0/1] fruit harvested flag (-> reserve refill)

% --- Organ allocation fractions ---
use_crop_alloc = isfield(phen_alloc, 'is_crop') && phen_alloc.is_crop == 1;
if use_crop_alloc
    f_leaf  = phen_alloc.crop_a_leaf;
    f_root  = phen_alloc.crop_a_root;
    f_stem  = phen_alloc.crop_a_stem;
    f_fruit = phen_alloc.crop_a_repr;
    if ~phen_alloc.crop_live || (isfield(phen_alloc,'crop_mature_flag') ...
            && phen_alloc.crop_mature_flag)
        f_leaf = 0; f_root = 0; f_stem = 0; f_fruit = 0;
    end
    fruit_active = f_fruit > 0;
elseif isfield(pftParams, 'use_fruittree_alloc') && pftParams.use_fruittree_alloc == 1
    % [13a 2026-06-04] FruitTree 4-phase GDD allocation (apple), mirroring
    % computeCarbonAllocation so the C-dominant organ N demand matches the C side.
    a_ft    = vegetation.carbon.fruitTreeAllocation(gdd_budburst, harvested, pftParams);
    f_leaf  = a_ft.leaf;
    f_root  = a_ft.root;
    f_stem  = a_ft.stem;
    f_fruit = a_ft.repr;
    fruit_active = (f_fruit > 0);                 % gates the D_N.fruit demand below
else
    % [CLM5 all-PFT allometry 2026-06-13] Mirror computeCarbonAllocation EXACTLY: natural woody PFTs use
    %   the nc f1/f3 ratios -> NPP shares (S=1+f1+f3: a_leaf=1/S, a_root=f1/S, a_stem=f3/S, no fruit);
    %   grass/crop (froot_leaf=stem_leaf=0) keep the legacy constant f_alloc with the DOY fruit window.
    froot_leaf_p = 0;  stem_leaf_p = 0;
    % Treat an unset (empty []) ratio as 0 (mirror computeCarbonAllocation): a PFT may
    % define froot_leaf but not stem_leaf (e.g. apple). Without ~isempty, 1+f1+[] -> []
    % -> 1/[] "dimensions must agree". [robustness 2026-06-28 Jizixing Zhu]
    if isfield(pftParams, 'froot_leaf') && ~isempty(pftParams.froot_leaf), froot_leaf_p = pftParams.froot_leaf; end
    if isfield(pftParams, 'stem_leaf')  && ~isempty(pftParams.stem_leaf),  stem_leaf_p  = pftParams.stem_leaf;  end
    if froot_leaf_p > 0 || stem_leaf_p > 0
        S_allom = 1 + froot_leaf_p + stem_leaf_p;
        f_leaf  = 1 / S_allom;                        % a_leaf
        f_root  = froot_leaf_p / S_allom;             % a_root = f1/S
        f_stem  = stem_leaf_p  / S_allom;             % a_stem = f3/S
        f_fruit = 0;
        fruit_active = false;                         % natural trees: no reproductive N demand
    else
        % Legacy constant f_alloc (grass/crop), DOY 150-290 fruit window.
        fruit_active = (pftParams.has_fruit == 1) && ...
                       (current_DOY >= 150) && (current_DOY <= 290);
        if fruit_active
            f_leaf  = pftParams.f_alloc_leaf;
            f_root  = pftParams.f_alloc_root;
            f_stem  = pftParams.f_alloc_stem;
            f_fruit = pftParams.f_alloc_fruit;
        else
            % Fruit window closed: fruit share -> stem (woody/storage).
            f_leaf  = pftParams.f_alloc_leaf;
            f_root  = pftParams.f_alloc_root;
            f_stem  = pftParams.f_alloc_stem + pftParams.f_alloc_fruit;
            f_fruit = 0;
        end
    end
end

% [F2 dormancy gate 2026-06-04] No leaf growth while dormant (deciduous): phen < 0.1
% -> zero leaf N demand, mirroring computeCarbonAllocation (C-dominant consistency).
% Leaf C/N stay in reserve for the spring flush.
if ~use_crop_alloc && phen < 0.1
    f_leaf = 0;
end

% [③ peaklai no-redirect 2026-06-12] Mirror the C-side EXACTLY (computeCarbonAllocation): at peak
% LAI stop leaf growth (f_leaf=0), NO redirect to root. Fine root is now bound to leaf (G_pot_root
% = froot_leaf*G_pot_leaf below), so f_leaf=0 auto-drops fine-root N demand at peak; the freed share
% flows to the FRUIT remainder. The old f_root += f_leaf inflated fine-root N (root doubled).
if ~use_crop_alloc && phen >= 0.1 && isfield(pftParams, 'laimx') && pftParams.laimx > 0 ...
        && isfield(pftParams, 'SLA0') && pftParams.SLA0 > 0 ...
        && isfield(phen_alloc, 'C_leaf_prev')
    if isfield(pftParams, 'SLAm') && pftParams.SLAm > 1e-9
        LAI_prev = (pftParams.SLA0 / pftParams.SLAm) * ...
                   (exp(phen_alloc.C_leaf_prev * pftParams.SLAm) - 1);
    else
        LAI_prev = pftParams.SLA0 * phen_alloc.C_leaf_prev;
    end
    if LAI_prev >= pftParams.laimx
        f_leaf = 0;
    end
end

% --- Potential growth C per organ [kg C/(m²·s)] ---
% [allocation restructure 2026-06-12] Mirror computeCarbonAllocation EXACTLY (C-dominant: organ N
%   demand = C allocation / CN_target). FruitTree: fine root = froot_leaf*leaf, coarse root =
%   croot_stem_alloc*stem, fruit = remainder. Non-FruitTree: legacy f_alloc fractions (no coarse root).
% [CLM5-exact 2026-06-13] mirror computeCarbonAllocation EXACTLY (C-dominant: N demand = C/CN_target).
%   FruitTree fork NutrientCompetitionFlexibleCN: organ = a_organ/denom, denom = 1 + a_stem*f2,
%   coarse root = f2*stem (src :563). afroot (Eq3) + arepr (Eq7), NOT froot_leaf-bind/remainder.
use_ft   = isfield(pftParams, 'use_fruittree_alloc') && pftParams.use_fruittree_alloc == 1;
f2_croot = pftParams.croot_stem_alloc;                       % f2 coarse-root:stem (citrus 0.08; natural woody 0.30; grass/crop 0)
denom = 1 + f_stem * f2_croot;                               % [CLM5 all-PFT 2026-06-13] UNIVERSAL (f2=0 -> denom=1); apple value unchanged
G_pot_leaf  = f_leaf * NPP / denom;
G_pot_stem  = f_stem * NPP / denom;
G_pot_root  = f_root * NPP / denom;                          % afroot (or legacy f_alloc; denom=1 non-FT)
G_pot_croot = f2_croot * G_pot_stem;                         % = f2*stem (src :563); 0 if f2=0 (non-woody)
G_pot_fruit = f_fruit * NPP / denom;                         % arepr (Eq7)

% --- Per-organ N demand [g N/(m²·s)] ---
% LPJ Eq. 12: D_N,i = G_pot,i / (C:N)_target,i
% Unit: [kg C/(m²·s)] × 1000 / [-] = [g N/(m²·s)]
ksmall = 1e-30;   % avoid div-zero on CN_target_*
D_N.leaf  = G_pot_leaf  * 1000 / max(pftParams.CN_L, ksmall);
D_N.root  = G_pot_root  * 1000 / max(pftParams.CN_target_root, ksmall);

if pftParams.CN_target_stem > 0
    CN_lw = pftParams.CN_target_stem;
    CN_dw = pftParams.CN_target_stem;
    flivewd = 1;
    if isfield(pftParams, 'CN_lw') && ~isempty(pftParams.CN_lw), CN_lw = pftParams.CN_lw; end
    if isfield(pftParams, 'CN_dw') && ~isempty(pftParams.CN_dw), CN_dw = pftParams.CN_dw; end
    if isfield(pftParams, 'flivewd') && ~isempty(pftParams.flivewd), flivewd = pftParams.flivewd; end

    % Live and dead woody pools use their own target C:N values; otherwise
    % the dead-wood pool receives too much N when CN_dw is high.
    D_N.livestem  = G_pot_stem  * flivewd       * 1000 / max(CN_lw, ksmall);
    D_N.deadstem  = G_pot_stem  * (1 - flivewd) * 1000 / max(CN_dw, ksmall);
    D_N.livecroot = G_pot_croot * flivewd       * 1000 / max(CN_lw, ksmall);
    D_N.deadcroot = G_pot_croot * (1 - flivewd) * 1000 / max(CN_dw, ksmall);
    D_N.stem      = D_N.livestem + D_N.deadstem;
    D_N.croot     = D_N.livecroot + D_N.deadcroot;
else
    D_N.livestem  = 0;
    D_N.deadstem  = 0;
    D_N.livecroot = 0;
    D_N.deadcroot = 0;
    D_N.stem      = 0;
    D_N.croot     = 0;
end

if fruit_active && pftParams.CN_target_fruit > 0
    D_N.fruit = G_pot_fruit * 1000 / pftParams.CN_target_fruit;
else
    D_N.fruit = 0;
end

D_N.total = D_N.leaf + D_N.root + D_N.stem + D_N.croot + D_N.fruit;

% --- Diagnostic ---
diag.G_pot_leaf   = G_pot_leaf;
diag.G_pot_root   = G_pot_root;
diag.G_pot_stem   = G_pot_stem;
diag.G_pot_croot  = G_pot_croot;
diag.G_pot_fruit  = G_pot_fruit;
diag.fruit_active = double(fruit_active);

end

% Change: 2026-09-30 Wenjie Liu -- annual WSC allocation is nitrogen-free;
% physiological maturity stops structural and grain growth N demand.
%   Phase C M2-B.1b (2026-05-19): new module. ClimaLand alignment — pure function,
%   no state mutation, called per KT in main loop AFTER computeNPP.
% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
