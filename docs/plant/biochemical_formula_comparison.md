# biochemical.m 逐行代码对照

**左边 = REF**（PhotosynthesisTemperatureEffects.m + Photofun_vec/LeafPhotosynthesis）  
**右边 = NEW**（biochemical.m）  
**标记：✅ 对应 | ⚠️ 差异 | ❌ 错误**

---

## A. 热适应参数（Kattge & Knorr 2007）

```
REF  PhotosynthesisTemperatureEffects.m L124-134          NEW  biochemical.m L144-149
─────────────────────────────────────────────────          ──────────────────────────────────────
T_growth_C = Tgrowth + Tfz;                                Tgrowth_clamped = max(11, min(35, Tgrowth));
T_acclim = max(11+Tfz, min(35+Tfz, T_growth_C));          % (Tgrowth 输入已是 Celsius，钳位等价) ✅
                                                           
R_JV = 2.59 - 0.035 * (T_acclim - Tfz);                   R_JV = 2.59 - 0.035 .* Tgrowth_clamped;    ✅
DeltaS_Vcmax = 668.39 - 1.07 * (T_acclim - Tfz);          DeltaS_Vcmax = 668.39 - 1.07 .* Tgrowth_clamped; ✅
DeltaS_Jmax  = 659.70 - 0.75 * (T_acclim - Tfz);          DeltaS_Jmax  = 659.70 - 0.75 .* Tgrowth_clamped; ✅
DeltaS_TPU = 490;                                          DeltaS_TPU = 490;                            ✅
DeltaS_Rd  = 490;  (在 L182 定义)                          DeltaS_Rd  = 490;                            ✅
```

---

## B. 25°C 参考参数

```
REF  PhotosynthesisTemperatureEffects.m L138-157          NEW  biochemical.m L153-159
─────────────────────────────────────────────────          ──────────────────────────────────────
Param_sun.Jmax25 = R_JV .* Vcmax25_sun_base;               Jmax25 = R_JV .* Vcmax25_base;              ✅
Param_sun.TPU25  = 0.167 .* Param_sun.Jmax25;              TPU25  = 0.167 .* Jmax25;                   ✅
% C3:                                                       
Param_sun.Rd25   = 0.015 .* Vcmax25_sun_base;              Rd25 = RdPerVcmax25 .* Vcmax25_base;        ✅ (Rdparam默认0.015)
% C4:                                                       
Param_sun.Rd25   = 0.025 .* Vcmax25_sun_base;              Rd25 = 0.025 .* Vcmax25_base;               ✅
```

---

## C. 温度响应函数定义

```
REF  PhotosynthesisTemperatureEffects.m L161-166          NEW  biochemical.m L161-162
─────────────────────────────────────────────────          ──────────────────────────────────────
f_Arr = @(DeltaHa, T)                                      f_Arr = @(Ha, Tk)
  exp((DeltaHa ./ (Tref * Rgas)) .* (1 - Tref ./ T));       exp((Ha ./ (Tref .* Rgas)) .* (1 - Tref ./ Tk)); ✅

f_Peak = @(DeltaHd, DeltaS, T)                            f_Peak = @(Hd, Ds, Tk)
  (1 + exp((Tref.*DeltaS - DeltaHd)./(Tref*Rgas))) ./       ((1 + exp((Tref .* Ds - Hd) ./ (Tref .* Rgas))) ./
  (1 + exp((T.*DeltaS - DeltaHd)./(Rgas.*T)));              (1 + exp((Tk .* Ds - Hd) ./ (Rgas .* Tk))));    ✅
```

---

## D. 动力学常数 Ha / Hd

```
REF  PhotosynthesisTemperatureEffects.m L170-172          NEW  biochemical.m L164-170
─────────────────────────────────────────────────          ──────────────────────────────────────
Ha_Kc = 79430;                                             % 内嵌在调用里:
Ha_Ko = 36380;                                             Kc = 404.9e-6.*p_bar.*f_Arr(79430, T)...    ✅
Ha_Gam = 37830;                                            Ko = 278.4e-3.*p_bar.*f_Arr(36380, T)...    ✅
Ha_Vc  = 72000;  Hd_Vc  = 200000;                         Gamma_star = 42.75e-6.*p_bar.*f_Arr(37830,T) ✅
Ha_J   = 50000;  Hd_J   = 200000;                         
Ha_TPU = 72000;  Hd_TPU = 200000;  ← 注意是200000!        Vcmax = ...f_Arr(72000,T).*f_Peak(200000,...) ✅
Ha_Rd  = 46390;  Hd_Rd  = 150650;                         Jmax  = ...f_Arr(50000,T).*f_Peak(200000,...) ✅
                                                           TPU   = ...f_Arr(72000,T).*f_Peak(200000,...) ✅ (已修正)
                                                           Rd    = ...f_Arr(46390,T).*f_Peak(150650,...) ✅
```

---

## E. Vcmax / Jmax / TPU / Rd 温度响应

```
REF  PhotosynthesisTemperatureEffects.m L211-214          NEW  biochemical.m L167-170
─────────────────────────────────────────────────          ──────────────────────────────────────
psnParams.sun.Vcmax = Param_sun.Vcmax25                    Vcmax = Vcmax25_base
  .* f_Arr(Ha_Vc, Tk)                                        .* f_Arr(72000, T)
  .* f_Peak(Hd_Vc, DeltaS_Vcmax, Tk)                        .* f_Peak(200000, DeltaS_Vcmax, T);         ✅
  .* f_water_biochem;                                      % f_water_biochem → sfactor 在下方乘入

psnParams.sun.Jmax  = Param_sun.Jmax25                     Jmax = Jmax25
  .* f_Arr(Ha_J, Tk)                                         .* f_Arr(50000, T)
  .* f_Peak(Hd_J, DeltaS_Jmax, Tk)                          .* f_Peak(200000, DeltaS_Jmax, T);          ✅
  .* f_water_biochem;

psnParams.sun.TPU   = Param_sun.TPU25                      TPU = TPU25
  .* f_Arr(Ha_TPU, Tk)                                       .* f_Arr(72000, T)
  .* f_Peak(Hd_TPU, DeltaS_TPU, Tk)                         .* f_Peak(200000, DeltaS_TPU, T);           ✅
  .* f_water_biochem;

psnParams.sun.Rd    = Param_sun.Rd25                       Rd = Rd25
  .* f_Arr(Ha_Rd, Tk)                                        .* f_Arr(46390, T)
  .* f_Peak(Hd_Rd, DeltaS_Rd, Tk)                           .* f_Peak(150650, DeltaS_Rd, T);            ✅
  .* f_water_biochem;                                      % ⚠️ Rd 不乘 sfactor（维持呼吸，胁迫不抑制）
```

---

## F. 水力胁迫乘入时机

```
REF  PhotosynthesisTemperatureEffects.m L191-199          NEW  biochemical.m L174-177
─────────────────────────────────────────────────          ──────────────────────────────────────
% 水分胁迫：阶段式 betaL                                   % 水力胁迫：VC(ψ_leaf)/gpmax
if betaL > 0.5                                             stressfactor = sfactor;
    f_water_biochem = 1.0;                                 Vcmax = Vcmax .* stressfactor;
else                                                       Jmax  = Jmax  .* stressfactor;
    f_water_biochem = max(0, betaL / 0.5);                 TPU   = TPU   .* stressfactor;
end                                                        % (Rd 不乘，见上)

% f_water_biochem 在 Vcmax/Jmax/TPU/Rd 计算时同步乘入     % sfactor 在四个参数算完后统一乘
```

---

## G. PSII 吸收光能 / 潜在电子传输

```
REF  Photofun_vec/LeafPhotosynthesis L172-183             NEW  biochemical.m L180-186
─────────────────────────────────────────────────          ──────────────────────────────────────
% C3 专用
f = 0.15;                                                  % po0 = Kp/(Kf+Kd+Kp) ≈ 0.80（T相关）
I_PSII = 0.5 .* (1 - f) .* PAR;                           Jpot = 0.5 .* po0 .* Q .* fPAR;
% = 0.5 * 0.85 * PAR = 0.425 * PAR                        % ≈ 0.5 * 0.80 * Q = 0.40 * Q
                                                           % ⚠️ 系数 0.425 vs ~0.40，差 ~6%

% Jmax 共限制二次方程
Theta_PSII = 0.7;                                          thetaJ = 0.7;                               ✅
a_quad = Theta_PSII;                                       
b_quad = -(I_PSII + Jmax);                                 Jact = sel_root(thetaJ,
c_quad = I_PSII .* Jmax;                                     -(Jpot + Jmax),
J = (-b_quad - sqrt(b_quad.^2-4*a_quad*c_quad))             Jpot .* Jmax,
    ./ (2 .* a_quad);  % 取小根                              -1);  % dsign=-1 = 小根              ✅
```

---

## H. Rubisco 限制（Vc）

```
REF  LeafPhotosynthesis L187-189                          NEW  biochemical.m L286-287
─────────────────────────────────────────────────          ──────────────────────────────────────
Km = Kc .* (1 + Oa ./ Ko);                                MM_consts = Kc .* (1 + O_bar ./ Ko);       ✅
Ac = Vcmax .* (Ci - Gstar) ./ (Ci + Km);                  Vc_loc = Vcmax .* (Ci_bar - Gamma_star)
                                                             ./ max(MM_consts + Ci_bar, 1e-12);        ✅
```

---

## I. 光限制（Ve）— ⚠️ 最大差异

```
REF  LeafPhotosynthesis L185                              NEW  biochemical.m L287-288
─────────────────────────────────────────────────          ──────────────────────────────────────
Al = (J ./ 4) .* (Ci - Gstar)                             cpe_loc = ((Ci_bar - Gamma_star)
     ./ (Ci + 2 .* Gstar);                                  ./ max(Ci_bar + 2.*Gamma_star, 1e-12))
                                                             .* effcon;        % effcon = 1/5
                                                           Ve_loc = Jact .* cpe_loc;
                                                           
% REF: 系数 1/4                    ⚠️⚠️                   % NEW: 系数 1/5（SCOPE 遗留）
% Ve_REF = J/4 * (Ci-Γ*)/(Ci+2Γ*)                        % Ve_NEW = Jact/5 * (Ci-Γ*)/(Ci+2Γ*)
% CLM5 标准是 1/4                                         % 结果：Ve 偏低 20%，需你决定是否修改
```

---

## J. TPU 限制（Vs）和共限制平滑

```
REF  LeafPhotosynthesis L192-209                          NEW  biochemical.m L290-292
─────────────────────────────────────────────────          ──────────────────────────────────────
Ae = 3 .* TPU;                                             Vs_loc = 3 .* TPU;                          ✅

% 第一步平滑：min(Ac, Al)                                 % 第一步平滑：min(Vc, Ve)
beta1 = 0.98;                                              V_loc = sel_root(atheta,           % atheta=0.8
b1 = -(Ac + Al);                                             -(Vc_loc + Ve_loc),
c1 = Ac .* Al;                                               Vc_loc .* Ve_loc,
Acol1 = (-b1 - sqrt(b1^2-4*beta1*c1)) / (2*beta1);         sign(Gamma_star - Ci_bar));
% ⚠️ beta1=0.98 vs atheta=0.8                             % ⚠️ 平滑因子 0.98 vs 0.8

% 第二步平滑：min(Acol1, Ae)                              % 第二步平滑：min(V, Vs)
beta2 = 0.95;                                              Ag_loc = sel_root(0.98,            % 硬编码0.98
b2 = -(Acol1 + Ae);                                          -(V_loc + Vs_loc),
c2 = Acol1 .* Ae;                                            V_loc .* Vs_loc, -1);
Agross = (-b2 - sqrt(b2^2-4*beta2*c2)) / (2*beta2);       % ⚠️ beta2=0.95 vs 0.98
```

---

## K. 气孔导度 gs（Medlyn）

```
REF  Solve_An_gs_Ci_Coupled.m L128-152（含gb修正）        NEW  biochemical.m L304-312 medlyn_gs()
─────────────────────────────────────────────────          ──────────────────────────────────────
% gb 包含在内（冠层整体求解）
Cs_est = Ca_pa_k - (An*1e-6 / gb_co2) * Pa_k;             % Cs 由 ebal 传入（已含 gb 修正）
D_leaf_kPa = (es_sat - ea_k) / 1000;                       VPD = max(0.05, (es_hPa - eb) * 0.1); % kPa ✅
D_leaf_kPa = max(0.05, D_leaf_kPa);                        
d = 1.6 * An / (Cs/Pa * 1e6);                             
% 二次方程含 gb：
b = 2*(g0*1e-6 + d) + (g1*d)^2 / (gb*D);                  % 简化叶级公式（无 gb 项）：
c = (g0*1e-6)^2 + [2*g0*1e-6 + d*(1-g1^2/D)] * d;        gs = g0 + 1.6*(1+g1/sqrt(VPD))*Aeff/Cs
gs_sol = (-b + sqrt(b^2-4c)) / 2;  % 大根                 gs = max(g0, gs);
                                                           % ⚠️ 架构差异：REF含gb，NEW不含
                                                           % 等价：ebal 的 heatfluxes 层处理 gb
```

---

## L. Ci 迭代

```
REF  Solve_An_gs_Ci_Coupled.m L87-170                    NEW  biochemical.m L205-220
─────────────────────────────────────────────────          ──────────────────────────────────────
for iter = 1:100                                           for iter_idx = 1:50                         ✅ (次数不同)
    % 调用 Photofun_vec 计算 An(Ci)                        % 直接调用 compute_assimilation(Ci*ppm2bar)
    An_iter = Photofun_vec(..., Ci_iter, ...);             [A,Ag,...] = compute_assimilation(Ci.*ppm2bar);

    % 更新 gs（含 gb 二次方程）                            % 更新 gs（简化）
    gs_sol = ...;                                          gs = medlyn_gs(Cs, A, VPD, g0, g1, type);
    gs_iter = 0.5*gs_iter + 0.5*gs_sol;  % 阻尼           
                                                           
    % 更新 Ci（含 gb 阻力）                                % 更新 Ci（Cs 已含 gb）
    Ci_new = Ca - An_mol*Pa*(1.4/gb + 1.6/gs);            Ci_next = Cs - 1.6*A/max(gs,1e-9);
    Ci_new = max(0.1, min(Ci_new, Ca));                    Ci_next = max(minCi*Cs, Ci_next);
                                                           Ci_next = min(Cs, Ci_next);
                                                           
    if abs(Ci_new - Ci_prev) < 0.01 [Pa]                  delta_pa = abs((Ci_next-Ci_ppm)*ppm2bar*1e5);
        break;                                             if all(delta_pa < 0.01 [Pa]): break;        ✅
    end
    Ci_iter = 0.5*Ci_iter + 0.5*Ci_new;  % 阻尼           Ci_ppm += 0.5*(Ci_next - Ci_ppm);           ✅
end
```

---

## M. 净光合 A

```
REF  LeafPhotosynthesis L212                              NEW  biochemical.m L301
─────────────────────────────────────────────────          ──────────────────────────────────────
An = Agross - Rd;                                          A_loc = Ag_loc - Rd;                        ✅
```

---

## 差异汇总（需你决定是否修改）

```
编号  位置              REF 值         NEW 当前值    影响估计
──────────────────────────────────────────────────────────────
D1   effcon（Ve系数）   1/4            1/5           Ve 偏低 20%，影响光限制 A
D2   第一步平滑 θ       beta1 = 0.98   atheta = 0.8  过渡区 A 差 1-3%
D3   第二步平滑 θ       beta2 = 0.95   0.98          TPU 限制区 A 差 <1%
D4   PSII 系数          0.425          ~0.40         Jpot 差 ~6%
D5   Rd 胁迫            随 f_water 缩放 不缩放        干旱极端条件下 Rd 偏高
D6   Ci 迭代次数        100            50            极端条件下可能未收敛
```

**D1 最优先**：effcon 从 1/5 改为 1/4 是最接近 REF 的修改，且有 CLM5 文献依据。
改之后 Ve 提升 25%，A 在光限制时约升 15-20%，需要重新确认 Vcmax25 标定值。
