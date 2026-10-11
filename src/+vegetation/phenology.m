classdef phenology
    % Plant.phenology  苹果树/通用植被物候 + LAI 动态模拟（EPIC/SWAT 框架）
    %
    % 实现来源：
    %   Yang et al. (2026) Agricultural Water Management 323, 110074
    %   SWAT 2012 / EPIC 作物生长模块 (Arnold et al., 2012)
    %
    % 调用示例：
    %   [CHU, HU, HUI, DRDS, RDS, HUF, dHUF] = ...
    %       Plant.phenology.CalculateHeatUnit(Tmax, Tmin, sowing_day, T_b, PHU, l1, l2);
    %
    %   [LAI_daily, f_min] = Plant.phenology.CalculateLAI( ...
    %       cfg.Veg, cfg.Grid, cfg.Time, ...       % txt 配置 / 网格 / 时间
    %       HU, CHU, dHUF, HUI, ...               % CalculateHeatUnit 输出
    %       Tmax, Tmin, ...                        % MeteoData
    %       state, SoilProperties, ...             % 状态变量 / 土壤性质
    %       ML_root);

    methods (Static)
        % ==================================================================
        function [CHU_record, HU, HUI_record, DRDS_list, RDS_list, HUF, dHUF] = ...
                CalculateHeatUnit(Tmax, Tmin, sowing_day, T_b, PHU, FRGW1, LAIMX1, FRGW2, LAIMX2, clockDates)
            % CalculateHeatUnit  逐日累积热量单位 + frLAImx 曲线
            %
            % 输入：
            %   Tmax, Tmin        – 日最高/最低气温 (°C)，等长列向量
            %   sowing_day        – 萌芽开始的 DOY 或数组下标
            %   T_b   (T_BASE)    – 基础温度 (°C)，Eq.1 中的 Tbase，论文取 5
            %   PHU               – 全生育期潜在热量单位总和 (°C·d)，Eq.2，论文取 4555
            %   FRGW1, LAIMX1     – 最优 LAI 曲线第1控制点：frPHU=0.1 时 frLAImx=0.15
            %   FRGW2, LAIMX2     – 最优 LAI 曲线第2控制点：frPHU=0.5 时 frLAImx=0.75
            %
            %   l1/l2 从两控制点反算（不再直接输入）：
            %     l2 = [ln(F1(1-L1)/L1) - ln(F2(1-L2)/L2)] / (F2-F1)
            %     l1 = ln(F1(1-L1)/L1) + l2*F1
            %
            % 输出：
            %   CHU_record   – 累积热量单位 CHU (°C·d)
            %   HU           – 逐日热量单位 (°C·d)
            %   HUI_record   – frPHU = CHU/PHU，即热量单位分数，Eq.3
            %   DRDS_list    – 逐日 frPHU 增量（= HU/PHU），与 HUI_record 互推
            %   RDS_list     – 逐日累积 frPHU（= HUI_record，保留接口兼容）
            %   HUF          – frLAImx，最大 LAI 分数 S 形曲线，Eq.6
            %   dHUF         – diff(HUF)，在 CalculateLAI 中用作 ΔfrLAImx

            if length(Tmax) ~= length(Tmin)
                error('Plant:phenology:DimMismatch', 'Tmax 和 Tmin 长度必须相同');
            end
            % Change: 2026-09-30 Wenjie Liu -- optional explicit calendar
            % drives recurring winter-crop seasons; legacy calls unchanged.
            if nargin>=10 && ~isempty(clockDates)
                [CHU_record,HU,HUI_record,DRDS_list,RDS_list,HUF,dHUF] = ...
                    vegetation.integrated.phenology.cropHeatUnits( ...
                    Tmax,Tmin,clockDates,sowing_day,T_b,PHU,FRGW1,LAIMX1,FRGW2,LAIMX2);
                return;
            end

            % 由两控制点反算 l1, l2（Eq.6 逆推）
            % 在 (FRGW1, LAIMX1) 和 (FRGW2, LAIMX2) 处展开：
            %   l1 - l2*F = ln(F*(1-L)/L)
            A1 = log(FRGW1 * (1 - LAIMX1) / LAIMX1);
            A2 = log(FRGW2 * (1 - LAIMX2) / LAIMX2);
            l2 = (A1 - A2) / (FRGW2 - FRGW1);
            l1 = A1 + l2 * FRGW1;

            total_days   = length(Tmax);
            HU           = zeros(total_days, 1);
            CHU_record   = zeros(total_days, 1);
            HUI_record   = zeros(total_days, 1);
            HUF          = zeros(total_days, 1);
            DRDS_list    = zeros(total_days, 1);
            RDS_list     = zeros(total_days, 1);
            CHU = 0;
            RDS = 0;

            for d = sowing_day:total_days
                Tav      = (Tmax(d) + Tmin(d)) / 2;
                % Eq.1：日热量单位
                HU(d)    = max(Tav - T_b, 0);
                % Eq.3 的日增量（= HU/PHU）
                DRDS     = max((Tav - T_b) / PHU, 0);
                RDS      = RDS + DRDS;
                DRDS_list(d) = DRDS;
                RDS_list(d)  = RDS;          % 等价于 CHU/PHU（保留接口）

                % 累积 CHU
                CHU         = CHU + HU(d);
                CHU_record(d) = CHU;

                % Eq.3：frPHU
                frPHU        = CHU / PHU;
                HUI_record(d) = frPHU;

                % Eq.6：frLAImx = frPHU / (frPHU + exp(l1 - l2*frPHU))
                HUF(d) = frPHU / (frPHU + exp(l1 - l2 * frPHU));

                if RDS >= 1.0
                    % 成熟后填平，防止 diff 产生虚假负跳变
                    CHU_record(d:end)  = CHU;
                    HUI_record(d:end)  = frPHU;
                    HUF(d:end)         = HUF(d);
                    RDS_list(d:end)    = RDS;
                    break;
                end
            end
            % diff(HUF) → ΔfrLAImx；末尾补零保持长度一致
            dHUF = [diff(HUF); 0];
        end

        % ==================================================================
        function [LAI_daily, f_min] = CalculateLAI( ...
                veg, grid, time, ...
                HU, CHU_record, dHUF, HUI_record, ...
                Tmax, Tmin, ...
                state, soil, ...
                ML_root, LAI_init)
            % CalculateLAI  逐日 LAI 动态模拟 + 各类胁迫因子
            % (注释与参数列表保持不变，在此省略以节省空间，请保留你原有的完整注释)

            % -----------------------------------------------------------------
            % [A] 从 veg 对象解包 txt 配置参数
            % -----------------------------------------------------------------
            GMHU       = veg.GMHU;
            XLAI       = min(veg.yr_cur / veg.yr_full, 1) * veg.LAI_max;
            T_b        = veg.T_BASE;
            T_op       = veg.T_OPT;
            DLAI       = veg.DLAI;
            ALAI_MIN   = veg.ALAI_MIN;
            C_f        = veg.C_f;
            sowing_day = veg.sowing_day;
            UN         = veg.UN;
            UNO        = veg.UNO;
            beta_d     = veg.beta_decay;
            H1 = veg.Feddes_H1; H2 = veg.Feddes_H2;
            H3 = veg.Feddes_H3; H4 = veg.Feddes_H4;

            % -----------------------------------------------------------------
            % [B] 逻辑控制开关
            % -----------------------------------------------------------------
            use_nut_stress   = false; if isfield(veg, 'use_nut_stress'),   use_nut_stress   = veg.use_nut_stress;   end
            use_wd_stress    = false; if isfield(veg, 'use_wd_stress'),    use_wd_stress    = veg.use_wd_stress;    end
            use_water_stress = false; if isfield(veg, 'use_water_stress'), use_water_stress = veg.use_water_stress; end
            is_deciduous     = true;  if isfield(veg, 'is_deciduous'),     is_deciduous     = veg.is_deciduous;     end

            % -----------------------------------------------------------------
            % [C] 从结构体解包 grid / time / state / soil 字段
            % -----------------------------------------------------------------
            NL        = grid.NL;
            DeltZ     = grid.DeltZ;
            DayNum    = time.JN(time.KT);
            Theta_LL  = state.Theta_LL;
            hh        = state.h;
            Theta_f   = soil.Field_mc;
            Theta_s   = soil.Theta_s;
            total_days = length(HU);

            % -----------------------------------------------------------------
            % [D] 变量初始化与预分配
            % -----------------------------------------------------------------
            LAI_daily  = zeros(total_days, 1);
            f_T        = zeros(total_days, 1);
            f_nut      = ones(total_days, 1);
            f_wd       = ones(total_days, 1);
            f_theta    = ones(total_days, 1);
            f_min      = zeros(total_days, 1);

            if nargin < 13 || isempty(LAI_init)
                LAI_init = ALAI_MIN;
            end

            LAI_daily(:) = LAI_init;
            LAI_prev = LAI_init;
            LAI_D    = 0;
            HUI_D    = DLAI;

            if DayNum < sowing_day
                return;
            end

            % =================================================================
            % 【核心亮点】：引入本年不可逆的"休眠锁"
            % 由于此函数通常按年调用（一次传进一整年的 Tmax/Tmin），
            % 可以在年内循环外部定义这把锁，每年调用时自动刷新为 false。
            % =================================================================
            has_entered_dormancy = false;

            for k = sowing_day:total_days
                Tav = (Tmax(k) + Tmin(k)) / 2;

                % ── 1. 温度胁迫 f_T ──────────────────────────────────────
                R_T = (Tav - T_b) / (T_op - T_b);
                if R_T >= 0 && R_T <= 2
                    f_T(k) = sin(pi/2 * R_T);
                else
                    f_T(k) = 0;
                end

                % ── 2. 养分胁迫 f_nut ─────────────────────────────────────
                if use_nut_stress
                    S_N       = 200 * (UN / UNO);
                    f_nut(k)  = S_N / (S_N + exp(4.065 - 0.0535 * S_N));
                end

                % ── 3. 通风胁迫 f_wd ───────────────────────────────────────
                if use_wd_stress
                    Theta_sum = 0; Theta_ff = 0; Theta_ss = 0;
                    for ML = NL:-1:ML_root
                        Theta_sum = Theta_sum + (Theta_LL(ML,1) + Theta_LL(ML,2)) / 2 * DeltZ(ML);
                        Theta_ff  = Theta_ff  + Theta_f(ML) * DeltZ(ML);
                        Theta_ss  = Theta_ss  + Theta_s(ML) * DeltZ(ML);
                    end
                    root_depth = sum(DeltZ(ML_root:NL));
                    S_AT = (((Theta_sum / root_depth) - (Theta_ff / root_depth)) / ...
                        ((Theta_ss / root_depth) - (Theta_ff / root_depth)) - C_f) / (1 - C_f);
                    if S_AT > 0
                        f_wd(k) = 1 - S_AT / (S_AT + exp(2.901 - 0.0387 * S_AT));
                    else
                        f_wd(k) = 1.0;
                    end
                end

                % ── 4. 水分胁迫 f_theta ────────────────────────────────────
                if use_water_stress
                    sum_alpha = 0;
                    for ML = 1:NL
                        alpha_L = 0; alpha_R = 0;
                        % 左节点
                        hL = hh(ML);
                        if hL > H2 && hL < H1, alpha_L = (hL - H1)/(H2 - H1);
                        elseif hL >= H3 && hL <= H2, alpha_L = 1;
                        elseif hL > H4 && hL < H3, alpha_L = (hL - H4)/(H3 - H4);
                        end
                        % 右节点
                        hR = hh(ML+1);
                        if hR > H2 && hR < H1, alpha_R = (hR - H1)/(H2 - H1);
                        elseif hR >= H3 && hR <= H2, alpha_R = 1;
                        elseif hR > H4 && hR < H3, alpha_R = (hR - H4)/(H3 - H4);
                        end

                        sum_alpha = sum_alpha + (alpha_L + alpha_R)/2 * DeltZ(ML);
                    end
                    f_theta(k) = sum_alpha / sum(DeltZ(1:NL));
                end

                % ── 5. 最小胁迫因子 ───────────────────────────────────────
                f_min(k) = min([f_T(k), f_nut(k), f_wd(k), f_theta(k)]);

                % ── 6. LAI 动态模拟 ───────────────────────────────────────

                % 6.1 检查是否触发不可逆的生物学休眠
                if is_deciduous && ~has_entered_dormancy
                    if k > 300 && (Tmin(k) < 0 || k >= 320)
                        has_entered_dormancy = true; % 一旦上锁，今年后续的循环都不再解锁
                    end
                end

                % 6.2 执行状态更新
                if has_entered_dormancy
                    % 已被死锁：保持光秃秃的树枝状态
                    LAI_daily(k) = ALAI_MIN;
                    LAI_prev     = ALAI_MIN;

                elseif CHU_record(k) <= GMHU
                    if k == 1, LAI_daily(k) = LAI_init; else, LAI_daily(k) = ALAI_MIN; end
                    LAI_prev = LAI_daily(k);

                elseif HUI_record(k) > DLAI
                    if LAI_D == 0
                        LAI_D = LAI_prev; HUI_D = HUI_record(k);
                    end
                    decay = max(0, (1.0 - HUI_record(k)) / (1.0 - HUI_D)) ^ beta_d;
                    LAI_daily(k) = max(ALAI_MIN, LAI_D * decay);
                    LAI_prev = LAI_daily(k);

                else
                    if k == sowing_day
                        LAI_daily(k) = LAI_init;
                    else
                        growth_res   = max(0, 1 - exp(5 * (LAI_prev - XLAI)));
                        LAI_daily(k) = LAI_prev + dHUF(k-1) * XLAI * growth_res * f_min(k);
                    end
                    LAI_prev = LAI_daily(k);
                end

                LAI_daily(k) = max(ALAI_MIN, LAI_daily(k));
            end
        end
    end % methods (Static)
end % classdef
