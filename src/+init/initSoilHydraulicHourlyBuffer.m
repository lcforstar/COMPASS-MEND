function acc = initSoilHydraulicHourlyBuffer(NN, NL, Nefc, Ntotal)
% initSoilHydraulicHourlyBuffer - 初始化土壤水文自适应步长小时级输出缓冲
%                                  Initialize hourly output buffer for soil hydraulic ATS
%
% 函数定义 / Function Definition:
%   acc = initSoilHydraulicHourlyBuffer(NN, NL, Nefc, Ntotal)
%
% 输入 / Inputs:
%   NN     - 节点数 / Number of nodes (scalar)
%   NL     - 土层数 / Number of layers (scalar)
%   Nefc   - 氮素方程开关 / Nitrogen equation switch (0=关闭, 1=开启)
%   Ntotal - 氮输入类型 / Nitrogen input type (1=总氮, 0=液相溶液浓度)
%
% 输出 / Outputs:
%   acc - 缓冲结构体 / Buffer struct
%         所有空间向量均为列向量 (NN×1 或 NL×1)，标量初始化为 0
%
% 使用示例 / Usage Example:
%   acc = init.initSoilHydraulicHourlyBuffer(NN, NL, Nefc, Ntotal);
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-13

% --- [1] ATS 时间控制 / Adaptive time stepping control ---
acc.output_interval  = 3600;   % 输出间隔 (s)
acc.next_output_time = 3600;   % 下一次输出的物理时刻 (s)
acc.out_idx          = 1;      % 输出行索引
acc.sum_dt_acc       = 0;      % 当前小时内已累积时间 (s)

% --- [2] 状态变量缓冲 (NN×1) / State variable buffer ---
acc.Theta_LLL = zeros(NN, 1);
acc.TTT       = zeros(NN, 1);
acc.P_ggg     = zeros(NN, 1);
acc.C_O2      = zeros(NN, 1);
acc.NNN       = zeros(NN, 1);
acc.NNno3     = zeros(NN, 1);
acc.NNurea    = zeros(NN, 1);
if Ntotal == 1
    acc.NNN_TOTAL   = zeros(NN, 1);
    acc.NNno3_TOTAL = zeros(NN, 1);
end

% 冠层/气象标量
acc.Tsh  = 0;  acc.Tbh  = 0;  acc.LAI  = 0;  acc.VPD  = 0;
acc.psil = 0;  acc.psir = 0;  acc.psis = 0;

% --- [3] 边界通量缓冲 (标量 / NL×1) / Boundary flux buffer ---
acc.WaterUptake  = zeros(NL, 1);
acc.QET          = 0;  acc.QEB          = 0;
acc.EVAP         = 0;  acc.Actual_Evap  = 0;  acc.Tp_t         = 0;
acc.Throughfall  = 0;  acc.Precip       = 0;
acc.Infiltration = 0;  acc.Actual_Inf   = 0;
acc.RS    = 0;  acc.QMT   = 0;  acc.QMB   = 0;
acc.cvBot  = 0;  acc.cvTop  = 0;  acc.cvBot1  = 0;  acc.cvTop1  = 0;
acc.cvSnit = 0;  acc.cvSmin = 0;  acc.cvSden  = 0;
acc.cvSinkS = 0; acc.cvSinkS1 = 0;
acc.Cvnit_N2o = 0;  acc.Cvdenit_N2o = 0;  acc.Cvdenit_N2 = 0;

% --- [4] 剖面通量与梯度缓冲 (NL×1) / Profile flux and gradient buffer ---
acc.QL     = zeros(NL,1);  acc.QL_h   = zeros(NL,1);
acc.QL_T   = zeros(NL,1);  acc.QL_A   = zeros(NL,1);
acc.QV     = zeros(NL,1);  acc.QV_h   = zeros(NL,1);
acc.QV_T   = zeros(NL,1);  acc.QVA    = zeros(NL,1);
acc.Qa     = zeros(NL,1);
acc.DPgDZ  = zeros(NL,1);  acc.DhDZ   = zeros(NL,1);  acc.DTDZ   = zeros(NL,1);
acc.KLhBAR = zeros(NL,1);  acc.KLTBAR = zeros(NL,1);  acc.DTDBAR = zeros(NL,1);
acc.Se     = zeros(NL,1);  acc.Sa     = zeros(NL,1);

% --- [5] 氮素通量缓冲 (NL×1，仅 Nefc=1 时分配) / Nitrogen flux buffer ---
if Nefc == 1
    acc.QNH4        = zeros(NL,1);  acc.QNH4_h      = zeros(NL,1);
    acc.QNH4_T      = zeros(NL,1);  acc.QNH4_A      = zeros(NL,1);
    acc.QNH4_conv   = zeros(NL,1);  acc.QNH4_conv_h = zeros(NL,1);
    acc.QNH4_conv_T = zeros(NL,1);  acc.QNH4_conv_A = zeros(NL,1);
    acc.QNH4_disp   = zeros(NL,1);  acc.QNH4_disp_h = zeros(NL,1);
    acc.QNH4_disp_T = zeros(NL,1);  acc.QNH4_disp_A = zeros(NL,1);
    acc.QNO3        = zeros(NL,1);  acc.QNO3_h      = zeros(NL,1);
    acc.QNO3_T      = zeros(NL,1);  acc.QNO3_A      = zeros(NL,1);
    acc.QNO3_conv   = zeros(NL,1);  acc.QNO3_conv_h = zeros(NL,1);
    acc.QNO3_conv_T = zeros(NL,1);  acc.QNO3_conv_A = zeros(NL,1);
    acc.QNO3_disp   = zeros(NL,1);  acc.QNO3_disp_h = zeros(NL,1);
    acc.QNO3_disp_T = zeros(NL,1);  acc.QNO3_disp_A = zeros(NL,1);
end

end
