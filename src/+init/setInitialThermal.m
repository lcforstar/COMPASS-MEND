function ThermalConductivity = setInitialThermal(SoilProperties, cfg)
% setInitialThermal - 从 SoilProperties 按需提取热学参数
%
% 输入:
%   SoilProperties - legacy SoilProperties 对象（已完成热学参数计算）
%   cfg            - legacy ModelSettings 对象
%
% 输出:
%   ThermalConductivity - 热学参数结构体，供 calculateSoilThermalProperites 使用
%
% lThermCond 对应方法:
%   1 = de Vries (1963)           适用: 矿质土，STEMMUS 默认
%   2 = Peters-Lidard / Johansen  适用: 冻融/通用矿质土
%   3 = Tian (2016) 简化 de Vries 适用: 高有机质/泥炭土
%   4 = Farouki 热容 + PL 热导    适用: 冻融 + 更准确热容

    % -----------------------------------------------------------------------
    % 公共字段（所有方法都需要）
    % -----------------------------------------------------------------------
    ThermalConductivity.HCAP  = SoilProperties.HCAP;   % [1×6] 各相热容
    ThermalConductivity.SF    = SoilProperties.SF;     % [1×6] 各相形状因子
    ThermalConductivity.TCA   = SoilProperties.TCA;
    ThermalConductivity.GA1   = SoilProperties.GA1;
    ThermalConductivity.GA2   = SoilProperties.GA2;

    ThermalConductivity.HCD   = SoilProperties.Heat_capacity_dry;
    ThermalConductivity.S1    = SoilProperties.S1;
    ThermalConductivity.S2    = SoilProperties.S2;
    ThermalConductivity.ZETA0 = SoilProperties.ZETA0;
    ThermalConductivity.CON0  = SoilProperties.CON0;
    ThermalConductivity.PS1   = SoilProperties.PS1;
    ThermalConductivity.PS2   = SoilProperties.PS2;
    ThermalConductivity.GB1   = SoilProperties.GB1;
    ThermalConductivity.GB2   = SoilProperties.GB2;

    % -----------------------------------------------------------------------
    % Tian (2016) 简化 de Vries（lThermCond == 3）
    % -----------------------------------------------------------------------
    if cfg.lThermCond == 3
        ThermalConductivity.TPS1  = SoilProperties.TPS1;
        ThermalConductivity.TPS2  = SoilProperties.TPS2;
        ThermalConductivity.TCON0 = SoilProperties.TCON0;
    end

    % -----------------------------------------------------------------------
    % Peters-Lidard / Johansen（lThermCond == 2 或 4）
    % -----------------------------------------------------------------------
    if cfg.lThermCond == 2 || cfg.lThermCond == 4
        ThermalConductivity.TCON_dry = SoilProperties.TCON_dry;
        ThermalConductivity.RHO_bulk = SoilProperties.Bulk_density;
        % TCON_s：Johansen(方法2) 几何平均 / Farouki(方法4) 线性加权，两者不同
        if cfg.lThermCond == 4
            ThermalConductivity.TCON_s = SoilProperties.TCON_s_Farouki;
        else
            ThermalConductivity.TCON_s = SoilProperties.TCON_s;
        end
    end

    % -----------------------------------------------------------------------
    % Farouki 热容（lThermCond == 4）
    % -----------------------------------------------------------------------
    if cfg.lThermCond == 4
        ThermalConductivity.FEHCAP = SoilProperties.FEHCAP;
    end

end
