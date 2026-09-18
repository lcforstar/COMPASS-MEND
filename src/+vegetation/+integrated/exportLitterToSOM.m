function litter_pack = exportLitterToSOM(NitrogenVariables, turn_diag, C_turn_diag, pftParams, KT)
%EXPORTLITTERTOSOM  Package per-KT C/N litter for SOM module (师弟 M3 接口).
%
%   每 KT 调一次, 把植物模块当 KT 产生的 litter (C + N + CN ratio) 打包给 SOM.
%
%   Returns:
%     litter_pack  struct  ready-to-send for soil.organic.SOMupdate(litter_pack, ...)
%       .C_leaf_litter    [g C/m²/KT]  叶凋落物 C
%       .N_leaf_litter    [g N/m²/KT]  叶凋落物 N
%       .C_root_litter    [g C/m²/KT]  根周转 C
%       .N_root_litter    [g N/m²/KT]  根周转 N
%       .CN_leaf_litter   [-]          叶 litter 品质 (= CN_L)
%       .CN_root_litter   [-]          根 litter 品质 (= CN_target_root)
%       .C_harvest_today  [g C/m²/KT]  harvest 出地 (不给 SOM, 但报告)
%       .N_harvest_today  [g N/m²/KT]
%       .KT               [-]          time step
%
%   注: harvest 出地 (商业果实带走), 不进 SOM, 但 export 接口里报告供诊断.
%
%   Reference:
%     cn_coupling §7b.8 (M2-B.8 师弟 SOM 接口规范)
%
%   Variable class (per compass-coding-standard §1):
%     NitrogenVariables, turn_diag, C_turn_diag — forcing slice (input)
%     pftParams — parameters
%     litter_pack — diagnostic (output, 给师弟 SOM 消费)
%
%   Usage:
%     M3 SOM 模块 (师弟做) 调:
%       [S_min_NH4_today, R_h_today, dCSOM] = soil.organic.SOMupdate( ...
%           litter_pack, SoilState, SOMpool_prev, dt);
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.8 S6)

    litter_pack = struct( ...
        'KT',              KT, ...
        'C_leaf_litter',   C_turn_diag.C_leaf_litter_today, ...
        'N_leaf_litter',   turn_diag.leaf_litter_today, ...
        'C_root_litter',   C_turn_diag.C_root_litter_today, ...
        'N_root_litter',   turn_diag.root_litter_today, ...
        'CN_leaf_litter',  pftParams.CN_L, ...
        'CN_root_litter',  pftParams.CN_target_root, ...
        'C_harvest_today', C_turn_diag.C_harvest_today, ...
        'N_harvest_today', turn_diag.harvest_today);
end
