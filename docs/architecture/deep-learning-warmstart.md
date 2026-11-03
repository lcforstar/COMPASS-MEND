# +solver/+warmstart

COMPASS 神经热启动属于数值求解辅助能力，现位于 `+solver/+warmstart`。五大物理域不依赖该可选子包；它也不拥有物理状态。

## 设计原则

1. **物理求解器是 ground truth**。深度学习模块只读取物理求解器的状态/矩阵，不修改求解器的核心逻辑。
2. **可关闭**。所有钩子默认 disable，通过环境变量或 `cfg.lDeepLearning` 开关激活。物理求解不会因为深度学习包不在而崩。
3. **MATLAB ↔ Python 严格对齐**。MATLAB 端导出快照、Python 端复刻同一公式，每次接入前先用 `compare_compass_water_residual_snapshot.py` 验证 `diff_max_abs < 1e-6`。

## 子模块

```
+solver/
  +warmstart/        ← 神经热启动 (PFEM 论文风格)
                       给 Picard 迭代提供更好的初值，加速收敛
                       当前状态：MATLAB 导出快照 + Python PINO 训练已通
                       待完成：MATLAB warm-start hook + 残差门控
```

### 未来可加

```
+emulator/           ← 整步代理模型 (替换 / 加速整个 KT 求解)
+super_resolution/   ← 时空超分 (从粗网格还原细网格)
+correction/         ← 后处理偏差校正
```

## 当前模块: +warmstart

| 文件 | 用途 |
|---|---|
| `calculateLinearWaterResidual.m` | 算 COMPASS 三对角水方程残差 `A*hh - RHS` |
| `saveWaterResidualSnapshotIfRequested.m` | 导出 C4/C4_a/RHS/hh 为 .mat 供 Python 训练 |

调用语法：
```matlab
r = solver.warmstart.calculateLinearWaterResidual(WaterMatrix, RHS, hh, cfg);
solver.warmstart.saveWaterResidualSnapshotIfRequested( ...
    WaterMatrix, RHS, SoilVariables, CHKWATER, cfg, KT, Delt_t, 'after_solve');
```

启用快照导出：
```powershell
$env:STEMMUS_WATER_RESIDUAL_SNAPSHOT_DIR='F:\COMPASSV1.0\data\water_residual_snapshots'
$env:STEMMUS_WATER_RESIDUAL_SNAPSHOT_STAGE='after_solve'
$env:STEMMUS_WATER_RESIDUAL_SNAPSHOT_MAX='200'
```

## Python 侧对应位置

- `tools/compass_richards_residual.py` — Python 端可微残差（PINO loss 核心）
- `tools/train_transolver_warmstart.py` — Transolver 训练入口
- `tools/compare_compass_water_residual_snapshot.py` — 对齐验证脚本
- `docs/neural_warmstart_pfem_to_compass.md` — 方法论 + 对照 PFEM 论文
