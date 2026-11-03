# COMPASS-MEND 1.0

长武冬小麦水、热、气与碳氮过程的 MATLAB 模拟基线。模型以 COMPASS Land System Model 为耦合框架，使用 MEND 描述微生物、酶及有机质碳氮过程。

*A MATLAB baseline for coupled soil water, energy and carbon–nitrogen dynamics in the Changwu winter-wheat field.*

本包对应用户确认的 `K1_mineral_0_20cm` 全年配置。新整理包已在Windows / MATLAB R2025a通过独立4小时原生检查，并通过Windows启动器1步运行；原全年对照与本次短步验证分别记录，详见[验证说明](docs/release/VALIDATION.md)。

## 快速运行

在项目根目录打开 MATLAB：

```matlab
model = run_site('changwu');
```

也可使用系统启动器：

```powershell
.\run.bat changwu
```

```bash
bash run.sh changwu
```

需要已安装且可正常使用的 MATLAB。启动器优先使用 `MATLAB_BIN` 指定的可执行文件，其次从 `PATH` 查找；本包没有声明未经验证的最低 MATLAB 版本。

原全年基线使用MATLAB R2025b；新包已在本地R2025a环境验证；不据此声明其他版本均已验证。核心模型和启动器不依赖Python；可选GUI需要Python与Tkinter，Pillow仅用于可选图像显示。

首次可用启动器限制为24步短跑：

```powershell
$env:COMPASS_MAX_KT = "24"
$env:COMPASS_LOG_WINDOW = "0"
.\run.bat changwu
```

```bash
COMPASS_MAX_KT=24 COMPASS_LOG_WINDOW=0 bash run.sh changwu
```

完整全年运行前取消 `COMPASS_MAX_KT`，或在未设置它的新终端启动。默认计算2000-06-15至2001-06-15，共8760个一小时步。`COMPASS_LOG_WINDOW=0` 关闭额外监控窗口。

结果写入 `data/output/changwu/`，日志位于其 `logs/` 子目录。每次启动在日志旁保存同名`.run.json`运行记录（launch sidecar），用于核对站点、步数上限、实际MATLAB进程PID、运行状态及最后完成KT；创建记录本身不表示计算完成或守恒通过。

## 当前科学基线

- 站点：长武冬小麦，单点一维土柱，物理柱深500 cm。
- 9月15日施入尿素120 kg N/ha及有机肥；本基线两者均按0–20 cm入土配置。120是氮量，不是尿素产品质量。
- 使用 `FAST10` 有机肥入口与 `B_ref` 初始有机氮分配。快入口将有机肥C和纯有机N各10%直接送入游离DOM，余量进入AOM，TAN单独处理。
- 保留根系水、氮吸收；关闭整株植物C/N、自动凋落物输入和指定根残体事件。
- MEND采用 `source_hourly`，DOM采用 `mend_dynamic`；矿肥N、有机肥有机N和有机肥TAN分别被动追踪。

11个MEND分层参数列在 `soil.csv`；`site.toml [microbe]` 保留35个全局/模式/初始化键。Qmax=1.5、Kba=6和Kdes_ref=0.006已由TOML迁至soil，18行原数值不变，旧TOML标量接口保留。详情见[接受基线说明](docs/release/ACCEPTED_BASELINE.md)和[当前吸附接口](docs/soil/biogeochemistry/MEND_LAYER_SORPTION_INTERFACE_20261031.md)。本次迁移回归与此前新包4小时/启动器检查分开记录，结果以本轮实际验证回执为准；未重新跑全年。

## 目录与数据

```text
COMPASS-MEND-1.0/
├── run_site.m, run.bat, run.sh
├── src/                         生产MATLAB源码
├── scripts/                     可移植启动与监控入口
├── data/sites/changwu/          接受的站点、强迫、初值及观测
├── data/model_library_nc/       随包静态模型库
├── data/global/                 可选区域数据
├── tests/, tools/               科学检查与诊断工具
└── docs/                        物理、接口及发布说明
```

输入日期、单位和跨年参照保留原始说明，见[数据说明](docs/DATA.md)。运行输出、verification历史、检查点、截图和工作缓存不进入源码版本控制。

## 来源与署名

COMPASS原始实现与历史署名予以保留；本分支MEND耦合及长武基线整理记录见源文件中的作者和Change条目。`Open` 描述模块化科学架构，不表示已授予统一开源许可。版权与第三方来源见[NOTICE.md](NOTICE.md)。

## 完整性与冒烟检查

```bash
python tools/verify_release.py
```

```matlab
addpath('tests');
receipt = verifyStandaloneRelease(pwd, [], 4);
```

参数总表位于 `docs/parameters/`。原生输出快速查看：

```matlab
addpath('PlotResults');
plot_profiles('data/output/changwu');
```

并行运行不同任务时，分别设置 `COMPASS_OUTPUT`，避免共享写目录。
运行结果、检查点和日志已由Git忽略；必需NetCDF输入库不会被忽略。
该目录尚未推送到GitHub，也未生成或绑定新的远程仓库。
