# 长武冬小麦叶片性状与光合证据检索

整理：Wenjie Liu，2026-09-30。范围为叶面积与叶干重的配对、叶 C/N、光合容量与叶角/光学性状；本次仅整理文献和来源，不修改模型、站点参数或观测。

用户随后明确本轮重点为水、碳、氮，不再审查SLA/Cdm。本报告保留为已完成的来源备查；下文叶性状缺口不再作为本轮待完成任务或模型验证前提。光合研究可作为碳交换的辅助证据，当前工作次序见[水碳氮检索清单](changwu_wheat_literature_data_inventory_20260930.md)。

**用户已有 Excel 实测 LAI，LAI 不列为缺失数据。** 这里寻找的叶面积是与同一叶片或同一批样品的干重、C/N 相配对的面积，用来检查 SLA、Cdm 和 Narea；其他站点的 LAI 不替换现有 Excel。

本次已取得两篇真正长武田间研究的期刊原文，可约束旗叶净光合速率、气孔导度和平均叶倾角；仍未取得长武同批叶片的面积＋干重＋碳分数，或可独立拟合 Vcmax/Jmax 的原始 A–Ci＋叶温/荧光记录。论文中出现“长武134”品种名并不证明试验站在长武；以下站点归属按原文 Methods 核定。

## 长武原文：可以使用的证据

| 文献与原文位置 | 试验、样品与条件 | 可用变量及数值 | 应用边界 |
|---|---|---|---|
| [Chen & Hao 2015, Photosynthetica 53:519–526](https://doi.org/10.1007/s11099-015-0147-9)，Methods pp.520–521；Results p.521；Figs.1–2 | 长武农生态站；2007/08、2008/09；8 品种，包括长武131/134和长旱58；随机区组3重复，N150、P₂O₅75 kg/ha。旗叶在拔节、开花、灌浆测量；09:30–11:00；CID-340；PPFD1000–1200 μmol photons m⁻² s⁻¹、RH60–80%；叶室温度未控制 | 拔节旗叶 Pn 的明确文字数值为乌枝9.19、长旱58 13.35 μmol CO₂ m⁻² leaf s⁻¹。gs、E 和另外两阶段 Pn 有图；逐点原始表未取得。原文称约9月25日播种、约6月30日收获 | 可以验证旗叶气体交换的量级及季节/品种差异；Pn 不是 Vcmax，旗叶 Pn 不是地面面积上的冠层 GPP。原文把环境 CO₂“浓度”360–400印成 μmol m⁻² s⁻¹，单位有误，保留疑点，不据此制造精确 CO₂ 条件或 Vcmax25 |
| [孙婴婴等2016，作物杂志(2):32–38](https://zwzz.chinacrops.org/CN/10.16035/j.issn.1001-7283.2016.02.006)，Methods p.33，Fig.3B p.36 | 真正长武站；2011/12；7 品种，包括小偃6、长武134、长旱58；100/250/350粒 m⁻²，3重复；播前N150、P₂O₅120 kg/ha，无追肥、无灌溉。花后3天旗叶Pn，LI-6400，09:30–11:30，每小区5次；LAI-2200在无直射光时测冠层性状 | 图3B提供平均叶倾角 MTA（°）；图中均值约43–57°，为目视判读而非原始逐点数据，未收入校准数值数组。Pn也有图。官方PDF已保存 | 可作为长武叶角先验的图形检验和后续数字化来源。MTA为冠层仪器估计，不能直接作为LIDFa或LIDFb；一个平均角不能唯一决定完整叶角分布。这里不提取其LAI来替换用户Excel |
| [Chen, Zhang & Deng 2019, The Crop Journal 7:89–100](https://doi.org/10.1016/j.cj.2018.10.004)，Methods p.90，Table2 p.95 | 长武站，长旱58；2011/12和2012/13，雨养；N120/N150/N195/N195＋覆盖处理。旗叶面积由叶长×叶宽×0.83计算 | 表2跨处理平均旗叶面积：2011/12为20.03 cm²，2012/13为11.03 cm²；这是单叶面积，不是LAI。收获2012-07-02、2013-06-26 | 可核对单叶尺度和器官面积；没有配对叶干重/%C，因此不能据此计算SLA或Cdm。两季收获日期是管理事实，不是本站模型所有年份必须固定的日期 |

第一篇的[期刊PDF](https://ps.ueb.cas.cz/pdfs/phs/2015/04/05.pdf)、第二篇的[期刊下载入口](http://zwzz.chinacrops.org/CN/article/downloadArticleFile.do?attachType=PDF&id=5409)、第三篇的[作者机构库出版版PDF](https://researchportal.murdoch.edu.au/view/pdfCoverPage?download=true&filePid=13136829620007891&instCode=61MUN_INST)已下载。原始PDF、文本及用于核读表/图的页面PNG位于 `outputs/changwu_wheat_literature_20260930/leaf_sources/`；JSON记录URL、字节数和SHA256。

叶氮还有一篇需要继续取得原文的候选：[Chen et al. 2019, Agronomy Journal 111:1218–1228，10.2134/agronj2018.04.0243](https://doi.org/10.2134/agronj2018.04.0243)。[作者所在UWA机构库](https://research-repository.uwa.edu.au/en/publications/nitrogen-vertical-distribution-differed-in-foliar-and-nonfoliar-o/)和Crossref核对了题名、作者、卷页及DOI；摘要说明它测了灌浆期叶与非叶器官的氮浓度和再转运。此次未获取可核读的出版原文，机构摘要本身不列试验站，故本报告不把搜索预览中的长武站点/处理数字作为已核结论，也不填入具体叶N浓度或Narea。

## SLA/Cdm缺口与可以落实的核对

当前合同中 `SLA0=0.035 m² leaf/g C`，`Cdm=0.012 g DW/cm² leaf`。如果它们描述同一批叶片，则

\[
\mathrm{LMA}=10^4 Cdm=120\ \mathrm{g\ DW/m^2},\qquad
f_C=\frac{1}{SLA_C\,(10^4 Cdm)}=0.238095\ \mathrm{g\ C/g\ DW}.
\]

这个组合需要实测叶C分数核对，不能用已有LAI单独判断两个参数谁应修改。[已核119条独立冬小麦配对记录，4TU DOI10.4121/21304917.v1](https://doi.org/10.4121/21304917.v1)给出的 `f_C=0.4198–0.4738`、`LMA=40.24–94.87 g DW/m²`、`SLA_C=0.02320–0.05796 m²/g C`、`Narea=1.076–3.033 g N/m² leaf`，是 Wageningen 2019/20控制环境/玻璃室的Thésée与Récital，**不是长武标定值**。完整原始工作簿、下载MD5核对及逐行换算保留在[既有叶氮与参数审计](winter_wheat_leaf_n_and_parameter_audit_20260930.md)。该样本内Cdm范围约0.004024–0.009487，不能把这个单一实验的上界当作所有冬小麦的生物学上界。

优先取得同一叶片/同一批样品的面积、恒重干重与%C；再增加%N、鲜重和叶龄/生育期。使用明确单位计算：

| 原始配对测量 | 推导量 | 用途 |
|---|---|---|
| 叶面积A（m²）、叶干重DW（g） | LMA=DW/A；Cdm=LMA/10000；SLA_DW=A/DW | 可直接检查光学干质量参数，但SLA_DW不是当前SLA_C |
| 同样品%C转换成fC（g C/g DW） | SLA_C=A/(DW×fC) | 检查当前按碳质量的SLA，不把干质量SLA直接代入 |
| 同样品%N转换成fN（g N/g DW） | Narea=DW×fN/A | 约束LUNA实际叶N输入；叶质量N浓度不能直接当Narea |
| 同样品鲜重FW（g） | Cw=(FW−DW)/A/10000（cm，水密度1 g/cm³） | 约束光学等效水厚度；必须同批、同面积且排除表面自由水 |

[CERN元数据16969](https://rs.cern.ac.cn/data/meta?id=16969)列有长武站及作物叶面积、叶干重、茎干重等指标，是取得配对干重的候选入口；**此次未下载原始记录，尚不清楚叶面积与干重能否按同一采样日期/地块/样品配对，也未证明提供%C或%N**。它解决的候选缺口是叶干重/生物量，不是用户LAI缺失。仅有这份元数据不等于已经完成SLA_C/Cdm验证。

## 明确为外站的原文与方法参考

| 原文 | 核定站点与可用内容 | 保留的限制 |
|---|---|---|
| [陈骎骎等2015，小偃60在不同水氮条件下的用水特性，10.13930/j.cnki.cjea.141429](https://www.ecoagri.ac.cn/cn/article/pdf/preview/10.13930/j.cnki.cjea.141429.pdf)，Methods pp.1254–1255，Tables4–5 p.1257 | **河北南皮，2013/14盆栽，小偃60**；50/65/80%田间持水量与0/0.10/0.20/0.29 g N/kg土壤组合。Table5 Pn均值8.34–13.80 μmol m⁻² s⁻¹，端点分别W0N0、W2N3；SPAD、叶水势另有表 | 非长武；SPAD为相对读数，未提供SPAD–Cab配对标定；不把SPAD转成Cab、Vcmax，也不把盆栽g N/kg换成站点施肥kg/ha |
| [贺鹏等2023，秸秆覆盖与施磷对西南冬小麦叶肉导度和CO₂同化效率的影响，10.20103/j.stxb.202301090057](https://www.sciengine.com/doi/pdf/8DFDAA81EE8B4D95BCC34346F4F80590)，Methods pp.9676–9677，Table2 p.9678，Figs.4–5 | **四川仁寿，2020/21、2021/22，川麦104，雨养**；有面积＋干重、LI-6800 A–Ci、gs、荧光变量J法gm及报告的Vcmax，适合作为验证方法参考 | 原文Table2报告LMA226–246 g/m²；Methods中LMA公式将LA写为cm²却直接标LMA为g/m²，PA公式的mg/g换算也有系数疑点。CO₂/补偿点单位亦有文字问题。保留原值与疑点，不擅自换单位或作为长武新默认；未取得原始A–Ci文件，未把图中容量当Vcmax25 |
| [武荣等2018，覆膜与补灌对冬小麦叶绿素和比叶面积的影响，节水灌溉(9):27–30](https://irrigate.whu.edu.cn/jsgg/CN/Y2018/V0/I9/27)，Methods p.28，Table1 p.29 | **陕西杨凌，2016/17，小偃22**；浸提测Chl a/b、叶面积与恒重干重，有SLA表 | 不是长武；表1的SLA单位确印cm²/g，数值12.82–52.18，与文中方法和常用叶片范围需再核实。文本PW2灌浆写12.92而表1写12.82。叶绿素仅标mg/g、鲜/干基准未明确，不直接转换为模型μg/cm² Cab；不猜测其表应为m²/kg |

IEEE论文 [Zhang et al., Retrieving Leaf and Canopy Water Content of Winter Wheat Using Vegetation Water Indices，10.1109/JSTARS.2017.2773625](https://doi.org/10.1109/JSTARS.2017.2773625) 的Crossref题名/作者/DOI已核。此次未取得可核读的出版方或作者机构原文；搜索到的非官方预览不作Cdm/Cw实测证据。未用其预览值填补本站Cw。

## 验证路线与尚未取得的量

1. 继续保留用户Excel实测LAI。消费该观测的数据模式可验证采样、插值和收获事件，不能用同一输入证明LAI预测独立准确；LAI预测误差需要预测模式或留出观测。另完成同批叶面积、干重和%C配对，核对SLA_C与Cdm的数学一致性；在原始样品缺失前保留当前参数及待核状态，不给出新默认。
2. 用已核长武旗叶Pn/gs作独立田间量级和阶段响应核对。模型需要输出与旗叶PPFD、温度、CO₂、湿度、单位和阶段一致的叶尺度量；其他年份/品种且缺叶温时，只能作为范围与趋势核对，不能声称已拟合当前2000–2001年本站冠层GPP或算出严格RMSE。
3. 如数字化孙等Fig.3B，应保留品种×密度、原图坐标、误差条、数字化误差及“图提取”标签；以预测平均倾角检查叶角分布。LIDFa/b的联合识别还需要叶倾角直方图或多角度冠层辐射/反射信息。
4. Vcmax/Jmax需要A–Ci原始点及叶温、光强、气压、湿度、CO₂控制条件；gm修正需要同步荧光/电子传递与Rd/补偿点方法。拟合时报告采用Ci还是Cc和是否归一到25℃。现有Pn端点、SPAD或N浓度本身不完成这项验证。

本次尚未取得**长武同批叶C分数与LMA、叶Narea、面积基叶绿素、Vcmax/Jmax原始曲线、gm原始记录、叶角分布直方图、鲜干重配对Cw与叶光谱**。这表示本次未获取可验证数据，不表示这些测量或论文不存在。已找到的本地/外站证据和访问缺口均列入 `outputs/changwu_wheat_literature_20260930/leaf_evidence.json`。

检索使用中英文“长武/Changwu＋冬小麦/wheat＋比叶面积/leaf dry weight/leaf nitrogen/Vcmax/Jmax/mesophyll conductance/叶绿素/叶角/水氮/喷灌”，先核期刊原文Methods与图表，再以Crossref/作者机构库核对题录。中文期刊DOI的Crossref404仅记录未收录，不能据此否定已取得的期刊原文。PDF表头与关键方法页已经渲染核读；本次没有开展模型拟合，也没有改写原始观测。
