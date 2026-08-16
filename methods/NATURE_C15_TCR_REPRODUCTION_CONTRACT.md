# Nature Climate Change C15–TCR 复现合同

状态：`FROZEN_PUBLIC_RECONSTRUCTION_V1`  
日期：2026-08-13

## 1. 唯一主线

本研究的风—雨灾害物理层以 Xi, Lin & Gori (2023), *Nature Climate Change* 的方法链为唯一主线：

> 同一个 C15 完整热带气旋风廓线用于风害模拟，并为 TCR 降雨模型提供所需风输入。

这意味着正式生产不采用“C15 风场 + ER11 降雨”的混合内核，也不继续修补当前自定义的 CLE–pyTCR 接口。论文中的物理方法与本文的全球道路风险问题相结合，但风廓线、降雨方程、环境输入和数值设置不得因工程方便而自行更改。

主文献链：

1. Xi, Lin & Gori (2023), *Increasing sequential tropical cyclone hazards along the US East and Gulf coasts*, **Nature Climate Change**, 13, 258–265. DOI: https://doi.org/10.1038/s41558-023-01595-7
2. Xi, Lin & Smith (2020), *Evaluation of a Physics-Based Tropical Cyclone Rainfall Model for Risk Assessment*, **Journal of Hydrometeorology**, 21, 2197–2218. DOI: https://doi.org/10.1175/JHM-D-20-0035.1
3. Lu et al. (2018), TCR 详细方程与风廓线敏感性。DOI: https://doi.org/10.1175/JAS-D-17-0264.1
4. Chavas, Lin & Emanuel (2015), C15 完整径向风廓线。DOI: https://doi.org/10.1175/JAS-D-15-0014.1
5. C15 官方代码档案（v1.0, CC0-1.0）：https://doi.org/10.4231/CZ4P-D448
6. Gori et al. (2022), Nature 研究采用的合成风暴尺度与 C15 链。DOI: https://doi.org/10.1038/s41558-021-01272-7
7. Chavas et al. (2016), 全球经验 outer-radius 分布。DOI: https://doi.org/10.1175/JCLI-D-15-0731.1

## 2. 一手来源已经明确的合同

### 2.1 Nature Climate Change 2023 明确规定

- 风害使用 C15 完整风廓线模型。
- 地表风叠加环境风修正，遵循 Lin & Chavas (2012)。
- 风害量为 10-minute sustained wind。
- 同一个 C15 模型用于准备 TCR 所需的风输入。
- 降雨使用 TCR；详细方程追溯至 Lu et al. (2018)，模拟设置追溯至 Xi et al. (2020)，环境场处理追溯至 Emanuel (2017)。

### 2.2 Xi et al. 2020 的历史事件验证设置

- 轨迹、强度和RMW由6-hourly线性插值到1-hourly；环境量先在每个六小时时次按TC中心做圆盘/环带平均，再将所得参数序列线性插值到1-hourly。
- TCR 水平网格为 `0.05° × 0.05°`。
- 风廓线采用 C15，而不是 pyTCR 1.2.3 的默认 ER11。
- 925-hPa specific humidity：以 TC 中心 200 km 内的空间平均值作为环境湿度输入。
- 200–850-hPa deep-layer wind shear：以 TC 中心 600–800 km 环带平均值作为环境切变输入。
- lower-troposphere height 为 `4000 m`。
- drag coefficient 由 surface roughness 推导，并遵循论文引用的方法。

这里的湿度和切变来自 NCEP 再分析，是 Xi et al. (2020) 对 1999–2018 年北大西洋历史事件进行模型验证时采用的输入链，不能直接当作 Nature 2023 合成事件的环境构造。

Xi et al. (2020) 正文印作 lower-troposphere height `4000 km rather than 1000 km`。根据 Lu et al. (2018) 对该量的定义、原设定 `1 km` 及文中“更符合理论”的表述，这显然是长度单位排印错误；合同记录其物理意图为 `4000 m (4 km)`，同时保留原文印刷事实，不把 `4000 km` 写进程序。

### 2.3 Nature 2023 合成事件的环境场与公开实现选择

- Nature 2023 明确将 TCR 详细方程追溯到 Lu et al. (2018)，数值设置追溯到 Xi et al. (2020)，但将合成事件所需环境场单独追溯到 Emanuel (2017)。
- 后续 Xi & Lin (2022) 对这条合成事件链的公开说明是：风暴中心近 `900 hPa` 饱和湿度由环境 `600 hPa` 温度及与 TC 强度相关的暖核温度构造。
- Kerry Emanuel官方MIT目录中的公开MATLAB包`v6.4`现已冻结。其`qs900b.m`给出可执行暖核公式，但变量名虽为`q900`，实际压力为`950 hPa`；`raingen.m`把合成切变构造为五倍的风暴平移相对850-hPa环境风，并含beta-drift修正；`pointwshortnq*.m`明确执行`w=min(w,7)`。
- 这套代码是Emanuel官方公开实现，但Nature 2023没有声明它就是生产版本，且其默认风廓线为ER11而非C15。本项目因此不声称源码同一，而采用[公开可审计重建合同](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)：官方C15生成唯一风廓线，CLIMADA core `v6.1.0` + Petals `v6.2.0`提供公开TCR数值骨架。
- 公开实现冻结为1 h、0.05°、`epsilon_p=0.9`、`H=4000 m`、随TC中心移动的300 km计算域、2 km径向差分、`w<=7 m s^-1`、`dM/dr>=10 m s^-1`和`Cd>=0.001`。Xi历史分支按Lu/Xi论文使用密度比`0.0012`，不沿用Petals常数`0.00117`。
- 合成事件仍不得套用Xi 2020历史NCEP链，也不得使用固定`q=0.01`；其湿度和切变必须来自该合成事件自己的CMIP6/Emanuel公开环境构造并写入manifest。
- Lin v1.1轨迹不含outer radius或逐时RMW。Gori et al. (2022)公开的方法是：每个TC从Chavas et al. (2016)经验lognormal分布抽取一次`r0`，整个生命周期保持不变；逐时再以轴对称强度、固定`r0`和科氏参数调用C15 `r0input`，由C15内生计算RMW。Nature 2023沿用了Gori的storm sets。因此Lin正式生产统一采用这条fixed-`r0`架构，撤出逐时Knaff/Nederhoff RMW。
- Gori未公开其精确lognormal参数、随机种子或抽样代码。本项目以Chavas et al. (2016)全球样本Table 1公开的`median/Q1/Q3=881.0/740.7/1054.4 km`冻结分位数匹配lognormal：`mu=6.78105762593618`、`sigma=0.261776756893889`；NumPy `PCG64`以seed `20260810`按冻结sample的`event_position=0..9999`各抽一次，生成后只读取不可变目录。详细合同见[Chavas 2016 fixed-r0公开重建合同](CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)。
- C15风幅度和`qs900b.m`实际950-hPa暖核湿度均使用Lin圆周/方位平均最大风`v_trks`；`vmax_trks`保留事件阈值/目录语义，不再作为经验RMW预测量。该选择是source-constrained public reconstruction，不是Gori私有随机序列的字节同一复刻；Knaff、项目outer-first尺度模型和残差分配只保留审计归档，不进入正式hazard production。

正式实现的逐项状态见 [C15–TCRM 来源与参数矩阵](C15_TCRM_SOURCE_PARAMETER_MATRIX.md)。冻结配置只接纳论文、所引原始方法、官方C15或公开TCR实现给出的内容；未公开的作者生产细节保留为主张边界，不再作为发送邮件或启动Irene的前置门槛。

## 3. 300 km计算域与明确撤回的错误口径

- Nature 论文中的 `250 km` 是逐点定义 hazard-producing day/event 的统计筛选距离，不是 TCR 的计算半径。
- Nature 正文和补充材料没有单独写出300 km，但Xi等人后续公开TCR应用明确说明：只模拟距TC中心300 km内的降雨（Xi et al. 2023, JAMC, DOI `10.1175/JAMC-D-22-0131.1`）；CLIMADA Petals公开实现同样将`max_dist_eye_km=300`设为默认计算域。因此本公开重建冻结**逐时随风暴中心移动的300 km计算域**。
- 300 km不是固定地理框、不是事后裁图，也不授权对C15外缘做taper；此前否定其公开来源的冲突口径废止。
- 道路栅格为 `0.1°` 不代表风雨物理计算也可直接使用 `0.1°`。正式顺序应为：

  `0.05° C15–TCR 计算 → 事件累计/极值 → 聚合或叠加到 0.1° 道路格网`

## 4. 禁止进入正式方法的自定义内容

以下内容全部退出科学生产链，不得被描述为 Nature/Xi 方法：

- pyTCR 1.2.3 默认 ER11 取代 C15；
- “C15 风 + ER11 雨”的混合内核；
- 第三方 `tcwindfields 0.2.0` 冒充官方 C15 原码；
- 自定义 `shared_cle_profile_provider.py`；
- `Vd = V`；
- C15 外缘手工置零、cosine taper、r0 clipping 或任意经验平滑；
- 为压低异常雨量而调整已冻结的`w<=7`、关闭时间导数或修改积分公式；
- 固定湿度常数替代逐时环境湿度；
- 以 250-hPa 风未经证据地替代论文规定的 200-hPa 风；
- 项目拟合的outer-first风暴尺度模型或事件尺度残差分配进入Lin正式风雨场；
- 逐时Knaff或Nederhoff RMW替代事件级fixed-`r0`，或在它们之间拼接/fallback；
- 直接在 0.1° 道路网格上运行 TCR；
- 把当前首事件 smoke/footprint 当作科学验证。

## 5. 当前代码与数据的处置

### 可继续使用

- Lin 100,000 条历史合成轨迹、事件权重和事件身份；
- 10,000 条固定 hazard sample；
- OSM 五级道路、0.1°道路长度格网和局部道路网络；
- 事件调度、NetCDF I/O、哈希、manifest 和 provenance 外壳；
- 原始环境场读取框架，但变量定义和适配器必须按 Nature/Xi 合同重写。

### 仅作失效诊断归档

- 当前 patched shared-CLE–pyTCR runtime；
- 当前自定义 provider；
- 首事件 `event_00000_complete_footprint.nc`；
- 由其得到的 `2929.97 mm` 累计雨量、`265.17 mm h-1` 雨强和对应道路暴露。

这些数值证明了旧接口的外缘时间导数问题，不是可用于论文的灾害结果。

## 6. 公开代码边界

- 官方 C15 代码档案已完成版本、许可、附件与函数级核对：PURR `v1.0`（2022-06-06，CC0-1.0）归档了 MATLAB bundle 和 Chia-Ying Lee 提供的 Python 版本。PURR 当前下载端点在本机存在 TLS 连接问题，但 PURR 官方落地页注明，该归档就是 Dan Chavas 实验室自 2020-06-23 公开的 primary version；因此已从两条 Purdue 官方原始附件 URL 的 Web Archive 保全副本取得完全对应的 2020-06-23 文件。
- `CLE15_windprofile_PUBLIC_2020-06-23.zip`：`1,455,158 bytes`；SHA-256 `5ebb6e7c253e927653b4515278ace60cb43d17504c720e2b6bc03ab2f1519dc2`。
- `CLE15_2020-06-23.py`：`22,762 bytes`；SHA-256 `6f1306fae71d0e772f17dbf67d5c8cfd94fa543dd122cbb390c6d50161325113`。
- MATLAB bundle 包含三种官方入口：已知 `r0`、已知 `rmax`、已知 `(rfit,Vfit)`。正式接入优先把原始 MATLAB bundle 作为数值基准；公开 Python 文件仍是 Python 2 语法，且 eye-adjustment 段被注释，不能未经验证地改写后冒充原始官方实现。
- 原始 ZIP 内 2015/2018 README 的“不得转发”文字与 2022 PURR 正式发布的 CC0-1.0 权利元数据存在历史冲突。项目 provenance 必须同时保留两者，不删除、不改写；正式许可声明以 DOI 存档元数据为准。
- Nature 论文的 Zenodo 仓库只公开物理风雨结果之后的概率模型、分析和绘图代码，不含完整 C15–TCRM 物理生成引擎：https://doi.org/10.5281/zenodo.7407013
- 尚未发现 Xi et al. (2020) 的完整可执行 C15–TCRM 源码仓库。
- Emanuel官方MIT目录的`scripts_ver6.4`已经逐字节冻结并审计。它公开了暖核湿度、合成切变、时间插值和`w<=7`等关键实现，但默认使用ER11，且Nature论文没有把该包声明为生产代码。原件与审计见`vendor/Emanuel_TCR/scripts_ver6.4/`和`methods/EMANUEL_SYNTHETIC_TCR_INPUT_CODE_AUDIT.md`。
- 本公开重建冻结CLIMADA core `v6.1.0`（commit `bb9595944fcf673152ea53e61efbafeb0b1e4406`）与Petals `v6.2.0`（commit `6ecd7af096f126df2da1023fbc5013765566d5e9`）。`tc_rainfield.py` SHA-256为`6f0bd30dc5532d907401a862f9d8b560c3feca6356c79ac7de444a52b315e062`。
- C15的含科氏风通过归一绝对角动量形状接入Petals；其无科氏伴随风使用公开代码派生关系，详见[C15–CLIMADA公开重建合同](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)。不使用`Vd=V`，也不另跑ER11或增加taper。
- 静态场使用CLIMADA正式package data：`c_drag_500.tif` SHA-256 `1c3f3b525f0c2a9e73f6fe6d3ba3caf7436699f8caeb3edb69548e04fd3f4a42`；`topography_land_360as.tif` SHA-256 `de8142fe9f50d0cfbd944884ee945bb355b09dfce2f214879e608f87ae0f0951`。它们可审计，但不是Xi作者静态场的字节同一副本。

因此，本项目执行的是 **method-faithful public reconstruction**，不能声称 **source-identical reproduction**。公开材料没有说明的作者生产细节不得被改写成“已确认一致”；本项目使用的替代公开实现、文件身份和差异必须在manifest与论文方法中直接披露。

## 7. 唯一执行顺序

1. **已完成：**取得并冻结 PURR 官方 C15 附件，记录版本、许可、逐文件 SHA256 和函数接口。
2. **已完成合同层：**从 Lu et al. (2018)、Xi et al. (2020)、Emanuel (2017) 与 Nature 2023 建立逐公式、逐参数、逐单位的来源矩阵，明确区分历史验证和合成事件。
3. **已完成默认例闭合：**以官方 MATLAB C15 bundle 为数值基准建立现代兼容层；官方 `r0/rmax/rfit` 三个默认示例共38个PDF曲线锚点的最大误差为`0.001609 m s⁻¹`，未改物理方程。该结论不扩张为全参数域等价，`eye_adj=1`仍禁止使用。
4. **已完成公开输入冻结：**Xi 2020 Irene案例所需的IBTrACS和NCEP 6-hourly `q925/u200/v200/u850/v850` 已逐文件保存、校验并记录SHA；两个CLIMADA官方静态GeoTIFF也已固定URL、字节数和SHA。
5. **已完成数值合同冻结：**以官方C15 + CLIMADA core/Petals为唯一实现，冻结1 h、0.05°、移动300 km、2 km差分、`epsilon_p=0.9`、`H=4 km`、`w<=7`、`dM/dr>=10`、`Cd>=0.001`及Xi历史密度比`0.0012`。
6. **已完成：**运行Irene全生命周期公开方法级重建，并与Stage IV完成216小时、陆地、逐时移动300 km共同支持对照；空间相关为0.7156，公开重建表现出明确空间信号和绝对量偏湿。Xi 2020 Table 1印刷时段与Fig.1图域/实际轨迹存在内部冲突，故不把该结果冒充Fig.1逐像元复刻。
7. **已完成：**冻结10,000-event fixed-`r0`不可变目录，把Lin生产接口从Knaff `rmaxinput`整体切换为官方C15 `r0input`；已用同一冻结内核重算当前首个合成事件。旧自定义接口的moving-`r0`外缘伪降雨带未再出现。合成环境输入未复用Irene NCEP序列。
8. **现在执行：**按 Gori et al. (2022) 对分析集合全算物理场。本项目的分析集合是约1000个模拟年里进入道路300 km的全部Lin事件；随机1万条只是方法试生产，不是终局样本。
9. 在近路全量历史hazard闭合前，不扩展未来情景、SFINCS、交通流或级联损失支线。

## 8. 完成标准

只有满足以下条件，才可以写“完成风雨场反演”：

- 风与雨由同一 C15 物理廓线驱动；
- 公式、单位、参数、软件提交和环境场均可追溯到一手论文或冻结的公开实现；
- 没有未披露的自定义 taper、截断、经验替代或混合内核；
- Irene历史案例完成公开方法级重建，并明确不宣称Xi作者静态输入/源码的字节同一；
- 当前合成首事件使用同一冻结实现成功生成可解释的风雨场；
- Lin合成事件的outer radius在事件内固定、逐时RMW由官方C15 `r0input`内生求得，且没有Knaff/Nederhoff拼接或fallback；
- 产物保存完整 provenance，且可区分论文原设定与本研究的全球应用层。
