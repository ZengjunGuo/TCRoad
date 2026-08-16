# C15–TCRM 来源与参数矩阵

状态：`SOURCE_CONTRACT_V3_FIXED_R0_PUBLIC_RECONSTRUCTION`  
日期：2026-08-13

本文件区分三层证据：论文明确项、官方/公开代码恢复项，以及仍未公开的作者生产细节。正式实现是`method-faithful public reconstruction`，不是`source-identical reproduction`；未公开项保留为主张边界，而不是邮件依赖或自由调参清单。

## 1. 风与降雨共同物理内核

| 项目 | 冻结口径 | 来源层级 | 当前状态 |
|---|---|---|---|
| 完整径向风廓线 | C15：ER11 内核与 E04 外区数值合并 | Chavas et al. 2015；官方代码 DOI `10.4231/CZ4P-D448` | 官方 MATLAB/Python 原件已冻结 |
| C15 输入 | 最大方位平均旋转风、径向尺度（`rmax`、`r0` 或 `(rfit,Vfit)`）及科氏参数；另含 `Cd`、`w_cool`、`Ck/Cd` 等官方参数 | 官方 C15 文档与函数 | 明确 |
| 风害 | C15 完整廓线；按 Lin & Chavas (2012) 加环境风修正；指标为 ten-minute sustained wind | Nature 2023 Methods | 明确到方法引用；实现仍需逐式恢复 |
| 降雨驱动风 | 与风害相同的 C15 用于准备 TCR 风输入 | Nature 2023 Methods | 明确；禁止 C15 风 + ER11 雨混搭 |
| TCR 雨率 | `P = epsilon_p (rho_air/rho_liquid) q_s w`，且仅在总上升速度 `w > 0` 时有雨 | Lu et al. 2018 Eq. 1/14 | 明确 |
| 上升速度组成 | 摩擦辐合、地形、涡旋随时间伸缩、深层切变和辐射冷却五项相加 | Lu et al. 2018 Eqs. 5–13 | 明确 |
| 绝对角动量 | `M=rV+0.5 f r^2` | Lu 2018理论定义 | Lu正文排印为量纲错误的`0.5 f V^2`；实现必须使用标准绝对角动量形式并记录勘误 |
| 数值雨强上限 | `w<=7 m s^-1` | Emanuel官方v6.4；CLIMADA Petals v6.2公开实现 | **公开重建已冻结**；不得声称Nature作者源码已逐字确认 |
| C15 外缘处理 | 论文及官方 C15 未给人为 taper、`r0` clipping 或经验尾部 | 一手论文与官方代码 | **禁止**自创 |
| 切变上升速度 | Lu Eq. 12给完整热力学形式；公开重建使用冻结的Petals v6.2实现 | Lu 2018；Emanuel v6.4；Petals `tc_rainfield.py` | **公开重建已冻结**；Nature未公开生产分支仍不声称源码同一 |
| C15风到TCR风量 | 官方C15给含科氏`V(r,t)`；无科氏伴随风保留同一归一绝对角动量形状：`mu=(rV+0.5fr²)/(rmVmax+0.5frm²)`，`Vd=rmVmax*mu/r`，`Vd(0)=0` | 官方C15角动量定义；Emanuel公开`windprofilem` companion构造 | **public-code-derived bridge已冻结**；不是`Vd=V`，不是ER11，也不声称Nature源码同一 |

## 2. Lu 2018 与 Xi 2020 明示参数

| 参数或设置 | 冻结值/定义 | 来源 | 说明 |
|---|---|---|---|
| `rho_air/rho_liquid` | Xi历史分支显式使用Lu/Xi论文的`0.0012`；Petals/v6.4常数`0.00117`只保留为实现差异记录 | 论文与公开代码 | **公开重建已冻结**；不得静默沿用Petals默认值 |
| precipitation efficiency `epsilon_p` | `0.9` | Lu 2018；Xi 2020 | 明确 |
| radiative cooling velocity `w_r` | `-0.005 m s^-1` | Lu 2018 | Xi 2020 除 `H_b` 和空间变 `Cd` 外沿用 Lu 设置 |
| topographic wind threshold `V_th` | `30 m s^-1` | Lu 2018 | Xi 2020 表述为其余参数沿用 Lu；进入实现前保留引用链 |
| lower-troposphere depth `H_b` | `4 km` (`4000 m`) | Xi 2020 | PDF 印作 `4000 km`，按物理定义和上下文记录为排印错误，禁止使用4000 km |
| 时间步 | 6-hourly轨迹插值到1-hourly；六小时时次先做环境空间平均，再将所得参数序列插值到1-hourly | Xi 2020 historical validation | 明确；不得颠倒空间平均与时间插值顺序 |
| TCR 网格 | `0.05° × 0.05°` | Xi 2020 | 明确；道路0.1°只用于事后叠加 |
| 地形、粗糙度 | CLIMADA package data：SRTM 0.1°地形与由ERA5 `forecast_surface_roughness`生成的0.25°drag | Petals v6.2公开实现与Data API metadata | **公开重建已冻结**；不是Xi作者ECMWF静态场的字节同一副本 |
| drag coefficient | 由 surface roughness 按 Feldmann et al. (2019) / Esau (2004) 计算 | Xi 2020 | `C_D'=[kappa/ln(500/z0)]^2`，`C_D=0.9 C_D'/(1+50 C_D')`，其中 `kappa=0.35`；不得用Lin强度模型常数Cd代替 |
| 物理计算半径 | 每时次仅计算距当时TC中心`<=300 km`的格点 | Xi et al. 2023 JAMC明确公开口径；Petals `max_dist_eye_km=300` | **公开重建已冻结**；移动域，不是固定bbox/事后裁图，也不是Nature的250 km统计距离 |
| 径向差分 | `2 km` | Petals `res_radial_m=2000` | **公开重建已冻结** |
| `dM/dr`下限 | `10 m s^-1` | Emanuel v6.4；Petals `np.fmax(10,dMdr)` | **公开重建已冻结** |
| drag下限 | `0.001` | Petals `min_c_drag=0.001` | **公开重建已冻结** |

### 2.1 公开实现与静态文件身份

| 对象 | 冻结身份 |
|---|---|
| CLIMADA core | `v6.1.0`, commit `bb9595944fcf673152ea53e61efbafeb0b1e4406` |
| CLIMADA Petals | `v6.2.0`, commit `6ecd7af096f126df2da1023fbc5013765566d5e9` |
| `tc_rainfield.py` | 75,745 bytes；SHA-256 `6f0bd30dc5532d907401a862f9d8b560c3feca6356c79ac7de444a52b315e062` |
| drag | <https://data.iac.ethz.ch/climada/cde33390-d2a2-4032-a26b-26ab285bcec0/c_drag_500.tif>；999,184 bytes；SHA-256 `1c3f3b525f0c2a9e73f6fe6d3ba3caf7436699f8caeb3edb69548e04fd3f4a42` |
| topography | <https://data.iac.ethz.ch/climada/bc5b8fc7-0e73-49d9-a73a-3fad22fdfff5/topography_land_360as.tif>；3,054,640 bytes；SHA-256 `de8142fe9f50d0cfbd944884ee945bb355b09dfce2f214879e608f87ae0f0951` |

## 3. 必须分开的两套环境输入

### A. Xi 2020 历史事件验证

| 输入 | 论文口径 | 数据源 | 当前服务器状态 |
|---|---|---|---|
| 湿度 | 六小时时次对`q925`做TC中心200 km圆盘平均，再将参数序列插值至1小时 | NCEP再分析 | **已冻结并执行**；37个六小时参数已提取并生成217个逐小时节点 |
| 深层切变 | 六小时时次对`(u200-u850, v200-v850)`做TC中心600–800 km环带平均，再将参数序列插值至1小时 | NCEP再分析 | **已冻结并执行**；37个六小时参数已提取并生成217个逐小时节点 |
| 轨迹/强度/RMW | IBTrACS 6-hourly→1-hourly；末端4个RMW缺测严格用Knaff & Zehr (2007) Eq. 6填补并保留estimated flag | IBTrACS；DOI `10.1175/WAF965.1` | **已冻结Irene原始文件**；不得用CLIMADA pressure-fit RMW估计器 |
| 地形/drag | CLIMADA `topography_land_360as` + `c_drag_500` | CLIMADA Data API公开package data | URL、字节数和SHA已冻结；**方法级公开替代，非Xi字节同一** |

### B. Nature 2023 合成事件

| 输入 | 已发表口径 | 状态 |
|---|---|---|
| 湿度 | Emanuel官方v6.4 `qs900b.m` 给出精确暖核实现，但实际目标为950 hPa：环境`T600`、圆周最大风（kt→m s⁻¹）、五次Newton迭代和`0.016V²`暖核项 | **公开候选已闭合**；Nature是否使用该版本仍是主张边界，不是邮件依赖 |
| 风切变 | Emanuel官方v6.4以`5×(storm translation−850-hPa flow)`构造，纬向分量扣除带半球符号的beta drift；不是直接200–850 hPa差 | **公开候选已闭合**；Nature是否使用该proxy仍是主张边界，不得与Irene历史切变混用 |
| 事件级outer radius与逐时RMW | 每个Lin事件从Chavas et al. (2016)全球经验lognormal分布抽取一次固定`r0`；公开重建以Table 1的`median/Q1/Q3=881.0/740.7/1054.4 km`冻结`mu=6.78105762593618`、`sigma=0.261776756893889`，用NumPy `PCG64`、seed `20260810`按`event_position=0..9999`各抽一次；逐时以`v_trks+r0+abs(f)`调用官方C15 `r0input`，由C15内生求RMW | **source-constrained公开闭合已冻结**；架构来自Gori et al. (2022)并由Xi, Lin & Gori (2023)沿用；不声称复现Gori私有参数或随机序列 |
| 时间与网格 | TCR setup follows Xi 2020：1-hourly、0.05° | 明确 |
| C15风输入 | 同一C15驱动风害与TCR | 明确 |

变量职责固定为：Lin的圆周/方位平均最大风`v_trks`进入C15 `r0input`风幅度与Emanuel `qs900b.m`实际950-hPa暖核湿度；`vmax_trks`只保留事件阈值和目录语义，不参与尺度闭合。正式目录不截断、不clip、不重抽，也不按盆地或强度改变draw；生产前逐事件核验`r0>302 km`，违反者标为`METHOD_DOMAIN_PENDING`而不改样本或权重。逐时Knaff/Nederhoff、项目outer-first模型及残差分配均不得进入Lin生产。

公开重建以Petals v6.2的径向/时间差分和静态场采样为单一数值骨架，冻结`deltar=2 km`、`dM/dr>=10`、`Cd>=0.001`和`w<=7`。C15 provider机械追加官方定义的精确端点`(r0,V=0)`；`r>r0`直接报错，不生成尾部、不taper。上述选择全部写入manifest，并明确不等同于Nature未公开生产源码。

## 4. 当前可复用与不可复用

### 可复用

- 100,000条Lin历史轨迹、事件权重和事件身份；
- 10,000条固定hazard sample；
- 轨迹的1小时位置与强度序列；
- OSM五级道路和0.1°道路长度格网；
- 调度、I/O、哈希和manifest外壳；
- 官方C15原始归档及其MATLAB参考输出。

### 不得进入正式物理结果

- 旧 `shared_cle_profile_provider` 和patched pyTCR；
- 固定 `q=0.01`；
- 旧runner中无版本/provenance的`translation-u850`实现；正式合成分支只能调用并记录冻结的Emanuel公开构造；
- 将 `u250` 改名为 `u200`；
- pyTCR内嵌但无来源元数据的地形/drag直接冒充Xi数据；
- 项目拟合的outer-first风暴尺度模型及事件尺度残差分配；二者仅保留审计归档，不进入首事件或后续10,000条正式hazard production；
- 自行删除或调整已冻结的`w<=7 m s^-1`；把移动300 km计算域误作固定bbox或事后裁图；手工外缘平滑；
- 旧首事件的2929.97 mm、265.17 mm h^-1及其道路暴露。

## 5. 下一步执行边界

1. **已完成：**官方MATLAB C15默认例与Python 3兼容层数值闭合。
2. **已完成公开输入：**Irene的IBTrACS和NCEP六小时动态输入、两个CLIMADA静态GeoTIFF及其SHA均已冻结。
3. **已完成公开实现冻结：**core v6.1.0 + Petals v6.2.0 + 官方C15与无科氏伴随风桥，参数详见[C15–CLIMADA公开重建合同](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)。
4. **Irene全生命周期公开重建与Stage IV共同支持对照已完成。**当前冻结10,000-event fixed-`r0`目录，并以官方C15 `r0input`重算首个合成事件；其环境输入不得复用Irene NCEP序列。fixed-`r0`首事件闭合并记录Irene偏湿解释边界后，才生产既定10,000条样本。

## 6. 一手来源

- Xi, Lin & Gori (2023), *Nature Climate Change*: <https://doi.org/10.1038/s41558-023-01595-7>
- Gori et al. (2022), event-fixed outer radius and C15 storm sets: <https://doi.org/10.1038/s41558-021-01272-7>
- Chavas et al. (2016), global observed outer-radius distribution: <https://doi.org/10.1175/JCLI-D-15-0731.1>
- Xi, Lin & Smith (2020), *Journal of Hydrometeorology*: <https://doi.org/10.1175/JHM-D-20-0035.1>
- Lu et al. (2018), *Journal of the Atmospheric Sciences*: <https://doi.org/10.1175/JAS-D-17-0264.1>
- Chavas, Lin & Emanuel (2015), *Journal of the Atmospheric Sciences*: <https://doi.org/10.1175/JAS-D-15-0014.1>
- Official C15 archive: <https://doi.org/10.4231/CZ4P-D448>
- Emanuel (2017), *PNAS*: <https://doi.org/10.1073/pnas.1716222114>
- Xi & Lin (2022), *Geophysical Research Letters*: <https://doi.org/10.1029/2022GL099196>
- Xi et al. (2023), 300 km TCR application domain: <https://doi.org/10.1175/JAMC-D-22-0131.1>
- Feldmann et al. (2019), TCR drag-coefficient update: <https://texmex.mit.edu/pub/emanuel/PAPERS/Feldmann_etal_2019.pdf>
- Esau (2004), surface-drag parameterization: <https://doi.org/10.5194/angeo-22-3353-2004>
- CLIMADA core v6.1.0: <https://github.com/CLIMADA-project/climada_python/tree/v6.1.0>
- CLIMADA Petals v6.2.0: <https://github.com/CLIMADA-project/climada_petals/tree/v6.2.0>
