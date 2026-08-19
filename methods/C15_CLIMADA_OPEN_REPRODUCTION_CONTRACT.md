# C15–CLIMADA 公开可审计重建合同

状态：`FROZEN_PUBLIC_RECONSTRUCTION_V2_FIXED_R0`  
日期：2026-08-13

## 1. 结论与主张边界

本项目不再等待未公开作者代码，也不发送方法询问信。唯一执行路线是：

> 以官方 C15 v1.0 生成统一径向风廓线，以 CLIMADA core/Petals 的公开 TCR 数值实现为骨架，按 Lu et al. (2018)、Xi et al. (2020) 和后续公开 TCR 应用冻结参数，先重建 Hurricane Irene 历史案例，再使用完全相同的冻结内核计算首个合成事件。

这一路线称为 **method-faithful public reconstruction（依据公开论文与公开代码的方法级忠实重建）**。它不是 Xi/Nature 未公开生产代码的 **source-identical reproduction**，也不得写成与作者原始静态场、插值器或二进制输出逐字节相同。

本合同只冻结一个候选，不设置 ER11、R-CLIPER、外缘平滑或其他并行生产分支。Irene 和首个合成事件闭合以前，不启动 10,000 条批量风雨场。

## 2. 冻结的软件与文件身份

| 组件 | 冻结版本 | Git commit / SHA-256 | 用途 |
|---|---|---|---|
| 官方 C15 | PURR v1.0，DOI `10.4231/CZ4P-D448` | 原件及逐文件哈希见 `vendor/CLE15/10.4231_CZ4P-D448/v1.0/PROVENANCE.md` | 统一轴对称完整风廓线 |
| CLIMADA core | `v6.1.0` | commit `bb9595944fcf673152ea53e61efbafeb0b1e4406` | 轨迹 SI 转换、距离和基础 TC 数值设施 |
| CLIMADA Petals | `v6.2.0` | commit `6ecd7af096f126df2da1023fbc5013765566d5e9` | 公开 TCR 数值骨架 |
| `tc_rainfield.py` | Petals `v6.2.0`，75,745 bytes | SHA-256 `6f0bd30dc5532d907401a862f9d8b560c3feca6356c79ac7de444a52b315e062` | 五项垂直速度、差分、静态场采样和降雨积分 |

Petals `v6.2.0` 与 core `v6.1.0` 是本合同冻结的兼容组合。运行 manifest 必须同时记录上述 tag、commit、文件哈希和环境依赖锁；只记录包名或浮动版本号不构成复现。

## 3. C15 接入 CLIMADA TCR 的唯一桥

### 3.1 含科氏风廓线

Irene历史基准仍以观测RMW（仅末端4个缺测值按第6.1节填补）调用官方C15 `rmaxinput`。Lin合成事件则在每个时次以`v_trks`、事件内固定的`r0`和`abs(f)`调用官方C15 `r0input`，并按冻结的Python返回顺序从第五个返回值读取该时次的\(r_m\)。两条分支只是在C15已知尺度输入上不同；TCR 的径向和时间差分均查询所得官方C15廓线，不得让CLIMADA默认ER11进入降雨分支。

官方 C15 使用绝对角动量

\[
M(r)=rV_{C15}(r)+\frac12fr^2,
\]

并以最大风半径处的值

\[
M_m=r_mV_{max}+\frac12fr_m^2
\]

定义归一形状

\[
\mu(r)=\frac{M(r)}{M_m}.
\]

### 3.2 CLIMADA 所需的无科氏伴随风

CLIMADA 的地形和切变项还需要一个保留同一径向角动量形状、但不含科氏项的伴随风。公开代码派生桥冻结为

\[
V_d(r)=\frac{r_mV_{max}\,\mu(r)}{r},\qquad r>0,
\]

并令 \(V_d(0)=0\)。这沿用官方 C15 的归一绝对角动量形状，并对应 Emanuel 公开 `windprofilem` 中以 \(r_mV_{max}\) 替代含科氏 \(M_m\) 的无科氏 companion 构造。它不是 `Vd=V`，也不是重新运行 ER11。

该桥属于 **public-code-derived bridge**：依据官方 C15 角动量定义和公开 TCR 数值代码恢复，可完全审计；没有证据声称它与 Nature 作者未公开桥接源码逐行相同。

C15 provider保持官方有限外半径定义：若官方函数返回数组未显式包含理论端点，则只机械追加精确点`(r0,V=0)`；任何`r>r0`查询直接报错。Lin正式目录须在运行前验证每个事件`r0>302 km`，以覆盖300 km支持域及2 km径向差分；违反者标为`METHOD_DOMAIN_PENDING`，不得重抽、clip、缩域、排除或改权重。不得外推尾部、将`r0`改写为计算域边界或增加taper。

### 3.3 明确禁止

- 不在 C15 外缘增加 cosine taper、平滑带、`r0` clipping 或经验置零；
- 不把 C15 风与 ER11 雨混合；
- 不把旧 `shared_cle_profile_provider` 或 `Vd=V` 带回生产；
- 不为压低雨量而调整 cap、关闭时间导数或修改五项上升速度公式。

## 4. 冻结的 TCR 数值设置

| 设置 | 冻结值 | 公开依据/实现位置 |
|---|---:|---|
| 轨迹与环境时间步 | `1 h` | Xi et al. (2020)；六小时量按其各自规则处理后线性插值 |
| 水平网格 | `0.05° × 0.05°` | Xi et al. (2020) |
| precipitation efficiency | `epsilon_p = 0.9` | Lu et al. (2018)；Xi et al. (2020)；必须显式覆盖 Petals 默认值 |
| lower-troposphere depth | `H = 4000 m` | Xi et al. (2020)；原文 `4000 km` 记录为排印错误 |
| 移动计算域 | 每一小时仅计算距当时 TC 中心 `<=300 km` 的格点 | Xi et al. (2023), JAMC, DOI `10.1175/JAMC-D-22-0131.1`；Petals `max_dist_eye_km=300` |
| 径向差分 | `2 km` | Petals `res_radial_m=2000` |
| 正上升速度上限 | `w <= 7 m s-1` | 公开 Emanuel v6.4 与 Petals `max_w_foreground=7.0` |
| 角动量径向梯度下限 | `dM/dr >= 10 m s-1` | 公开 Emanuel v6.4 与 Petals `np.fmax(10, dMdr)` |
| drag coefficient 下限 | `Cd >= 0.001` | Petals `min_c_drag=0.001`；静态场梯度按公开实现处理 |
| radiative cooling | `-0.005 m s-1` | Lu et al. (2018) |
| Xi 历史分支密度比 | `rho_air/rho_liquid = 0.0012` | Lu/Xi 论文合同；显式覆盖 Petals 的 `0.00117` |

`300 km` 是随风暴中心逐时移动的物理计算/输出域：任何格点只在该时次位于中心 300 km 内时参与计算。它不是固定 bbox，不是 Nature 统计事件定义中的 `250 km`，也不是对已算全场做事后裁图。C15 本身仍按官方方程求解；计算域边界不授权对 C15 廓线做 taper。

这里选用 `w<=7`、`dM/dr>=10`、2 km 差分和 `Cd>=0.001`，是因为它们在公开 TCR 代码链和独立 CLIMADA 实现中有明确、可执行的同向证据。它们是冻结的公开实现规格，不得表述为 Nature 2023 未公开生产代码已被逐字确认。

## 5. 冻结的公开静态场

为了不再猜测 Xi et al. (2020) 未公开的 ECMWF 文件身份，本重建使用 CLIMADA Petals 为 TCR 正式提供的 package data，并按原始字节冻结：

| 字段 | CLIMADA 身份与来源 | 字节数 | SHA-256 |
|---|---|---:|---|
| drag coefficient | `c_drag_500` v1；ERA5 `forecast_surface_roughness` 按 Feldmann et al. (2019) Eqs. 7–8 转换；[官方文件](https://data.iac.ethz.ch/climada/cde33390-d2a2-4032-a26b-26ab285bcec0/c_drag_500.tif) | 999,184 | `1c3f3b525f0c2a9e73f6fe6d3ba3caf7436699f8caeb3edb69548e04fd3f4a42` |
| topography | `topography_land_360as` v1；SRTM 以 GDAL average 聚合至 0.1°、海洋为 0；[官方文件](https://data.iac.ethz.ch/climada/bc5b8fc7-0e73-49d9-a73a-3fad22fdfff5/topography_land_360as.tif) | 3,054,640 | `de8142fe9f50d0cfbd944884ee945bb355b09dfce2f214879e608f87ae0f0951` |

这两个字段是 CLIMADA 的公开、带来源说明、可哈希静态输入，因此优于无 provenance 的 pyTCR 内嵌文件。它们**不等于** Xi et al. (2020) 作者所用 ECMWF 0.25°地形/粗糙度文件的字节同一替代：drag 数据家族为 ERA5，topography 为 SRTM，分辨率也不完全相同。Irene 输出因此只能称公开方法级重建，不能称 Xi Fig. 1 的逐像元复制。

## 6. 两套环境输入仍严格分开

### 6.1 Irene 历史基准

- IBTrACS `2011233N15301` 的位置、最大风和 RMW：六小时标准时次先按既定来源处理，再线性插值到 1 小时；37个标准时次中末端4个RMW缺测，只用Knaff & Zehr (2007) Eq. 6（DOI `10.1175/WAF965.1`）填补并记录`rmw_estimated`/来源标志，不使用CLIMADA pressure-fit估计器；
- NCEP R1 `q925`：每个六小时时次先取距 TC 中心 200 km 圆盘平均，再把所得标量序列插值至 1 小时；
- NCEP R1 shear：每个六小时时次先计算并平均 `u200-u850`、`v200-v850` 于 600–800 km 环带，再插值至 1 小时；
- 使用本合同第 5 节的 CLIMADA 静态场；
- 降雨换算使用 Lu/Xi 论文的密度比 `0.0012`。

当前冻结数据与当前IBTrACS release支持的是Irene全生命周期方法重建。Xi et al. (2020) Table 1印作`0000 UTC 21 Aug–0000 UTC 24 Aug 2011`，但该时段风暴位置与Fig.1(c–d)约33–39°N、81–74°W的NC/VA陆地图域及正文landfall语境不相容。这一项必须标记为`PUBLISHED-TABLE CONFLICT`，不能擅自替作者改日期。Xi et al. (2020)的整体评估使用中心600 km范围；本重建的300 km移动域来自后续Xi et al. (2023)公开运行协议与CLIMADA实现。因此Irene只用于同一C15–TCR方法族的历史量级与空间结构对照，不把完整生命周期累计或Table 1印刷时段冒充Xi Fig.1的精确窗口，也不声称逐像元复制。

Irene正式运行与Stage IV对照已经完成。公开重建输出包含217个逐小时节点（首节点按公开实现为0）和216个可与Stage IV累计结束时刻一一对应的小时。对照使用33–39°N、81–73°W近似论文视觉框、Natural Earth陆地、Stage IV到原始0.05°TCR质心的固定球面最近邻，以及逐时距中心不超过300 km的共同支持。6,881个陆地质心均有Stage IV 216/216小时完整覆盖；TCR与Stage IV累计的空间相关系数为0.7156，平均值为175.09与101.78 mm，平均偏差为+73.31 mm。该结果用于冻结“空间结构有历史信号、绝对量偏湿”的解释边界，不构成Xi Fig.1精确复刻，也不授权调参或加入经验修正。

### 6.2 合成事件

合成事件不得复用 Irene 的 NCEP 环境序列。其近 900/950 hPa 湿度、环境风和切变必须由合成事件自身对应的 CMIP6/Emanuel 公共环境构造产生，并在 runner manifest 中逐项记录。C15–CLIMADA 物理内核、网格、时间步、移动 300 km 域和静态场保持与 Irene 一致；环境输入来源不同。

Lin v1.1保存的合成轨迹既没有outer radius，也没有逐时RMW。Gori et al. (2022)公开的方法是：每个TC从Chavas et al. (2016)经验lognormal分布抽取一次outer radius，并在生命周期内保持不变，再用强度、该outer radius和科氏参数通过C15完整风廓线模型内生求RMW；Xi, Lin & Gori (2023)沿用了Gori storm sets。因此首事件和后续10,000条Lin合成事件统一采用这一fixed-`r0`架构。

Gori未公开其精确分布参数、RNG、seed或event-to-draw映射。本项目冻结的公开重建按Chavas et al. (2016)全球样本Table 1的`median/Q1/Q3=881.0/740.7/1054.4 km`匹配lognormal：`mu=6.78105762593618`、`sigma=0.261776756893889`；以NumPy `PCG64`和seed `20260810`按冻结sample的`event_position=0..9999`各抽一次，生成后只读取不可变目录，不截断、不clip、不重抽、不按盆地或强度改变draw。完整目录与失败语义见[Chavas 2016 fixed-r0公开重建合同](CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)。

每个时次以`v_trks + r0_event + abs(f)`调用官方C15 `r0input`，由其内生给出RMW和完整廓线；`v_trks`同时进入Emanuel `qs900b.m`实际950-hPa暖核湿度。`vmax_trks`只保留事件阈值和目录语义，不参与尺度闭合。逐时Knaff/Nederhoff RMW、项目outer-first尺度模型和事件尺度残差分配均只保留审计归档，不得拼接、补洞、fallback或进入正式hazard production。

### 6.3 合成事件二维风场与道路共同评价域

第一个Lin合成事件的二维风场使用与TCR完全相同的官方C15逐时径向廓线。切向投影的半球符号冻结为 CLIMADA core `v6.1.0` `tctrack_to_si` 的整轨多数表决：若 `count(lat<0) > count(lat>0)` 则 `latsign=-1`（南半球、顺时针），否则 `latsign=+1`（北半球或平局，含赤道节点不计入任一侧）。整条轨迹只用这一个符号；科氏参数仍按该时次 `sin(lat)`。不按过赤道拆轨、丢弃或改权重。北半球逆时针、南半球顺时针；再按Lin & Chavas (2012, DOI `10.1029/2011JD017126`)叠加空间均一的地表背景风：平移矢量乘`0.55`并按同一 `latsign` 气旋向旋转`20°`。不增加径向背景风衰减、外缘taper、经验尾部、入流角或额外`0.85`地表折减。

Lin v1.1对`v_trks`只标为maximum azimuthal wind，现有论文、正式文件attrs和冻结源码均未给出1-min或10-min平均时段。因此正式输出命名为`model-native near-surface wind`，不乘`0.893`，也不冒充Nature文中的10-minute sustained wind。该时段语义不会阻止无阈值的暴露叠加，但在以后引入风损伤函数时必须与函数的风速口径统一。

首事件风场只在与TCR相同的逐时移动`<=300 km`共同评价域上物化，以便风—雨—道路使用完全相同的支持域。这里的300 km不是C15物理截断、不是完整`r0`风足迹，也不是Nature统计模型的250 km POI事件筛选距离；域外保持未评价，而非零风。

道路叠加使用正式0.1°五级机动车道路长度网格。只抽取与0.05°事件质心精确重合的道路格中心，周期统一经度后逐格保存五类道路长度、事件最大模型原生近地面风和事件累计雨量；不做空间插值、域边外推、灾害阈值、脆弱性函数、失效规则或损失估算。该产物是描述性联合暴露底图，不是影响或风险结果。

## 7. 唯一执行顺序与验收

1. **已完成：**冻结并哈希 CLIMADA core、Petals、`tc_rainfield.py` 和两个静态 GeoTIFF；
2. **已完成：**实现 C15 含科氏廓线与第 3.2 节无科氏伴随风桥，完成官方 C15 示例回归；
3. **已完成：**运行 Irene 全生命周期历史基准，保存逐小时 rain rate、事件累计、最大24 h累计和完整 manifest；
4. **已完成：**核对单位、时间积分和移动300 km掩膜，排除旧接口的外缘伪雨带，并完成Stage IV共同支持历史对照；
5. **已完成审计基准：**旧Knaff-`rmaxinput`首事件已生成降雨、二维风场和道路叠加，用于验证数值链可运行；该结果已被fixed-`r0`生产接口取代，不再作为正式hazard结果；
6. **正在执行：**冻结10,000-event fixed-`r0`不可变目录，并以官方C15 `r0input`重算首事件风—雨—道路闭环；
7. **下一步：**fixed-`r0`首事件通过数值与方法核验后，结合Irene偏湿边界冻结批量输出语义，进入既定10,000条历史hazard sample。

“完成风雨场反演”至少要求：Irene 与首个合成事件均由同一 C15 驱动风和 TCR；无 ER11 混搭或外缘 taper；所有输入、版本、参数和单位可追溯；结果明确标注为公开方法级重建而非作者源码逐字复现。

## 8. 来源

- Chavas, Lin & Emanuel (2015), C15: <https://doi.org/10.1175/JAS-D-15-0014.1>
- Official C15 archive: <https://doi.org/10.4231/CZ4P-D448>
- Lu et al. (2018), TCR equations: <https://doi.org/10.1175/JAS-D-17-0264.1>
- Xi, Lin & Smith (2020), historical validation: <https://doi.org/10.1175/JHM-D-20-0035.1>
- Xi, Lin & Gori (2023), Nature method chain: <https://doi.org/10.1038/s41558-023-01595-7>
- Gori et al. (2022), event-fixed outer radius and C15 storm sets: <https://doi.org/10.1038/s41558-021-01272-7>
- Chavas et al. (2016), observed global outer-radius distribution: <https://doi.org/10.1175/JCLI-D-15-0731.1>
- Xi et al. (2023), explicit 300 km TCR application domain: <https://doi.org/10.1175/JAMC-D-22-0131.1>
- Knaff & Zehr (2007), Irene-only missing-RMW fill: <https://doi.org/10.1175/WAF965.1>
- Lin & Chavas (2012), surface background-wind correction: <https://doi.org/10.1029/2011JD017126>
- CLIMADA core `v6.1.0`: <https://github.com/CLIMADA-project/climada_python/tree/v6.1.0>
- CLIMADA Petals `v6.2.0`: <https://github.com/CLIMADA-project/climada_petals/tree/v6.2.0>
