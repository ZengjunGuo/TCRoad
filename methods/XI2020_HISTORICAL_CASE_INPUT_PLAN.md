# Xi et al. (2020) 历史案例最小输入合同：Hurricane Irene (2011)

**用途：**用一个公开、可追溯的历史事件检查官方C15–CLIMADA TCR实现是否形成Xi et al. (2020)方法链的公开忠实重建。它不是合成气候生产输入，也不把Xi的历史NCEP环境量移植到Nature Climate Change 2023的合成事件链。

状态词只有四种：

- **VERIFIED**：论文或数据机构的一手来源已经明确，可以冻结和下载；
- **PUBLIC-RECONSTRUCTION DIFFERENCE**：作者原始身份未公开；本项目使用已冻结且完整披露的公开实现，因此可执行但不能声称字节/源码同一；
- **PUBLISHED-TABLE CONFLICT**：论文表格印刷值与论文图域/实际轨迹互不相容；必须原样记录，不能擅自替作者改成一个“合理”日期；
- **FORBIDDEN**：不得为了跑通而自行替换或猜选。

## 1. 案例与时段

选择 **Hurricane Irene (2011)**，因为 Xi et al. (2020) Fig. 1(c–d) 直接给出了 Irene 的 C15–TCRM 与 Stage IV 对照，它因而比任意另选风暴更接近论文已有的可核结果。

| 项目 | 冻结值 | 状态 |
|---|---:|---|
| IBTrACS SID | `2011233N15301` | **VERIFIED**（当前 v04r01） |
| ATCF/NHC ID | `AL092011` | **VERIFIED** |
| 盆地/名称 | North Atlantic / IRENE | **VERIFIED** |
| 完整六小时时段 | `2011-08-21 00:00 UTC` 至 `2011-08-30 00:00 UTC`，首尾均含，共 37 个 NCEP 时次 | **VERIFIED**（当前 v04r01 事件边界） |
| 模型时间步 | 6-hourly轨迹线性插值至1-hourly；六小时时次先提取环境参数，再将参数序列线性插值至1-hourly | **VERIFIED**（Xi 2020） |
| Table 1印刷的selected period | `0000 UTC 21 Aug–0000 UTC 24 Aug 2011` | **PUBLISHED-TABLE CONFLICT**：该时段Irene仍位于加勒比/巴哈马一带，与Fig.1(c–d)约33–39°N、81–74°W的NC/VA陆地区域和正文landfall语境不相容 |
| Fig. 1 精确降雨累计窗口 | 不能由上述冲突表格值或完整事件边界替代；本次主基准使用完整Irene生命周期，保留后续独立冻结美国影响窗口的可能，但不冒充Fig.1逐值复刻 | **PUBLIC-RECONSTRUCTION DIFFERENCE** |
| TCR计算网格与时间步 | `0.05° × 0.05°`，`1 h` | **VERIFIED** |
| 计算域 | 每小时距当时TC中心`<=300 km` | **PUBLIC-RECONSTRUCTION DIFFERENCE**：来自Xi et al. (2023) JAMC公开TCR协议与CLIMADA Petals实现；Xi et al. (2020)历史评估关注的是中心600 km内降雨 |

本合同下载并运行完整Irene生命周期的最小环境包；它避免遗漏TCR所需的时间演化，但不扩展到2011全年或1999–2018全时段。Table 1的印刷日期必须在provenance中原样记录为论文内部冲突，不得自行改写成8月27–29日并称其为作者原值。

## 2. IBTrACS：轨迹、强度与 RMW

### 2.1 官方来源

- 数据集：NOAA/NCEI **IBTrACS v04r01**, DOI [`10.25921/82ty-9e16`](https://doi.org/10.25921/82ty-9e16)
- 官方说明与引用要求：<https://www.ncei.noaa.gov/products/international-best-track-archive>
- 当前 North Atlantic NetCDF（生产只取这一份原始盆地文件）：
  <https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/netcdf/IBTrACS.NA.v04r01.nc>
- CSV 只作人工核对备选，不与 NetCDF 重复进入生产：
  <https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv>
- 官方字段文档：<https://www.ncei.noaa.gov/sites/g/files/anmtlf171/files/2025-04/IBTrACS_v04r01_column_documentation.pdf>

IBTrACS v04r01 是持续更新产品；因此版本字符串本身不足以固定数据，必须同时固定下载日期、原始字节数和 SHA-256。

### 2.2 事件提取合同

只保留 `sid == "2011233N15301"`，并读取：

| 类别 | 原始字段 | 用途/单位合同 |
|---|---|---|
| 标识与审计 | `sid`, `season`, `name`, `iso_time`, `usa_atcf_id`, `iflag` | 事件、时间与插值来源审计 |
| 位置 | `usa_lat`, `usa_lon` | NHC/USA 同源位置，degree north/east |
| 强度 | `usa_wind` | 1-min sustained wind，原始 `kt`；模型副本乘 `0.514444` 转 `m s-1` |
| RMW | `usa_rmw` | 原始 `n mi`；模型副本乘 `1.852` 转 `km` |
| 辅助 QC | `usa_pres`, `usa_status` | 不作为必需动力输入，但保留用于来源与生命周期核对 |

约束：

1. 北大西洋位置、风速和 RMW 必须采用同一 `usa_*`/NHC 来源，不能把 IBTrACS 跨机构平均位置与 `usa_wind`、`usa_rmw` 混合。
2. 模型六小时时序只取 `00/06/12/18 UTC` 的原始同化时次；IBTrACS 中 3-hourly 插值行和非标准 landfall special rows 不直接加入六小时输入。`iflag` 必须随子集保存，以证明没有把插值值当原始观测。
3. 当前37个Irene标准时次中有33个`usa_rmw`观测值，末端4个时次缺测。只对这4个缺测时次使用Knaff & Zehr (2007), Eq. (6), DOI [`10.1175/WAF965.1`](https://doi.org/10.1175/WAF965.1)：`RMW_km = 66.785 - 0.09102*Vmax_kt + 1.0619*(lat_deg-25)`。Irene位于北半球，`lat_deg`按北纬正值代入。保存原值、估计值及`rmw_source=KZ07_EQ6`/`rmw_estimated=true`；不得重估已有RMW，也不得调用CLIMADA基于压强拟合的RMW估计器。
4. 原始 `kt`、`n mi` 字段与 SI 派生字段并存，不覆盖原始值。

**版本边界：**Xi et al. (2020) 没有公开其下载时使用的IBTrACS release/文件哈希。当前v04r01形成公开数据的**方法级忠实重建**，不能声称与作者2020年输入逐字节相同。

## 3. NOAA PSL NCEP/NCAR Reanalysis 1：最小环境场

Xi et al. (2020) 的历史链需要：

- 在每个六小时时次，以对应TC中心为圆心，将 `shum` at 925 hPa 在 `r <= 200 km` 圆盘内平均；
- 在每个六小时时次，先算 `(u200-u850, v200-v850)`，再在 `600 <= r <= 800 km` 环带内平均；
- 将上述六小时环境参数序列线性插值到1-hourly；TC位置、强度和RMW则由各自六小时时序独立线性插值到1-hourly。

这一区分直接来自 Xi et al. (2020) 的表述：文中先定义从NCEP场提取的圆盘/环带平均环境参数，随后说明“6-hourly environmental data”插值至1小时。因空间平均与时间插值不对易，正式复现不得改成“先插值整幅NCEP场，再逐小时做空间平均”。

### 3.1 官方文件、变量和服务

官方目录：<https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/pressure/catalog.html>

| 年文件 | 变量/单位 | 层 | 官方 OPeNDAP 数据集端点 | 最小下载方式 |
|---|---|---|---|---|
| `shum.2011.nc` | `shum`, `kg kg-1` | 925 hPa | <https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/pressure/shum.2011.nc> | NCSS 单层、单事件 bbox/time |
| `uwnd.2011.nc` | `uwnd`, `m s-1` | 200, 850 hPa | <https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/pressure/uwnd.2011.nc> | 两个单层 NCSS 子集 |
| `vwnd.2011.nc` | `vwnd`, `m s-1` | 200, 850 hPa | <https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/pressure/vwnd.2011.nc> | 两个单层 NCSS 子集 |

官方 OPeNDAP 元数据显示：时间为 6-hourly；水平网格为 2.5°，纬度 `90` 至 `-90`，经度 `0` 至 `357.5`；`shum` 文件包含 925 hPa，`uwnd`/`vwnd` 文件包含 850 和 200 hPa。不要误用 NCEP/DOE Reanalysis 2、daily mean、sigma-level 或 pressure-level climatology。

### 3.2 最小固定包络

包络由 Irene 全轨迹和论文平均半径向外扩展，再向外吸附至 NCEP 2.5°格点：

| 输入 | 南—北 | 西—东（NCEP 0–360 经度） | 说明 |
|---|---|---|---|
| `shum925` | `12.5–55.0°N` | `277.5–305.0°E` | 轨迹外扩 200 km 后的固定下载包络 |
| `u/v` 200、850 hPa | `7.5–62.5°N` | `272.5–312.5°E` | 轨迹外扩 800 km 后的固定下载包络 |

这些矩形只是避免下载全年全球场的**提取包络**，绝不能代替逐时的 200 km 圆盘或 600–800 km 环带。掩膜距离必须用球面大圆距离，不能用经纬度欧氏距离。

### 3.3 可冻结的 NCSS 请求

NCSS 根端点：

- `https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/shum.2011.nc`
- `https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/uwnd.2011.nc`
- `https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/vwnd.2011.nc`

精确请求如下；每个 URL 和返回文件分别进入 manifest。

**q925**

```text
https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/shum.2011.nc?var=shum&north=55&west=277.5&east=305&south=12.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=925&accept=netcdf4
```

**u200**

```text
https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/uwnd.2011.nc?var=uwnd&north=62.5&west=272.5&east=312.5&south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=200&accept=netcdf4
```

**u850**

```text
https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/uwnd.2011.nc?var=uwnd&north=62.5&west=272.5&east=312.5&south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=850&accept=netcdf4
```

**v200**

```text
https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/vwnd.2011.nc?var=vwnd&north=62.5&west=272.5&east=312.5&south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=200&accept=netcdf4
```

**v850**

```text
https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure/vwnd.2011.nc?var=vwnd&north=62.5&west=272.5&east=312.5&south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=850&accept=netcdf4
```

推荐文件名：

```text
ncep_r1_shum925_irene2011_20110821T00_20110830T00_bbox12p5_55_277p5_305.nc
ncep_r1_uwnd200_irene2011_20110821T00_20110830T00_bbox7p5_62p5_272p5_312p5.nc
ncep_r1_uwnd850_irene2011_20110821T00_20110830T00_bbox7p5_62p5_272p5_312p5.nc
ncep_r1_vwnd200_irene2011_20110821T00_20110830T00_bbox7p5_62p5_272p5_312p5.nc
ncep_r1_vwnd850_irene2011_20110821T00_20110830T00_bbox7p5_62p5_272p5_312p5.nc
```

## 4. Xi原始静态场身份与本项目的公开重建选择

Xi et al. (2020)引用Dee et al. (2011)，其作者原始静态数据家族可追溯到 **ERA-Interim**。论文只写“ECMWF 0.25° topographic height and surface roughness”，没有给出`paramId`、analysis/forecast类型、日期、原始文件、插值算子或哈希。因此，作者究竟使用surface geopotential、`paramId=173`还是`paramId=244`，以及海洋编码和重网格方式，仍属 **PUBLIC-RECONSTRUCTION DIFFERENCE**。

这一区别限制的是`source-identical`主张，不再阻止透明的公开方法级重建。禁止按“看起来合理”选择一个ERA-Interim字段后称其为作者原文件，也禁止使用无provenance的pyTCR静态场。

### 4.1 冻结的CLIMADA公开静态场

Irene公开重建只使用CLIMADA Petals正式提供的package data，并冻结其原始字节：

| 字段 | 官方URL | 字节数 | SHA-256 |
|---|---|---:|---|
| `c_drag_500.tif`，0.25°，ERA5 `forecast_surface_roughness`按Feldmann Eqs. 7–8转换 | <https://data.iac.ethz.ch/climada/cde33390-d2a2-4032-a26b-26ab285bcec0/c_drag_500.tif> | 999,184 | `1c3f3b525f0c2a9e73f6fe6d3ba3caf7436699f8caeb3edb69548e04fd3f4a42` |
| `topography_land_360as.tif`，0.1°，SRTM average、海洋为0 | <https://data.iac.ethz.ch/climada/bc5b8fc7-0e73-49d9-a73a-3fad22fdfff5/topography_land_360as.tif> | 3,054,640 | `de8142fe9f50d0cfbd944884ee945bb355b09dfce2f214879e608f87ae0f0951` |

这两个字段有公开metadata、稳定URL和可复核哈希，是本项目的正式公开重建输入。其数据家族和分辨率与Xi原始描述不完全相同，所以manifest与论文必须明确标记`PUBLIC-RECONSTRUCTION DIFFERENCE`，不能称为Xi静态输入的字节同一副本。

## 5. Feldmann/Esau 粗糙度到拖曳系数

Xi et al. (2020) 明确沿用 Feldmann et al. (2019) 的做法；可执行公式应忠实采用 Feldmann 的 Eqs. (7–8)，而不是另写一般 10 m 对数律：

\[
C_D^{0}=\left[\frac{\kappa}{\ln(500/z_0)}\right]^2,
\qquad
C_D=0.9\,\frac{C_D^{0}}{1+50C_D^{0}},
\qquad \kappa=0.35,
\]

其中 `z0` 以米计，`500` 也是米；所得 `Cd` 无量纲。Feldmann 将此描述为“loosely following Esau (2004)”，所以实现的直接规范是 Feldmann Eqs. (7–8)，不能用 Esau 论文中其他边界层参数化替换。

来源：

- Feldmann et al. (2019), DOI [`10.1175/JAMC-D-19-0011.1`](https://doi.org/10.1175/JAMC-D-19-0011.1)，作者公开 PDF：<https://texmex.mit.edu/pub/emanuel/PAPERS/Feldmann_etal_2019.pdf>
- Esau (2004), DOI [`10.5194/angeo-22-3353-2004`](https://doi.org/10.5194/angeo-22-3353-2004)

Xi作者对`z0<=0`、海洋值和缺测值的具体处理仍未公开。本项目不重新处理原始`z0`，而直接读取上节冻结的CLIMADA成品drag GeoTIFF，并按Petals公开规则使用`Cd>=0.001`；因此不存在本地猜选阈值。该选择不被写成Xi源码同一。

## 6. 原始文件、SHA-256 与许可清单

每个下载对象都执行同一冻结协议：

1. HTTP 返回体原样写入 `raw/`；在解码、裁剪、重命名变量或经度转换前计算 SHA-256。
2. manifest 至少记录：完整 URL（含 NCSS query）、UTC 获取时间、HTTP status、`ETag`/`Last-Modified`/`Content-Length`（若服务器返回）、字节数、SHA-256、数据集 ID/版本、变量/单位/dtype/dimensions、坐标方向与范围、层、时间范围和 bbox。
3. 对 IBTrACS 同时记录：原始 NA 盆地文件 SHA、Irene 事件派生表 SHA、提取脚本 Git commit/SHA、筛选条件、单位转换和 RMW 缺测标志。
4. 本地与服务器传输后分别重算 SHA-256，完全一致才允许进入计算；不一致即停止。
5. 任何持续更新 URL 都不能仅靠 URL/文件名复现，必须由 manifest 中的 SHA-256 锁定具体字节。
6. NOAA NCSS对完全相同的查询可能生成数值相同但HDF5容器字节不同的文件。manifest因此同时保存原始文件SHA-256和解码后变量身份/掩膜/数值的规范哈希；前者固定本次输入字节，后者判断重复获取的科学内容是否一致。

许可与引用：

| 数据 | 使用边界 | 论文引用/致谢 |
|---|---|---|
| IBTrACS v04r01 | NOAA/NCEI 公开数据；持续更新，必须固定 access date 和 SHA | Knapp et al. (2010) BAMS + Gahtan et al. (2024) 数据集 DOI `10.25921/82ty-9e16`，注明 subset/access date |
| NCEP/NCAR R1 via NOAA PSL | PSL 公开数据政策下无使用限制；保留来源信息 | Kalnay et al. (1996), *BAMS* 77, 437–471，并致谢 NOAA PSL；参考 <https://psl.noaa.gov/data/atmoswrit/reference.html> |
| CLIMADA package data | CC BY 4.0；保留Data API身份、URL、SHA和转换说明 | CLIMADA core/Petals；drag同时引用Feldmann et al. (2019)，topography注明SRTM-derived |
| ERA-Interim | 只作为Xi原始文献来源边界，不作为本次运行文件 | Dee et al. (2011)；不得暗示本次静态场与作者字节同一 |

## 7. Go/No-Go

### 已完成

- IBTrACS NA v04r01 原文件已经冻结，并唯一定位Irene；
- 上述五个NCEP NCSS子集已经冻结；
- schema、坐标、37个时次、单位、原始文件SHA和解码变量规范哈希已经完成；
- CLIMADA core v6.1.0、Petals v6.2.0、`tc_rainfield.py`及两个静态GeoTIFF已经冻结；
- Irene末端4个RMW缺测的Knaff–Zehr Eq. 6填补规则已经冻结。
- Irene全生命周期C15–TCR公开方法级重建已经完成；Stage IV对照使用216个一小时累计、Natural Earth陆地、固定球面最近邻和逐时移动300 km共同支持，6,881个陆地格点全部有完整观测覆盖。

正式服务器目录：`data/observations/xi2020_validation/irene2011_public_inputs_v1`。

### 仍然不能声称

- 不能把已完成的公开方法级重建称为Xi 2020作者生产代码或Fig.1精确窗口的复刻；
- 不能称静态场、作者代码或Fig.1累计窗口与Xi逐字节/逐像元相同；
- 不能把空间相关0.7156且绝对量偏湿的历史结果包装成“完全验证”或`source-identical reproduction`。

### 唯一下一步

不发送邮件。Irene历史重建与Stage IV共同支持对照已经完成；现在按[C15–CLIMADA公开重建合同](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)和[Chavas 2016 fixed-r0公开重建合同](CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)重算首个Lin合成事件。物理内核保持1 h、0.05°、移动300 km、`epsilon_p=0.9`、`H=4 km`、2 km径向差分、`w<=7`、`dM/dr>=10`和`Cd>=0.001`；合成湿度、切变和轨迹必须来自该事件自身的公开环境构造，不得复用Irene NCEP序列。每个Lin事件按Chavas et al. (2016)全球Table 1分位数匹配lognormal，以NumPy `PCG64`、seed `20260810`只抽一次生命周期固定`r0`；逐时用`v_trks+r0+abs(f)`调用官方C15 `r0input`内生求RMW。`vmax_trks`只保留事件阈值/目录语义；Knaff/Nederhoff、项目outer-first尺度模型及残差分配不进入正式hazard production。

## 8. 方法来源

- Xi, Lin & Smith (2020), *Evaluation of a Physics-Based Tropical Cyclone Rainfall Model for Risk Assessment*, DOI: <https://doi.org/10.1175/JHM-D-20-0035.1>
- Xi, Lin & Gori (2023), *Increasing sequential tropical cyclone hazards along the US East and Gulf coasts*, DOI: <https://doi.org/10.1038/s41558-023-01595-7>
- Dee et al. (2011), ERA-Interim, DOI: <https://doi.org/10.1002/qj.828>
- Kalnay et al. (1996), NCEP/NCAR Reanalysis, DOI: <https://doi.org/10.1175/1520-0477(1996)077%3C0437:TNYRP%3E2.0.CO;2>
- Knaff & Zehr (2007), RMW Eq. 6, DOI: <https://doi.org/10.1175/WAF965.1>
- Gori et al. (2022), event-fixed outer radius and C15 storm sets, DOI: <https://doi.org/10.1038/s41558-021-01272-7>
- Chavas et al. (2016), global observed outer-radius distribution, DOI: <https://doi.org/10.1175/JCLI-D-15-0731.1>
- Xi et al. (2023), moving 300 km TCR application domain, DOI: <https://doi.org/10.1175/JAMC-D-22-0131.1>
