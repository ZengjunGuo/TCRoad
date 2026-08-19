# 历史窗口风致道路资产损失：用了哪些风、算了哪些路、怎么算的

状态：对象级账本已冻结方法并开跑历史窗口。  
合同：[`WIND_ASSET_IMPACT_CONTRACT.md`](WIND_ASSET_IMPACT_CONTRACT.md)  
代码：`code/road_wind_object_join.py`（采风）+ `code/road_wind_asset_impact.py`（定价核）

---

## 人话

台风风把路边的树刮倒，要清；特别大的阵风还可能把桥刮垮。这一步只算这两笔 **2025 美元**，算在每一条 OSM 机动车路上，不算车流、不算水淹、不算未来气候。

路面本身 **不按风速乘损坏比重建**。隧道两笔都是 0。

---

## 用了哪些历史风速

只用已经收口的历史风雨场：

- 运行：`lin_road_domain_300km_v1`
- 气候窗：Lin 合成轨迹，MPI-ESM1-2-LR `r1i1p1f1`，**1995–2014**
- 文件：每个事件一个 compact NetCDF
- 字段：**`event_maximum_near_surface_wind_speed`**（官方 C15 加 Lin–Chavas 背景后的 **模型原生近地面风**，不是 10 分钟持续风、也不改名）
- 参加算风的事件：**99,234** 条（有 compact + overlap）
- **不参加** 的 8 条 `METHOD_DOMAIN_PENDING`（没有足迹、没有风、没有美元，样本权重不动）：

  `11902, 11944, 12357, 50194, 62311, 68925, 72126, 86977`

- 未来 8 个 SSP 窗口的轨迹 **不用**。

经度按周期折到 `[-180, 180)`，再贴到 compact 的 **0.05°** 格子上。格子里有这个事件的足迹，就用该格的事件最大 C15；没有足迹就是“这场台风没扫到这条路”，风速记为缺失，美元为 0。不做插值，不用 0.1° 道路密度网格当美元网格。

Koks / Virot 的树折阈值是 **3 秒阵风 151 km/h**。本库 C15 按 Harper 陆上 \(G_{3s/10min}=1.66\) 换成：

\[
151 / 1.66 = 91.0\ \mathrm{km\,h^{-1}} = 25.3\ \mathrm{m\,s^{-1}}.
\]

主文不用未换算的 151 km/h 直接卡 C15。

---

## 算了哪些路

- OSM 快照 **`planet-260803`（2026-08-03）** 上已经赋过 2025 重建成本的机动车 way
- 代表点经纬度与提取表中点相同（赋值 CSV 本身没存 lon/lat，从 `planet-260803_motor_ways` 按 `way_id` 接上）
- 接受的机动车对象约 **1.108 亿条**（`110,822,264` accepted ways）
- **不是** 0.1° 美元网格
- 隧道：清树、塌桥都是 0
- 桥：可以有塌桥费，没有清树费

---

## 具体怎么算

树密 \(N\)（茎 km⁻²，DBH ≥ 10 cm）来自 Crowther et al. (2015) 生物群系 WGS84 GeoTIFF（Yale EliScholar `Revision_01`），取代表点所在像元。sha256 `1812e5cbb17f91f3a1dfc3033e9cc402bc557ad6ed3827c84ce1fc3f8f05c338`。

\[
P(N)=\begin{cases}
0 & N\le 0\\
\min(N,\,10^{4})/10^{4} & N>0
\end{cases}
\]

清树费（每条非桥、非隧道的接受 way，每个有足迹的事件）：

\[
L_{\mathrm{cleanup}}
=
L_{\mathrm{km}}
\times 4979
\times \mathbf{1}\{V_{\mathrm{C15}}\ge 25.3\,\mathrm{m\,s^{-1}}\}
\times P(N)
\]

\(4979\) 美元/km 是 Escobedo (2009) 中等清理体积 × 28.25 美元/m³，按同一套美国 GDP 平减指数从 2005 胀到 2025。不用 Koks 补充材料里假设的 5k–50k 美元/km 档。

塌桥费（每座 OSM 桥，每个有足迹的事件）：

- 阵风换到 C15 后必须 **大于** `gmtra` 该类 \(V^{*}\)（primary 218.4 / secondary 203.3 / other 188.3 km/h C15）
- **并且** 用这 20 年 compact 在该桥上的事件最大风速算经验重现期，必须 **严于** `gmtra` 设计重现期（按世界银行收入组和道路类）
- 两关都过：损失 = 该桥已冻结的 GIRI 2025 重建成本；否则 0
- 8 条 pending 不进入这座桥的峰值序列

资产损失 = 清树 + 塌桥。与用户（车流）损失分开记账。

---

## 入口

```bash
python3 code/road_wind_object_join.py score-event ways.csv \
  --compact EVENT.nc --output impact.csv
python3 code/road_wind_object_join.py score-historical \
  --valued-dir runs/valuation/planet-260803_valued \
  --extract-dir runs/valuation/planet-260803_motor_ways \
  --compact-dir runs/hazard_production/lin_road_domain_300km_v1/compact_hazard_footprint \
  --trees data/impact/raw/crowther_extract/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif \
  --output-dir runs/impact/historical_wind_asset_v1
```
