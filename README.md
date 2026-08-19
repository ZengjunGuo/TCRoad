# TCRoad

全球热带气旋对机动车道路的影响：先算当前气候，再扩到同一 GCM 的四个 SSP。

详细进度（从生产机读取，不以本页记忆为准）：[RESEARCH_STATUS.md](RESEARCH_STATUS.md)

大文件（planet PBF、Lin 轨迹、事件 NetCDF、Crowther GeoTIFF）在服务器  
`/mnt/sdb_test/tang/zengjun/TC_Road_Risk`。本仓库放方法、代码、造价账本、测试和定稿图。

---

## 这条研究怎么走（已冻结）

1. **危险性：** Nature / Gori / Xi 公开重建。Lin v1.1 轨迹，MPI-ESM1-2-LR `r1i1p1f1`。官方 C15 风 + CLIMADA TCR 雨，移动 300 km / 0.05° 盘。历史 1995–2014 近路样本 **99,242** 条里，**99,234** 条风雨场已收口；8 条 `METHOD_DOMAIN_PENDING` 留在样本里，不重抽、不改权重。
2. **暴露：** OSM `planet-260803`（2026-08-03）全部机动车道，对象级，不是 0.1° 美元网格。
3. **重建成本：** 长度 × 2025 美元/km。发展中国家用 ROCKS Actual new-build；美日澳欧中用国家账本。全球合计约 **$49.230 T**（含地方路）/ **$35.805 T**（不含）。合同：[ROAD_ASSET_VALUATION_CONTRACT.md](methods/ROAD_ASSET_VALUATION_CONTRACT.md)
4. **风致资产损失（正在对象级生产）：** 路面不按风速乘损坏比重建。清树用 Escobedo 中等价 × Crowther 树密；桥要阵风和设计重现期两关都过才记整桥重置。合同：[WIND_ASSET_IMPACT_CONTRACT.md](methods/WIND_ASSET_IMPACT_CONTRACT.md)，算法说明：[HISTORICAL_WIND_ASSET_LEDGER.md](methods/HISTORICAL_WIND_ASSET_LEDGER.md)
5. **流量（方法已定，未生产）：** WorldOD / GlODGen 通勤 OD 再分配到同一 OSM 图，**不用观测 AADT**。[ROAD_FLOW_ASSIGNMENT_CONTRACT.md](methods/ROAD_FLOW_ASSIGNMENT_CONTRACT.md)
6. **雨变水深（方法已定，未跑）：** SFINCS 类二维淹没。雨毫米不是水深，没有水深不算洪水美元。[RAIN_TO_INUNDATION_CONTRACT.md](methods/RAIN_TO_INUNDATION_CONTRACT.md)

资产损失和用户损失分开报，不混成一笔。未来 8 窗 Lin 环境/轨迹已经有了，**尚未**做未来 C15–TCR。

---

## 定稿图

### 1. 当前气候合成台风活动与道路层级

加权 TC 通行频率与五级全球路网的空间重叠；Gulf Coast、Pearl River Delta、Bengal Delta 三处真实路网放大。描述性底图，不是损失图。

- [PNG](figures/global_tc_roads_baseline/global_tc_roads_baseline.png)
- [PDF](figures/global_tc_roads_baseline/global_tc_roads_baseline.pdf)
- [作图代码](figures/global_tc_roads_baseline/plot_global_tc_roads_figure.py)

### 2. 全球机动车道路 2025 重建成本

国家尺度美元/km，以及孟加拉三角洲、珠江口、美国湾岸三处对象级放大。造价图，不是损失图。

- [PNG](figures/road_replacement_cost/road_replacement_cost.png)
- [PDF](figures/road_replacement_cost/road_replacement_cost.pdf)
- [作图代码](figures/road_replacement_cost/plot_road_replacement_cost_figure.py)

---

## 现在在做什么

历史风雨场已经按合同收口。对象级风致资产账（清树 + 塌桥）已在 2026-08-19 算完 99,234 场；8 条 pending 为 0。加权清树约 **$5.389 B**，加权塌桥约 **$5.419 M**（2025 美元；未加权事件加总不是年均）。流量、SFINCS、未来情景风雨场都还没开。

---

## 方法合同（完整列表）

- [C15–CLIMADA 公开重建](methods/C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)
- [Nature C15–TCR](methods/NATURE_C15_TCR_REPRODUCTION_CONTRACT.md)
- [Chavas 2016 事件固定 r0](methods/CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)
- [道路重建成本](methods/ROAD_ASSET_VALUATION_CONTRACT.md)
- [风致资产损失](methods/WIND_ASSET_IMPACT_CONTRACT.md)
- [流量赋值](methods/ROAD_FLOW_ASSIGNMENT_CONTRACT.md)
- [降水转淹没](methods/RAIN_TO_INUNDATION_CONTRACT.md)

Irene 2011 是同一 C15–TCR 内核的历史方法对照（偏湿，不是 Xi 图逐像元复刻）。对照摘要在 [validation/irene2011](validation/irene2011) 与 [code/compare_irene_stage4.py](code/compare_irene_stage4.py)。
