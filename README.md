# TCRoad

研究说明与 **2026-08-16 生产进度**（从服务器当场读取，不是旧 README 记忆）：

[RESEARCH_STATUS.md](RESEARCH_STATUS.md)

方法合同：道路重建成本、[流量（WorldOD，不用 AADT）](methods/ROAD_FLOW_ASSIGNMENT_CONTRACT.md)、[降水转淹没（SFINCS，今夜不跑）](methods/RAIN_TO_INUNDATION_CONTRACT.md)。

本目录保存可以复现的方法、代码、造价账本和测试。大型气候、轨迹、OSM 与风雨场数据继续保存在服务器：

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk`

## 当前正式成果

### 全球合成热带气旋活动与道路层级图

- [PNG 预览](figures/global_tc_roads_baseline/global_tc_roads_baseline.png)
- [矢量 PDF](figures/global_tc_roads_baseline/global_tc_roads_baseline.pdf)
- [绘图脚本](figures/global_tc_roads_baseline/plot_global_tc_roads_figure.py)
- [数据与输出摘要](figures/global_tc_roads_baseline/global_tc_roads_baseline_data_summary.json)
- [图件 QA](figures/global_tc_roads_baseline/global_tc_roads_baseline_qa.txt)

图件范围：MPI-ESM1-2-LR、r1i1p1f1、historical、1995–2014。它展示加权合成 TC 通行频率与五级全球道路系统的空间重叠，并包含 Gulf Coast、Pearl River Delta 和 Bengal Delta 三个真实道路网络放大图。

这是一张当前气候基线的描述性地图，不是道路损失、网络中断或多 GCM 集合结论。

正式服务器代码提交：`d3d2296`

| 文件 | SHA-256 |
|---|---|
| PDF | `86822fddc4a69418c34f621e47e13adfcc16fa217bf981ee7cd3e5a65fd35f44` |
| PNG | `04e62f4dbba9a97e4bfd4b4719d83f9f3abc0f7f69201773d8c9c4f5d9a59f05` |
| 绘图脚本 | `f038f0b2f7ae3d578bfbe19c6ad8e609556726b5b945e54fdbb62269d7a2c200` |

## 研究当前状态

- 已完成一个 GCM（MPI-ESM1-2-LR）、一个成员的1995–2014当前气候基线：20个环境年、每年5000条接受轨迹，共100,000条。
- 已形成按 IBTrACS 1995–2014校准的加权事件目录，全球年均频率为84.9 storms/year。
- 已固定抽取10,000条 hazard production sample。既有 outer-first 风暴尺度模型及事件残差分配产物保留作方法审计，但不进入正式 hazard production。
- 正式 10,000 事件 fixed-`r0` 目录已经冻结：PCG64 `standard_normal` 显式变换，存到最近 1 mm；最小 `r0=308.736422 km`，没有任何事件落入 302 km 诊断域以内。
- Lin合成事件的尺度闭合已按Nature/Gori原方法架构纠正：每个事件从Chavas et al. (2016)全球经验lognormal分布抽取一次固定outer radius，逐时RMW由官方C15 `r0input`内生求得。Knaff/Nederhoff和项目outer-first模型均退出正式生产链。
- 旧 shared CLE–pyTCR 已完成一次首事件接口诊断，但该自定义耦合不属于当前冻结的 C15–TCR 公开重建路线，结果已退出科学生产链。
- 全球 OSM 五级道路0.1°长度栅格与三个局部真实道路网络已经完成。
- 同一 GCM 的四个 SSP 原始 CMIP6 输入已经下载，但未来 Lin 环境场和未来轨迹尚未生成。
- 已完成 Hurricane Irene 2011 的官方 C15 + CLIMADA TCR 全生命周期公开方法级重建，并完成与 NCEP Stage IV 216 小时降雨的同窗、同陆地、逐时移动 300 km 共同支持对照。
- Irene 对照在6,881个陆地格点上的空间相关系数为0.716；模型与Stage IV平均累计分别为175.09和101.78 mm。它证明公开重建已产生有意义的历史空间结构，同时记录了绝对量偏湿的解释边界；不表述为Xi Fig. 1逐像元复刻或“完全验证”。
- 已完成Knaff版本首事件的94小时风—雨—道路数值闭环及独立QA；它证明了TCR、二维风场和道路叠加代码可运行，但尺度接口现已被更忠实的fixed-`r0`方法替代，因此该组数值降为回归/审计基准，不再作为正式hazard结果引用。
- 已完成fixed-`r0`首事件重算：`r0=1000.708405 km`，逐时RMW 62.58–256.31 km，最大小时雨 91.91 mm h⁻¹（位于当时RMW环），最大累计 959.29 mm，最大近地面风 30.47 m s⁻¹。外缘伪雨带未再出现。Lin `v_trks`的平均时段仍未公开，因此风场继续保持model-native near-surface语义，不擅自乘时段换算系数，也不冒充10-min持续风。
- 官方 C15 v1.0 原件已经冻结；Python 3兼容层已对照官方三张MATLAB示例图的38个数值锚点闭合，三类默认入口最大误差为0.001609 m s⁻¹。
- Xi et al. (2020) Hurricane Irene历史案例的IBTrACS与五个NCEP最小公开输入已经在服务器冻结。
- 唯一公开重建实现已经冻结：官方C15 v1.0 + CLIMADA core v6.1.0 + Petals v6.2.0；以Petals公开TCR数值骨架接入C15，不使用默认ER11，不增加外缘taper。
- Irene静态场改用CLIMADA正式提供且可哈希的`c_drag_500`与`topography_land_360as`。这使历史基准可执行、可审计，但它们不是Xi作者未公开ECMWF文件的字节同一副本，因此结果口径是`method-faithful public reconstruction`，不是`source-identical reproduction`。

### Hurricane Irene 2011 历史重建与 Stage IV 对照

- [原始 C15–TCR 降雨场](validation/irene2011/irene_c15_tcr_raw_rainfall.nc)
- [Irene 运行 manifest](validation/irene2011/irene_c15_tcr_run.manifest.json)
- [Stage IV 共同支持对照 NetCDF](validation/irene2011/irene_stage4_c15_tcr_approx_fig1_bbox_comparison.nc)
- [Stage IV 对照摘要](validation/irene2011/irene_stage4_c15_tcr_approx_fig1_bbox_summary.json)
- [对照生成脚本](code/compare_irene_stage4.py)

正式服务器提交：Irene重建`6134dff`；Stage IV对照`3ba7b42`。对照严格使用2011-08-21 01 UTC至8月30日00 UTC的216个一小时累计、Natural Earth陆地掩膜、固定球面最近邻匹配及逐时300 km共同支持；没有设置pass/fail阈值。

## 下一步唯一优先任务

Irene历史重建和fixed-`r0`首事件都已完成；当前唯一任务是把同一冻结内核参数化为可重启的10,000条历史hazard sample批处理。方法冻结文件：

- [C15–CLIMADA公开重建合同](methods/C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)
- [Nature C15–TCR 复现合同](methods/NATURE_C15_TCR_REPRODUCTION_CONTRACT.md)
- [Chavas 2016 fixed-r0公开重建合同](methods/CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)

执行顺序只有三步：

1. **已完成方法与历史案例：**官方C15、无科氏伴随风桥、CLIMADA core/Petals、两个静态GeoTIFF与参数身份已经固定；Irene全生命周期与Stage IV共同支持对照已经完成。
2. **已完成首事件重算：**事件级`r0`、官方C15 `r0input`逐时RMW、TCR降雨、二维风场和道路交叠已经闭合；解释时必须同时保留Irene偏湿边界。
3. **近路全量已开跑：**1°道路占用地图经保守膨胀后，100,000条里有99,242条进入道路300 km分析域（权重合计84.28/年，接近全球84.9/年）。C15–TCR 48路可重启分片正在生产：`runs/hazard_production/lin_road_domain_300km_v1`。3条`r0<302 km`按合同留在样本中，标为`METHOD_DOMAIN_PENDING`，不重抽、不删权。此前不扩展未来情景、SFINCS或交通级联支线。

官方源代码与方法边界：

- [官方C15 v1.0原件与provenance](vendor/CLE15/10.4231_CZ4P-D448/v1.0/PROVENANCE.md)
- [官方C15 Python 3兼容层与数值测试](adapters/c15_python3/README.md)
- [C15–CLIMADA公开重建合同](methods/C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)
- [C15–TCRM来源与参数矩阵](methods/C15_TCRM_SOURCE_PARAMETER_MATRIX.md)
- [Xi 2020 Irene历史案例输入合同](methods/XI2020_HISTORICAL_CASE_INPUT_PLAN.md)
- [Emanuel官方合成事件TCR代码审计](methods/EMANUEL_SYNTHETIC_TCR_INPUT_CODE_AUDIT.md)
- [本地与服务器合成环境代码审计](methods/LOCAL_SYNTHETIC_ENVIRONMENT_CODE_AUDIT.md)

服务器已冻结的Irene公开输入：

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/data/observations/xi2020_validation/irene2011_public_inputs_v1`

对应下载器提交：`211b208`；正式manifest SHA-256为`214fff675f7747bcb30b180befc155135ecd8e5888657ddd7c02aa83898c323d`。下载器同时记录原始文件SHA与解码变量规范哈希；本次六个输入的解码内容与上一正式包逐项完全一致，五个NOAA NCSS返回文件仅HDF5容器字节发生变化，未被误判为科学数据变化。上一正式包已可逆归档在服务器同一数据层的`archive/`目录。

旧首事件产物仅作失效方法诊断归档，其中2929.97 mm累计降雨和265.17 mm h⁻¹雨强不得用于论文、道路暴露或风险结果。在历史 hazard–road 链闭合前，不批量生成未来轨迹、不启动 SFINCS、不开展交通流或级联损失支线。
