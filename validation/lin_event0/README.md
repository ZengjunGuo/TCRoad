# Lin synthetic event 0: superseded Knaff-interface regression archive

状态：`SUPERSEDED_BY_CHAVAS2016_FIXED_R0`。本目录只保留旧接口的数值回归与审计证据；下列结果不得作为Lin正式hazard production、fixed-`r0`首事件结果或Nature/Gori方法复现引用。

事件：`stream0000-year1995-track000002`  
气候输入：MPI-ESM1-2-LR / r1i1p1f1 / historical / 1995–2014  
历史方法口径：官方 C15 v1.0 + CLIMADA core v6.1.0 + Petals v6.2.0；94 个原生 1 h 节点；0.05°移动 300 km 计算域；所有时次统一使用 Knaff et al. (2015) Eq. (1) RMW和官方C15 `rmaxinput`。该闭合现已被事件级固定`r0`与官方C15 `r0input`整体取代。

这是旧Knaff接口下首个Lin合成事件的风—雨—道路闭环。它证明TCR、二维风场和道路叠加数值链能够运行，但不证明现行尺度方法。现行生产按[Chavas 2016 fixed-r0公开重建合同](../../methods/CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md)：每个事件只抽一次生命周期固定`r0`，逐时以`v_trks+r0+abs(f)`调用官方C15 `r0input`内生求RMW。下列数值只用于回归比较；Irene–Stage IV历史对照的绝对量偏湿边界仍须用于解释重算结果。

## 核验结果

- 网格点：48,858。
- 最大逐小时雨率：91.2433 mm h⁻¹。
- 最大事件累计：907.4122 mm。
- 最大 24 h 累计：907.4122 mm。
- 全部数组有限且非负；300 km 之外无正降雨。
- 最大逐小时雨率距当时中心109.43 km；未触发TCR垂直速度上限。
- 最大累计格位于墨西哥高地，约23个连续有雨小时构成全部累计；不是94小时机械重复累加。
- 最大模型原生近地面风：30.4665 m s⁻¹；峰值距中心62.92 km，位于逐时RMW环附近。
- 风场平均时段在Lin `v_trks`中未记录；没有应用0.893换算，也不标为10-min持续风。
- 共同移动300 km评价域内匹配12,222个0.1°道路格中心，其中6,245格含道路；五级道路总长669,542.60 km。
- 道路叠加没有灾害阈值、脆弱性函数、失效规则或货币损失。

## 文件与SHA-256

| 文件 | SHA-256 |
|---|---|
| `lin_event0_c15_tcr_raw_rainfall.nc` | `f2abd54e10ae47f52f7011b8c48be9a49a4c7f5369e123074319b5a2a9b134cf` |
| `lin_event0_c15_tcr_run.manifest.json` | `65aef2d81c986a0206feed146123c44e666ed910a3dd91ad4c0564d79b4dbcf8` |
| `lin_event0_moving_300km_grid.nc` | `d219e760c8135df397b2f310cf6a7f60724bf91a70d0f32068b49c4dc31afbb3` |
| `lin_event0_one_hourly_climada_track.nc` | `4122de095b65fd24465bac4b2565015a3907f9a6f41b04fea80f65cf87d5cb23` |
| `lin_event0_public_inputs_prepare_only.manifest.json` | `fc2304751631c3d4960e75777c80a7f20533b10c79cf86fd0f992d549b28b1e3` |
| `wind/lin_event0_c15_lin_chavas_windfield.nc` | `141efe2487eb5e237265eea51197f6562c78c98ad4406382adbf68e3bff36332` |
| `wind/lin_event0_c15_windfield.manifest.json` | `8f894cbf5091c70b6d39c2f442ce847d878d6b74b2b29930b0a0e1316255bc8a` |
| `road_overlap/lin_event0_road_grid_joint_exposure.nc` | `9332c85191417ac7b02d36360881d81a7323527ae9bf13c48f801f13e33c3fc6` |
| `road_overlap/lin_event0_road_class_joint_exposure_summary.json` | `d6489f930ba298664d70caae7f3ae5b972d40e31c31cd16c4fd258d954c4c11b` |
| `road_overlap/lin_event0_road_class_joint_exposure_summary.csv` | `5449a3fcb08dcf92986242c6d1fce57c00a539ec55d22facc206177fd6dfae31` |
| `road_overlap/lin_event0_road_overlap.manifest.json` | `9f17e2eeeec19cf24a6c573f61349064f5e1cc1fa3376fa1c717a12f2cb356d2` |

服务器审计归档产物：

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/runs/hazard_production/lin_event0_c15_climada_public_reconstruction_v1`
