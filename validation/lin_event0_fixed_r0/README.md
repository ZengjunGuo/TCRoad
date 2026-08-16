# Lin synthetic event 0: Gori/Chavas fixed-r0 public reconstruction

状态：`FIRST_EVENT_CLOSED_FIXED_R0_V1`。这是现行正式方法下的首个合成事件风—雨—道路闭环，取代 Knaff `rmaxinput` 版本。旧目录 `../lin_event0/` 只作回归审计。

事件：`stream0000-year1995-track000002`  
气候输入：MPI-ESM1-2-LR / r1i1p1f1 / historical / 1995–2014  
方法：官方 C15 v1.0 `r0input` + CLIMADA core v6.1.0 + Petals v6.2.0；事件级固定 `r0 = 1000.708405 km`（Chavas 2016 全球 lognormal，PCG64 seed `20260810`，存到最近 1 mm）；94 个原生 1 h 节点；0.05° 移动 300 km 计算域。

## 核验结果

- 网格点：48,858。
- 逐时 RMW：62.58–256.31 km，全部有限且 `RMW < r0`。
- 最大逐小时雨率：91.9089 mm h⁻¹，出现在距当时中心 128.8 km，与当时 RMW 128.3 km 重合；不是外缘伪雨带。
- 最大事件累计：959.2935 mm；最大 24 h 累计：957.2990 mm。
- 全部数组有限且非负。
- 最大模型原生近地面风：30.4716 m s⁻¹。风场平均时段在 Lin `v_trks` 中未记录，没有应用 0.893 换算，也不标为 10-min 持续风。
- 共同移动 300 km 评价域内匹配 12,222 个 0.1° 道路格中心。道路叠加没有灾害阈值、脆弱性函数、失效规则或货币损失。

与已作废接口比较：自定义 CLE–pyTCR 曾给出 2929.97 mm / 265 mm h⁻¹ 外缘伪雨，已退出生产。Knaff 版同事件为 907.41 mm / 91.24 mm h⁻¹ / 30.47 m s⁻¹，证明数值链可跑，但尺度方法不正确。现行结果与 Knaff 版同量级，峰值雨率落在 RMW 环，符合 C15–TCR 结构。Irene–Stage IV 对照的绝对量偏湿边界仍然适用，不得把 959 mm 当成观测真值。

## 文件与 SHA-256

| 文件 | SHA-256 |
|---|---|
| `lin_event0_public_inputs_prepare_only.nc` | `45e16d3cbcffbe4d65d336852444eaef66aaac8d720f2b588d6922a7d0941937` |
| `lin_event0_public_inputs_prepare_only.manifest.json` | `8f5d5adeba92c00e606877bc4340bc378a530751c5c1c111ada20327326027ee` |
| `lin_event0_c15_tcr_raw_rainfall.nc` | `e5bae280e683293fb6c2996009d1b949387bddee5fe81b6ccab3bc50799371a5` |
| `lin_event0_c15_tcr_run.manifest.json` | `2e3c6a02ba656ee5c41680ec622db3ff6232e5cef21b3811e698c1571ca785b7` |
| `lin_event0_moving_300km_grid.nc` | `d219e760c8135df397b2f310cf6a7f60724bf91a70d0f32068b49c4dc31afbb3` |
| `lin_event0_one_hourly_climada_track.nc` | `0bb3afc8976250f327af2da904e563a7a1f41f87b91f1f578592bd8e7d6d7c75` |
| `wind/lin_event0_c15_lin_chavas_windfield.nc` | `b99ea35845c725bbf15759b1d08fd31d053e8add4650f3bdd1e61d2a7193d47f` |
| `wind/lin_event0_c15_windfield.manifest.json` | `16854c7f9ccc12ce27511ce120c6c40ff16cff9ce359771ada8071abc0fbbb08` |
| `road_overlap/lin_event0_road_grid_joint_exposure.nc` | `3b481edd6eaaf5b73b2d3c47f261782ebe8c56941d68c2c30d9d0d3a9604fdd0` |
| `road_overlap/lin_event0_road_class_joint_exposure_summary.json` | `01536b9a903ed8f361d8d6240e115d952505636e35d97067e5ae27871469bdd6` |
| `road_overlap/lin_event0_road_class_joint_exposure_summary.csv` | `5f1e79a08f81d9e5e4b58eea5303447291c11f8cc4819090dc6d349b74a92579` |
| `road_overlap/lin_event0_road_overlap.manifest.json` | `342a668bbd4d33327b0d3873150f6b59ce19d049413130714e92fab721619919` |
| `fixed_r0_catalogue_v1.manifest.json` | `23bd306536df38367a75ea6e7b34437d2d9bd77f24cc52647f53eeebf29ffd10` |

服务器正式产物：

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/runs/hazard_production/lin_event0_c15_fixed_r0_public_reconstruction_v1`

正式 10,000 事件 `r0` 目录：

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/data/lin/samples/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/hazard_pilot_10000/fixed_r0_catalogue_v1.nc`  
NetCDF SHA-256 `14e633ad5f7a9ac3c0b5676bdb45b4b3a89b1a8e2f7ed64ddea1e25c7e9daf6d`
