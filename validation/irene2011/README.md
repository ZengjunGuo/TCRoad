# Hurricane Irene 2011 C15–TCR public reconstruction

本目录保存 Hurricane Irene 2011 的公开方法级 C15–TCR 重建与 NCEP Stage IV 对照产物。它们用于验证统一 C15–TCR 数值链的时间、空间和量级语义，不声称逐字节复现 Xi et al. (2020) 的未公开生产代码或 Fig. 1 精确累计窗口。

## 文件

- `irene_c15_tcr_raw_rainfall.nc`：217个逐小时节点、事件累计和最大24小时累计的原始C15–TCR降雨场。
- `irene_c15_tcr_run.manifest.json`：输入、代码、参数、版本和输出哈希。
- `irene_stage4_c15_tcr_approx_fig1_bbox_comparison.nc`：Stage IV与TCR在共同支持上的配对figure data。
- `irene_stage4_c15_tcr_approx_fig1_bbox_summary.json`：描述性统计与全部输入、软件、输出provenance。

## 固定比较口径

- 216个一小时累计结束时刻：2011-08-21 01 UTC至2011-08-30 00 UTC；
- 33–39°N、81–73°W近似论文视觉框；
- Natural Earth 1:110m陆地掩膜；
- Stage IV原生格点到TCR 0.05°质心的固定球面最近邻，无插值和平滑；
- 每小时距当时TC中心不超过300 km的共同物理支持；
- 无pass/fail阈值，无观测偏差校正。

6,881个正式陆地比较格点均有Stage IV 216/216小时完整覆盖。TCR与Stage IV累计的空间相关系数为0.7156，平均累计分别为175.09和101.78 mm，平均偏差为+73.31 mm，第95百分位为494.22和264.44 mm，最大值为668.09和532.06 mm。

科学解释：统一公开重建已捕捉到明确的历史空间结构，但绝对累计偏湿。该边界必须随首个合成事件和后续批量结果一并保留；它不授权调参、平滑RMW或加入自创的降雨修正。

## SHA-256

- `irene_c15_tcr_raw_rainfall.nc`: `54b59c3f73c4155aba2b53baec157f04c8e33edd6ffdff6cb9da9882c5920643`
- `irene_c15_tcr_run.manifest.json`: `f86603eedb40dd59a0b50d57928e568b09db15c4471dda32e93222a7ae290695`
- `irene_stage4_c15_tcr_approx_fig1_bbox_comparison.nc`: `47a5ed197182cc3e471f7ca09044797d5e8c9e74e97e5e942fb927f93ea63910`
- `irene_stage4_c15_tcr_approx_fig1_bbox_summary.json`: `a8f560b0ccfd5305bbab70aceb3ed371989f0b8609d95426799f74b0af054800`

正式服务器代码提交：`3ba7b426974c3cab936a2b84f0a30f41343f37dc`。
