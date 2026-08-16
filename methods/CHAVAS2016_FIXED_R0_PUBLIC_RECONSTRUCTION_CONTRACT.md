# Chavas 2016 fixed-r0 公开重建合同

状态：`FROZEN_FOR_IMPLEMENTATION_V1`  
日期：2026-08-13

## 1. 方法选择

Lin 合成事件的正式尺度闭合采用 Gori et al. (2022) 及 Xi, Lin & Gori
(2023) 的方法架构：每个热带气旋抽取一个事件级外半径 `r0`，该值在整个
生命周期内保持不变；每个时次再以 Lin 的方位平均最大风 `v_trks`、固定
`r0` 和风暴中心科氏参数绝对值调用官方 C15 `r0input`，由 C15 内生求得
逐时 RMW 与完整径向风廓线。

这条路径整体取代逐时 Knaff/Nederhoff RMW。二者不得拼接、补洞或作为
fallback；项目既有 outer-first 尺度模型也不进入正式 hazard production。

## 2. 公开分布

Gori et al. (2022) 只说明 `r0` 来自 Chavas et al. (2016) 的经验
lognormal 分布，未公开其拟合参数、随机种子、随机数生成器、global/NA
选择或 event-to-draw 映射。因此本项目不声称恢复作者私有 storm set，而
冻结以下“公开汇总统计约束的透明再实现”：

```text
ln(r0 / km) ~ Normal(mu, sigma^2)
mu    = ln(881.0) = 6.78105762593618
sigma = [ln(1054.4) - ln(740.7)] / [2 Phi^-1(0.75)]
      = 0.261776756893889
```

参数精确匹配 Chavas et al. (2016) Table 1 全球样本 `N=578` 的中位数、
第一和第三四分位数 `881.0/740.7/1054.4 km`。该分布隐含均值
`911.71 km`、标准差 `242.81 km`，与论文报告的 `909.4 km` 和
`248.5 km` 分别相差约 `+0.25%` 和 `-2.29%`。

正式目录采用全球 pooled 分布，不以北大西洋分布替代，也不按盆地拼接。
原因是目标目录为全球事件，而 Gori 未公开盆地条件抽样规则；北印度洋样本
仅 `N=6`，南大西洋无单独观测样本。南大西洋事件因此明确属于全球 pooled
先验的空间外推，不另造一套经验关系。

## 3. 随机与目录合同

- 随机数生成器：NumPy `PCG64`。
- seed：`20260810`，沿用已冻结 10,000-event sample 的项目种子。
- 抽法：`Z = Generator(PCG64(seed)).standard_normal(10000)`，再
  `r0_km = exp(mu + sigma * Z)`。禁止调用 `Generator.lognormal`：该助手函数
  跨 NumPy 版本不保证同一 seed 得到同一字节序列。`standard_normal` 序列在
  本机 NumPy 2.3.4 与服务器 NumPy 2.2.6 上按位一致，SHA-256
  `c503c52f7f84a93f90b3d85fc8dddf5eb2e244652cfdbf6b2f99b83229abdf10`。
- 存储量子：将每个 `r0` 四舍五入到最近 1 mm（ties to even）。这只消除
  `exp` 的跨平台 libm 舍入，不截断、不 clip、不重抽。量化后
  `outer_radius_m` 序列 SHA-256
  `be647f18bb05dafd0bf54e1343ad1e52a4e6fa7f39921faced5a22869a8eef71`。
- draw 顺序：冻结 sample 中的 `event_position = 0..9999`。
- 每个 event 只抽一次；运行、分片、并发和重启时只读取冻结目录，不再次抽样。
- 不截断、不 clip、不重抽、不按强度或盆地改变 draw，不修改事件权重。
- NetCDF 必须保存 `event_position`、`event_id`、`outer_radius_m`；manifest 必须
  保存 sample SHA、分布参数、RNG/seed、标准正态序列哈希、量化后目录哈希、
  NumPy 版本、目录及脚本 SHA、实际 `<302 km` 计数和完整分位数摘要。

seed 是实现层复现元数据，不是新的物理参数，也不支持“与 Gori 原随机序列
一致”的主张。

## 4. C15 接口

每个时次调用官方：

```text
ER11E04_nondim_r0input(v_trks, r0_event, abs(f), ...)
```

冻结的 Python 适配器返回顺序为：

```text
rr, VV, rmerge, Vmerge, rmax
```

因此 `rmax` 必须从第五个返回值读取。官方数值网格若止于 `r0` 内侧，可追加
C15 定义点 `(r0, 0)`；不得追加尾部、taper 或负风外推。`vmax_trks` 不参与
尺度闭合；它继续用于事件阈值和目录语义。C15 强度与 Emanuel `qs900b`
暖核湿度均使用 `v_trks`。

## 5. 300-km 域与失败语义

TCR 仍在逐时随中心移动的 300-km 支持域运行。不得把支持域缩成
`min(300 km, r0)`，因为平移风、环境风和地形机制并不随 C15 旋转涡旋在
`r0` 消失。

公开 C15 只定义 `0 <= r <= r0`，而当前公开 TCR 桥的无科氏伴随风在
`r0` 外没有作者公开实现。因此生产前必须逐项验证冻结目录满足：

```text
r0_event > 302 km
```

其中额外 2 km 覆盖 TCR 径向差分。若实际 draw 违反该条件，事件进入
`METHOD_DOMAIN_PENDING`：不重抽、不 clip、不缩域、不排除、不改权重。

所有逐时 C15 解还必须有限、RMW 为正且 `RMW < r0`。任何失败均触发审计，
不得静默切换回 Knaff、Nederhoff 或 ER11。

## 6. 证据与主张边界

- Gori et al. (2022), *Nature Climate Change*：每个 TC 从经验 lognormal
  抽取 outer radius、生命周期固定，并以强度和 outer radius 通过完整风廓线
  模型估算 RMW。DOI: https://doi.org/10.1038/s41558-021-01272-7
- Xi, Lin & Gori (2023), *Nature Climate Change*：沿用 Gori storm sets，
  C15 同时服务风害和 TCR 风输入。DOI: https://doi.org/10.1038/s41558-023-01595-7
- Chavas et al. (2016)：全球 `r0` 分布及 Table 1 汇总统计。DOI:
  https://doi.org/10.1175/JCLI-D-15-0731.1
- 官方 QSCAT-R 数据集：https://verif.rap.ucar.edu/tcdata/quikscat/dataset/
- 官方 C15 v1.0：https://doi.org/10.4231/CZ4P-D448

准确表述为：`source-constrained, method-faithful public reconstruction`。
禁止表述为：Gori/Xi 私有生成代码、精确 MLE、原随机序列或 bitwise storm-set
复现。
