# ChatGPT 回答收割状态 — 2026-07-31（摸底 + 开工）

> 由 2026-07-31 下午 session 摸底。范围：chatgpt-answers/ 与 drops/ 下
> 与 2.1–2.7 相关的回答（P3.2 的回答由 dm 窗口 P3.2 session 负责，不在此列）。

## 总览

2.x 相关的 Q55xx/Q59xx/Q60xx 回答共 19 条。已整合进 Lean 的：5 条；
部分整合的：8 条；完全未整合的：6 条。**已整合 ≠ 缺口已关**——
部分整合的回答对应的 Lean 定理常仍以假设形式存在（详见下）。

## 完全未整合

| Q | 主题 | 对应 Lean 缺口 |
|---|---|---|
| Q6047 | P2.4 Q⁻ 六积分证书 | `alternatingQuadraticEulerTerm24_hasSum`（**正在做**，见 Problem24QuadraticAlt.lean） |
| Q6052 | P2.7 有限规范化证书（偏分式+Bernoulli） | `hBarnes : zudilinCombinedError = zudilinBarnesErrorIntegral27`（Problem27Barnes.lean:1075 仍为假设） |
| Q6024 | P2.4 辅助 Euler 积分证书（文件截断，仅 235 字节） | 无法判定对应物 |
| Q6014 | P2.4 乘积 FTC 路径（F=log(1+x)·Li2(1−x)） | 结论已另路证明（coreCrossIntegral24），路径未采纳 |
| Q6028 | P2.4 两个线性 weight-four 和的证书 | **已另路证明**（cubicLinear/alternatingCubicLinear 的 hasSum @3800/@6470）——证书本身未入库 |
| Q6033 | P2.4 积分证书（∫Li2·log(1+x)/x 链） | 目标值（2Li4(1/2)+L⁴/12−π²L²/12+(7/4)LZ3−π⁴/60）无匹配定理，疑似对应 Q+ 或某成分，未核实 |

## 部分整合（证书输出层仍缺）

| Q | 主题 | 已进 Lean | 仍缺 |
|---|---|---|---|
| Q6045 | P2.2 有限矩鞍证书 | Stein 层（rivoalWeightStein22 等，Problem22Concentration.lean） | 矩界 (A)/(B) 与 `tendsto_saddle_mean22`（Problem22.lean:1534 的 RivoalHarmonicConcentrationClaim22 仍为假设） |
| Q5972 | P2.4 五个 level-2 Euler sums | T3±、Shift 三条 HasSum（3047/3053/3059） | Q± 两条 |
| Q6035 | P2.7 Barnes 几何界 | 实变量层（zudilinMidpointValue_le27 等） | Barnes 等式本身（hBarnes 假设） |
| Q6044 | Q6035 审计 | 实现建议（真导数定义）已采纳 | 两个 gap 均未解决（打包进 hBarnes） |
| Q6039 | P2.4 交替三次 Euler sum（F） | F 闭式已证（alternatingOrdinaryHarmonicCubicTerm24_hasSum @6422） | 无（K 的 D/E 分解证书未独立形式化，非缺口） |
| Q6037 | P2.6 实变量证书链 | **P2.6 已整体完成**（problem26 @Problem26Assembly.lean:24） | 无（走 halfMap 路线，非证书链） |
| Q6038 | P2.4 半区间对数积分 | quarticCoreIntegral24/quarticCoreHalfIntegral24 已证 | 无 |
| Q6029 | P2.4 非循环证书（R=∫Li2/(2−t)） | 未核实（task-2 中断） | 未核实 |

## P2.4 主定理 `problem24Statement` 的剩余缺口（9 个假设中 4 个未证）

1. `HasSum quadraticEulerTerm24 quadraticEulerValue24`（Q+）— 无现成完整证书
2. `HasSum alternatingQuadraticEulerTerm24 alternatingQuadraticEulerValue24`（Q⁻）— Q6047 证书，**进行中**
3. `HasSum leshchinerWeightFourTerm24 ((7/8)·π⁴/90)` — 无证书
4. `HasSum bbbWeightFourTerm24 (π⁴/90)` — 无证书

已证 5 个：cubicLinear(@3800)、alternatingCubicLinear(@6470)、
pairedAlternatingHarmonic(@871)、quadraticLinear(@3047)、alternatingQuadraticLinear(@3053)。

## 当前工作

- 工作副本：本地 clone（APFS clone + rsync；canonical checkout 对本 session 只读）
- 新文件：`lean/RamanujanChallenge/Problem24QuadraticAlt.lean` — Q6047 六积分证书
  的对象定义（W0/H1/H2/I10..I22/bridgeValue）+ 目标定理占位
- 依赖：6 个线性定理 public 可用；K=quarticPlusIntegral24、B3=halfLogCubeOneSubIntegral24
  闭式已证（private，值可复制）；J3two(∫Li3/(2−t)) 与 E(∫log(1−t)log(1−t/2)/t) **不存在**，需新建
- 注意：Problem24Euler.lean 几乎全 private，跨文件只能引用 6 个公开 HasSum 定理
