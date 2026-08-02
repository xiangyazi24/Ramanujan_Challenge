# CRON DOCTRINE — 战役三最后一击：family 相容性定理（2026-08-01）

**主目标（一句话）**：证出（或严格判死）"固定尺度免费/增长尺度缺失"的那个 family 相容性定理——
使 (R)+(C) ⟹ GARQI-1′ ⟹ 一阶强度闭合的链条通电；cron 侧主攻数学重构与最小实例，
life 侧主攻层论（B 包二参数转移层），对表在 CAMPAIGN_MAP §7。

## Avenues（按初始把握排序）

### (a) 两点函数重构【主攻】
核心观察：N_h(r) 不是"无结构增长度数家族"——它是**同一个固定 rank-2 Apéry crystal
在两点 r, r+h 之间平行输运的矩阵元**。W(r,h):=N_h(r) 是**固定二元 holonomic 系统**的解：
- h 方向：w_{h+1}=A_h(r)w_h−B_h(r)w_{h−1}（transfer 递推，已知）
- r 方向：continuant 底部递推 N_h(r) ~ N_{h−1}(r+1)（待机器验证精确形式）
Q6463 §2.6(iv) 的反对（"非固定二元多项式的特化"）是对的，但对象是**固定二元
holonomic 系统的解族**——比任意族强得多的结构。此外 Lucas/Gessel digit 结构 =
h 方向的 Frobenius 压缩机制（已证定理）。
**攻击序列**：A1 机器验证 r-方向递推（2D holonomic 结构坐实）→ A2 成文"固定 crystal
两点等分布"精确命题 + holonomic-family BV 命题 → A3 SOL 派单（qK）：固定对象两点框架
能否绕过增长导手障碍；Frobenius/digit 压缩能否给出跨 h 相容性 → A4 与 life 二参数层对表。
**终止条件**：命题被证 / 被反例杀死 / 被归约到一个已知 open 问题（写明是哪个）。

### (b) Colored BDH 方差正面强攻
Q6491 次优路线。非对角 = 双返回约束下两转移态的联合分布。
**攻击序列**：B1 先测经验：跨 (h₁,h₂) 根集相关性（N̂_{h1}, N̂_{h2} mod p 根的联合统计，
机器直算）——若经验独立（相关 →0），方差路线的对角主导有数据支撑 → B2 成文精确方差和 +
对角计算 → B3 SOL 派单（qL）：带数据的 off-diagonal 估计攻坚。
**终止条件**：非对角有界性被证 / 被数据反驳（发现强相关结构 → 该结构本身入账）/
归约到 (a) 的同一命题。

### (c) 小 h 色饱和 + polylog 桥
(R) 无条件对 h ≤ polylog 可用（固定 h Chebotarev + 平凡求和）。缺口在 h 大。
**攻击**：C1 精确写出"h ≤ log²X 无条件 + h > log²X 需色"的分界账；C2 检查大 h 端
有没有被忽略的免费结构（如 h 接近 p−1 时 N_h 与 N_{p−1−...} 的反射/互补恒等式——
镜像之外还有没有第二条代数压缩）。
**终止条件**：分界账成文入账（这是记账 avenue，必成）。

## Fallbacks
- 全部死：把三张终局判决 + 归约图写成"family 相容性定理 = 唯一剩余公开问题"的
  正式 open-problem 文档（这本身是战役三的合法交付物——问题定形即胜利的一半）。
- qI（M_h 不可约性）落地若给出全 h 定理 → 直接强化 (a)/(b) 的 κ=O(1) 输入。

## 协同边界
- life：Q6462 + marked-dec 在飞；二参数层归 B 包（他们的工地）。我不碰层构造，
  只供数学重构 + 数据 + 最小实例。
- 单写者：本文件 + CRON_RUN_LOG.md + 账本附录归 cron。

## SEGMENT 5 DOCTRINE (2026-08-01 16:1x — 18:00 交卷冲刺)

Goal: 破 3/2 能量墙。墙谱已定格（AT.38）：[BDH-LAG] ≤ [SAME-LAG-L2] ≤ [ZERO-TAIL-2] ≤ [MESO-TOTAL]=W1；组合路线三重判死 ⟹ 纯算术。

Avenues (ranked):
(a) [主攻] [BDH-LAG]/[SAME-LAG-L2] 多路并击：
    (a1) 逐点反集中 [PT-ANTICONC]: R_h ≤ N^{1/4−δ} ⟹ 全部（新蒸馏面，经验 max=8 恒定）— cron2
    (a2) shifted-gcd/结式 ρ_h(T) 散布 — cron9 首答 PARTIAL（deg=d², ord₀=d, 缺失命题=固定素数结式反集中）
    (a3) [ZERO-TAIL-2] 切片不交双计数（跨 lag 结式 + 全对互素引理）— cron6 + cron7(全对互素 descent)
    (a4) fiber/color 同 lag 交换律 + 分层图 trace — cron3
    (a5) 镜像骨架 κ_h：偶 h 中心根 ⟹ κ_h ≥ 1_{h even} + R_h 奇偶律（新发现候选，机验中）— cron10 + 本地 sympy
    (a6) 文献移植（Sawin-Shusterman/Katz/连分数家族）— cron5
    终结条件: 任一 W 面证明（成功）/ 各向量书面死亡证书 + 更弱残差定名。
(b) [CRIT-2H] all-h: Q6697 共轭引理完整证明（cron1 在飞）→ K_∞ 分支表 → h₀≤40 ⟹ 有限计算闭链。
(c) 证书护航: 度 4/5/6 PTE 归约 all-h（cron4 在飞，四元组等前三幂和分拆判死路线）; A_h/W_h 归 life 侧。
(d) codex gpt-5.6-sol ULTRA 档（爸爸指令试新档）单挑 breakwall-2（spec: CODEX_SPEC_CRON_breakwall2.md）。
敌审记录: [NO-RIGHT-2-3] + 反模型 Thm5.1 两票 CONFIRMED（Q6709/Q6713，快 tab 档）。
Poisson 判决（Q6714）: max=8 恒定拒绝 iid Poisson @ 7e-5 ⟹ 存在精确机制（奇偶律候选）。
17:4x 收口: 全部 commit+push，账本 AT.40 起 append。
