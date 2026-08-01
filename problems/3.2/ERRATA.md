# ERRATA — Problem 3.2 (常驻勘误表)

制度: 每个形式化/证明项目一张常驻表 (同 UNDERSTANDING 地位)。
分类: [F]=我方 formalization/writeup 错误; [P]=外部论文/答案自身错误;
[P-intent]=字面 bug 但按意图可修。随发现随编。

| # | 日期 | 类 | 位置 | 错误 | 修正 | 状态 |
|---|------|----|------|------|------|------|
| 1 | 08-01 | [P] | cron q15/Q6375 主张 | "FFK⟹每层Mellin零点O(1)可证" 混淆 exact(ℓ-adic) 与 mod-𝔭 消失 | 拒收; 正确口径=[GAP-CARTIER]/[GAP-DCM] | 闭 (账本 P) |
| 2 | 08-01 | [P] | Q6375 rank 归属 | "F = rank-3 ₃F₂(1/3,2/3,1;1,1)" 字面识别 | Clausen 平方: ₂F₁²=₃F₂, Franel=rank2, Apéry=Sym²=rank3; [FALSE-3F2] 定理级 | 闭 (§56) |
| 3 | 08-01 | [F] | cron 附录 [CONJ-MASS-1] | "均值恰1 ⟹ mass formula 级结构" | 均值→1 = 均匀模型一阶预测, 零张力; 降级+功效分析(任何规模无功效) | 闭 (S.5) |
| 4 | 08-01 | [P] | cron q6329 初报 | 中点形式 = level-6 K3 在 CM 点 | 正确: level-8 order-4 刚性 CY (η(2z)⁴η(4z)⁴, S₄(Γ₀(8))); proof.tex 本就 level-8 无需改 | 闭 (§46) |
| 5 | 08-01 | [F] | bounded-object 定理初稿 (life) | 层归属 "φ_*Sym²E cond 20 + A_- cond 11" 二层对 | 投影公式: 两全推前字面差恒零; 正确基=本征降 A_± (Artin cond 9/11, Sym²E=11); 三源汇合(q18+codex-fm NEGATIVE-PAIR+cron codex 数值) | 闭 (§62, 已改文) |
| 6 | 08-01 | [F] | 早期口径 "truth is O(1)" | 混淆 o(p) 与一致 O(1) | 分裂: [GAP-LT-MELLIN] (o(p), sup 缓慢无界=正常随机) vs [+] 版(一致界=额外刚性) | 闭 (§54) |
| 7 | 08-01 | [P] | ChatGPT 四次编造 | 幻影 Q32_paper tex/Theorem E(1)/commit 831e147; 假 commits eda609dd/732d85a1/2336ac53 | 全部 git cat-file 证伪; 纪律=每条 repo 引用必验 | 闭 (§前期) |
| 8 | 08-01 | [F] | "唯一剩余目标=零密度" 口径 | 被 referee 判 REFORMULATION 非定量归约 | 收录口径: 归约精确但不降难度; 提交文档按此措辞 | 闭 (§58) |
| 9 | 08-01 | [F] | `CODEX_SPEC_laststand_cocycle.md` / `CODEX_SPEC_laststand_norun.md` | 把单解的 two-time companion action 错套到 `xi_r=(b_r:c_r)`，故 fixed-point、2-run eigenvector、3-run curve 全部靶错 | 正确状态是 fundamental frame；return 当且仅当 `(T_{r,d})_{21}=0`（固定 Borel hit） | 闭 (`CODEX_LASTSTAND_COCYCLE_report.md` §1) |
| 10 | 08-01 | [F] | laststand `[FR_eta]` 口径 | 把 shell incidence `sum_d |Z_d|`、dyadic union、primitive count 称为等价；且字面 `L -> infinity` 未写 slow/nonempty-window 限制 | 分开记 `S_B` 与 `U_B`，仅有 `U_B <= S_B`；能量链使用 `L=N^{o(1)}`, `D=o(N)` | 闭 (`CODEX_LASTSTAND_COCYCLE_report.md` §1.5) |
| 11 | 08-01 | [P] | Q6730 “honest poles remain nonzero mod p” | 整数非零误推工作素数下非零；`N_3(-3)=584=8*73`，且 `gcd_F73(N_3,N_4)=X+3` | fixed-p correlation 必须先 saturation；不能直接用特征零 pole/resultant | 闭 (`CODEX_LASTSTAND_COCYCLE_report.md` §3.4) |
| 12 | 08-01 | [F] | AT.71 / 初版 `CODEX_MARGINS_report.md` | 只检查每个有限高度的 far minimum 为正，就误报成随 `h` 一致的 FAR 下界，首版写 `TWO-REGIME CONFIRMED` | 固定 `k=1` band 实测按 `h^-1.088` 衰减；`h=52` 另有宏观 inner/outer 近碰撞；判决改为 `TWO-REGIME REFUTED` | 闭 (`bfc611c`) |

## [F] 2026-08-02 — §134 [NO-RUN] 推导整体作废 (codex 驳倒, 活反例)

- **错误**: §134/DOCTRINE LS-b/CTX_LASTSTAND 声称 ξ_{r+1}=M_r ξ_r (ξ_r=(b_r:c_r))
  及其推论 (2-run ⟹ ξ_r=v(r,d); 3-run ⟹ F(r,d)=0 ⟹ 3-run 总数 O(D))。
- **错因**: 射影对象混淆。伴随阵 M_n 传播**单解**二时刻状态 (u_{n−1},u_n)
  (即 b-ratio walk η_r=(b_{r−1}:b_r)), 不传播**双解行** (b_n:c_n)。
- **反例** (CODEX_NORUN_report.md, 机器验证): p=997, r=248,249∈Z_182, c_248≠0,
  ξ_248 key=758 ≠ v(248,182) key=798。§134 的"实测 2-run 吻合"从未真正逐点
  验证 (脚本 placeholder 未跑) — 流程教训已提取。
- **幸存部分**: (a) F=−24d²G 符号证书对 **ratio-walk collision 族**有效 (另一族,
  非能量对象); (b) [GAP-ONE-CASORATIAN] W_r=6/(r+1)³ 被 codex 独立证明;
  §135 strip-分离/Heilbronn 归约链不依赖 run-rigidity, 完好。
- **新正果** (codex, 同报告): U_r=[u_{r−1};u_r] 可逆 (det=W≠0) ⟹
  **collision ⟺ (T_{r,d})_{21}=0** — 碰撞条件从头就是轨道无关多项式恒等
  (=continuant N_d(r)); 正确的 2-run 问题 = p | Res_x(N_d(x),N_d(x+1))。
- 波及: Q6733/Q6735/Q6742 等快答中所有引用 run-rigidity 之处作废;
  CTX_LASTSTAND.txt 已修正; q63 发出的 "NO-RUN certified" 状态更新作废。

## [F] 2026-08-02 — §135 "FR_η ⟸ 单变量 Weyl 幂省" 归约命题作废 (自审判死)

- **错误**: §135/§137 声称 |σ(t,s)|≤p^{1−δ} (任意/门槛 δ) ⟹ [FR_η]。
- **判死推导** (自算, ledger §138): strip 展开 (Dirichlet/Fejér) 是恒等式;
  max-route 需 (1/p)·p·(Σ|ĝ|)·p^{2−2δ} ≤ D^{2−η} ⟺ **δ ≥ 1/2+η/4** — 超过
  平方根; 而 Weyl 真值恰 δ=1/2 (实测 4√p) ⟹ max-Weyl 类输入原则性不足。
  矩路线 (Parseval over t) 精确回代 T 自身 — 循环。合成重纤维测试 (§137)
  是同一漏洞的实例化。
- **幸存**: Weyl 实测校准 (两相位皆随机地板); 双线性分离**恒等式**本身;
  调和玩具 strip 靶文献空缺的校准; 新死亡证书 [DEAD-MAX-WEYL-CLASS]:
  任何只用单变量 max-Weyl (乃至平方根最优) 的方法差 η/4 幂, 不可闭合
  [FR_η] — 缺的是 (值,位差) 联合二维相消, 即对象本身。
- 波及: §135 归约表述、q73/q74/q75/q76 brief 中 "would imply" 措辞、
  TG 报告 (已发更正)。q76/Q6778 audit 落地后对表本判。

## [制度] 2026-08-02 — 验收门制度 (因本日五次改口而立)

**背景**: 2026-08-02 life 战役单日五次公开改口 ([NO-RUN] / §135 Weyl 归约 /
∞-纤维洞 / [FR_η] 的 S_B-U_B 口径 / [AVG-ZERO] 余量口径)。其中至少两次
**手上已有验证手段却未执行**——[NO-RUN] 的验证脚本里逐点核对是 placeholder,
未跑即入账并对外报告, 后被 codex 用活反例 (p=997,r=248,d=182) 推翻。
根因不是数学难度, 是把**探索期的每次摆动当作进展发布**。

**制度 (自即日起, 本项目所有 lane 适用)**:
1. 任何陈述进入账本"已证"栏之前, 必须在 `LASTSTAND_VERIFY.py` 里有**一个对应
   gate**, 且该 gate 实跑通过。没有 gate 的陈述一律标 provisional, 不得引用。
2. gate 必须与被验陈述**独立**: 用独立的 ground truth (如精确二项式和 vs 递推),
   多参数 (多个 p / 多个 X), 不接受单点验证。
3. **判死也要证据**——反例或完整簿记, 不接受"看起来不行"。
4. 对外报告 (TG) 只在 gate 全绿后发; 探索期的 bank-then-reverse 留在 commit
   与账本里, 不发布。
5. 验证脚本中出现 placeholder / 跳过分支 = 该 gate 视为**未跑**。
