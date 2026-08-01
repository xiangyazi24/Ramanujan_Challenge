# DOCTRINE — 最后一搏: fully unconditional (2026-08-02, life, fresh eyes)

## 主目标 (一句话)

无条件证明能量猜想的任一充分靶: 对某个无界 L(p), 证 [MESO-TOTAL] S_D≪N
(D=√N·L) 或 (4.9) P_D≪N/(L²logD) 或 [ANTI-COEXIST] — 即破 3/2; 乃至 E1。

## Fresh-eyes 核心重述 (自推, 本 session 独立导出)

Δ_{r,h}≡0 ⟺ ξ_{r+h}=ξ_r, 其中 ξ_r=(b_r:c_r)∈P¹(F_p) 是射影 Apéry 轨道
(Casoratian 恒等式)。⟹ collision 关系是**等价关系**, collision graph = 纤维团
(union of cliques)。S_D/P_D/Q_D 全由纤维位置集决定。ξ 由非自治 Möbius 走动
ξ_{r+1}=M_r·ξ_r 驱动 (三项递推转移阵)。
- 纤维内 k-AP (公差 d) ⟺ 连续同长转移积 T_{r,d}, T_{r+d,d}, … 共享不动点 —
  强代数刚性候选, 对表已证 self-twist rigidity / Sp-full / honest poles。
- 预算恒等: C_g ≤ 3(g−1) 全局 (deg N_g); Σ_{g≤D} 3(g−1)≈(3/2)D²=(3/2)NL²,
  距 [MESO-TOTAL] 只差 L² 因子 — 要的只是族平均省 L², 不是逐个 poly 的根数控制。
- Strike2 padded word 证明: 词级+团级输入**不可能**推出 (4.8)/(4.9)。必须新增
  Apéry-特定算术。未入词级约束表的已证算术: apparition 四分律 mod 24, parity law,
  quarter-value eta-乘子, addition law (双线性), [PER-H-WEIL-4H-1], self-twist
  rigidity, honest triple poles。**其中哪条杀 padded word = 缺失机制的名字。**

## Avenues (fresh 编号 LS-a…, 按初始赢面排序)

### LS-a (修订版, §7 check-existing 后): [FR_η] 一阶幂省 + 未审计算术
Strike2 §7 已审计九条 banked 输入: **无一**能排除 padding (各有精确死因), 缺失
输入被命名为两形态: (1) **[FR_η] P_D≪D^{2−η}** (4.11, 任意 η>0 即破 3/2, 等价
dyadic 版 (4.13) |∪_{D/2<d≤D}Z_d|≪D^{2−η}); (2) 共存 P_{D,1}·Q_D≪N² (7.1)。
§7 表**未含**: quarter-value eta-乘子, self-twist rigidity, honest triple poles/
四分律 mod 24 — 三条未审计。任务: 补审这三条 vs padding + 直攻 [FR_η]。
- 关键新对偶 (本 session 自推, 待验证): **Δ_{r,·} 作为 h 的函数是 Apéry 递推
  的解** (两解线性组合), 故 h-方向每行有 continuant/restart 型 2/3 界 — 与
  d-方向 deg 界 3(d−1) 形成 transpose 双界; 单用任一侧都不够 (N·D^{2/3}≫D²),
  问题=能否 dyadic-block 双计数插值。
- 终止 (成功): [FR_η] 任一 η>0 或 (7.1) 证明骨架。
- 终止 (失败): 三条补审全不咬 + transpose 插值判死 (写明为何) ⟹ 转 LS-b 主攻。

### LS-b: 【已判死 2026-08-02, 见 ERRATA】转移积不动点/run 刚性
**整条 avenue 建立在 ξ_{r+1}=M_rξ_r 之上, 该式为假 (射影对象混淆, codex 活反例
p=997/r=248/d=182)。** 正确重述: collision ⟺ (T_{r,d})₂₁=0 (轨道无关);
2-run@d ⟺ p|Res_x(N_d(x),N_d(x+1))。终判: 此 avenue 关闭, 主攻转 LS-e。

### LS-e: strip-分离 ⟹ 单变量 Weyl 幂省 (Heilbronn 类) 【§135, 现主攻】
[FR_η] ⟸ |Σ_r e_p(t·φ_r+sr)| ≤ p^{1−δ} (φ=调和玩具 H_r 或 Apéry 卡 u_r)。
玩具障碍: 标准 Stepanov 死于极点增长 (Q6761); 反制在飞: mod-p² 提升/Γ_p
分布关系 (q75), Stepanov-on-coincidence-variety (q74, 需 addition law 低次
闭包定理; exp3 已判碰撞簇无镜像线外低次分量)。codex ultra 班并行主攻。

### (原 LS-b 正文, 已作废存档)
自推链: r,r+1∈Z_d ⟺ ξ_{r+1+d}∥M_{r+d}ξ_r 且 ξ_{r+1}=M_rξ_r ⟺ ξ_r 是
M_r^{−1}M_{r+d} 的特征向量 (M_s=显式 Apéry 转移阵)。⟹ Z_d 内长度-k run 强制
特征向量链 = k 条有界次代数条件。padded word 的 Z_{d_j}∩block=**整区间** —
run-rigidity 若证 "Z_d 无长 run" 即杀该构造; 残余敌人=scattered pairing (块内
变 gap 配对), 需二阶矩/gap-多重性反集中补刀。破 3/2 ⟸ 排除"近满配对 D-块"
占正比例位置 (本 session 块局部重述)。数值先行: 实测 Z_d 最长 run / 纤维最长
窗口 AP / 满配对块统计 (p→10⁵)。
- 终止 (成功): "Z_d run ≤C" + "满配对块稀疏" ⟹ [FR_η] 或直接 (7.1)。
- 终止 (失败): 特征向量链可满足 (显式见证) + scattered 不可排 — 记录, 转 LS-c。

### LS-c: d≤32 精确 census 账 (handoff 活口2, 没人算过)
Sp-full h≤32 无条件 ⟹ 该段的等分布/矩输入。算清: 这段对 [MESO-S1-2/3] /
near-wall K_p(H,D) / [TR6] 水平线各买到什么。若 bounded-d 段+尾部预算的组合
给出任何 S₁(D)≪N^{2/3} 型链, 立即接 [TRIANGLE-DICTIONARY]⟹MESO-PAIR。
- 终止 (成功): 明确的非平凡 census 定理 + 链条。
- 终止 (失败): 账算清后 O(1)/无关 — 记账, 关此活口。

### LS-d: [MESO-TOTAL] 最小机制逆向工程
从目标反推: S_D≪N 等价于纤维位置集的"D-窗口平均邻居数 O(1)"。已有: restart
Z(H)≤1+(3/2)H^{2/3}, 窗口 2/3 律, 团结构。问: 加什么最弱的种子陈述即闭合?
(候选: 纤维位置集的 additive energy 控制 / gap 分布反集中)。把种子陈述精确化后
回喂 LS-a/b 检验其可证性。
- 终止: 种子陈述清单 + 各自可证性判决。

### Fallbacks
- 全部失败 ⟹ 死亡证书补全 + 五路刻画 program 文档升级 (已是可发表认识, §133)。
- 任何 avenue 中途发现新活口 ⟹ 追加 LS-e… 不请示。

## 并行算力配置

- codex 双开 (gpt-5.6-sol, ultra 档若存在否则 max/high): #1=LS-a 机械测试班,
  #2=LS-b 刚性推导+数值班。
- ChatGPT 9-tab 饱和: 空闲 tab 全部填 LS-a/b/c/d 的 grounded 问题 (q57+);
  在飞 life3/5/8 落地即收。backlog 常备 ≥3 题。
- 我 (Fable 主循环): LS-b 独立推导 (verify-don't-transcribe), 统筹, 审计, 落账。

## 终止条件 (全局)

- 硬停: 破 3/2 达成 (任一充分靶无条件闭合) / 全 avenue 带终判耗尽 / 外部资源。
- 报告: avenue 边界粗粒度 TG; 未破就说未破。
