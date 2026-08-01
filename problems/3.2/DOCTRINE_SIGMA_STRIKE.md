# DOCTRINE — σ-ladder strike (campaign 3, life 窗口, 2026-08-01 午后)

## 主目标 (一句话)

证明三靶之一的任何 σ>1 进展 = 新定理交卷; σ=2 (E1) = 全胜。
权威量词: proof.tex §17 (FABLE_SECTION_orbit_energy.tex) — ATR / GPRV / GAP-2D-SQRT,
并列充分非等价 (勿写等价)。

## 本战役新杠杆 (session 起点自推, 待审计)

**Second-moment 归约 (绕过逐 t 的 2D-SQRT):** 由加法正交 (Parseval),

    Σ_{t≠0} |B_t(H)|² = p·N_coinc(H) − (#I_H)²,
    N_coinc(H) := #{((r,h),(r',h')) ∈ I_H² : Δ_{r,h} = Δ_{r',h'} in F_p}.

若 **[COINC]**: N_coinc(H) ≤ (#I_H)²/p + C·#I_H, 则 Cauchy–Schwarz 直接给

    Σ_{h≤H} C_h ≪ H + (pH)^{1/2}   (与 2D-SQRT 的推论同强度, 无需逐 t).

代入 σ-梯 (Prop oe-sigma-ladder): 在 Δ=√p 处 Q_p(Δ)≪Δ^{3/2} ⟹ σ=3/2,
**F_p ≪ p^{4/3}** — 破 3/2 纪录即新定理。

**结构账:** Δ_{r,h} = U_{h−1}(r+1)/∏_{j=1}^h (r+j)³ (Casoratian 推导, 待机器验证)。
N_coinc 逐对 (h,h') 是差曲线 X_{h,h'}: Δ_{r,h}=Δ_{r',h'} 的点计数;
逐对 Weil 误差 Σ ~ H³√p, 仅在 H≤p^{1/4} 可控 — 但平凡界 H² 在 H<p^{1/3} 更优,
所以逐对 Weil **不够**; 需族内 Frobenius 误差相消 (√H² 级) → 推到 H=√p。
这独立重导出战役图 §7 "全前沿单点收敛 = family 相容性定理" — 交叉验证通过。

## 战况更新 (12:05, 自己推导+数值双验)

1. **单尺度修正 (keystone)**: σ-梯证明只用一个 Δ=平衡点。σ=3/2 ⟹ 真靶只有一个:
   **[COINC @ H=p^{2/3}]**: N_coinc(H) ≤ (#S_H)²/p + K·#S_H。
   代入梯子 (Q≤ΣC_h≤p^{5/6}): F ≪ 3p²/Δ + 2√2·p^{23/12}/Δ|_{Δ=p^{2/3}} ≪ p^{4/3}。✔
2. **曲线族路线封顶**: 逐对 Weil 死 (H≤p^{1/4} < 平凡界交叉点 p^{1/3});
   即使族内完美平方根相消也只到 H≤√p, 而 H≤√p 的 [COINC] 只复现 3/2 (Δ≤√p
   时梯子给 p^{3/2}, 无增益)。**破 3/2 的活路全在大 gap 区 h>√p (deg N_h>√p,
   Weil 死区)** — 与战役图"family 相容性定理"单点收敛第三次交叉验证。
   ⟹ avenue (b) Dwork/crystal 升为主攻实现器。
3. **循环性消解**: zero-zero 项 (ΣC_h)² 在逐曲线点计数内自动包含 (0维交点
   C_h·C_h'≤9hh' 逐对吸收), 不需 a priori C_h 界; Fourier 也可拆掉——
   centered Cauchy 直接给 ΣC_h ≤ #S/p + √(N_coinc−#S²/p)。
4. **数值 go/no-go 通过** (campaign3_scripts/coinc_gonogo.py):
   D1 恒等式 300/300 ×5 素数; K_emp(p, H=p^{0.66}) = 1.24/1.08/1.08/1.08/1.05
   (p=499..7919), **有界且随 p 递减**; n₀=ΣC_h 远低于 H+√(pH)。
   K_emp ≈ 1.05 与 GPRV 归一方差 1.08 疑似同一常数 (待审计确认联系)。

## Avenues (按初始 promise 排序)

### (a) COINC 曲线族攻坚 [主攻]
1. **(a1)** 机器验证 Δ_{r,h} 恒等式 + N_coinc 的 Parseval 账 (sympy/F_p 直算)。
2. **(a2)** 差曲线 X_{h,h'} 的分量控制/绝对不可约性 (generic (h,h')), 工具:
   更新恒等式 N_{m+g+1}=N_{m+1}N_g(shift)−(s+m+1)⁶N_mN_{g−1}(shift) + 相邻互素
   + 反射定理强制分量的显式剥离 (κ_h 镜像碰撞 → 结构化主项, 不是误差)。
3. **(a3)** 族平均 Frobenius 误差: 离散增亏格族 {X_{h,h'}}_{h,h'≤H} 的
   Σ误差 ≪ pH 机制 — 无现成定理 (Katz vertical 是固定亏格); 攻击点:
   (i) 把 (r,r',h,h') 整体当高维簇切片计数 (h 是整数移位非代数变量, 需换参);
   (ii) 更新恒等式给出的 h-方向自相似 → 转移算子/母函数;
   (iii) Dwork 迹公式一次吃掉全族 (接 avenue b)。
   终止条件: [COINC] 在 H=p^{β} (β>1/3) 成立 ⟹ σ>1 新定理 = 成功;
   或三条子攻击各有书面死因 = 失败判词。

### (b) 秩四 mixed cocycle / Dwork 迹 [奖品线]
素材: Q6517_formulas.txt (L(C)=t 非齐次 ⟹ Y=(C,θC,θ²C,t) rank-4, θY=G(t)Y,
Frobenius 半线性 θΦ=GΦ−pΦG(t^p); 深度层塔 p³c_{p+v}≡b_v+2pD_v+p²E_v (p³))。
攻击: 碰撞计数 = 定 r 的 fibre 结构 → p-adic 迹公式把 Σ_r 一次算完,
正是 (a3)(iii) 的实现器。终止条件: 迹公式给出 N_coinc 主项 + 可控误差,
或明确写出为何 mixed(非齐次)结构阻断标准 Dwork 理论。

### (c) GPRV 分层方差 [cron 并行线, 不重复主攻]
cron 侧 qS-gprv/qU-bilinear/qR-sigma 在飞。我方只做: 落地答案的对表 + 差异化
补角 (薄例外族 −51 分裂律分层)。不自起炉灶。

### (d) Fallback: vdC 变体 + RADIAL-BDH 谱对象
(a)(b) 均有终止判词后才启用。素材在账本 §64–110。

## 死亡证书 (勿重访)

proof.tex §17 + DOSSIER + 账本 §64–110: jet 线性杠杆 / St⊠St / 仿射 Parseval /
Stepanov=Bézout 无盈余 / 结式只控 Q / 九证书 (CAMPAIGN_MAP §2)。
⚠️ 注意: "仿射 Parseval" 死因需先核对是否覆盖本战役 second-moment 账 —
(a1) 第一步就做这个对表, 若已判死则 (a) 降级, 直上 (b)。

## 协同

- ChatGPT life1–10 全空 → 立即饱和 (a1 审计 / a2 不可约 / a3 族平均 / b Dwork /
  qFS2 重发)。cron 十连发在飞, 黑名单纪律照旧。
- Codex 双开 (high+max): 数值 N_coinc 主项验证 + 差曲线分量异常扫描; spec 文件化。
- 落地即收: 每个 answer 先补发再处理。

## 终止条件 (全局)

- 成功: 任一 σ>1 定理端到端书面证明 + 机器数值双验 → 入 proof.tex 新小节。
- 失败: (a)(b)(d) 各三次具体尝试 + 书面死因 → 并入死亡证书, 交 RUN_LOG。
