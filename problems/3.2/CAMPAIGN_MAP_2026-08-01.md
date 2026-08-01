# P3.2 战役总图 — 2026-08-01 (收网固定版 / 下一场接战文件)

一句话: 点态 P3.2 已被**精确等价**为一个有界几何对象的定义特征零计数问题;
八条现成武器全部判死(各有证书); 交付定理群 + 死亡证书 + 施工图。
新 session 从本文件 + FABLE_NOTES §43–62 + CRON_FRESH_EYES 附录 K–AA 接战。

## 1. 定理群 (全部已入 proof.tex, 144 页构建绿, 全部机器认证)

| 定理 | 内容 | 认证 |
|------|------|------|
| strong reflection | y_{p−1−n}=y_n 全解 germ | strong_reflection.py 80/80 |
| structure | A_p=S_p² / qS_p², S_p=固定√F/√(F/q)截断 | sp_spotcheck.py, p<150 |
| tower | b_{p+r} 五层展开+master corner b_p−5≡−7p²H⁽²⁾ (p⁵) | 63/63+31/31+21/21 |
| apparition | lc=(−2\|p), 支反射, 四分律 mod 24 | 165/165+92/92, 双路证明 |
| self-twist rigidity | conjugate-twist={1,δ}, order≤4; Sym²无扭; 3+1仅对角 | q6372_local_system_check.py (40阶+指标) |
| **bounded object** | b_r = 本征降 A_± (rank3, Artin cond 9/11) 的 Mellin 迹归约; Beauville-IV 椭圆(Γ₁(6)) Fricke w₆ 商链 | referee 升级条件已兑现; isogeny 2148点; a_p(E_u)=H_p(u) 1084纤维 |
| **quarter value** | τ_{(p−1)/4} ≡ [Q^p] 权-3 eta-乘子形式 (level-6 四重覆盖) | OOS 5/5 (p≤1153) 双机复跑 |
| parity law | \|Z_p\| ≡ 1[p\|b_{(p−1)/2}] (2); 例外={η形式非常规素数} | 三层: 定理+模形式直算+10⁶扫描 |
| counting lemma | 低阶字符律无害; 前沿=高阶 Mellin 对角 | — |
| Zudilin 数据点 | {11,3137} 纪录 2×10⁴→10⁶ 零新增 (50×) | 双源 78,462 零点逐位一致 |

实证大厦 (双源, p<10⁶): Poisson(1/2) 律 λ=0.4998 χ²=3.43; mean→1; max=12~log p;
α_emp=0; AQI k=2 七尺度平稳; F₄≡0; maxH=3@2×10⁶。

## 2. 九张死亡证书 (每条有精确死因, 勿重访)

1. **Katz 等分布** — 只管特征零归一迹, F_p 归约/整除不可及 (Q6377/Q6405/Q6415)。
2. **T-adic/eigencurve** — 四票: 解析变量=加法/wild/单tame分量内; tame=连通分量, halo 逐分量 (Q6415/Q6411+2)。
3. **Stickelberger 赋值** — slope-0 块=ρ+1 单位项其和≡b_ρ; 零=单位抵消, 赋值盲 (Q6371)。
4. **Galois 轨道范数** — 阈值错配: 整除价 p, 共轭尺寸 p^{3/2}; 轨道离婚(exact全轨道 vs mod-𝔭 分解群) (Q6405/CODEX_LT/q18)。
5. **HGM/Hasse 度数** — 字符概形零维长度 p−1, Hasse 截面低次代表=Apéry 和自身 (Q6393)。
6. **mass formula** — Fourier 对角化后一切改写=|Z_p| 精确别名 (Q6417, nullity=零计数 5/5)。
7. **黑盒有界导手定理** — 常值层反例 (a_p≡3 ⟹ p−2 例外扭): 原理性不可能, 必须用 Apéry 特定结构 (Q6417)。
8. **色散 k=2** — 零检测器无有界导手版本: deg≥|S_p|−1≫p^{1/3} 无条件, SL₂ 万有检测器秩 Ω(p), 指标插值恰 p−1 次 (Q6420)。
9. **符号因子化 (q23)** — 行结式≠值循环阵 (Fourier 对偶不相似, p=17 见证: R_row=2 单位 vs R_correct=0); Chebyshev/CFVZ 因子化真实但在错误行列式上; 值侧奇偶分解=精确 Fourier 别名; 宽-3 递推控制系数不控制 Smith 零度 (Q6429, 我方 49/49 复验 CRON_q6429_verify.py)。

记录: |Z_p| ≤ 3p^{2/3} (无条件, 递归/continuant); N_p(c) ≤ 8p^{3/4}; E(p) ≪ p^{5/3}。

## 3. 存活纲领 = 下一场战役的四个方向

- **[GAP-2] Tannakian 群 + Katz 27.1 水平钩子** (= q18 Package B 正式施工单):
  Katz Thm 27.1 允许有限域序列——天然跨素数定理, 缺输入: A_± 的公共算术/几何
  Tannakian 群 + 相容整实现。经验地基已交付(自扭排除三源: qD 定理 + codex-fm
  单值性 + cron codex ord≤30 数值唯平凡)。常值层反例说明这是必经之路。
- **定义特征相位模** (q22 判决书语言): p^{o(1)} 需要非张量机制的全新
  "phase modulus"——q18 Package E, 唯一直接触及 |Z_p| 的 package, 无文献。
- **坏对角逆定理** (Q6339 口径 + q18 逆定理例外表 10 项, Apéry 二次伴侣=
  第一个 graph 案例): 双层 Mellin 逆定理 M₁(χ)=M₂(χ) 多 χ ⟹ 同层差自扭。
- ~~q23 符号因子化支线~~ — 已判死 (证书 9, Q6429/附录 AB); 复活三条件恰好
  映射到上面三个 Package: (a) W_p 有界次关系→B; (b) cross-r 小矩阵→E;
  (c) 保零击自治化定理 (无来源, 挂起)。

辅助资产: (DRS)/(RLL) 接口引理(α<1/2 关 k=2); Tests A–G 数值任务单(cron
对相关账本续跑); 中点零可剥离小定理; Fricke w₆ 结构(矩常数 4)。

## 4. 文件索引

- 主文: proof.tex (§15 value-distribution, §16 apparition tower, 144pp)。
- 账本: research/working_notes/FABLE_NOTES_energy_bootstrap.md §0–62 (life);
  CRON_FRESH_EYES_pointwise.md 附录 K–AA (cron)。
- 勘误: ERRATA.md (本目录, 常驻)。
- 审计原文: research/working_notes/GPT_Q63xx_*.md; chatgpt-answers/。
- Codex 报告: CODEX_FRANEL_MELLIN.md, CODEX_LT_MELLIN.md, CODEX_HORIZONTAL_GK.md,
  CODEX_JACOBSTHAL_DEEP.md。
- 验证脚本: research/scripts/ (~40 个, 全部当日复跑绿)。
- 施工图: q18 Package A–F (cron 附录 Z); F 已由 quarter-value 定理交付。

## 5. 未结尾声 (新 session 先查收)

- qB2 (Q64xx, Shimura lift 重发, life7) — 判 quarter-value 定理的 Shimura
  对应面, 落地对表 Q6360。
- Q6416 (轨道乘积审计, life3 悬案) — 若落地, 对表 CODEX_LT_MELLIN §1–2。
- cron: q16 (Q6412, cron7) / q19 (Q6414, cron4) 两发 (q23=Q6429 已收, 附录 AB);
  Q6375 已被 tab 自身正式撤回 (Q6438/附录 AC, 仅 sympy 已验部分可引);
  Tests A/F 数值。
- 收网提交文档框架 — 等爸爸定 (TG 已问)。
