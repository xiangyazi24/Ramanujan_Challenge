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
