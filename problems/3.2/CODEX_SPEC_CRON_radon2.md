# CODEX SPEC — Radon 谱精化实验协议（P3.2，Q6556 §7-8）

## 任务
读本目录 `chatgpt-answers/Q6556.md` 的 §7（central-binomial 阴性对照）与 §8（Apéry 实验协议），照其设计实现并运行。设计冲突处以 Q6556 为准。已有基础设施：`CRON_radon_spectrum.py`（Apéry b,c 序列、hist+双 FFT 骨架，可复用）。

## 交付物（本目录）
1. `CRON_radon2.py` — 实现：(a) §7.3 central-binomial 对照（F_bin 谱，预言 η=1 线性增长 + Dirichlet ridge）；(b) §8 Apéry 协议全套：cyclic table 口径、raw transform 的 max/median/99分位/超阈值点数、zero-push (0.1) 诊断量 B_p 与 B_p^max、§8.3 exceptional-locus 测试（最大值坐标是否落固定直线族）。素数扫到计算可承受（至少 101..3001 数个）。
2. `CRON_RADON2_REPORT.md` — 结果 + 按 Q6556 的判读（Apéry vs binomial 对照是关键图；zero-push B_p=O(1) 是 GPRV 主诊断）。

## 纪律
小素数先验证正确性再放大；REPORT 只写脚本真实输出；完成后 git add+commit（不 push）。
