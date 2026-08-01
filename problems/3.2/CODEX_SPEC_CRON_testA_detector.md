# CODEX SPEC — Test A: exact zero-detector complexity (Q6420 numerical dispatch)

## Context
Apery numbers b_r (1, 5, 73, 1445, ...), recurrence (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}. For prime p let T_p(r) = b_r mod p for 0 <= r <= p-2. Z_p = {r: T_p(r)=0}. The Q6420 dispersion verdict predicts: the minimal polynomial detector of the zero set has degree |S_p| - 1 with |S_p| ~ (1 - e^{-1/2}) p ≈ 0.3935 p — i.e. LINEAR in p, killing any fixed-degree FD hypothesis. This test makes that prediction data.

## Task
Write a single self-contained C program (or C + a thin python driver) `CRON_testA_detector.c` in THIS directory (~/repos/Ramanujan_Challenge/problems/3.2/) that, for EVERY prime p in [500, 4000]:
1. computes the row T_p(r), 0<=r<=p-2, by the recurrence mod p (modular inverse via Fermat);
2. computes the value set S_p = Image(T_p), the multiplicities N_p(v), and collision energy E_p = sum_v N_p(v)^2;
3. verifies the Cauchy–Schwarz bound |S_p| >= (p-1)^2 / E_p exactly;
4. builds the minimal detector Q_{p,min}: the unique polynomial over F_p of degree |S_p \ {0}| supported on the value set, with Q(0)=1 and Q(v)=0 for every v in S_p, v != 0 (Lagrange on the value set — do NOT interpolate on all of F_p); verifies Q_{p,min}(T_p(r)) = 1_{T_p(r)=0} for ALL r (exact, every r, no sampling);
5. records: p, |Z_p|, |S_p|, deg Q_{p,min} = |S_p|-1 (assert 0 in S_p or handle the empty-zero case where deg = |S_p| and the indicator is identically 0 — then Q_min is the zero function; record it as such), ratios |S_p|/p and E_p/p^{5/3}.

## Deliverables
- `CRON_testA_detector.c` (+ optional driver), compiled with plain cc -O2, no external libs.
- `CODEX_TESTA_report.md` in this directory: table of summary statistics (mean/min/max of |S_p|/p vs 0.3935 benchmark with a proper birthday-model prediction derivation: values hit ~ p(1-(1-1/p)^{p-1}) -> (1-e^{-1})p for a uniform random map — note the verdict's 1-e^{-1/2} figure came from a different normalization; CHECK WHICH ONE THE DATA MATCHES and say so plainly), E_p/p^{5/3} profile (does E_p actually grow like p^2/(something)? for a uniform random map E_p ~ 2p; report the true growth exponent by regression), deg Q_min / p profile, and the verdict line: "fixed-degree detectors are off by factor X at p=4000".
- Progress prints every ~50 primes (long-run discipline).

## Hard constraints
- Pure integer arithmetic, no floats in the exact parts.
- Do NOT touch CRON_FRESH_EYES_pointwise.md, ERRATA.md, CAMPAIGN_MAP*, chatgpt-answers/, or any existing CRON_*.py. New files only.
- Do NOT dispatch anything to the ChatGPT bridge (no ask-gpt.py, no scripts/). No network.
- Do NOT git commit; leave files for the session owner to verify and commit.

## Acceptance
- Program runs clean end-to-end for all primes in range; internal asserts (3),(4) never fire.
- Report states the measured |S_p|/p limit constant to 3 decimals with its correct random-map explanation.
- Spot-check hook: for p=17 the program must report Z_17={3,13} and for p=13,29 empty Z (known ground truth).
If stuck > 30 min on any sub-item, deliver what runs + a stall note in the report.
