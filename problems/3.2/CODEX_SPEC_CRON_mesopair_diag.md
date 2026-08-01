# CODEX SPEC: MESO-PAIR five-quantity diagnostic (saturated resultant vs regular roots)

## Mission

Implement the diagnostic computation recommended by the MESO-PAIR attack (Q6567 §11, archived at `chatgpt-answers/Q6567.md` in this directory — READ IT FIRST, especially §3, §4, §8, §11). The point: raw resultant-divisibility statistics are polluted by cut-edge closure and are NOT diagnostic; the decisive objects are the saturated/shifted quantities in the balanced region.

## Definitions (positional gap gauge)

Continuants over Z:
```
K_0(X)=1, K_1(X)=P(X), P(X)=34X^3+51X^2+27X+5
K_{m+1}(X)=P(X+m)K_m(X)-(X+m)^6 K_{m-1}(X)
N_h(r)=K_{h-1}(r+1)
```
Apery numbers b_m: (m+1)^3 b_{m+1} = P(m) b_m - m^3 b_{m-1}, b_0=1, b_1=5.

Orbit: pi(r) = projective Apery state mod p. Existing scripts in this directory (e.g. `CRON_radon_spectrum.py`, `CRON_b1_crosscorr.py`, `CRON_avggcd.py`) contain verified orbit computation (`orbit_keys`-style) and the verified h-direction recurrence. REUSE, do not rewrite from scratch; cross-check any new implementation against them on one prime.

## Task 0 — calibration (must pass before anything else)

1. Machine-verify the exact adjacent-resultant formula for m<=8 over Z (exact integer arithmetic):
```
A_m := Res_X(K_m, K_{m+1}) = (-1)^(m(m+1)/2) * prod_{j=1..m} ((j!)^3 b_j)^6
```
2. Machine-verify the renewal factorization for a few (m,g):
```
Res_X(K_m, K_{m+g+1}) = A_m * Res_X(K_m(X), K_g(X+m+1))
```
3. Reproduce the pollution example: b_3=1445=5*17^2; verify that for p=17, all 15 same-base resultants Res(K_m,K_k), 3<=m<k<=8, are divisible by 17, with common root X=-3 mod 17.

Record PASS/FAIL for each in the report. If any FAILS, stop and report the discrepancy precisely (do not silently continue).

## Task 1 — five quantities per prime

For primes p in {211, 499, 1009, 3001} (extend to 10007 if runtime permits) and D = floor(sqrt(p)*log(p)) (also record subcritical D=floor(sqrt(p)) for comparison), compute over the pair triangle {(a,g): a,g>=1, a+g<=D}:

1. RAW: #{(a,g) : p | Res(N_a(X), N_g(X+a)) over Z reduced mod p} — computed mod p via gcd degree: gcd_{F_p}(N_a(X) mod p, N_g(X+a) mod p) nontrivial. (For speed work directly in F_p[X]; the Z-resultant itself is not needed.)
2. DESATURATED: same count after removing cut-edge roots: exclude common roots X=x0 with x0 in the singular clock set {0,-1,...,-(a+g)} (mod p shifts per the gauge — read Q6567 §4 for the exact cut-edge locus and implement its removal).
3. SHIFTED pairs: #{(a,g): deg gcd_{F_p}(N_a(X), N_g(X+a)) > 0} after desaturation.
4. RHO: sum over (a,g) of rho_p(a,g) = #{r regular, r+a+g <= p-2 : N_a(r)=0 and N_g(r+a)=0 mod p} — the TRUE regular two-edge chain count. Cross-check identity: sum over full triangle of rho = Q_D computed independently from the orbit (Q_D = sum_r C(d_D(r),2)). This identity MUST match exactly per prime — it is the main correctness gate.
5. WEIGHTS: sum of deg gcd_{F_p}(N_a, N_g(X+a)) (desaturated), reported separately for axis region min(a,g)<=G and balanced region min(a,g)>G, with G = floor((p/ (24*D**(2/3)))**0.5) per Q6567 (8.4).

## Task 2 — stratification

Stratify quantities 3-5 by: min(a,g); max(a,g); dyadic ratio a/g; parity (a mod 2, g mod 2); rows with leading-coefficient degree drop; low-Apery-zero rows (does p divide some b_j, j<=D? list them); small-field classes (a=2 or g=2: the -51 layer).

## Task 3 — verdict table

Per prime: axis vs balanced share of each quantity; whether the balanced region shows any excess over the Poisson-scale expectation (E[rho per pair] ~ small; the whole balanced sum should be O(p) if MESO-BALANCED-PAIR holds — report the observed balanced sums vs p); any hidden rows/diagonals after saturation (the decisive test per Q6567 §11).

## Deliverables (all in this directory)

- `CRON_mesopair_diag.py` — the script, with per-step progress prints (every <=10s of runtime; long loops MUST print progress).
- `CODEX_MESOPAIR_DIAG_report.md` — calibration PASS/FAIL, the five-quantity table per prime, strata tables, verdict.
- `mesopair_diag_results.json` — raw numbers.

## Rules

- Exact arithmetic for Task 0 (sympy/flint integers). F_p polynomial arithmetic for Tasks 1-2 (sympy GF or python-flint; flint strongly preferred for speed).
- The rho/Q_D identity check is the correctness gate: if it fails, debug until it passes; do not report results with a failing gate.
- No effort cap. Work until the acceptance criteria are met: Task 0 all PASS, Task 1 five quantities for >=4 primes with gate passing, Tasks 2-3 tables complete.
