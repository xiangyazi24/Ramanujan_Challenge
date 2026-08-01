# CODEX SPEC: supercritical Q_D quantiles + arithmetic stratification (orbit side)

## Mission

Close the owed numerical item (Q6554 blind spot 2): supercritical Q_D statistics beyond band means — max, high quantiles, and arithmetic strata, across many primes, versus a mirror-random baseline. Band means are known to be blind to thin structured lag families ("moving resonance" h=h(p)); this computation is designed to see them.

## Objects

Orbit pi(r) = projective Apery state mod p, r=0..p-2 (nonwrapping). Existing verified orbit code is in this directory: `CRON_radon_spectrum.py` (orbit_keys), `CRON_b1_crosscorr.py` (also has the mirror-random baseline generator `mirror_random_key`), `CRON_stratify_t34.py`. REUSE these; cross-check on one prime before scaling.

```
I(r,h)=1 iff r+h<=p-2 and pi(r)=pi(r+h)
R_h   = sum_r I(r,h)
d_D(r)= sum_{h<=D} I(r,h)
Q_D   = sum_r C(d_D(r),2)
C_p(a,g) = #{r: r+a+g<=p-2, pi(r)=pi(r+a), pi(r+a)=pi(r+a+g)}
```
Identity (correctness gate, check per prime): Q_D = sum_{a,g>=1, a+g<=D} C_p(a,g). Also global gate: 6*Q_{p-2} = T3 = sum_v m_v(m_v-1)(m_v-2). Known calibration values p=10007: Q_D = 0/10/124/1089 at D=100/316/1000/3162.

## Task 1 — supercritical sweep

Primes: all primes in [3000,4200] (144 primes) + {10007, 30011}. Scales: D = floor(sqrt(p)*L) for L in {log p, p^0.1, p^0.2} (cap D <= p-2). For each (p, L):

- Q_D, S_D, M_D=max_r d_D(r).
- Per-lag profile R_h for h<=D: max_h R_h, 95%/99% quantiles, and argmax lags.
- Per-pair C_p(a,g): total (=Q_D via gate), max over pairs, 99% quantile, argmax pairs.

## Task 2 — strata

For each (p,L), split Q_D mass and R_h mass by lag strata:
(i) mirror-forced (even h, and specifically the forced root layer); (ii) axis strips min(a,g) <= G = floor((p/(24*D**(2/3)))**0.5); (iii) h | p-1; (iv) h | p+1; (v) parity classes; (vi) h prime vs smooth (smallest prime factor > sqrt(h) vs not); (vii) the h=2 (-51) layer.

## Task 3 — baseline

Mirror-random baseline (reuse `mirror_random_key`): 3 replicas per prime for the [3000,4200] band, 5 replicas for 10007/30011. Same statistics. Report data/baseline ratios for every statistic; flag any statistic where data exceeds all replicas AND exceeds baseline mean by >3 baseline-sd, per prime and pooled across primes (count of flagged primes per stratum with binomial p-value).

## Task 4 — report

`CODEX_QDSUP_report.md`: calibration gates (PASS/FAIL with numbers), tables per L-scale (pooled over band + the two large primes separately), strata tables, flagged anomalies with their (p, h or (a,g)) coordinates, and a one-paragraph verdict: does any thin lag family carry Q_D mass above the mirror-random level at supercritical scales, or is the picture pure Poisson+mirror everywhere?

## Deliverables

- `CRON_qdsup.py` (progress prints every <=10s of runtime — long loops MUST print).
- `CODEX_QDSUP_report.md`, `qdsup_results.json`.

## Rules

- numpy vectorized where possible; the band has 144 primes x 3 scales — budget accordingly, print per-prime timing.
- Gates must pass (Q_D pair identity per prime; T3 identity on at least 3 primes; the p=10007 calibration row) before Tasks 2-4 are reported.
- No effort cap; run to completion of all four tasks.
