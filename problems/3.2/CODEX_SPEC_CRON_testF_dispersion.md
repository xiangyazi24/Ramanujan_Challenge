# CODEX SPEC — Test F: dyadic gap profile, exact indicator vs fixed-degree surrogates (Q6420 numerical dispatch)

## Context
Apery numbers b_r mod p: T_p(r), 0<=r<=p-2, zero set Z_p. The Q6420 k=2 dispersion verdict (file: chatgpt-answers/Q6420.md, sections 2.4, 8, 9) predicts a NUMERICAL SIGNATURE: in the cross-prime pair statistic, FIXED-DEGREE polynomial surrogates of the zero indicator (T_p, T_p^2, T_p^3, centered) DO disperse (square-root cancellation across dyadic gap bins), while the EXACT zero indicator does NOT. Confirming this signature is the empirical fingerprint of the "detector degree must grow" obstruction.

⚠️ The markdown export of Q6420.md LOST the display equations (blank lines where formulas were). Your first task is to RECONSTRUCT the intended definitions from the surrounding prose of sections 2.4 (fourfold Fourier formula), 9 Test E (direct cyclic-centered pair C_{p,q} and its Fourier reconstruction — use the stated equality direct==reconstruction as your self-consistency check), and 9 Test F. Read the file. If a reconstruction remains genuinely ambiguous, FALL BACK to this explicit definition and FLAG the deviation in the report:
  For primes p<q in the window, d=q-p, define for a centered row function f_p(r)=F(T_p(r)) - mean_r F(T_p(r)):
  C_{p,q}(F) = sum_{r=0}^{p-2} f_p(r) f_q(r)   (common index range, i.e. the aligned-index pair sum; also compute the shifted variant sum_r f_p(r) f_q(r+delta) maximized/profiled over |delta|<=d if section 2.4's prose indicates a shift structure).
  Surrogates: F(x)=x, x^2, x^3 (as integers in [0,p) reduced, centered). Exact indicator: F(x)=1_{x=0} (centered by |Z_p|/(p-1)).

## Task
1. Reconstruct/settle the pair statistic as above. Primes p,q in [500,2000], pairs binned dyadically by D=q-p: D in [2^k, 2^{k+1}).
2. For each bin and each F in {indicator, T, T^2, T^3}: record bin count, sum of C_{p,q}, RMS of C_{p,q}, and the normalized ratio RMS / sqrt(expected variance under independence) (compute the independence benchmark exactly from single-prime moments of f_p, not by simulation).
3. Midpoint/reflection removal: recompute the indicator statistics after deleting the midpoint r=(p-1)/2 contribution and symmetrizing out the reflection pairing r <-> p-1-r (project f_p onto the reflection-even component and use only one representative per orbit). Report how much of any excess correlation is explained by these forced components.
4. Near-resonance bookkeeping if the fourfold Fourier reconstruction succeeds: contribution of modes |u| <= X/D, v ≈ -u per bin, else skip with a flag.
5. Verdict table: per bin, does each surrogate show RMS ratio ~ 1 (dispersion) while the exact indicator shows ratio significantly > 1 or a drift? State the signature verdict plainly: CONFIRMED / NOT CONFIRMED / MIXED, with numbers.

## Deliverables
- `CRON_testF_dispersion.c` (heavy loops in C, -O2, no external libs) + optional thin python driver `CRON_testF_driver.py` for the statistics/report stage.
- `CODEX_TESTF_report.md` in this directory: reconstruction notes (what you recovered from the prose, what you fell back on), the bin tables, the verdict.
- Progress prints every ~30s of runtime.

## Hard constraints
- Do NOT touch CRON_FRESH_EYES_pointwise.md, ERRATA.md, CAMPAIGN_MAP*, chatgpt-answers/* (read-only), or existing CRON_* files. New files only.
- No ChatGPT bridge dispatch (no ask-gpt.py), no network.
- No git commit; leave for session owner.

## Acceptance
- Exact ground truth honored: p=17 zeros {3,13}; parity law |Z_p| odd iff p | b_{(p-1)/2} — assert on a sample of 20 primes.
- Independence benchmark derived in closed form in the report (show the two-line variance computation).
- The final verdict table distinguishes surrogate vs indicator behavior per dyadic bin with explicit numbers.
If the full pair loop is too slow for [500,2000], shrink q-range but keep >= 10^4 pairs and say so.
