# UNDERSTANDING — Problem 3.2, state at 2026-08-02 (life lane, last-stand run)

Written at session close. Everything here is either (A) gate-verified, (B) explicitly
provisional, or (C) a death certificate. Nothing else is banked. Next session /
Codex: read this file, then `DOCTRINE_LASTSTAND_2026-08-02.md` (route board), then
`ERRATA.md` (the two corrections and the process rule).

## 0. The target — do not drift off it

**MAIN THEOREM: `gcd(d_n a_n, d_n b_n) = e^{o(n)}`**, i.e. `log G_n = o(n)`, where
`d_n = lcm(1..n)^3` and `a_n, b_n` are Apéry's sequences for ζ(3).

The energy line (breaking the 3/2 exponent) is a **means**, not the goal. Earlier
sessions of this campaign were aimed at the energy exponent; that framing was
narrowed today after the main statement was re-read from `notes.md`.

Banked reduction chain: `v_p(G_n) <= 6` for `sqrt(n) < p <= n`; then
`log G_n <= 6 M(n) + O(n^{2/3} log n)` with
`M(n) = sum_{sqrt n < p <= n} log p * |R_p(n)|`,
`R_p(n) = {r < n/2 : r == n mod p, p | b_r}`, `|R_p(n)| <= 1 + n/(2p)`.

## 1. Gate-verified today (`LASTSTAND_VERIFY.py`, 8 gates, exit 0)

Run `python3 LASTSTAND_VERIFY.py`. Every claim below has one independent gate
(independent ground truth, multiple parameters). No gate = not banked.

1. **Digit criterion.** For p prime, `p | b_n` <=> some base-p digit of n lies in
   `Z_p = {d < p : p | b_d}`. (From Gessel/Lucas `b_{mp+r} = b_m b_r mod p`
   iterated, plus primality.) Gate: 630 (p,n) pairs over 7 primes, ground truth =
   exact binomial sum (independent of the recurrence), 0 mismatches.
2. **Top-window collapse.** For `p in (n/2, n]`: `p | b_{n-p} <=> p | b_n`.
   Gate: 855 (n,p) pairs, 0 mismatches.
3. **Averaging identity.** `sum_{n<=X} #{p in (n/2,n] : p|b_n} = #{(p,r) : r in Z_p,
   r>=1, p+r<=X}`. Gate: exact equality at X = 150/300/450/600 (20=20, 37=37,
   61=61, 82=82).
4. **Parity law.** `|Z_p|` is odd <=> `(p-1)/2 in Z_p`. Gate: all p<4000,
   even=546, odd=2, 0 violations. (Mechanism: reflection pairing r <-> p-1-r.)
5. **Reflection.** `b_{p-1-r} = b_r mod p`. Gate: 5 primes checked pointwise.
6. **`sum_{p<=X} log p |Z_p| = Theta(X)`** — measured ratio 0.954/0.907/0.887 at
   X=1000/2000/4000. This **refutes** the tempting sufficient condition
   "`sum_p log p |Z_p| = o(X)`": its hypothesis is false.
7. **Threshold calibration.** `theta(X) ~ X` (0.990/0.997 at X=10^4/10^5), so a
   pointwise bound `|Z_p| <= C log p` yields only `Theta(X)`, not `o(X)`.
8. **Determinant-bilinear scale.** `mean_t |B(t)|` is far below the trivial bound
   (p=1009: 204 vs 29295; p=2003: 403 vs 84172), i.e. the square-root regime.

Separately measured (stable over a 10x range, `scratchpad_laststand/zp_max.py`):
`mean |Z_p| = 1.0140` over 4201 primes below 40000, `sum |Z_p| / pi(X) = 1.014`,
max `|Z_p| = 8`. The zero **pairs** follow Poisson(1/2) — the pure random floor.

## 2. The three live routes (ranked by my judgement at close; codex was asked to re-rank)

**LS-f — `[AVG-ZERO]`: `sum_{p<=X} |Z_p| = o(X)`.** Cleanest target in the campaign.
By items 1+3 above this is EXACTLY the average-over-n version of the main theorem's
top window. Truth `~ pi(X) = X/log X`; unconditional record only `|Z_p| <= 3p^{2/3}`
giving `X^{5/3}/log X`; nothing proved in between. Threshold is delicate: needs
`mean |Z_p| = o(log p)` (item 7), so any `p^alpha` bound is useless. Confirmed
**strictly weaker** than the full conjecture (not equivalent). Equivalent arithmetic
form: `sum_{r<X} #{primes p in (r,X] dividing b_r}` — a cross-prime large-prime-factor
count. The named missing input: **Apéry-specific Cartier rigidity** (bounded
zero-producing state mass in the p-kernel, on average over p).

**LS-g — digit/sieve, pointwise.** Target `#{p in (n/2,n] : (n mod p) in Z_p} =
o(n/log n)` for EVERY n. Structure: one prescribed residue per prime (top window:
`n - p`, affine in p) tested against a sparse set of average size ~1. The gap is a
**quenched diagonal anti-concentration** — ruling out alignment between the affine
clock and the prime-indexed zero sets. This is a **cross-prime correlation** problem,
orthogonal to the whole fixed-p campaign. Note: averaging over n is easy and equals
LS-f; the difficulty is entirely in "every n".

**LS-e' — chart-free determinant/symplectic bilinear (energy gateway).**
`T = sum_{shell} 1[det(u_r, u_{r+d}) = 0]`, error term
`B(t) = sum_{r,d} W(d) e_p(t det(u_r,u_{r+d}))`. This formulation includes the
infinity fibre automatically, repairing the hole that any affine-chart route has
(and that hole is `|Z_p|`, i.e. the main theorem's own object — so chart routes are
conditional on what we are trying to prove). Caveat of record: it counts the
incidence sum `S_B`, which is **stronger** than the real target `U_B` (union);
only `U_B <= S_B`.

## 3. Provisional — NOT banked, do not use as a premise

**§139 completion bookkeeping.** The claim that the completion route lands at
`N log p` (hence not better than trivial) and that "any `N^{eta/2}` power saving
closes `[FR_eta]`". Written out line by line in ledger §146 with the single unproved
hypothesis isolated (square-root cancellation `|F(t,xi)| <= Cp` of the complete
two-variable sum, which has numerical support from cron AS.1 but no proof). I redid
this derivation once and got the same answer, but it has had no independent check.
It is the item most likely to move.

## 4. Death certificates from today — do not re-walk

- **`[NO-RUN]` REFUTED.** The claim `xi_{r+1} = M_r xi_r` for the two-solution row
  is false; the companion matrix propagates the two-time state of ONE solution.
  Live counterexample p=997, r=248, d=182. Correct collision criterion:
  `(T_{r,d})_{21} = 0` (fixed-Borel), orbit-free.
- **`[DEAD-MAX-WEYL-CLASS]`.** Single-variable max-Weyl input cannot close
  `[FR_eta]`: the max route needs `delta >= 1/2 + eta/4` (beyond square root, and
  the truth is exactly square root — measured 3.5–4.8 sqrt(p) for both the harmonic
  and the Apéry phase), and the Parseval-in-t route is an identity (zero information).
- **`[PRIVATE-SINGLETON]`** (codex cocycle lane). A construction with ALL
  pair-correlations vanishing yet union mass `~N/2 = D^{2-o(1)}`: every
  incidence/intersection/cascade/restart combinatorial input is provably insufficient.
- **Chebotarev/fixed-point shortcut to `[AVG-ZERO]`: dead.** The `|Z_p|` distribution
  is the pure random floor (mean 0.945→1.014, var = 2×mean, pairs Poisson(1/2));
  there is no fixed finite group action whose Frobenius fixed-point count it could be.
- **Christol/automatic shortcut to `[AVG-ZERO]`: dead** (two independent tabs).
  Algebraicity does NOT imply zero-coefficient sparsity — `1/(1-t^m)` has a tiny
  automaton and nearly all-zero first p coefficients; no p-uniform automaton size.
- **`sum_p log p |Z_p| = o(X)` as a sufficient condition: hypothesis is false** (gate 6).
- **Gessel conversion is not a free lunch**: it returns the original statement rather
  than producing a local Frobenius invariant.
- **Characteristic-zero nonvanishing has the wrong quantifier order at a fixed prime**:
  `N_3(-3) = 584 = 8*73` vanishes at p=73, and `gcd_{F_73}(N_3,N_4) = X+3`. Do not
  use honest-poles / char-0 resultant coprimality as fixed-p inputs.
- Earlier (still valid): padded-word theorem; the twelve audited banked inputs; the
  nine certificates in `CAMPAIGN_MAP_2026-08-01.md` §2; low-genus collapse;
  char-0 H-algebraization; size arguments (7x slack).

## 5. Process rule instituted today (see ERRATA.md)

Five public retractions in one day; root cause was announcing at exploration
checkpoints rather than verified-stable ones, and in the worst case a **placeholder
branch in my own verifier** let `[NO-RUN]` through unchecked. Rule now: a statement
enters the "proved" column only with a corresponding gate in `LASTSTAND_VERIFY.py`
that actually ran and passed; gates must use independent ground truth and multiple
parameters; death verdicts need a counterexample or complete bookkeeping; a
placeholder in a verifier means that gate did not run.

## 6. In flight at close (not my results — do not bank without checking)

- **codex ultra in tmux window `rc`**, executing `CODEX_SPEC_MAINTHEOREM_laststand.md`
  with the corrections sent afterwards (S_B/U_B, scale regime, char-0 quantifier,
  R1/R2/R3 closed, plus the digit-criterion and `[AVG-ZERO]` results). Its task list
  included `CODEX_MAINTHM_report.md` / `CODEX_MAINTHM_verify.py` and an
  `RC_BREAKTHROUGH_report.md`. **Read its report before trusting anything it wrote,
  and check its verifier for placeholder branches.**
- ChatGPT q77 (the amplification crux for the determinant route) timed out on
  delivery; the tab was likely still running. Do not re-dispatch blindly.
- Answers already landed and archived under `chatgpt-answers/`: Q6798, Q6801, Q6804,
  Q6809, Q6816, Q6817 (plus the earlier wave Q6730–Q6762). All are fast-tier
  triage — re-derive before citing.

## 7. Files

- Ledger: `research/working_notes/FABLE_NOTES_energy_bootstrap.md` §139–146 (today).
- Route board: `DOCTRINE_LASTSTAND_2026-08-02.md`.
- Errata + process rule: `ERRATA.md`.
- Gates: `LASTSTAND_VERIFY.py`; supporting scripts `scratchpad_laststand/`.
- Codex reports today: `CODEX_LASTSTAND_COCYCLE_report.md` (+ verifier),
  `CODEX_NORUN_report.md` (the refutation), `CODEX_MAINTHM_*` (in progress at close).
- Run log: `research/working_notes/Q32_RUN_LOG.md`.
