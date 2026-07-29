# Q3.2 FRONTIER ATTACK — the second-moment / large-sieve reframing (2026-07-22, automode)

**Standing order (Xiang, 2026-07-22):** pursue fully-unconditional pointwise P3.2, do NOT concede.
The disproof has been argued untenable (life window) ⟹ the conjecture is TRUE; the job is to PROVE it.
Keep a precise record of what we tried, what we found, and *why* each wall stands (or dissolves).
Prior campaign declared a "seven-fold barrier" + terminus; the handoff explicitly warns those walls may be
**mischaracterized** (5+ premature terminus declarations last campaign, each dissolved under push).

---

## 0. The reframing (why this is not a re-run of the closed families)

The closed families (R9–R14) were argued at PROSE level and — critically — the (SA) formulation restricted
to a SHORT interval `n∈[2X, X·log^A X]` (only small slopes `q=⌊n/p⌋`). This document puts the pointwise
target on an **explicit second-moment / large-sieve footing over the HONEST interval** `n∈[X,4X²)` and asks,
quantitatively, exactly where the power-saving is lost. This is a sharper object than the prose barrier.

**Verified ground truth (Python, 2026-07-22, `scratchpad/verify_lucas.py`):**
- **(L) Lucas 2-digit:** for `n=qp+r`, `q=⌊n/p⌋<p`, `b_n ≡ b_q·b_r (mod p)` — **1300/1300 checks, 0 fail.**
  ⟹ for middle prime `p∈(X,2X]` and `n<p²`, `p|b_n ⟺ b_q≡0` OR `b_r≡0`; the `b_q≡0` channel is O(polylog)
  primes total (q small, `b_q` fixed), so up to polylog, `p|G_n ⟺ (n mod p)∈Z_p`. (Recovers support law.)
- **(Z) zero-set structure:** `mean|Z_p|≈1.019` (doubled-Poisson E≈1); `Z_p` is EMPTY or a union of
  **reflection pairs `{z, p−1−z}`** (verified 109/109) plus possibly the center `(p−1)/2`; `max|Z_p|=6`, p<600.

## 1. The exact second-moment target

For dyadic block `p∈(X,2X]`, `K(n)=#{p∈(X,2X]: (n mod p)∈Z_p}`. Pointwise P3.2 follows from
> **(2M)** `N_coll(X) = Σ_{n∈[X,4X²)} C(K(n),2) ≪ X^{2−δ}` for some `δ>0`.
because a single term bounds the max: `K(n)² ≤ 2N_coll+K(n) ≪ X^{2−δ}` for EVERY `n` ⟹ `K(n)≪X^{1−δ/2}=o(π(X))`.

**Additive-character expansion.** With `Ẑ_p(a)=Σ_{z∈Z_p} e(−az/p)`:
```
K(n) − M(p-block mean) = Σ_{p∈(X,2X]} (1/p) Σ_{a≠0 (p)} e(an/p) Ẑ_p(a),   M = Σ_p |Z_p|/p ≈ log 2.
Σ_{n∈I} |K(n)−M|²  =  DIAG (p1=p2,a1=a2)  +  OFFDIAG.
DIAG = Σ_p (1/p²) Σ_{a≠0} |Ẑ_p(a)|² · |I|  ≈ |I|·Σ_p |Z_p|/p ≈ |I| log 2   [= Var ≈ mean, POISSON-consistent]
OFFDIAG = Σ_{p1≠p2} (1/p1p2) Σ_{a1,a2≠0} Ẑ_{p1}(a1) conj(Ẑ_{p2}(a2)) · \hat{1_I}(a1/p1 − a2/p2).
```

**The quantitative wall (computed).** Using only the TRIVIAL bound `|Ẑ_p(a)| ≤ |Z_p| = O(1)`:
`OFFDIAG ≪ Σ_{p1≠p2}(1/X²)·Σ_{a1,a2} min(|I|, 1/‖a1/p1−a2/p2‖)`. The inner double sum ≈ `X²·log X` per
prime pair (m=a1p2−a2p1 ranges over |m|≤X², each hit O(1) times, Σ_m min(|I|,p1p2/|m|)≈X²logX with |I|=X²);
times `(X/log X)²` prime pairs `×(1/X²)` ⟹ **OFFDIAG ≪ X²/log X.** A single log, **NO power saving.**

**Conclusion — the wall, stated exactly:** the trivial Fourier bound gives `N_coll ≍ X²/log X`; the power
saving `X^{2−δ}` is EQUIVALENT to genuine **cross-prime phase cancellation** in the bilinear form
`Σ_{p1,p2} Ẑ_{p1}(a1)\overline{Ẑ_{p2}(a2)}·(phase)`. This IS "Chebotarev without a quotient", now pinned to a
concrete bilinear exponential-sum inequality rather than prose.

## 2. New attack vectors on the cross-prime cancellation (this campaign)

- **(F1) Reflection leverage.** `Z_p` = union of pairs `{z,p−1−z}` ⟹ `Ẑ_p(a)` is REAL up to the fixed phase
  `e(−a(p−1)/2·2/... )`: `Ẑ_p(a)=e(−a(p−1)/p·?)·Σ_{pairs}2cos(2πa(2z−(p−1))/(2p))`. Does the forced realness +
  the specific cosine phases create cancellation the trivial bound misses? [→ Fable (A)]
- **(F2) Hypergeometric/Weil bound on `Ẑ_p(a)`.** `b_r mod p` = trace of a hypergeometric motive / `₄F₃` over
  F_p (Ahlgren–Ono/Kilbourn, weight-4 level-8). "b_r≡0" = vanishing locus of an explicit character sum; try to
  estimate `Ẑ_p(a)=Σ_{r:b_r≡0}e(−ar/p)` by completing against the full trace-function sum. Trace-function /
  Fouvry–Kowalski–Michel machinery for `Σ_{zeros} e(−ar/p)`? Obstruction = zero-detection conductor. [→ both]
- **(F3) Bilinear / Type-II via Lucas multiplicativity.** `b_d` is base-p-digit-multiplicative (Lucas). On the
  q=1 arc, N_coll's core is `#{d<n/2: (n−d) prime, (n−d)|b_d}`. Split `b_d` via a Vaughan/Heath-Brown identity
  into bilinear pieces exploiting digit-multiplicativity → Type-II cancellation. New because Lucas gives a
  genuine multiplicative structure a bilinear sieve can bite. [→ ChatGPT (1)]
- **(F4) Interval-regime correction.** Prior (SA) sat in `|I|=X^{1+o(1)}≪Q²=X²` (large sieve loses X^{1−o(1)}).
  The honest interval `[X,4X²)` is BALANCED `|I|≈Q²`. Is the difficulty genuinely still only the q=1 arc, or
  does the balanced-regime large sieve over ALL slopes give something the short-interval prose missed? [→ Fable (C)]
- **(F5) NEW avenue (not in DOCTRINE): Stepanov-style direct auxiliary polynomial.** Stepanov counts points on
  curves mod p WITHOUT a Galois group (elementary Weil). Vertical `|Z_p|` is already O(1); the horizontal
  decorrelation is the target. Check whether the degree-≥p/C ("entropy-positive") claim for Z_p's defining
  locus is PROVED or asserted (handoff: walls may be mischaracterized) — if Z_p were bounded-genus, Stepanov +
  large sieve would close it. VERIFY THE DEGREE CLAIM NUMERICALLY FIRST.

## 3. Terminal conditions
- Success: a proof of `N_coll ≪ X^{2−δ}` (any of F1–F5), OR a different sufficient condition.
- Proof-of-failure (per vector): a written, mechanism-level reason the vector cannot beat `X²/log X`, checked
  against the ACTUAL structure (reflection + hypergeometric + honest interval), not inherited prose.
- Do NOT accept the inherited seven-fold terminus without re-deriving the wall on THIS concrete bilinear object.

## 4. SELF-CORRECTION (2026-07-22, measured) — the interval is FORCED short; the battlefield is exact

Two numerical experiments (`scratchpad/measure_ncoll.py`, `measure_short2.py`) corrected the §1 framing:

**(i) The full interval `[X,4X²)` is USELESS.** Measured full-interval `N_coll ≈ X²/(2log²X) = X^{2−o(1)}`
(X=2600: predicted 54600, measured 67993). No power saving ⟹ `K(n)²≤2N_coll` gives only `K(n)≤X/logX` =
TRIVIAL. The full second moment is dominated by irrelevant LARGE-SLOPE collisions spread over a huge n-range.
**My morning "honest interval escapes the bad large-sieve regime" idea was WRONG.** Recorded as a dead end.

**(ii) The forced target is the SHORT interval, and it IS provable-in-principle (huge margin).** For a FIXED n,
decompose its middle primes `p∈(√n,n]` by dyadic block. BULK blocks `X∈[√n, n/L]` (slopes `q∈[L,√n]`) are
handled by the TRIVIAL bound: `Σ_{X≤n/L} π(2X)log(2X) ≈ Σ_X 2X ≈ 4n/L = o(n)`. Only TOP blocks `X∈[n/L,n]`
(slopes `q≤L`) need a nontrivial bound — exactly the small-slope hard core (Fable §5.27 / Q4928 q=1). There
`n∈[2X, XL]` IS in range, so the pair-rigidity `sup_n K_X(n) ≤ 1+√(2·N_coll(short))` applies. Measured
short-interval `N_coll(J=[2X,X(logX)²]) ≈ X^{1+o(1)}` (X=8000: 1932, ≈X/2); gives `K_X(n)≤X^{1/2+o(1)}=o(π(X))`,
top-block contribution `Σ_{top}X^{1/2}logX = o(n)`. **So the whole pointwise problem = prove `N_coll(short) ≪
X^{2−δ}` (true value ≈ X^{1+o(1)}, δ-margin ≈ 1).**

**(iii) The wall, now EXACT.** Short interval `|J|=X^{1+o(1)}`, moduli `p~X`, so `|J| < p1p2 ≈ X²`: the CRT
pair-residues are UNDER-sampled (interval shorter than the pair-modulus). Trivial Fourier bound `|Ẑ_p(a)|=O(1)`
gives `N_coll ≍ X²/log` — a FULL POWER above the truth `X^{1+o(1)}`. Capturing ANY power saving `δ>0` needs
cross-prime cancellation. Since `|J|>p` (interval longer than a SINGLE modulus, `|J|/p≈(logX)²`), the
per-prime sum COMPLETES to a full character sum mod p — so the natural tool is a **Weil/Deligne bound on the
complete exponential sum attached to `Z_p`** (trace-function machinery), NOT the large sieve. This is F2, and it
is the only surviving hope besides F3 (Type-II via Lucas). F4/F5 are dead (recorded above).

**(iv) F5 dead:** measured `max|Z_p|` grows (6@p<600 → 8@p=3727), sizes ~Poisson-in-pairs — `Z_p` is NOT
bounded-genus F_p-points; Stepanov's bounded-degree auxiliary polynomial does not apply.

## RUN NOTES (append as we go)
- 2026-07-22 open: reframing written; ground truth (L)+(Z) verified; Fable+ChatGPT dual oracle dispatched.
- 2026-07-22 self-correction (§4): full interval useless; forced short-interval target N_coll≈X^{1+o(1)},
  wall = complete-sum Weil cancellation for Z_p (F2) in the |J|<p² under-sampled regime. F4,F5 dead. Live: F2,F3.
- Next self-driven: numerically probe Ẑ_p(a) statistics (does the complete sum attached to Z_p exhibit Weil
  √p-cancellation / bounded conductor?) — the decisive test for whether F2 can beat trivial. Then synthesize
  with Fable/ChatGPT oracle returns.

## 5. POTENTIAL CRACK (2026-07-22, measured) — the complete sums have BOUNDED-CONDUCTOR √p cancellation

`scratchpad/probe_tracefn.py`: the two-variable complete sum on the graph {(r,b_r): r∈F_p},
`S(t,a)=Σ_{r∈F_p} e((t·b_r − a·r)/p)`, satisfies (over all t≠0, a):
```
 p     max|S|/√p    mean     p99    p99.9
 101     3.02       0.888    2.13    2.63
 211     2.99       0.888    2.09    2.58
 401     3.32       0.886    2.15    2.65
 809     3.87       0.886    2.13    2.61
 1601    3.63       0.886    2.15    2.62
```
**`max|S|/√p` is BOUNDED (~3–4), does NOT grow with p; mean = 0.886 = E|complex Gaussian| — the values
`S/√p` are distributed like standard complex Gaussians (Sato–Tate/CUE).** So `S(t,a)` has genuine Weil
√p-cancellation with bounded sup-norm — i.e. `r↦b_r mod p` behaves like a **bounded-conductor trace function**
at the level of the additive-character complete sum. This is what FKM (Fouvry–Kowalski–Michel) trace-function
machinery needs as input, and it **appears to contradict** the R14.4 "entropy-positive / deg Φ_p=4(p−1) / no
bounded complexity" reading — at least, the additive complete sum sees NO such obstruction; it cancels cleanly.

**Caveat (do not overclaim).** Bounded sup-norm of the complete sum is NECESSARY for the trace-function route
but not sufficient: (i) a merely-pseudorandom `b_r` also gives √p cancellation (random-sum max/√p ≈ √(4log p)≈5
at p=1601 — our 3.6 is slightly MORE rigid than random, mild evidence beyond pseudorandomness); (ii) bounding
`S(t,a)` per prime does NOT by itself bound the CROSS-PRIME bilinear `Σ_{p1,p2}Ẑ_{p1}(a1)\overline{Ẑ_{p2}(a2)}·
(phase)` — that needs an FKM "sum/bilinear of trace functions across moduli" theorem, and whether `b_r`'s sheaf
is geometrically generic enough (large monodromy, no special multiplicative relations across p) is the open
point. R14's claim was that no cross-p structure exists to EXPLOIT for largeness; here we want it for
DECORRELATION (independence ⟹ collisions rare), the opposite use.

**Decisive next tests:** (a) push p→~4000 to confirm `max|S|/√p` stays bounded (⟹ bounded geometric
conductor, not slow growth); (b) exhibit or refute an actual SHEAF: is `r↦b_r mod p` the trace function of a
constructible ℓ-adic sheaf of bounded conductor on 𝔸¹? (the Apéry local system is the symmetric square of a
weight-2 / a rank-4 K3 motive — the sheaf should be its (r-parametrized) specialization). If yes, FKM gives the
cross-prime bilinear bound and F2 CLOSES the pointwise problem. This is the live crack to drive.

## 6. F2 KILLED by Fable oracle (2026-07-22) — the √p noise floor vs the O(1) target (decisive)

Fable's correction (verified — it is a scale argument, airtight): the √p cancellation of `S(t,a)` is REAL but
USELESS. `Ẑ_p(a)=(1/p)Σ_{t≠0}S(t,a)`; summing `p` terms each `≤C√p` gives only `Ẑ_p=O(√p)` — WORSE than the
trivial `|Ẑ_p|≤|Z_p|=O(1)`. **Completion has an irreducible √p noise floor; a set of size O(1) is invisible
below it. No Weil/trace-function/Deligne input on the complete sum can EVER beat the trivial bound on `Ẑ_p`.**
The cross-prime cancellation must come from POSITIONAL correlations of `z_p` across `p`, never from per-p
Fourier size. My §5 "crack" was a correct fact with a wrong use — F2 (per-prime algebraic cancellation) is
DEAD. (Also: `r↦b_r mod p` is p-adically analytic / Mahler-degree ~p, not a bounded-degree rational function
of the geometric variable `r`; `r` sits in the character-DUAL slot, so Deligne is for the wrong variable.)

## 7. THE WALL, re-pointed exactly by Fable (2026-07-22) — horizontal equidistribution of z_p

- **The reduction confirmed**, and the target is `N_coll(short) ≪ X^{2−δ}` (short interval; the FULL-interval
  pair count is `≍X²/log²X` provably by elementary CRT and CANNOT go lower — so only the short interval is a
  live target; second moments cap at trivial on the full interval). Reflection buys only a factor 2 (it
  re-centers 1_{Z_p} about the half-integer (p−1)/2, matching the kernel frequency — cannot cancel).
- **The adversarial-repositioning obstruction (clean, concrete).** Model `Z_p={n₀ mod p, p−1−(n₀ mod p)}` for a
  fixed `n₀`: `|Z_p|=2`, full reflection symmetry, identical trivial Fourier bounds — yet `K(n₀)=π(X)`-many and
  `N_coll` saturates. So ANY argument using only cardinality + symmetry + positional-agnostic per-p bounds
  proves a false statement for this model and cannot exist. New math must use the POSITION of `z_p`.
- **The exact missing inequality (Fable):** `Σ_{p∈(X,2X]} α_p · e(c·z_p·\overline{M_p}/p) ≪ π(X)^{1−δ}`
  uniformly in `c,M` (`z_p∈Z_p`, `\overline{M_p}=M^{−1} mod p`, `|α_p|≤1`) — horizontal (p-varying)
  equidistribution WITH POWER SAVING for the Apéry zero position `z_p`. With `z_p` an explicit algebraic
  function of `p` this is DFI/Kloosterman/dispersion (hard-but-done). With `z_p` defined by `b_{z_p}≡0 (p)` it
  needs a reciprocity organizing a p-varying congruence whose modulus exceeds the **weight-4 Hasse rigidity
  threshold**: weight-2 has `|a_p|<2√p<p` so `a_p≡0(p)⟹a_p=0` (supersingular geometry — rigid); **weight-4 has
  `|a_p|≤2p^{3/2}≫p`, so `a_p≡0(p)` forces nothing — rigidity evaporates.** Same wall as Lehmer `τ(p)≢0` /
  weight-≥4 non-ordinariness; untouched by GRH (fixed-modulus, not horizontal).

## 8. LIVE CRACK to drive (Fable's blind spot) — the OFF-CENTER zeros

Fable (and the whole prior campaign) reduces the wall to the CENTRAL digit `z=(p−1)/2 ⟺ a_p` (weight-4).
But `Z_p` is empty or a reflection pair at a GENERIC position; the central zero (`a_p≡0 mod p`) is RARE. The
bulk of `Z_p` is OFF-CENTER: `b_z≡0 mod p` with `z≠(p−1)/2`. **What governs the off-center zero position?**
If the off-center `z_p` reduces to `a_p` it inherits the weight-4 wall; if it has DIFFERENT (more tractable, or
provably-equidistributed) structure — e.g. a reciprocity `b_z mod p ↔ (·) mod z`, or a modular parametrization
of the off-center vanishing — that is the quotient-free-locating handle. NEXT: measure the distribution and
arithmetic of the off-center `z_p` (is `z_p/p` equidistributed? any p-varying functional relation? reciprocity?).
Also await ChatGPT q2 (automatic-sequence-over-primes route: Müllner/Drmota — orthogonality of automatic
sequences to multiplicative functions along primes — which Fable did NOT address; F3 still live).

## 9. CONSOLIDATED STATE after round 1 (2026-07-22) — wall held under push, but sharpened; routes status

Every route funnels to ONE wall: **horizontal (p-varying) equidistribution with power saving of the Apéry
zero position z_p** (`Σ_p α_p e(c z_p \overline{M}/p) ≪ π(X)^{1−δ}`). Corroboration this session:
- Large-deviation / higher-moment route: `#{n∈J_short: K(n)≥k}` for k≥2 needs the CRT-lift of the k-tuple of
  zeros to equidistribute in the short interval (modulus X^k > |J|) — SAME wall. So max_n K(n) small ⟺ same.
- Anatomy route (q=1 arc): `K(n)=#{d<n/2: (n−d) prime, (n−d)|b_d}` — the "diagonal coincidence" prime = n−d
  is an additive/collision condition; density bound `Σ_p|Z_p|≪π(X)` is an AVERAGE that does NOT stop one n₀
  from grabbing many (the adversarial alignment is not excluded by any provable Apéry constraint found). SAME.

**Routes status:** F1 (reflection) = factor-2 only, DEAD. F2 (per-prime algebraic/Weil cancellation) = DEAD
(√p floor vs O(1)). F4 (interval regime) = DEAD (full interval provably X²/log², short interval forced). F5
(Stepanov bounded-genus) = DEAD (Z_p entropy-positive). **F3 (digit/transfer-matrix bilinear) = LIVE.**

**F3 SHARPENED — the transfer-matrix bilinear split (fresh sub-avenue, in neither framework).** `b_z mod p` =
a fixed entry of the transfer-matrix product `M(1)M(2)···M(z) mod p` (M(j) polynomial-entry, from the Apéry
recurrence). Split at a midpoint m: `M(1)···M(z) = A·B`, `A=M(1)···M(m)`, `B=M(m+1)···M(z)` — `b_z` is
BILINEAR in (entries of A, entries of B). On the q=1 arc `z=n−p`, summing the vanishing condition over p is a
genuine Type-II bilinear form in the transfer-matrix entries — a structure the automatic-along-primes machinery
(single digit z<p has no self-similarity) and the trace-function machinery (r in the dual slot) both MISS. This
is the concrete object to test for cancellation (dispersion/DFI-style), and the one place the wall might be
mischaracterized. Also awaiting ChatGPT q2 (automatic-over-primes lit) for an independent read.

**NOT a terminus declaration.** The wall held under this round's push and is now pinned to an exact bilinear /
horizontal-equidistribution inequality tied to the weight-4 Hasse threshold — sharper than the prior prose. Two
of my own framings were corrected mid-round (full-interval; F2 usefulness). F3-transfer-matrix + the two
pending ChatGPT reads are the live continuation; the conjecture is TRUE (z_p equidistribution measured).

## 10. F3a transfer-matrix bilinear — VERIFIED structure, but the concrete re-collapse (2026-07-22)

VERIFIED (`scratchpad/transfer_matrix.py`, 202/202): `b_z mod p` = a fixed entry of `Π_{j=1}^{z-1}M(j)·[5,1]`,
`M(j)=[[P(j)inv((j+1)³), −j³inv((j+1)³)],[1,0]] mod p`. Midpoint split `Π = B·A`, `A=Π_{1..m}`, `B=Π_{m+1..z-1}`:
```
b_z ≡ B₀₀·u + B₀₁·w (mod p),   u=5A₀₀+A₀₁, w=5A₁₀+A₁₁.
```
So `b_z≡0 mod p ⟺ (B₀₀,B₀₁)·(u,w)≡0` — a genuine BILINEAR condition (upper half × lower half).

**Why it re-collapses (concrete, not hand-waving).** On the q=1 arc `z=n−p`: (i) fix small `m` ⟹ `A(1..m) mod p`
is the reduction of a FIXED rational matrix, so `(u:w)` is an essentially FIXED direction mod p — the split
isolates a fixed line and ALL the difficulty stays in the long product `B(p)` (= the full horizontal problem,
no gain); (ii) balanced `m~z/2` ⟹ BOTH halves `A,B` are prefix/suffix of the SAME product `Π_{1..n-p}`, hence
BOTH locked to the single parameter `z=n−p` — they do NOT range independently, so there is no Type-II bilinear
variable. The ONLY genuinely independent bilinear pair in the whole problem is `(p₁,p₂)` in
`N_coll=Σ_n[p₁|b_{n mod p₁}][p₂|b_{n mod p₂}]` — and decoupling THAT pair IS the off-diagonal cancellation =
the wall. So the transfer-matrix bilinear is real but provides no new independent variable; F3a funnels to the
same horizontal-equidistribution wall. (Awaiting Fable #2 independent verdict.)

## 11. Fable #2 verdict (2026-07-22) — F3 refuted, barrier STRENGTHENED, + a new positive deliverable

**F3 is not a crack; both branches funnel to the horizontal wall with named collapse points:**
- **F3a dies at Fourier completion.** The matrix split reorganizes the COMPUTATION of z_p but cannot decouple
  z_p from p; after completing the short-interval sum in n, what survives is exactly Σ_p e(−c z_p q̄_p/p) (the
  horizontal Weyl sum). z_p has Christol degree ~p^C (Adamczewski–Bell 2013, diagonals of rational functions),
  so Kloosterman completion lands on unbounded-complexity varieties, outside Weil/Deligne. Floor calibrations:
  even z_p = root of a FIXED cubic congruence to prime moduli is open (DFI 1995 does only quadratics via Salié);
  the 1×1 case n! mod p (perfect bilinear split) already defeats all machinery (Erdős/Garaev–Luca–Shparlinski).
- **F3b dies by single-digit triviality + repositioning-invariance.** Automatic content = digit-multiplicativity
  = "output = product of single-digit values", INVARIANT under repositioning the single-digit table ⟹ the
  adversarial model is realizable within Lucas-type families ⟹ any theorem using (automatic + |Z_p| + reflection)
  is false. This STRENGTHENS the barrier to cover all automatic-sequence methods (Müllner/DMR/Konieczny).

**NEW structure (passes the adversarial filter) — the Katz–Mellin identification.** Via Gross–Koblitz,
`b_z mod p` = a finite-field ₄F₃ with character parameter ω^z = the MELLIN transform of a FIXED bounded-
complexity trace function (Straub's rational function generating the diagonal) at χ_z. Explains both regimes:
z=(p−1)/2 → quadratic char → completes to the 8.4.a.a motive (Ahlgren–Ono); off-center → char parameter
z/(p−1) moves with p → incomplete. Consequences:
- **PROVABLE VERTICAL PROGRAM (new bankable deliverable):** Katz's Mellin–Sato–Tate convolution machinery
  (Ann. Math. Studies 180) applies to rank-bounded hypergeometric families with varying character twist ⟹
  should UNCONDITIONALLY prove the vertical statistics of z↦b_z mod p: mean|Z_p|=1, Poisson zero counts,
  equidistribution of z_p/p — turning the measured χ²=4.0 / Weyl≈0.03 into THEOREMS. Strongest rigorous support
  for P3.2 short of P3.2. [FLAGGED: verify the Gross–Koblitz elementary-factor bookkeeping numerically for a few
  (p,z); check Katz applicability — rank, wild ramification at the character-degenerate locus.]
- **EXACT horizontal funnel:** P3.2 = the INVERSE Mellin problem (which character vanishes, tracked across p) =
  horizontal Sato–Tate (Katz conjecture: vertical known, horizontal open in every instance; only Fouvry–Michel
  averages). The missing object, named precisely: a theta-Waldspurger-type correspondence for TRUNCATED weight-4
  hypergeometric vanishing (analogue of Lang–Trotter a_p=0 ⟺ supersingular/CM ≪X^{3/4}, Fouvry–Murty, via
  weight-3/2 theta descent). Off-center INCOMPLETENESS blocks attaching a bounded-conductor motive ⟹ genuinely
  new mathematics, not an application of existing theorems.

**Fourth-idea inventory (master filter: must (α) distinguish b from repositioned b AND (β) bind z_p across p;
every structure has α XOR β):** large sieve (size-only, dead); fixed-integer binding (= 2^n mod p class,
hopeless); p-adic/crystalline via Morita Γ_p + Weierstrass (per-p rigidity, proves mult-1 & |Z_p| bounds, but
per-p not cross-p — funnels); reciprocity z↔p (the only solution-shaped idea = the nonexistent theta descent).

## 12. Central anchor VERIFIED (2026-07-22) — Ahlgren–Ono b_{(p−1)/2} ≡ a_p mod p, 548/548

`scratchpad/verify_ahlgren_ono.py`: computed `a_p` as the q-expansion coefficient of `η(2τ)^4 η(4τ)^4`
(the weight-4 level-8 newform 8.4.a.a) and `b_{(p−1)/2} mod p` via the recurrence, for all primes 5≤p<4000:
**548/548 exact match, 0 fails.** Confirms the CENTRAL-digit end of Fable #2's Katz–Mellin identification (the
central character is quadratic ⟹ the sum completes to the 8.4.a.a motive). The off-center ₄F₃ (general ω^z)
end still needs the Gross–Koblitz elementary-factor bookkeeping check before the full Katz–Mellin vertical
program is quotable as rigorous — but the anchor is solid.

## 13. ChatGPT Q_A + Q_B (2026-07-22) — TRIPLE independent convergence on the wall

**Q_A (horizontal-equidistribution feasibility, Pro).** Pointwise `K_X(n)≪X^{1−δ}` is OPEN; bottleneck is
STRUCTURAL: need a fixed-complexity horizontal description of the Apéry zero detector. Key refinements:
- **Correction to my framing:** the unweighted Weyl sum `Σ_p e(cz_p/p)` is NOT equivalent to `K_X(n)`; the exact
  detector carries the reciprocal twist `e(−hn/p)` and needs a cross-prime DISPERSION / small-CRT theorem.
  Arbitrary-weight (|α_p|≤1) horizontal cancellation is FALSE (adversary picks the conjugate).
- **DFI (Duke–Friedlander–Iwaniec, roots of quadratic congruences to prime moduli, Ann.Math.141 1995) = closest
  MODEL**, but every load-bearing ingredient absent: Apéry zero has no fixed bounded-degree algebraic
  parametrization; interpolation degree ~p; r is a Mellin-character index of the K3 Hasse function, not a point
  on a fixed curve. [Converges with Fable's Christol-p^C + Mellin.]
- **FKM (Fouvry–Kowalski–Michel, "Algebraic trace functions over the primes", Duke 163 2014, Thm 1.5) = wrong
  quantifiers**: fixes ONE field F_q, varies args inside; Apéry changes field to F_p per outer prime. Needs a
  compatible system over ℤ with Frob trace = e(cz_p/p): none known. [Converges with Fable.]
- **Large sieve** (Montgomery–Vaughan): power-saving for ALMOST ALL m<X², NOT one prescribed n; the adversarial
  model (z_p=n₀ mod p) not excludable. [Independently reconstructs the adversarial obstruction.]
- **Exact next-theorem menu (none known):** (1) fixed algebraic correspondence F(R,P)∈ℤ[R,P] encoding R∈Z_p;
  (2) bounded-conductor zero detector (phase = Frobenius trace); (3) direct horizontal dispersion; (4) **higher
  CRT non-clustering (3rd/4th factorial moment)** — "(3)/(4) would solve the needed part directly." [(4) echoes
  Fable #1's m≥3.]

**Q_B (automatic-sequences-along-primes, Pro).** DEAD OUTRIGHT. The prime in Müllner (Duke 2017) /
Lemańczyk–Müllner / Mauduit–Rivat / Byszewski–Konieczny–Müllner (Gowers) is the INPUT INDEX of ONE fixed-base
sequence; here p is the varying CHARACTERISTIC/base and n mod p is a one-symbol lookup. **Decisive no-go
(independently reconstructs Fable's repositioning-invariance):** for ANY S_p, the digit-product weight
w_p(d)=0 iff d∈S_p has exact Lucas structure + two states + singleton zero; setting S_p={n mod p} gives
`#{p: W_p(n mod p)=0}=π(2X)−π(X)`. So no per-prime-automaticity + bounded-state + Lucas + singleton theorem
gives a horizontal bound. Apéry's cross-prime coherence is from the integral recurrence/modular family, which
automaticity cannot access. Same bottleneck: bounded-conductor Chebotarev parametrization OR small-CRT dispersion.

**TRIPLE CONVERGENCE (Fable#1, Fable#2, ChatGPT Q_A, Q_B):** the wall is real and correctly characterized; two
independent reconstructions of the adversarial-repositioning obstruction; the missing input is a fixed-complexity
horizontal description of the moving Apéry zero (= horizontal Sato–Tate / theta-descent), open in every instance.
The wall is NOT mischaracterized. NEW positive deliverable stands: the Katz–Mellin VERTICAL program (central
anchor 548/548 verified). Awaiting Q_C (different sufficient condition), Q_D (Katz vertical), Q_E (theta-descent),
Q_F (calibration).

## 14. ChatGPT Q_C (2026-07-22) — candidates a/b/c dead, but TWO genuinely-different open reformulations

**(a) anatomy of b_n DEAD** — top-range Lucas makes `p|b_n ⟺ p|b_r`, so b_n's medium-prime anatomy is not
independent; size bound gives O(n/log n) (wrong scale). **Decisive no-go: the central binomial `C(2n,n)` —
coefficient of the ALGEBRAIC G-function `(1−4z)^{−1/2}` — has `~n/(6 log n)` medium prime divisors (Legendre,
v_p=1 for n/2<p<2n/3).** So NO theorem based only on integral / P-recursive / G-function / period / exp-bounded
can give `o(n/log n)`; any anatomy theorem needs special Apéry structure. (Luca–Shparlinski Apéry results are
lower-bound / density-one / below the p~n scale.)

**(b) radical/valuation/abc DEAD** — G-function theorems control DENOMINATORS (vacuous, b_n∈ℤ); fixed-S doesn't
cover a moving ~n/log n prime set; abc proves LARGE radical (wrong direction); uniform v_p bound irrelevant
(obstruction = # distinct primes). **(c) reciprocity r↔p DEAD** — no law swaps index and characteristic; all
identities (Lucas, reflection, finite-field hypergeometric/Mellin, central supercongruence) are same-prime.

**TWO genuinely-different OPEN reformulations (avoid the collision moment entirely) — record as live avenues:**
- **(MTG) Moving-target gcd sum.** `∀Q fixed: Σ_{q≤Q} Σ_{r≡n (q), r<n/(q+1)} log gcd(b_r, (n−r)/q) = o_Q(n)` ⟹
  pointwise P3.2. Does NOT mention Z_p / collisions / prime pairs. Connects to Corvaja–Zannier / Grieve–Wang
  (Trans.AMS 2020) / Xiao (Math.Z. 2024) moving-target gcd — but those need torus / C-finite (constant-coeff)
  structure with an exp-polynomial rep `b_r=Σ P_i(r)α_i^r`, which Apéry (P-recursive) LACKS. OPEN; needs a
  Subspace-Theorem / Diophantine-divisibility analogue for 2nd-order P-recursive sequences vs a moving LINEAR
  target. **Different technology class (Diophantine geometry), not analytic NT / trace functions.**
  Equivalent: a subexponential-height integer certificate `C_{n,Q}≠0` divisible by all bad primes,
  `log|C_{n,Q}|=o_Q(n)` (continuant certificates don't compress to subexp height).
- **(Diophantine) Sharp irrationality along the Apéry subsequence.** `|ζ(3)−a_n/b_n| ≥ q_n^{−2−o(1)}`
  (q_n=b_n/G_n) ⟹ `G_n=e^{o(n)}`, avoiding prime distribution entirely. = irrationality EXPONENT 2 along the
  Apéry subsequence; at least as deep (Rhin–Viola μ(ζ(3))≤5.514 gives only exponential). ORTHOGONAL route.

**Synthesis note:** the problem sits at the intersection where BOTH analytic NT (horizontal Sato–Tate /
trace functions, §13) AND Diophantine geometry (moving-target gcd / Subspace, MTG) fall just short — linked by
Fable #1's critical adelic ledger (subspace surplus = (λ−3)n = the ζ(3) irrationality margin, already spent).
MTG and sharp-irrationality are the two live NON-collision avenues; both open, both a documented research target.

## 15. ChatGPT Q_E (2026-07-22) — confirms wall, CORRECTS Fable #2 twice, tempers Katz vertical

**Confirms:** no theta-lift / Waldspurger / automorphic object locates the off-center Apéry zero across p;
central index (quadratic char → 8.4.a.a) is a fixed-parameter exception; general ω^{z_p} order grows with p,
no compatible motive. Genuinely open, needs new mathematics. Off-center is NOT a finite-branch set (matches the
measured z_p/p uniformity — "computations do not suggest a finite-branch description").

**CORRECTION 1 to Fable #2 (incompleteness ≠ per-prime conductor growth).** At each single prime, the Kummer
sheaf L_χ (χ=ω_p^{−z}) has rank 1, tame, BOUNDED geometric conductor INDEPENDENT of ord(χ); its cohomology has
bounded dimension. So the off-center incompleteness does NOT make the local conductor grow like p. The true
obstruction is HORIZONTAL COMPATIBILITY: char order varies with p, coefficient-field degree grows, monodromy
eigenvalues are roots of unity of unbounded order, and z_p is selected only after reducing mod the same p — no
fixed compatible family across primes. (Sharper than Fable's "incompleteness blocks the motive".)

**CORRECTION 2 to Fable #2 (supersingular attribution).** The x^{3/4} supersingular bound is GRH-conditional
effective Chebotarev/Lang–Trotter (Zywina), NOT Elkies's CM infinitude (which gives only infinitude, no power
saving). The decisive elliptic feature: a_p=0 is a FIXED-TRACE condition in ONE fixed 2-dim compatible system
+ Deuring/class-number after averaging. Apéry has NEITHER. So "unconditional theta descent gives horizontal
power saving" was mis-stated; even the elliptic power saving is conditional.

**TEMPERS the Katz–Mellin vertical program (Fable #2's §11 deliverable).** Katz's Convolution & Equidistribution
FIXES the sheaf and equidistributes normalized Frobenius as the EXTENSION FIELD grows, against CONTINUOUS class
functions. "Divisibility by the characteristic" (exact b_z≡0 mod p) is NOT a continuous class function of the
compact monodromy group — it is a local-limit/Lang–Trotter event. So Katz does NOT directly deliver the exact
zero-COUNT statistics (mean|Z_p|=1, Poisson); those are harder exact-vanishing events. The vertical program is
weaker than Fable #2 claimed (equidistribution of NORMALIZED traces yes; exact zero counts no). [Q_D pending for
the precise boundary of what Katz gives.]

**Viability (3 breakthroughs, none known):** finite rational-branch description of z_p/(p−1) (numerics say NO);
a special-cycle/class-number identity (Shimura variety Z(D)); or a horizontal residual-Mellin theorem (new
p-adic analogue of Katz). Net: confirms the wall, corrects the mechanism, tempers the positive deliverable.

## 16. ChatGPT Q_F (2026-07-22) — DEFINITIVE calibration: state of the art is the TRIVIAL bound

**The strongest unconditional POINTWISE (every-n) bound for `ω_{(X,2X]}(b_n)` is still the TRIVIAL
`O(X/log X)`.** No published theorem gives even `X/(log X)^{1+δ}` (one extra log) for Apéry. **Even
`ω_{(X,2X]}(b_n) ≪ X/log²X` would cross a real frontier — the first pointwise horizontal non-concentration
estimate of this kind for the sequence.**
- **C(2n,n) obstruction (2nd independent confirmation, cf. Q_C):** central binomial = diagonal of 1/(1−x−y),
  hypergeometric/holonomic/integer/exp-growth, yet `v_p=1` for n<p<2n gives `π(2n)−π(n)~n/log n` medium prime
  divisors. NO theorem from holonomic/diagonal/hypergeometric/automatic properties alone gives `o(π(X))`.
- **Lucas-sequence calibration benchmark:** rank-of-apparition (`p|U_n ⟺ ρ_p|n`, `ρ_p|p±1`) gives a genuine
  pointwise `Ω_X(n) ≤ X^{1/2+o(1)}` uniformly for X>√n (elementary, §6 of the answer). The power saving comes
  from `ρ_p|p±1` + the divisibility law — NOT from holonomicity. **This is exactly the fixed-quotient/rank law
  Apéry provably LACKS** — the benchmark pins the missing mechanism.
- Stewart (largest prime factor), Sanna (rank distribution), Luca–Shparlinski (density-one lower bounds),
  Bugeaud–Evertse (fixed S), BCZ/Fuchs (needs an independent 2nd sequence containing the Apéry prime factors —
  none known; Wronskian gives only coprimality) — all adjacent, differently quantified, none gives one log of
  saving.

**SIX-WAY CONVERGENCE (Fable#1, Fable#2, ChatGPT Q_A/B/C/E/F).** Pointwise P3.2 is beyond ALL current
unconditional technology — not by a log exponent, but by the total absence of a horizontal non-concentration
mechanism. The one mechanism that delivers such bounds (rank of apparition / fixed quotient) is provably absent
for Apéry (Integrality Criterion; SL₂ difference-Galois; measured entropy-positive Poisson zeros). The trivial
bound is the state of the art; even one log of saving is an open research frontier.

## 17. ChatGPT Q_G (2026-07-22) — higher CRT moments: sharper restatement of the wall, NOT an escape

**Verdict: higher factorial moments are a real sufficient condition but a methodological mirage with present
inputs.** Exact: `T_k(X)=Σ_{n<X²}(K(n))_k = Σ_{k-tuples}Σ_{z_i∈Z_{p_i}} 1[CRT lift < X²]`. Random scale
`≈X²·A_X^k` (A_X=Σ|Z_p|/p≈1/logX); CRT-uniqueness bound `T_k ≪ R_X^k ≈ X^k/log^k X`. **The gap is X^{k−2} and
GROWS with k.** "Product modulus > interval" only turns many-lifts into 0/1-lift; it does NOT show the unique
lift lands in [0,X²) with prob X^{2−k} — that IS the missing horizontal equidistribution. Adversarial alignment
`Z_p^adv={n₀ mod p, …}` saturates T_k for EVERY k (K(n₀)=π(X)); higher moments only DIAGNOSE alignment, don't
exclude it. No precedent supplies it: Gallagher conditional-Poisson (needs Hardy–Littlewood k-tuples = the hard
input), Maynard–Tao (needs Bombieri–Vinogradov, doesn't manufacture joint equidist from cardinality),
Matomäki–Radziwiłł (almost-all not every-n), EKR (star = aligned = EKR-extremal, wrong way). Large sieve WORSENS
(modulus ~X^k, spacing X^{2k}≫N).

**Quantified target (the round's most concrete forward handle):** with expected density A_X≈1/logX, **D_3 or D_4
suffices** (K(n)≪X^{2/3}/logX resp. X^{1/2}/logX); under only the proved |Z_p|≪p^{2/3} row bound need k≥7. And
a WEAKER power-saving form is enough: `T_k ≪ X^{k−η}/log^k X` (ANY η>0) ⟹ `sup_n K(n) ≪ X^{1−η/k}/logX`. **Most
realistic research target: any aggregate power saving in #{k-tuples whose CRT representative is anomalously
small} — but it must use genuine cross-prime Apéry arithmetic (continuants / resultants / K3).** Connects to the
MTG integer-certificate idea (Q_C).

## 18. ROUND-1 FINAL SYNTHESIS (2026-07-22) — SEVEN-way convergence; wall confirmed+sharpened; new deliverables

**Seven independent deep analyses** (Fable#1, Fable#2, ChatGPT Q_A/B/C/E/F/G; Q_D failed delivery) CONVERGE:
pointwise P3.2 is beyond ALL current unconditional technology. The state of the art is the TRIVIAL bound
`O(X/log X)`; even `X/log²X` is an open frontier (Q_F). Every route — second moment / large sieve (F4), per-prime
Weil/trace-function on Z_p (F2, √p floor), Stepanov (F5), automatic-over-primes (F3b), transfer-matrix bilinear
(F3a, Fourier completion), DFI/FKM (wrong quantifiers / no fixed-complexity realization), anatomy of b_n (C(2n,n)
no-go, twice), radical/abc (wrong direction), reciprocity r↔p (none), higher CRT moments (gap grows with k) —
FUNNELS to the same positional wall. The missing mechanism is a horizontal non-concentration / rank-of-apparition
analogue for the moving Apéry zero, = horizontal Sato–Tate for truncated weight-4 hypergeometric vanishing,
provably absent (no fixed quotient; entropy-positive Poisson zeros; SL₂ difference-Galois). The wall is NOT
mischaracterized — it is now characterized to unprecedented precision, and the adversarial-repositioning
obstruction was independently reconstructed THREE times (Fable, Q_A large sieve, Q_B digit-product, Q_G every-moment).

**This round's genuinely NEW products (beyond confirming the wall):**
1. Explicit second-moment/bilinear reframing pinning the wall to a concrete exponential-sum inequality (not prose).
2. F2 killed rigorously (√p noise floor vs O(1) target) — a clean mechanism, not prose.
3. **Katz–Mellin identification** (b_z mod p = Mellin transform of the Straub trace function at χ=ω^z) → a VERTICAL
   Poisson program (TEMPERED by Q_E: normalized-trace equidistribution yes; exact zero counts are Lang–Trotter
   events, harder). Central anchor Ahlgren–Ono `b_{(p−1)/2}≡a_p` VERIFIED 548/548; multiplicity-1 verified (1
   small-prime exception p=17). A bankable structural result even short of P3.2.
4. **Three concrete forward research targets, in DIFFERENT technology classes** (the value for future rounds):
   (a) **MTG** moving-target gcd sum for a P-recursive sequence vs a moving linear target — Diophantine/Subspace
       class (Grieve–Wang/Xiao need C-finite; Apéry P-recursive = open). [Q_C]
   (b) **Sharp irrationality** exponent 2 along the Apéry subsequence — orthogonal, ≥ as deep. [Q_C]
   (c) **Aggregate power-saving small-CRT discrepancy `D_k`** (any η>0) using Apéry continuant/resultant/K3
       arithmetic — "most realistic". [Q_G]
   (d) horizontal residual-Mellin theorem / special-cycle (Shimura Z(D)) identity. [Q_E]
5. Two self-corrections (full-interval useless; F2 usefulness) + two corrections OF Fable#2 by ChatGPT
   (incompleteness≠per-prime-conductor-growth; supersingular x^{3/4} is GRH-conditional) — verify-don't-transcribe
   working in both directions.

**Standing verdict:** conjecture TRUE (numerics + Katz–Mellin vertical support); pointwise P3.2 needs genuinely
new mathematics, now named to a point. The productive continuation is NOT more collision/automatic/trace-function
attacks (all certified funneling) but the four forward targets above — especially MTG (a) and D_k (c), which are
in the Diophantine/continuant technology class the analytic-NT campaign never fully engaged.

## 19. Round 2 open — MTG DEFLATED numerically before full campaign (2026-07-22)

Before committing to an MTG campaign, verified its q=1 core numerically (`scratchpad/verify_mtg.py`):
`W1(n) = Σ_{r<n/2} log gcd(b_r, n−r)`, using `gcd(b_r,m)=gcd(b_r mod m, m)`:
```
 n     W1(n)    W1/n
 200   394.2    1.97
 400   925.2    2.31
 800  2136.4    2.67
1500  4476.3    2.98
3000 10023.4    3.34   (GROWING, not →0)
```
**MTG-as-stated is FALSE: `T_Q(n)` is NOT o(n) — it is ~3n and growing.** Reason: `gcd(b_r, n−r)` includes ALL
prime factors of `n−r`, and SMALL primes p divide b_r with positive density, contributing
`Σ_r Σ_p (c_p/p)log p ~ O(1)·(n/2)` = LINEAR mass. ChatGPT Q_C's derivation (2.2)–(2.3) only proved
`W(n) ≤ T_Q(n)` (a lower bound) + "remaining primes ≤ n/(Q+1)"; it never checked `T_Q(n)` CAN be o(n). It cannot:
the small-prime part dominates linearly. So MTG is a valid inequality with a USELESS (never-small) upper object.
**The useful (large-prime, p>√n) part of the gcd IS exactly the original middle-prime W(n)** (o(n)); the
small-prime part is a linear over-count that must be excluded, and once excluded MTG = the original problem. The
CERTIFICATE version (subexp-height integer divisible by the bad LARGE primes) = the R14 construction/eliminant
gate, already CLOSED (D³ tax → height ≥λn, universal divisibility, zero localizing content).

**Verdict (pending Fable confirmation): MTG likely RE-COLLAPSES** — the naive full-gcd form over-counts (false),
and the large-prime-restricted form is the original W(n) / the closed R14 certificate gate. This is a correction
to §18's forward-target (a): MTG is NOT the promising different-class escape it appeared. The remaining
genuinely-open forward targets are the sharp-irrationality route (§18 (b), Diophantine, ≥ as deep) and the D_k
small-CRT discrepancy (§18 (c)) — both also hard, both linked to the same ledger. Fable dispatched for the exact
re-collapse step (MTG vs critical adelic ledger / 3n cap / D³ tax).

## 20. Fable MTG verdict — MTG CLOSED (rigorous), frontier UNCHANGED (2026-07-22)

Fable confirms the numerical deflation with an UNCONDITIONAL disproof + a verbatim re-collapse:

**(1) MTG false as stated (rigorous, unconditional).** Mod 5 (Gessel Lucas): `5|b_r ⟺ r` has a base-5 digit in
`{1,3}` (b_1=5≡0, b_3=1445≡0). For q=1, every `r≡n (mod 5), r<n/2` has `5|(n−r)`; among these `n/10` values,
those with `5∤b_r` (all base-5 digits avoid {1,3}) number `O(n^{log3/log5})=O(n^{0.683})`. So
`T_1(n) ≥ (n/10 − O(n^{0.683}))log5 ≈ 0.1609n` for EVERY large n. Full small-prime mass `≍ n log log n`.
[Confirms measured W1/n≈2–3.3 with a mechanism.]

**Repair re-collapses VERBATIM.** Restrict to `p>n/(Q+1)`: low digit `r_0=n mod p` is CONSTANT along the
progression `r≡n (mod pq)`, so by Gessel `p|b_r ⟺ p|b_{n mod p}`; thus `T'_Q(n)=Θ_Q(Σ_{p|b_{n mod p}}log p)+
O_Q(log n)` — an EQUIVALENCE, not a new sufficient condition. **Structural reason moving-target gives ZERO
leverage:** in Corvaja–Zannier the moving gcd produces DISTINCT quasi-independent Diophantine events; Apéry's
Lucas law FREEZES the low digit, so the `n/(2p)` events "p|b_r" on the progression are the SAME event
`p|b_{n mod p}` (perfectly correlated copies). MTG is a REARRANGEMENT of the original prime-indicator vector.
"The same self-similarity that gives the reduction annihilates every moving-target average."

**(2) No G-function Subspace analogue.** Rank blowup (mult. span of {b_n} not f.g.; f.g. ⟹ exp-polynomial ⟹
C-finite, contradiction; ζ(3) op = Calabi–Yau AESZ #1, Zariski-dense symplectic monodromy, no torus quotient =
3n cap). L¹-vs-L∞: Chudnovsky D³ / André–Bombieri are GLOBAL averages (proved the density version), provably
blind at a single (p, n mod p) → re-enters the density-vs-pointwise = quotient-free-locating frontier.

**(3) No o(n) certificate.** Residue vector ~n/2 bits, no known cross-p correlation; universal construction must
vanish on the whole zero locus → height ≥ λ·scale (R13/R14). Hankel wrong scale; continuant/resultant reproduce
Lucas (relate, don't compress); Beukers per-prime interpolations don't glue below Σlog p. Wronskian
`a_{s+1}b_s−a_s b_{s+1}=6/(s+1)³` isolates bad digits per p-line (p>s+1 can't divide both b_s,b_{s+1}) but
doesn't correlate across p — thins constants, not the n-scale. Only cross-p special-digit correlation =
Ahlgren–Ono apex (exhausted).

**(4) MTG joins the closed families.** Frontier UNCHANGED = quotient-free locating (turn the L¹ density/
André-size ledger into an L∞ statement at one (p, n mod p) = any cross-p correlation of `b_{n mod p} mod p`
beyond the Ahlgren–Ono apex). Genuinely new mathematics ("a Subspace theorem for G-operators", or a cross-p
residue-correlation theorem); no literature instance. **MTG neither reaches nor reformulates it.**

## 21. CAMPAIGN STATE (2026-07-22, Round 1+2) — EXHAUSTIVE closure, all avenues terminal
Every avenue now has a mechanism-level terminal verdict: collision 2nd moment (short-interval forced, adversary
saturates), per-prime Weil on Z_p (√p floor), Stepanov (entropy-positive), automatic-over-primes (single-digit +
repositioning-invariance), transfer-matrix bilinear (Fourier completion), DFI/FKM (wrong quantifiers / no fixed-
complexity realization), anatomy (C(2n,n) no-go), radical/abc (wrong direction), reciprocity (none), higher CRT
moments (gap grows with k), MTG moving-target gcd (unconditionally false + verbatim re-collapse), G-function
Subspace (rank blowup + L¹/L∞ blindness), integer certificate (entropy/D³ tax). Frontier = quotient-free locating
/ horizontal Sato–Tate for truncated weight-4 hypergeometric vanishing = a cross-p correlation of `b_{n mod p}
mod p`. Conjecture TRUE. Deliverables: (i) sharpened status resolution (7-way convergence, weight-4 Hasse
threshold, C(2n,n) calibration); (ii) Katz–Mellin VERTICAL Poisson program (bankable structural result, anchor
548/548, mult-1 verified, tempered by Q_E).

**Next move — build it ourselves (遇山开山, per Xiang 2026-07-22; see [[feedback_yushan_kaishan]]).** The
frontier is a CONCRETE construction target, not a stopping point: build the horizontal Mellin count for the
explicit Straub sheaf `K` (the sheaf whose trace function gives `b_z mod p`). Round 3 plan in
`Q32_KatzMellin_vertical.md` §3–4: Phase A = pin `G_geom` of `K` (transport the known Sp₄ monodromy of AESZ #1
to the Mellin sheaf; makes the vertical program rigorous), Phase B = construct the horizontal residual-Mellin
count directly from `G_geom` + a p-averaged trace estimate for OUR convolution sheaf (prove the specialized case
for `K` from its computable monodromy rather than seeking an off-the-shelf theorem). Attack directly.
