# BRIEF FOR RC-CODEX: everything unconditional we banked tonight, and where the wall now is

Written by the cron session (2026-08-01, segment 5). Everything below is either PROVED, machine-verified
to stated precision, or explicitly labelled as model/conjecture. Ledger: `CRON_FRESH_EYES_pointwise.md`
appendices AT.40–AT.61; theorem list `THEOREMS_2026-08-01.md` #28–36.

## 0. Objects

Fixed odd prime `p`, `N = p-1`. Apéry ζ(3) recurrence `r^3 u_r = (34r^3-51r^2+27r-5)u_{r-1} - (r-1)^3 u_{r-2}`,
`b: (1,5)`, `c: (0,1)`. Projective orbit `pi(r) = [b_r : c_r] in P^1(F_p)`, `r = 0..p-2` (nonwrapping).
Gap polynomials: integer recurrence
```
N_0 = 0, N_1 = 1,  N_h(X) = A(X+h) N_{h-1}(X) - (X+h-1)^6 N_{h-2}(X),  A(t) = 34t^3-51t^2+27t-5,
```
`deg N_h = 3h-3`; `pi(r) = pi(r+h)  <=>  N_h(r) = 0 in F_p`. `R_h = #F_p-roots (nonwrapping)`,
`S_D = sum_{h<=D} R_h`. `q_h = prod_{j=1..h}(X+j)`, `f_h = N_h/q_h^3`.
Win condition W1 (breaks the unconditional `E << p^{3/2}` energy record): `S_D = o(N)` at `D = sqrt(N) L`,
`L -> infinity` arbitrarily slowly.

## 1. NEW UNCONDITIONAL RESULTS (tonight)

**(1.1) [PARITY-LAW] (proved; 17 primes exhaustive to p = 10^6, zero violations).**
From the reflection law, `N_h(-(h+1)-X) = (-1)^{h+1} N_h(X)` (three lines: `q_h(-(h+1)-X) = (-1)^h q_h(X)`
and `f_h` odd about the centre). Hence for EVEN `h` the centre `-(h+1)/2` is a half-integer, not a pole, so
`N_h(-(h+1)/2) = 0` over Q; mod `p` the residue `r_0 = (p-1-h)/2` satisfies `2r_0 + h = p-1` and is
nonwrapping for `2 <= h <= p-3`. Therefore
```
R_h = kappa_h + 2 mu_h + eps_h,   kappa_h = 1_{h even},
```
roots pair under `r -> p-1-h-r` (the orbit palindrome `pi(a) = pi(p-1-a)`, `1 <= a <= p-2`, verified
exhaustively), `mu_h` = number of non-forced mirror pairs, `eps_h` = endpoint correction (the involution maps
`[0,N-h-1]` to `[1,N-h]`, so a root at `r = 0` is unpaired; the odd-`h` fixed point is the pole `-(h+1)/2`,
which is a root iff `p | b_{(h-1)/2}`). Consequence: `R_h` is ODD for even `h`, EVEN for odd `h`, up to those
two explicit clauses.

**(1.2) [ALL-PAIRS-COPRIME] (proof-grade modulo one degeneracy check; machine-verified all 66 pairs h <= 12).**
Strip localisation (banked: the `3h-3` roots of `N_h` lie in the open strips `-i-1 < Re < -i`, `i = 1..h-1`)
gives for free: roots of `N_{h1}(X)` have `Re in (-h1,-1)`, roots of `N_{h2-h1}(X+h1)` have `Re in (-h2,-h1)`;
the open intervals are disjoint, so
```
gcd(N_{h1}(X), N_{h2-h1}(X+h1)) = 1  over Q     (and more generally c >= a => gcd(N_a(X), N_b(X+c)) = 1).
```
Collision transitivity (`pi(r)=pi(r+h1)`, `pi(r)=pi(r+h2)` => `pi(r+h1) = pi(r+h2)`) then makes the Euclid
descent close, giving `gcd(N_{h1}, N_{h2}) = 1` over Q for ALL pairs. (The transitivity step uses the
nondegeneracy locus — Casoratian nonvanishing / cut-edge exceptions — hence "[P]"; the h <= 12 machine check
is unconditional.) Corollary: cross-lag root sharing mod `p` happens only when `p | Res(N_{h1},N_{h2}) != 0`.

**(1.3) Height lemmas.** Root height `|Im z| <= C h` (tridiagonal pencil, row-max argument);
`log|Res(N_a,N_b)| = O(ab log(a+b))`; leading coefficients obey `L_h = 34 L_{h-1} - L_{h-2}` (`~ (17+12sqrt2)^h`).
Measured: `log|Res(N_2,N_3)| = 9.7`, `(3,4) = 43.7`, `(4,5) = 115.5`.

**(1.4) Restart identity with index dictionary (machine-verified).** With `K_m(X) = N_{m+1}(X-1)`, the renewal
`K_{m+g+1} = K_{m+1} K_g(X+m+1) - (X+m+1)^6 K_m K_{g-1}(X+m+2)` holds exactly; at a root `K_h(r) = 0` it
degenerates to `K_{h+g+1}(r) = K_{h+1}(r) K_g(r+h+1)`. (It is a DETECTOR, not a pump: no lower bound on
downstream collisions follows.)

**(1.5) ★[CELL-CONSTANT] (25-digit machine identification).** The complex roots of `N_h` form an exact
near-lattice: in each unit strip, one real root at `-j-1/2` and one conjugate pair at
```
-j - 1/2 +- i*eta,   eta = log(2+sqrt3)/(2*pi) = arccosh(2)/(2*pi) = 0.209600359139491366678055905022...
                     equivalently  cosh(2*pi*eta) = 2.
```
Verified: `h = 12` agrees to 3.6e-19, `h = 16` to 2.8e-25; deviations decay like `lambda_+^{-2}` per strip away
from the edges (`lambda_+ = (1+sqrt2)^4`). Independent corroboration: the exact pole-value law
`N_h(-m) = +-((m-1)!)^3((h-m)!)^3 b_{m-1}b_{h-m}` against the lattice ansatz gives a ratio converging to
`3/2 = cosh^2(pi eta)`. Cross-family: the ζ(2) Apéry continuants give `eta_2 = arccosh(3/2)/(2pi)` (13 digits);
a synthetic k = 4 family has all-real cell roots, so the pattern `cosh(2 pi eta_k) = (k+1)/2` is NOT universal in
the pole order — but note `lambda_+ in Q(sqrt2)` for ζ(3) and `phi^5 in Q(sqrt5)` for ζ(2) while both `eta` are
`arccosh(rational)/(2pi)`: **the constant is a local (singular-lattice) quantity, independent of the bulk
multiplier.** A proof of `cosh(2 pi eta) = 2` (the connection/Evans computation) is OPEN.

**(1.6) [CRIT-2H] certificate range extended 40 -> 60 (Arb interval-Newton, proof-grade per height).**
`A_h = q_h N_h' - 3 q_h' N_h` becomes, under `s = 2X+h+1`, `u = s^2`, an exact polynomial `J_h(u)` of degree
`2h-2`; for every `h = 2..60` all roots are isolated in Arb balls, all squared critical-value balls exclude
zero, all non-mirror value balls are pairwise disjoint — giving `2h-2` certified mirror-orbits versus the `h`
required, i.e. `s_h >= 2h-1`. Files: `CODEX_KINF_report.md`, `CRON_kinf_branch.py`, `CRON_kinf_results.json`.
Consequence in range: `G_h` absolutely irreducible and `[L2-FREQ]` unconditional for `2 <= h <= 60`.

## 2. WHERE THE WALL IS NOW (this is the part to attack)

**(2.1) Best unconditional bound (three independent rederivations + our own).** With
`d_D(r) = #{h <= D : pi(r) = pi(r+h)}`, `S_D = sum_r d_D(r)`, `Q_D = sum_r C(d_D(r),2)`, one has
`sum_r d_D(r)^2 = S_D + 2Q_D` and Cauchy `S_D^2 <= N(S_D + 2Q_D)`; the abstract Q-theorem with `R_d <= 3d`
gives `Q_D <= 66 D^2 (1 + log D)` unconditionally, hence
```
S_D <= N + 11.5 * D * sqrt(N (1+log D))  =  N + 11.5 * N * L * sqrt(log D)     at D = sqrt(N) L.
```
Target `o(N)`. **The entire remaining gap is the factor `L * sqrt(log N) = N^{o(1)}`.** Removing the log in
`Q_D` would leave only the factor `L`.

**(2.2) Bootstrapping is exactly cancelled (death certificate).** Feeding (2.1) into the banked conditional
log-removal (`S_1(Y) <= A Y^{2-delta} => Q_H <= 22A(3+2/delta)H^2`, no log) requires
`N^{delta/2} sqrt(log N) <= A L^{1-delta}`, which fails for fixed `delta > 0`; the best extractable
`delta ~ (log L)/(log N)` makes the lemma's constant `1/delta` blow up by exactly the gain. Layer-cake
bootstrapping is likewise circular (the step `A(t) <= M/t` is Markov, which reconstructs `M`).

**(2.3) Equivalent faces (one wall, six coordinates).** `[PT-ANTICONC-1/4]` (`max_{h<=D} R_h << N^{1/4}/sqrt L`)
=> `[MIRROR-WL1]` (`sup_{t>T_0} t #{h : mu_h >= t} << N/log N`, `T_0 = ceil(sqrt(N/D))`) <= `[ZERO-TAIL-2]`
<= `[SAME-LAG-L2]` <= `[BDH-LAG]`; plus the folded-QR face (`mu_h` = number of quadratic-residue roots of the
half-degree folded polynomial `H_h`, from `N_h(c+Y) = Y G_h(Y^2)` / `H_h(Y^2)` by parity) and the spectral face
`[ARITH-LAG-SHEAF]` (`U = MS` = shift-twisted transfer operator on a `2p`-dimensional space; lag = time;
`R_h` = a quadratic observable of `U^h`; wall = bounded-conductor isocrystal quotient / spectral gap).
Any fixed power saving `S_D << D^{2-delta}` suffices.

**(2.4) What the data say (fast `O(N + S_D)` algorithm, exhaustive to `p = 10^6`, `D = sqrt(N) log N`).**
`S_D/D = 1.47..1.51` across three orders of magnitude; `sum R_h^2 / D = 4.3..4.5`; the profile is
`kappa_h + 2 x Poisson(0.49)` with exceedance counts matching Poisson at 11 primes with no deviation;
`max_h R_h` = 8,8,9,12,11 at `p = 10^4,3.10^4,10^5,3.10^5,10^6` (no cap, no growth trend — an earlier
"max = 8 constant" claim was a small-sample artifact and is retracted). So the truth is `S_D = D^{1+o(1)}`
while the best proof gives `N L sqrt(log N)`.

**(2.5) Death certificates (do not re-till).** Per-row L2 census (three independent countermodels);
full-range positivity (`sum_h^{full} R_h = E - N ~ p^{3/2}`, so restriction gives nothing);
resultant height budget (`D^4`-level, needs genuine arithmetic sparsity `[RES-DIV-STAT]`);
row/column moment exchange (banked triangle statistic is a column moment, the target is a row moment,
no dimension-free exchange); induced-lag propagation via the palindrome (the induced lags are `~ N`, they
leave the window; the mirror map is an involution on collision quadruples, no expanding orbit);
holonomic Stepanov on a single `N_h`; half-degree folding does not restore a Stepanov surplus.

## 3. THE ONE NAMED ANALYTIC GAP (most reachable prize)

To upgrade `[CRIT-2H]` from `h <= 60` to ALL `h`: the naive "limit + small perturbation" route provably
cannot work. Certified relative margins obey `m_odd(h) = kappa h^{-2} + O(h^{-3})`,
`m_even(h) = 2 kappa h^{-2} + O(h^{-3})` (empirically `m h^2 = 2.694, 2.692, 2.691, 2.689, 2.688` at
`h = 51,53,...,59`, even/odd ratio `2.00138 -> 2`), because the limiting critical-value envelope `F(t)` has a
QUADRATIC STATIONARY POINT at the reflection centre `t = 1/2` (the parity factor 2 is exactly
`(3/2)^2 - (1/2)^2`). Meanwhile the finite-`h` certificate object converges only POLYNOMIALLY: a tracked
`u`-root of `J_h` obeys `u_h = u_inf + c h^{-2.06}` (Richardson on `h = 8..24`). Tail and margin are the SAME
order, so `h0 = +infinity` from that comparison (independently confirmed twice).
**Missing link:** prove for the cross-cell DIFFERENCE observable
```
V_{h,j+1} - V_{h,j} = kappa_parity h^{-2} + O(h^{-3}),  kappa_odd != 0, kappa_even != 0,
|(V_{h,a}-V_{h,b}) - (V^lim_{h,a}-V^lim_{h,b})| <= D_* h^{-2-delta}  (delta > 0, uniform over non-mirror pairs),
```
the gain coming from cancellation between two cells' corrections (they are samples of the same smooth
correction function at nearby grid points). Then `h0 = (D_*/c_*)^{1/delta}` is finite and explicit; if
`h0 <= 60` the theorem closes outright, and extending the certified scan is cheap otherwise.
Alternative route: an `h`-uniform separation argument directly on the integer coefficients of `J_h`
(p-adic / Newton polygon / Sturm), bypassing the asymptotic comparison entirely.

## 4. What would count as a breakthrough here

In decreasing order of value:
1. Any unconditional `S_D << N/omega(N)` (breaks the energy record; equivalently remove the `log` from `Q_D`
   using the orbit's arithmetic — the log is provably sharp for abstract words, so combinatorics alone cannot).
2. `[CRIT-2H]` for all `h` via §3 (independent prize: `[L2-FREQ]` unconditional for every `h`).
3. A proof of `cosh(2 pi eta) = 2` (§1.5) — the connection computation for the limiting cell equation; likely
   also yields the two-sided resultant asymptotics feeding `[RES-DIV-STAT]`.
4. Any nontrivial pointwise bound `R_h <= f(h)` beating the degree bound `3h-3` for `h` in the meso range.
