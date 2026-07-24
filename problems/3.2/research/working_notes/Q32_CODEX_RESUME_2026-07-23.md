# Ramanujan Challenge Problem 3.2 — Codex resume log (2026-07-23)

## Logging doctrine

Record every serious route, including failed ones.  For each route keep:

1. the exact target and quantifiers;
2. why the route looked plausible;
3. the strongest computation or theorem actually checked;
4. the precise failure mechanism (not merely "open" or "did not work");
5. what was learned and what new input would justify reopening the route.

The goal remains the fully unconditional pointwise statement
`gcd(d_n a_n, d_n b_n) = exp(o(n))` for every `n`.

## 0. Inherited state checked on restart

Primary local records read:

- `Q3.2_density_theorem.md` (the long hostile-audited campaign);
- `Q32_RUN_LOG.md`;
- `Q32_MAIN_THEOREM_writeup.md`;
- `Q32_KatzMellin_vertical.md`, especially §§40–48;
- ChatGPT answers Q304, Q309, Q313, Q315, Q316.

Established and reusable:

- Wronskian control makes the small-prime contribution `o(n)`.
- The correct two-channel support law for `sqrt(n) < p <= n`,
  `n = qp+r`, is
  `p | G_n <=> a_q b_r = 0 (mod p)`.
- The companion `a_q` channel is pointwise controllable by height plus the
  short interval for primes with `floor(n/p)=q`.
- The lower `b_r` channel is the true pointwise obstruction.
- The Apéry zero fiber admits the unconditional continuant bound
  `|Z_b(p)| << p^(2/3)`, enough for a density theorem.
- The corrected continuants are
  `N_0=0`, `N_1=1`,
  `N_{h+1}(x)=P(x+h)N_h(x)-(x+h)^6N_{h-1}(x)`,
  `P(x)=34x^3+51x^2+27x+5`, with `deg N_h=3(h-1)`.

### 0.1 Possible regression in the later Q32 note (under hostile audit)

`Q32_MAIN_THEOREM_writeup.md` reverted to the false one-channel law
`p | G_n <=> b_{n mod p}=0`.  The earlier hostile audit records explicit
counterexamples `13 | G_26,G_27` while `b_0=1,b_1=5`.  The same later note
uses a continuant normalization with cubes rather than the audited sixth
power.

Consequences not yet accepted:

- the claimed `X^(2/3+o(1))` exceptional-set theorem as written;
- the claimed equivalence of full pointwise Problem 3.2 with non-ordinary
  primes of the single fixed form `8.4.a.a`.

These may have correct repaired versions after restoring the companion-channel
argument, but every dependency is being rechecked.  ChatGPT Q335 is the
independent hostile audit.

### 0.2 The latest continuant-first-moment lead

For ordinary tested primes, numerics gave
`S_p(H)=sum_{h<=H} #roots_Fp(N_h) ~= 2.5 H`, and the close-pair inequality
would turn a critical-scale estimate into a `sqrt(p)` zero-fiber bound.
Q313/Q316 correctly identified the missing input as a first-moment/singleton
bound; pairwise GCD/Cassini information controls repeated roots but not
singletons.

New quantifier concern found on restart: for each fixed `H`, primes splitting
the compositum of the splitting fields of `N_2,...,N_H` may force many roots
simultaneously.  This appears to refute an all-`(p,H)` estimate
`S_p(H) << H`, although it need not refute a critical-scale statement with
`H` comparable to a power of `p`.  The radical degrees, separability, and
quantifier order must be checked before this is called a counterexample.
ChatGPT Q337 is auditing exactly this issue.

## 1. ChatGPT round A (dm1–dm7)

All questions are non-load-bearing: local derivation and verification continue
while the web long-thinks run.  Answers are accepted only after checking their
equations against the audited normalization.

| Q | Angle | Why this angle | Status |
|---|---|---|---|
| Q335 | Hostile audit of support law and fixed-form reduction | Later notes conflict with an explicit counterexample and the earlier audited law | returned; checked against local counterexamples |
| Q336 | Direct single-`n` middle-prime attack | Full pointwise needs cross-prime coherence; avoids spending another round on average zero counts alone | dispatched |
| Q337 | Chebotarev/splitting-prime quantifier audit | The proposed uniform continuant first moment may be false even though generic-prime numerics are clean | returned; uniform all-`(p,H)` theorem rigorously refuted |
| Q338 | Closed form / orthogonal-polynomial identification of `N_h` | A special-function formula could add structure unavailable to generic sheaf/GCD arguments | returned; exact Casoratian/J-fraction found, positive-orthogonality route refuted |
| Q340 | Rational gauge / projective-cocycle expansion | A special factorization could supply genuine deterministic mixing rather than merely Zariski density | returned; canonical Borel-coset form found, standard gauge/expansion route blocked |
| Q341 | `8.4.a.a` versus the moving hypergeometric family | Tests whether the late fixed-form reduction silently replaced a moving residue by the central fiber | dispatched |
| Q342 | Broad independent completion attempt | Gives one tab freedom to look for a missed global denominator/valuation/resultant identity while respecting all dead-route verdicts | dispatched to degraded/disconnected `dm3`; do not re-fire until liveness is resolved |

## 2. Reopen conditions for inherited failed routes

- **Cassini / shifted GCD.** Reopen only if a new identity controls singleton
  roots or the first moment, not merely pairwise intersections.
- **Full-cycle transfer.** Reopen only with a normalization retaining prefix
  information; trace/determinant of the cleared full cycle are nilpotent and
  forget the statistic.
- **Generic Katz/FKM monodromy.** Reopen only after a bounded-complexity joint
  `(x,h)` object is constructed; fixed-`h` bounded-conductor results have the
  wrong quantifiers.
- **Weight-4 fixed-form ordinarity.** Reopen for the full problem only if the
  parameter mismatch is repaired.  A theorem about one central fiber cannot
  silently replace the moving residue `r=n mod p`.
- **Height / abc / generic prime-factor bounds.** Reopen only with a genuinely
  Apéry-specific cross-prime invariant; size bounds alone permit the adversarial
  section concentrated at one integer.

## 3. Q335 foundation audit and immediate local checks

### 3.1 Q335 verdict

Q335 independently recovered the two-channel support law

`p | G_n <=> a_q b_r = 0 (mod p)`, where `n=qp+r` and `sqrt(n)<p<=n`.

It confirmed:

- the one-channel law is valid only for `p>n/2` (`q=1`, `a_1=6`);
- the companion channel has total pointwise logarithmic weight
  `O(n^(2/3))`;
- the zero-fiber theorem `|Z_b(p)| << p^(2/3)` is unaffected;
- the density-one theorem is repairable by counting the lower channel and
  adding the companion channel separately;
- the single fixed modular form `8.4.a.a` records only the central residue
  `r=(p-1)/2`, not the moving residue `r=n mod p`;
- a central-residue contribution at fixed `n` satisfies
  `2n+1=(2q+1)p`, hence its total logarithmic weight is only `O(log n)`.

Therefore the exact remaining pointwise target is

`L(n) = sum_{sqrt(n)<p<=n, p | b_(n mod p)} log p = o(n)`.

Equivalently, after discarding `p <= n/(f(n) log n)` by prime counting, it is
enough to prove this uniformly in the slowly growing quotient range
`q=floor(n/p)<f(n) log n`.

### 3.2 Explicit support-law counterexamples

Exact recurrence arithmetic gives

- `a_2=351/4`, so `a_2=0 (mod 13)`;
- for `n=26=2*13+0`, `b_0=1`, yet `v_13(G_26)=1`;
- for `n=27=2*13+1`, `b_1=5`, yet `v_13(G_27)=1`.

These examples refute every proof step that uses
`p | G_n <=> p | b_(n mod p)` throughout the entire range `p>sqrt(n)`.

### 3.3 Newly isolated proof gap in the support law

Q335 points out that the current note's Wronskian relation

`D_q b_(q-1) - D_(q-1) b_q = 6/q^3 (mod p)`

does not by itself imply `D_q=a_q (mod p)` when `b_(q-1)=0 (mod p)`.
The support law is numerically and independently supported, but the written
derivation must use one of:

1. a direct proof of the full block congruence
   `p^3 a_(qp+r) = a_q b_r (mod p)`;
2. an explicit companion-sum identity;
3. a proof that `D_q` satisfies the complete Apéry recurrence in `q`.

Until one of these is written out, the block-Frobenius step remains a
foundation lemma to repair rather than a black box.

The earlier audited note does contain the needed independent route: the
explicit Apéry companion sum `(A)`, two Kummer valuation tables, and a
termwise proof of

`p^3 a_(qp) = a_q (mod p)`.

Its valuation cases were re-read on restart.  In particular, non-multiples
of `p` in the inner sum vanish after multiplication by `p^3`; multiples of
`p` reduce to the corresponding companion summands at index `q`; outer
indices carrying a base-`p` carry have at least one residual factor of `p`.
This is the correct way to close the Wronskian gap, subject to turning that
audit note into the final proof.

A fresh exact-rational computation checked the full congruence

`p^3 a_(qp+r) = a_q b_r (mod p)`

for all `1<=q<p`, `0<=r<p`, and
`p in {7,11,13,17,19,23,29,31,37,41,43}`: 7,948 identities, through
index 1,848.  This is verification, not a replacement for the Kummer proof.
The same run reconfirmed `v_13(G_26)=v_13(G_27)=1`.

### 3.4 Fixed-form route: exact failure and reopen condition

The fixed form controls only `r=(p-1)/2`.  At a fixed `n`, such a prime must
divide `2n+1`, so this slice is already harmless without any density theorem.
The route fails because it has the wrong parameter, not because the desired
ordinarity theorem is merely difficult.

Reopen only if a modular/motivic construction retains a free residue
parameter and produces a uniform theorem for the moving family
`r=n mod p`.

## 4. Q337: the ambient continuant first moment has false quantifiers

### 4.1 Exact radical-degree input

The tridiagonal splitting at `x=-a` gives, for `1<=a<=h`,

`N_h(-a)=(-1)^(a-1) ((a-1)!(h-a)!)^3 b_(a-1)b_(h-a)`.

All `b_j` are positive.  Hence the values at `-1,-2,...,-h` alternate in
sign, so `N_h` has a distinct real root in every intervening interval.
Consequently

`deg rad(N_h) >= h-1`

without any irreducibility or squarefreeness conjecture.  For even `h`, after
removing the forced reflection factor `2x+h+1`, the primitive radical still
has degree at least `h-2`.

### 4.2 Unconditional Chebotarev obstruction

Fix `H` and take the compositum of the splitting fields of
`rad(N_2),...,rad(N_H)`.  Infinitely many good primes split completely in
this fixed finite extension.  Such a prime can be chosen arbitrarily large,
in particular `p>H^2`, while

`S_p(H)=sum_(h<=H) #roots_Fp(N_h) >= sum_(h<=H)(h-1) >> H^2`.

Therefore no absolute estimate `S_p(H)<<H` can hold for every pair
`H<=sqrt(p)`.  In fact no `p`-independent `o(H^2)` bound has those
quantifiers.  Deleting the forced central factors does not repair it.

The earlier generic-prime numerics were not wrong; they sampled away from a
thin but rigorous complete-splitting obstruction.

### 4.3 What survives

The counterexample fixes `H` first and may require a prime astronomically
larger than `H`.  It does not refute a prescribed critical-scale statement

`S_p(floor(sqrt p)) << sqrt(p) polylog(p)`.

Generic effective Chebotarev cannot prove this repair: it bounds the least
split prime from above or counts split primes beyond a discriminant-dependent
range, whereas the desired result would exclude unusually small completely
split primes in a special growing tower.

More importantly, the block proof only uses actual returns of the
distinguished Apéry solution,

`A_p(H)=#{(x,h): b_x=b_(x+h)=0 (mod p), 1<=h<=H}`,

not every ambient root of every `N_h`.  Fixed-field splitting creates
potential continuant roots but does not create Apéry zeros.  Q351 asks dm4
for a critical-scale bound on `A_p(H)` and, explicitly, for a bridge (if any)
from such a bound to the cross-prime all-`n` target.  Without that bridge,
even a square-root zero-fiber theorem improves density estimates only.

### 4.4 Operational failure recorded

An attempted `ask-gpt.py --help` call was interpreted by the positional
client as a real question and submitted as Q349 to dm1.  Cause: the client
has no option parser; every argument after the script name is prompt text.
No mathematical conclusion is attached to Q349.  Reopen condition: none;
future invocations use source inspection for usage and pass only a prompt.

## Q338 quick audit: exact continuant structure, but no critical-scale bound

Q338 supplied and checked several exact identities:

- `N_h` is the Apéry Dirichlet continuant and a normalized Casoratian of the
  balanced `4F3` solution and its companion;
- for `m>=0`,
  `N_h(m)=((m+h)!/m!)^3/6 * (a_(m+h)b_m-b_(m+h)a_m)`, equivalently a positive
  ratio-difference sum;
- the transfer addition law, Dodgson identity, matching expansion, and an
  adjacent resultant express the existing gap certificates exactly;
- the centered parity factor for even `h` is forced, while a useful positive
  Wilson/Racah orthogonality identification is ruled out already by `N_3`'
  s four nonreal roots.

The exact formulas recover the `p^(2/3)` zero-fiber argument.  They do not
control isolated roots: the resultant only controls repeated/adjacent roots
and itself factors through earlier Apéry numbers.  A fixed-`h` splitting
field also gives complete splitting at infinitely many primes.  No checked
theorem gives the required aggregate bound when `h` grows with `p`.

Route status: **closed at the current technology boundary**.  Reopen only
with a genuine two-variable cancellation theorem for the nonautonomous
transfer word, or a resultant-dispersion identity that constrains singleton
roots and does not reintroduce the original factors `b_j`.

## Q340 quick audit: projective cocycle normalization obstruction

Q340 found an exact `SL_2` normalization
`B(t)=t^(-3)A(t)` with `det B(t)=1`, the Bruhat form
`B(t)=u(P(t)) a(t^3) w`, and the relative-Borel identity
`B(t)B(s)^(-1) in B_upper`.  A rational diagonal gauge flattening the hopping
would require `g(t)g(t+1)=t^3`; divisor propagation proves no rational solution
in characteristic zero or odd characteristic.  The weaker projective equation
`g(t)g(t+1)=t^6` has only a degree-Theta(p), everywhere-singular finite-field
solution, so it cannot give a bounded-complexity finite-field orbit theorem.

The sliding identity
`M_h(x+1)=B(x+h)M_h(x)B(x+1)^(-1)` is deterministic endpoint transport,
not independent convolution.  Thus group generation and sum-product/expander
theorems do not apply to the prescribed word.  This route currently yields
no two-prime dispersion or pointwise `K(n)` bound.

Reopen only with a prefix-sensitive invariant for the sliding word, or a
cross-prime identity that survives the Borel-coset collapse.

Operational note: after Q338/Q340 completed, two follow-ups (Q364/Q365) were
submitted while their client processes still reported the channels busy; the
allocator queued both on dm1 rather than selecting dm5/dm6.  This is a routing
failure, not a mathematical result.  They remain useful prompts and will be
read when they land; no duplicate submissions will be made until the queue
state clears.

## Q340 detailed audit: exact cocycle normal form, and why generic expansion does not apply

Q340 returned the canonical determinant-one normalization

`B(t)=t^(-3) [[P(t),-t^6],[1,0]]
     =u(P(t)) a(t^3) w`,

with projective dynamics

`phi_t(z)=P(t)-t^6/z`.

For any nonzero `s,t`,

`B(t)B(s)^(-1)`

lies in the upper Borel; explicitly its diagonal is
`((t/s)^3,(s/t)^3)`.  Thus the whole one-step family is a single
Borel coset.  The matrices may generate `SL_2`, but the actual word is a
prescribed arithmetic progression of point masses, not a random walk.
Helfgott/Bourgain--Gamburd expansion therefore supplies no prefix mixing.

The exact affine comparison is

`phi_t(z)=(t/s)^6 phi_s(z)+P(t)-(t/s)^6 P(s)`.

This correction is recorded because the text extraction of Q340 displayed
only the multiplicative term; direct substitution shows that an additive
term is essential.

### 5.1 Gauge obstruction

A diagonal `SL_2` gauge making both off-diagonal entries constant would
require

`g(t)g(t+1)=t^3`.

The valuation equations around the odd `p`-cycle force `2v_0(g)=3`, so no
rational solution exists over either `F_p(t)` or characteristic zero.
A projective gauge only requires `g(t)g(t+1)=t^6`, but every solution has
degree `Theta(p)` and a zero or pole at every finite-field parameter.  It is
not an invertible bounded-complexity conjugacy of the actual orbit.

This closes the natural constant-off-diagonal/Jacobi flattening.  A
non-diagonal gauge is not formally excluded, but it must genuinely mix the
two flags rather than hide the same scalar coboundary.

### 5.2 What was learned and what is still missing

The Borel-coset identity explains why Zariski density and generic sum-product
theorems were giving no theorem with the required quantifiers.  It does not
produce a cross-prime invariant, and even a critical-scale mixing result
would initially improve only a per-prime zero-fiber estimate.

Q359 asks whether reflection can pair a prefix with its suffix while
retaining the lost Borel coordinate and thereby produce a small global
certificate for the top-half condition `p=n-r | b_r`.  The route is to be
declared dead if the factorization reduces to the already-forced reflection
root or if every certificate has exponential height.

### 5.3 Local probe: lifting reflection modulo `p^2`

Potential target: the exact mod-`p` symmetry

`b_(p-1-r)=b_r (mod p)`

might strengthen at a zero `p|b_r`, thereby giving a second independent
congruence.  This would make the reflected pair useful rather than
tautological.

Fresh exact computation for every prime `7<=p<=499` found 47 reflection
pairs with a zero in the half interval.  Apart from the central identity
`p=11,r=(p-1)/2`, none satisfied

`b_(p-1-r)=b_r (mod p^2)`.

The first-order defects

`(b_(p-1-r)-b_r)/p (mod p)`

were nonzero and varied without a simple proportionality to the companion
value `a_r`, to `2r+1`, or to their quotient in the tested data.

Failure mechanism: at a zero, mod-`p` reflection supplies the same condition
twice, but it does not automatically lift one extra `p`-adic digit.  Thus it
does not create the hoped-for second certificate.  Reopen only if an explicit
first-order reflection formula exposes a new global invariant; mere
mod-`p^2` symmetry is empirically false.

## Bounded-order Mellin pruning

For `n=qp+r`, the lower-digit Mellin character satisfies

`omega_p^r = omega_p^(n-q)`,

because `r=n-qp=n-q-q(p-1)`.  Put

`d_p=ord(omega_p^r)=(p-1)/gcd(p-1,n-q)`.

If `d_p=d`, then

`p-1 | d(n-q)`.

For fixed `q,d`, every candidate `p` is therefore one plus a divisor of
`d(n-q)`, independently of whether the Apéry value vanishes.  Using
`tau(ab)<=tau(a)tau(b)`, the standard uniform divisor bound, and
`sum_(d<=D)tau(d)<<D log D`, one obtains

`sum_(q<=Q, d_p<=D) log p <= Q D n^(o(1))`.

After the established reduction to a slowly growing quotient range `q<=Q`,
choosing, for example, `D=n^(1-epsilon)/Q` makes this `o(n)`.  Hence the
true residual may be restricted to almost-maximal-order characters:

`d_p > n^(1-epsilon)/Q`,

or equivalently

`gcd(p-1,n-q) < Q n^epsilon`.

This rigorously generalizes the harmlessness of the quadratic central slice:
all bounded-order slices, and in fact a large power range of orders, are
negligible at fixed `n`.  It does not solve the problem because almost all
candidate primes may remain in the small-kernel regime.  Q362 asks whether
the cyclotomic/hypergeometric realization can exploit near-primitivity or
whether a zero at one residue prime fails to amplify across conjugates.

Quantifier correction: after the quotient reduction one only knows
`p > n/(Q log n)`.  Thus `gcd(p-1,n-q)<Q n^epsilon` implies the weaker but
correct lower bound
`d_p > n^(1-epsilon)/(Q^2 log n)`, up to harmless endpoint constants; the
stronger displayed `n^(1-epsilon)/Q` should not be quoted literally.

An exact integer recurrence sweep for `n<=500` checked the congruence
`r-(n-q)=0 (mod p-1)` for all 269 detected lower-channel bad pairs.  With the
illustrative split `Q=floor(n^(1/3))` and `D=floor(n^(4/5)/Q)`, the low-order
part had maximum observed weight `8.21` (at `n=257`); this is only a sanity
check, but it confirms the divisor parametrization and the intended pruning.

## Q341 fixed-slice audit

Q341 confirms the exact scope of Ahlgren--Ono: the weight-four form
`8.4.a.a` detects the unique reflection-fixed central coefficient and hence
the parity/unpaired part of `Z_b(p)`.  It does not detect the off-center
reflection pairs, which occur even at ordinary primes.  The moving object is
the Mellin/Kummer character family; its order and coefficient field vary
with `p`, so fixed-form Chebotarev, ordinary-reduction, K3 special-cycle,
and Deligne--Katz equidistribution theorems have the wrong quantifiers.

The central remainder slice is already harmless at fixed `n` because its
primes divide `2n+1`; the central quotient slice has at most one candidate.
Thus this route fails by parameter mismatch, not by an unsolved estimate
that would unlock Problem 3.2.  Reopen only with a horizontal theorem for
the moving Kummer character, or a bounded-complexity correspondence turning
every coefficient zero into a fixed-fiber condition.

## Q338 detailed closed form and its limit

### 7.1 Bankable identities

Let

`B(z) = _4F_3(-z,-z,z+1,z+1;1,1,1;1)`

be the meromorphic Apéry solution, and let `A(z)` be the companion normalized
by

`A(z+1)B(z)-A(z)B(z+1)=6/(z+1)^3`.

Q338 derived the exact Casoratian formula

`N_h(x) = [prod_(j=1)^h (x+j)^3 / 6]
          [A(x+h)B(x)-B(x+h)A(x)]`.

At an integer `m>=0`,

`N_h(m)=((m+h)!/m!)^3
        [a_(m+h)b_m-b_(m+h)a_m]/6`.

Equivalent descriptions are:

- the Dirichlet determinant of the symmetric tridiagonal Apéry difference
  operator;
- the numerator of its finite J-fraction;
- the path-matching continuant;
- the equality locus of two Apéry solution ratios.

The exact addition and Dodgson identities, and the adjacent resultant,
are also explicit.  In particular,

`Res_x(N_h,N_(h+1))
 = (-1)^(h(h-1)/2) prod_(j=1)^(h-1) [(j!)^3 b_j]^6`.

This explains precisely why adjacent common roots are controlled by earlier
Apéry zeros.

### 7.2 Orthogonal-polynomial recognition fails

The family is not a Wilson/Racah/dual-Hahn family with a fixed spectral
variable.  Matching the Wilson recurrence would require the spectral
parameter to depend on the recurrence index.  More decisively, centered
`N_3` has four nonreal roots, so no useful positive-measure orthogonality or
real interlacing is available.

The gap generating function is D-finite but satisfies a sixth-order
equation; it does not put `(h,x)` on a fixed bounded-degree algebraic curve.
Fixed-`h` irreducibility is also useless for an all-prime root bound because
complete-splitting primes exist.

### 7.3 Exact failure mechanism and reopen condition

The Casoratian, Dodgson, and resultant identities detect repeated or
clustered roots.  They are compatible with a population of isolated roots,
which is the empirically dominant case.  Consequently they recover the
`p^(2/3)` close-gap theorem but give no aggregate
`sum_(h<=H)#roots(N_h)=o(H^2)` at the critical scale.

Reopen this route only if one can turn the exact Casoratian into either:

1. a theorem dispersing isolated roots across `h`; or
2. a cross-prime norm/certificate at the fixed integer `n`.

Recognition inside the Askey scheme, by itself, is closed.

### 7.4 Local verification of the negative-integer factorization

Exact SymPy expansion checked
`N_h(-a)=(-1)^(a-1)((a-1)!(h-a)!)^3 b_(a-1)b_(h-a)` for all
`2<=h<=8` and `1<=a<=h` (35 identities), with degrees
`0,3,6,9,12,15,18,21`.  This confirms the Chebotarev radical-degree input
and the continuant normalization; it does not alter the singleton-root
obstruction.

## 9. Q351: exact same-prime target and the unavoidable cross-prime bridge

Q351 sharpened the surviving continuant route.  Each regular root of `N_h`
has a unique projective solution direction (a color); the Apéry zeros are one
distinguished color class.  Thus Q337's complete-splitting obstruction attacks
the ambient first moment but not the Apéry color.

The credible same-prime targets are either:

- a direction-energy bound at `H=sqrt(p)`, or
- a common-root bound for pairs of short continuants.

Neither is currently proved; Cassini and Christoffel--Darboux only exclude
local adjacent/triangular degeneracies, while isolated returns remain.

Q351 also gives a decisive logical audit: even a square-root bound for every
`|Z_p|` or a bound on the Apéry short-return count does **not** imply the
pointwise moving-residue theorem.  A reflection-invariant abstract set with
one prescribed residue per prime can have tiny same-prime incidence while
every prime hits one fixed integer.  A valid completion additionally needs
either a centered two-prime variance power saving or a five-prime CRT
factorial-moment estimate.  Ordinary large sieve bounds retain the full local
mass and do not supply this saving.

Q370 is now asking dm4 for the most concrete restricted version of this
cross-prime theorem.

### 9.1 Near-primitive numerical obstruction

An exact recurrence sweep for primes `p<1000` found 163 off-center Apéry zeros;
70 (43.0%) have `gcd(p-1,r)=1`, i.e. genuinely primitive Mellin order.  Thus
the low-order character pruning is a real theorem but cannot be extended by a
simple claim that primitive characters never vanish.

## 10. Q342: forced-factor obstruction to the natural global resultant

Q342 audited the factorial-scaled determinant identities and the reflected
zero-pair route.  For a noncentral zero `t` and its reflection, set
`h=p-1-2t`.  The gap certificate always has

`N_h(t) = (2t+h+1) C_h(t) = p C_h(t)`.

This is an exact central linear factor, independent of whether `b_t=0`; hence
the resultant eliminating `N_h(t)` and `p=2t+h+1` is identically zero.  The
certificate cannot couple distinct bad primes.

Dividing by the factor does not recover support information.  The quotient is
an explicit linear combination of the normalized endpoint layers
`beta_t=b_t/p` and `beta_(p-1-t)=b_(p-1-t)/p`, so it is invisible modulo `p`.
Q342's exact example is

`p=17, t=3, h=10, b_3=1445=5*17^2`,
`N_10(3)=34 (mod 17^2)`, hence `C_10(3)=2 (mod 17)`.

My independent recurrence computation reproduced the residue sequence
`(0,1,215,110,58,11,6,186,130,287,34)` modulo `17^2` and the same quotient.
Thus a genuine reflected Apéry zero pair does not force a second `p` factor.

The route is closed unless one proves a new horizontal theorem on the
normalized `p`-adic layers, or finds a second independent certificate not
sharing the reflection factor.  Q376 is asking whether higher `p`-adic
supercongruences can supply such a certificate.

## 11. Q359: reflection folding gives a new fixed-gcd formulation

Q359 retained the prefix at the full `SL_2` matrix level.  The exact reflected
step identity telescopes, but the prefix cancels by `M^T w M=w`; it yields
only the known palindromy and no second scalar condition.  This closes the
matrix-reflection rescue sharply.

It does give a useful non-tautological reformulation.  If `n=qp+r` and the
lower channel has `p|b_r`, fold `r` to `j=min(r,p-1-r)`.  Then

`p=(n-j)/q` or `p=(n+1+j)/(q+1)`,
`j<=(n-q)/(2q+1)`, and `p|b_j`.

Therefore every lower-channel bad prime divides both `b_n` (by Lucas) and
the lcm of the first third of the sequence.  Define

`C_n^fold=gcd(b_n, lcm(b_0,...,b_floor((n-1)/3)))`.

The estimate `log C_n^fold=o(n)` would be sufficient for the pointwise
lower-channel theorem, without a cross-prime collision moment.  No theorem
currently proves it: the lcm is a new moving-target object, and generic
P-recursive gcd results do not match Apéry's non-C-finite recurrence.

The exact matrix algebra explains why the fold is the strongest consequence
of reflection alone; it cannot bound this gcd.  Q378 is testing this new
integer-gcd route.

### 11.1 Numerical sanity check

For `n<=120`, exact integer recurrence and lcm/gcd arithmetic gave
`log C_n^fold/log b_n <= 0.1182` (maximum at `n=5`); for the last ten indices
the gcd logarithms were between `0` and `10.27`, while `log b_n` was
`383`--`414`.  This supports the conjectured sublinear common part but is
not a proof.

## 12. Q362: near-primitive cyclotomic norm audit

Q362 confirms the order split and closes the proposed norm-amplification
rescue.  For `g=gcd(p-1,n-q)` and `d=(p-1)/g`, the prime `p` splits
completely in the relevant cyclotomic/trace field.  The decomposition group
at the chosen prime is trivial, so one simple Apéry zero contributes exactly
one factor `p` to the minimal norm, not `p^(phi(d))`.  Conjugation moves the
character and the prime above `p` together.

The Weil height grows proportionally to the trace-field degree, cancelling
any putative cyclotomic gain.  Gross--Koblitz/Stickelberger valuations have
minimum zero for every character; `p|b_r` is cancellation among unit terms,
not an order-dependent slope.  The numerical primitive-zero sweep is
consistent: 70 of 163 zeros for `p<1000` had primitive order.

Thus the low-order character pruning is theorem-level progress, but the
near-primitive remainder still requires a fixed carrier, horizontal
non-ordinary-prime estimate, or cross-prime CRT theorem.  Q382 is now
attacking the folded `q=1` gcd directly.

At the earlier checkpoint Q336, Q365, Q368, Q370, Q376, and Q378 remained in
the bridge queue/long-think state; no conclusions were inferred from their
absence of output.  Q364 and Q362 have since landed and were audited below.

## 13. Q364: Casoratian certificates are exact but non-compressive

Q364 gives two exact certificates for a lower-channel prime in a fixed
quotient slice: consecutive Casoratian values `Lambda_m(n)` and
`Lambda_(m+1)(n)`, and candidate-specific endpoint continuants.  The
two-value map is invertible modulo every relevant `p` because its determinant
is the Apéry Wronskian unit.  Thus simultaneous divisibility is only a change
of coordinates for the original row `(d_n a_n,d_n b_n)`.

The endpoint certificate evaluates modulo `n-r` to `(r!)^3 b_r`; it therefore
reproduces the original moving condition exactly.  Fixed-q certificates have
height `Theta(n)` for bounded `q`, and multiplying endpoint certificates has
height `Theta(n^2 log n)`.  Exact resultants insert the earlier `b_j` factors
rather than removing them.  Q387 is testing whether three or more shifted
Casoratians can evade this height/coordinate obstruction.

Operationally, Q336's connector failure is now confirmed in bridge status
(`dm2: failed=1, processing=0`).  It was not silently re-fired; Q387 is a
new, narrower multi-shift certificate question sent after confirmation.

## 14. Local q=1 mod-p^2 Lucas-defect probe

After reflection folding, both q=1 branches have `j<=n/3` and `p>2j`:

`n=p+j`, or `n=2p-1-j`.

For every zero `p|b_j` with `p<1000` and `0<=j<=(p-1)/2`, exact integer
recurrence arithmetic tested the two natural Lucas defects

`b_(p+j)-5b_j`, and `b_(2p-1-j)-5b_j`.

There were 82 folded zero pairs.  In each branch, 81 defects had exact
`p`-adic valuation 1 and only the central `p=11,j=5` case had valuation 2.
In particular `p=17,j=3` has `v_p(b_j)=2` but both Lucas defects have
valuation 1.

Failure mechanism: a lower-channel zero does not generically promote the
one-digit Lucas congruence to a second `p`-adic digit.  Thus the most natural
mod-`p^2` difference certificate is unavailable even when `b_j` itself has a
double zero.  Reopen only with a different normalized defect or a theorem
using several independent shifts; fixed-depth supercongruence would improve
only a constant factor and cannot by itself yield `o(n)`.

## 15. Local rank-two collapse for multi-shift Casoratians

For fixed `n`, define

`Lambda_m(n)=d_n(a_n b_m-b_n a_m)`.

As a sequence in `m`, `Lambda_m(n)` satisfies the same second-order Apéry
recurrence, because it is a fixed linear combination of `a_m` and `b_m`.
Consequently the entire multi-shift family has rank two:

- any two consecutive values are coordinates for the original row
  `(d_n a_n,d_n b_n)` through the Wronskian-unit change of basis;
- every third and later value is a recurrence combination of those two;
- all exterior minors of order at least three vanish identically.

For coefficients `c_m`,

`sum c_m Lambda_m(n)
 = d_n[a_n sum c_m b_m-b_n sum c_m a_m]`.

Thus a nonzero multi-shift certificate merely chooses another vector in the
same two-dimensional solution space.  Choosing coefficients that exactly
eliminate one column via the Wronskian recovers a scalar multiple of
`d_n a_n` or `d_n b_n`; choosing a recurrence relation gives zero.

This is a structural obstruction to gaining independent congruences merely
by adding more `m` values.  A possible reopening would need nonlinear
operations with a provably smaller height, or data from a genuinely new
solution/motive rather than further shifts of the same rank-two system.

## 16. Primitive-divisor literature check for the folded gcd

The folded certificate

`gcd(b_n, lcm_(j<=n/3)b_j)`

is the old-prime part of `b_n`: it retains primes already appearing in the
first third of the sequence.  A targeted literature search found primitive
divisor/Zsigmondy theorems for:

- powers, Lucas/Lehmer and other C-finite divisibility sequences;
- elliptic divisibility sequences (often using their strong-divisibility
  or rank-of-apparition structure);
- arithmetic dynamical orbits, frequently conditional on Vojta.

No theorem found treats Apéry numbers or a general nonautonomous
P-recursive sequence strongly enough to bound this old-prime part.  The
available holonomic representation theorems concern definability/closure,
not prime apparition.  This matches the structural diagnosis: Apéry lacks a
fixed quotient/rank-of-apparition law.

Reopen the primitive-divisor route only with an Apéry-specific apparition
principle or a new Zsigmondy theorem for second-order P-recursive sequences
with polynomial coefficients; merely proving the existence of one new prime
per term would still be far weaker than showing the entire old part has
subexponential weight.

## 17. ChatGPT workhorse outage and lost follow-up batch

The follow-up batch Q365, Q368, Q370, Q376, Q378, Q382, and Q387 produced no
answer files.  This was an infrastructure failure, not negative mathematical
evidence.

Read-only inspection of the `ask-gpt-git` bridge established:

- the bridge server remained alive continuously and its browser connectors
  continued to heartbeat;
- the last successful Notion delivery in `runs.log` was at 05:33;
- from 06:06 onward, completed `ask-gpt.py` client runs returned
  `PENDING:pending` with zero bytes and then exited at their deadlines;
- Q365, Q368, Q370, and Q376 later failed after the bridge's 2400-second
  processing timeout;
- Q387 failed after five dispatch retries because the composer remained busy
  for 60 seconds each time;
- Q378 and Q382 remained stuck in `processing`, with no answer;
- the same failure pattern affected unrelated queues, so it was not caused by
  the Problem 3.2 prompts.

Diagnosis: the browser extensions were alive enough to heartbeat and sometimes
claim work, but the ChatGPT-tab automation/composer/result-delivery layer was
not completing jobs.  Since the original polling clients had already timed
out, no process remained to write a late answer to `/tmp/gpt_Q*.md`.

Operational lesson: do not treat a connector heartbeat as proof that a tab can
dispatch and retrieve an answer.  Before resubmitting this batch, first run one
controlled end-to-end smoke test on a known idle `dm` tab, verify a nonempty
git-drop/Notion result, then resubmit only the failed question IDs.  Do not
clear the global bridge queue: it contains unrelated work.

## 18. Q378: the folded gcd target is too broad

Q378 eventually landed in Notion after its original polling process had exited.
Its strongest unconditional result is the following old-prime-window bound.
For

`R_n(J)=prod {p>sqrt(n): p|b_n and p|b_j for some j<=J}`,

distinctness of the primes gives

`log R_n(J) <= sum_(j<=J) log b_j
             = (log(17+12 sqrt(2))/2) J^2 + O(J log J)`.

Consequently every window `J=o(sqrt(n))` contributes `o(n)`.  This is a real
partial theorem, but folding needs `J` of order `n` (already `n/3` for `q=1`),
where this estimate becomes useless.

The exact obstruction is Lucas propagation.  If `n=qp+j`, then

`b_n = b_q b_j (mod p)`;

on the reflected branch the palindromic congruence gives the same implication.
Thus `p|b_j` already forces `p|b_n`.  Membership in
`gcd(b_n,b_j)` is not a second independent divisibility condition.  The full
folded gcd also contains unrelated old factors and prime powers, so it is
strictly broader than the support needed for Problem 3.2.

The sharper remaining target is therefore a moving-prime divisor estimate for
the two exact linear forms

`p=(n-j)/q` and `p=(n+1+j)/(q+1)`,

with `p` prime and `p|b_j`, summed over the permitted `j` and fixed small `q`.
This is an Apéry-specific prime-divisor sieve problem; recurrence gcd,
continuant, primitive-divisor, and generic C-finite Subspace-Theorem arguments
do not match it.

## 19. Polling-script failure and repair

Q378 also supplied a concrete counterexample to the old polling policy:
connector delivery can land after the bridge client has stopped waiting.
Inspection of `scripts/ask-gpt.py` found three interacting defects:

1. the overall caller deadline was 90 minutes, while a task may spend substantial
   time pending before extended thinking begins;
2. bridge status `failed` was treated as terminal even though the independent
   Notion connector could still create the requested page later;
3. exact Notion matches older than one hour were discarded despite the random
   `drop_id`, making later recovery unnecessarily fragile.

The active workspace polling script was changed to:

- use a six-hour default overall deadline;
- continue polling an expected connector for four hours after bridge failure;
- inspect browser capture before stopping on a terminal bridge status;
- accept an exact `Q<number> <drop_id>` Notion page regardless of age.

`python3 -m py_compile` and `git diff --check` both passed.  This repairs future
long-running calls.  It does not retroactively revive already-exited pollers, so
the remaining Q365/Q368/Q370/Q376/Q382/Q387 pages must still be recovered by
exact-title search if they land.

## 20. Recovered late answers Q365, Q368, Q370, Q376, and Q382

Exact-title Notion search recovered five pages that the dead pollers had missed:

- `/tmp/gpt_Q365.md`: the Borel/sliding cocycle is characteristic-local and
  time-synchronous.  Distinct `p,q>sqrt(n)` necessarily use distinct residues
  `n mod p` and `n mod q`; universal transfer identities act componentwise over
  `F_p x F_q` and do not synchronize those two orbit times.
- `/tmp/gpt_Q368.md`: the moving-character order split is a genuine
  unconditional reduction.  With
  `d=(p-1)/gcd(p-1,n-q)`, the total logarithmic weight for `q<=Q,d<=D` is
  `<< Q D n^{o(1)}`.  Taking `D=n^{1-epsilon}/Q` makes it sublinear.  The
  residual characters have small kernel, but cyclotomic conjugation gives no
  amplification because `p` splits completely in `Q(zeta_d)`.
- `/tmp/gpt_Q370.md`: a power-saving centered second moment is incompatible
  with the Bernoulli diagonal.  Even a natural optimal variance does not give
  the required pointwise bound.  Conditional on a square-root same-prime zero
  bound, the first robust fixed-order bridge is a five-prime factorial moment.
- `/tmp/gpt_Q376.md`: after removing the forced reflected-continuant factor,
  the quotient is exactly a unit multiple of Gessel's derivative defect:
  `C_h(t)b_{t+1} = -Q_h(t) dot(b)_t (mod p)`.  Hence further division is
  possible exactly on the exceptional mod-`p^2` reflection digits; subsequent
  quotients expose unconstrained normalized beta layers.  This closes the
  formal first-order p-adic quotient route.
- `/tmp/gpt_Q382.md`: for `q=1`,
  `rad_(n/2,n](C_n^fold)=rad_(n/2,n](b_n)`.  Thus the folded gcd adds no
  independent condition at the support level.  However the folded indices
  `j<=n/log^2(n)` contribute only `O(n/log n)` by direct candidate counting,
  so the hard `q=1` primes lie strictly inside
  `(n/2+n/(2log^2 n), n-n/log^2 n)`.

Q387 had not landed when the exact-title recovery was run.

Q410 later exhausted its dispatch retries because the dm2 composer remained
busy; it never produced a Notion page or browser capture.  The threshold
optimization question was resubmitted once, as Q416/dm2, after the channel
returned idle.  The other Q409/Q411--Q415 tasks remained in processing.

## 21. New kernel identity after folding

Combining the order split with reflection produces an exact simplification not
used in Q368:

`gcd(p-1,n-q)=gcd(p-1,j)`,

where `j=min(r,p-1-r)` and `n=qp+r`.  Indeed:

- direct branch: `n-q=q(p-1)+j`;
- reflected branch: `n-q=(q+1)(p-1)-j`.

Writing

`p-1=g d`, `j=g a`, `gcd(a,d)=1`,

therefore gives:

- direct: `n-q=g(qd+a)`;
- reflected: `n-q=g((q+1)d-a)`.

For `q=1` these become the fixed-linear relations

`n-1=g(a+d)` and `n-1=g(2d-a)`.

An exhaustive exact check over all `n<5000` and all prime
`sqrt(n)<p<=n` verified 1,726,830 tuples, including both branches.

This parameterization combines the proved small-order saving with a divisor
condition `g|n-q`.  It does not itself bound the residual zeros `p|b_(ga)`, but
it is a sharper starting point: the remaining near-primitive obstruction can be
studied at fixed `q`, fixed divisor `g`, and one linear relation between the
coprime character parameters `a,d`.

## 22. Follow-up batch after polling repair

Seven new tasks were dispatched with the repaired six-hour poller:

- Q409/dm1: exploit the folded kernel parameterization for fixed `q,g`;
- Q410/dm2: optimize all quotient/order/edge thresholds and count the residual;
- Q411/dm3: attack the `q=1` reflected relation `n-1=g(2d-a)`;
- Q412/dm4: attack the direct fixed-sum relation `n-1=g(a+d)`;
- Q413/dm5: seek a batch resultant/determinant compression across all `j`;
- Q414/dm6: develop or rule out a diagonal large-divisor sieve;
- Q415/dm7: seek a bounded-degree binomial/finite-field rewrite for the linear
  moduli.

## 23. Completion-ack deadlock and durable late recovery

The first Q409--Q415 dispatch remained pending even though every `dm` connector
reported healthy.  The decisive observation was that the channels had stopped
calling `/api/pending`: connector-health reports were still arriving, but the
task-dispatch loop considered the tabs busy with the old Q365--Q382 work.

The missing transition was `/api/gitdrop-done`.  The original pollers had
exited before the Notion pages appeared, so they never acknowledged those
successful deliveries to the bridge.  After exact-title recovery, completion
acks were sent for Q365, Q368, Q370, Q376, Q378, and Q382.  Within seconds all
seven new tasks changed from `pending` to `processing`.  This experimentally
identifies the deadlock:

`late Notion success -> dead poller -> no completion ack -> extension remains busy
-> no /api/pending -> later tasks never dispatch`.

The polling script was therefore strengthened again.  Every future submission
now spawns a detached low-frequency late watcher.  It:

- survives the foreground caller;
- searches by the exact random `Q<number> <drop_id>` title for up to 24 hours;
- waits for the required `ANSWER Q<number> <drop_id>` first line, rather than
  treating a newly searchable but still-empty page as complete;
- writes `/tmp/gpt_Q<number>.md`;
- sends `/api/gitdrop-done`, releasing the extension's busy state;
- nudges the owning tmux window once when it recovers a previously absent file.

The foreground poller now also requires the answer echo line before declaring a
Notion page content-ready.  Syntax and diff checks pass.

## 24. Initial data for the folded kernel parameterization

Exact modular recurrence computation for all 668 primes `3<=p<=5000` found 305
folded first-half zeros.  Their kernel factors `g=gcd(p-1,j)` began:

`g=1:123`, `g=2:85`, `g=3:19`, `g=4:18`, `g=6:9`.

Thus small `g` is the generic empirical regime, not an exceptional set:
approximately 60.3% of the zeros had `g<=p^0.1`, 78.7% had `g<=p^0.2`, and
87.2% had `g<=p^0.3`.  The new parameterization is algebraically sharper but
candidate counting cannot discard its small-kernel range.

The same q=1 data already exhibits genuine cross-prime pile-ups after the
direct/reflected maps.  For example `n=321` receives three distinct reflected
hits:

- `(p,j,g,a,d)=(179,36,2,18,89)`;
- `(193,64,64,1,3)`;
- `(211,100,10,10,21)`.

Other tested `n` receive two distinct near-primitive hits with `g=1`.  Therefore
neither fixed `n`, the two linear fold relations, nor coprimality `g=1` forces
uniqueness.  Reopen the parameter-counting route only with an Apéry-specific
restriction on the zero condition `p|b_(ga)`.

## 25. Exact q=1 primorial-gcd reformulation

For `n/2<p<=n`, write `n=p+j`.  Gessel--Lucas gives

`b_n = b_1 b_j = 5b_j (mod p)`.

Since `p>5` in the asymptotic range, the q=1 lower-channel condition is
therefore exactly

`p|b_j <=> p|b_n`.

Consequently the direct q=1 contribution is the fixed-term quantity

`S_1(n)=log gcd(b_n, prod_(n/2<p<=n) p)`.

This is sharper than the folded-old-prime formulation: primes larger than
`n` which divide `b_n` are irrelevant, and primes at most `n/2` belong to
higher quotient slices.  A temporary numerical probe mistakenly summed all
prime factors `p>n/2`, including `p>n`; after imposing the required upper
bound `p<=n`, the apparent linear counterexample vanished.  For `10<=n<=45`,
the only top-half hits found were

`(n,p)=(16,11),(20,17),(27,19),(29,19),(30,17),(39,31)`.

This computation is only a sanity check.  The proof problem is now visibly a
small-prime-part theorem for one exponentially growing Apéry integer:

`log gcd(b_n,n#)=o(n)`

would be more than enough, while the required top-half version is weaker.
No general height argument proves either statement: `log b_n=Theta(n)` and
the primorial interval itself has logarithm `Theta(n)`.

## 26. 2025 truncation-square theorem: new route under audit

Caruso--Fürnsinn--Vargas-Montoya--Zudilin prove for

`A_p(t)=sum_(j=0)^(p-1) b_j t^j in F_p[t]`

that

- `A_p=B_p^2` for `p=1,5,7,11 (mod 24)`;
- `A_p=(t^2-34t+1)B_p^2` for `p=13,17,19,23 (mod 24)`.

This is an exact theorem, not the character-side Hasse-polynomial heuristic.
It is relevant because q=1 asks for a coefficient zero of `A_p` in its first
half.  The first obstruction is also exact.  After choosing the constant term
of `B_p`, its coefficients below half the degree are forced recursively to be
the reductions of one of two fixed formal power series:

`sqrt(f_b(t))`, or `sqrt(f_b(t)/(t^2-34t+1))`.

For example

`sqrt(f_b(t))=1+(5/2)t+(267/8)t^2+...`.

Thus `b_j=0 (mod p)` becomes the vanishing of one convolution coefficient.
It does not imply that either square-root coefficient vanishes.  Moreover,
because `f_b` has a nonzero finite value at its dominant square-root
singularity, taking a formal square root is not expected to halve the
exponential coefficient height; the singularity, rather than the value, is
halved only in a nonlinear local expansion.  This route remains open only if
the polynomial degree constraint on `B_p`, together with reciprocity, creates
an additional first-half relation not already equivalent to reflection.

Seven new ChatGPT audits Q448--Q454 were dispatched to test the primorial-gcd,
truncation-square, Franel substitution, batch determinant, companion
`p^2`, horizontal large-sieve, and adversarial-threshold routes.  The older
Q409--Q416 late watchers remain active, so resubmission did not discard any
late Notion delivery.

## 27. Low-jet obstruction to the truncation-square and Franel routes

The first-half range used after folding is exactly the range in which the
square factorization has no access to its polynomial endpoint.

Let `R` be any ring in which 2 is invertible.  The map

`1+tR[t]/(t^(m+1)) -> 1+tR[t]/(t^(m+1)),  V |-> V^2`

is a triangular bijection: if `U=1+sum u_k t^k`, its unique square root with
constant term 1 is determined by

`2v_k=u_k-sum_(1<=i<k)v_i v_(k-i)`.

Therefore, for `m<(p-1)/2`, the low coefficients of `B_p` are merely the
unique formal square root of the low coefficients of `A_p`.  In the second
congruence class the same statement applies to

`A_p/(t^2-34t+1)`,

whose denominator has constant term 1.  A zero `b_j=0` is one convolution
equation among the root coefficients; it does not force an additional
root-coefficient zero or divisibility condition.  The termination and
reciprocity of `B_p` occur at index about `p/2`, outside this low-jet
calculation.

The modular substitution

`f_b(t)=(1+x)h(x)^2`, `t=x(1-8x)/(1+x)`

has the same obstruction.  The change of variable has linear coefficient 1
and is an automorphism of every truncated formal-jet ring; multiplication by
`1+x` is by a unit; taking the square root is the bijection above.  Hence the
Franel representation is also an invertible nonlinear change of coordinates
on the relevant low jet.  By itself it cannot turn `b_j=0` into a second
independent congruence.

Reopen either route only if the degree-`~p/2` endpoint of `B_p` can be coupled
to the moving relation `p=n-j` with a certificate of subexponential height.
Any argument using only coefficients through degree `j<p/2` is formally
noncompressive.

## 28. Full lower channel is equivalent to a primorial gcd

The fixed-term reformulation is not limited to q=1.  Define

`P_b(n)=sum_(p<=n, p|b_n) log p
       =log gcd(rad(b_n),prod_(p<=n)p)`.

For `p>sqrt(n)`, write `n=qp+r`.  Lucas gives

`p|b_n <=> p|b_q b_r`.

The primes with `p|b_r` contribute exactly the remaining lower-channel sum
`L(n)`.  The primes with `p|b_q` are the already controlled companion channel,
whose total logarithmic weight is `O(n^(2/3))`.  The primes `p<=sqrt(n)` have
total logarithmic weight at most `theta(sqrt(n))=O(sqrt(n))`.  Hence

`P_b(n)=L(n)+O(n^(2/3))`.

In particular, after the existing support-law and companion-channel work,
the pointwise core of Problem 3.2 is equivalent to the single fixed-integer
statement

`log gcd(rad(b_n),n#)=o(n)`.

This removes the moving residue from the statement, but not from known proof
technology.  The fixed carrier `b_n` has optimal index but exponential height
`exp(Theta(n))`, so its size gives only `O(n)`.  Lucas amplification preserves
a zero only at indices congruent to `n (mod p)`; synchronizing this for all
bad primes requires a shift divisible by their product (or by a comparable
primorial), recreating exponential CRT height.  This explains why the
fixed-carrier reformulation is exact but does not by itself close the last
linear factor.

## 29. Exact modular scan through n=200,000

`projects/Q-series-and-Chan-s-work/q32_scan.cpp` computes every single-digit
Apéry zero by the audited recurrence modulo each prime, then inverts the
incidence relation

`n=qp+j`, `q<p`, `p|b_j`.

Thus it evaluates the lower-channel mass without factoring the enormous
integer `b_n`.  The scan through `n=200,000` found:

- 18,126 zero pairs `(p,j)` in total;
- at most 6 simultaneous lower-channel primes at any one `n`;
- at most 3 simultaneous q=1 primes at any one `n`;
- worst dyadic lower mass ratios
  `0.12598` up to 100,
  `0.02999` up to 1,000,
  `0.005047` up to 10,000,
  `0.0009247` on `(50,000,100,000]`, and
  `0.0005225` on `(100,000,200,000]`;
- the last worst point is `n=100,766`, with six lower-channel hits.

For q=1 the corresponding last two dyadic maxima were `0.0005086` and
`0.0003204`, and the maximum hit count remained three.

The behavior matches a sparse Poisson-type model: q=1 has expected mass much
smaller than the whole quotient range, while all `sqrt(n)<p<=n` have a
constant expected hit count.  It gives no proof, but it rules out an early
linear pile-up and shows that the needed theorem is vastly weaker than the
observed behavior.  A putative absolute bound of 6 or 3 is not inferred:
random maxima also grow extremely slowly.

## 30. Q452: exact horizontal Fourier remainder

Q452 confirms that the 2025 truncation-square/Galois theorem has no direct
horizontal coefficient content.  It also gives an exact Fourier formulation
of the q=1 problem.  Put

`Z_p={r in F_p:b_r=0}`,

`U_p(h)=sum_(r in Z_p) exp(-2 pi i h r/p)`.

Then Fourier inversion gives

`K_1(n)
 =sum_(n/2<p<=n)|Z_p|/p
  +sum_(n/2<p<=n) p^(-1) sum_(h!=0)U_p(h)e_p(hn)`.

The proved `|Z_p|<<p^(2/3)` makes the first term
`O(n^(2/3)/log n)`.  Thus the pointwise obstruction is now entirely the
cross-prime oscillatory remainder.  Any uniform estimate

`sup_n |sum_(p~n) p^(-1)sum_(h!=0)U_p(h)e_p(hn)|
 << n/log^(1+delta)n`

would close q=1.

Parseval alone cannot supply this.  The frequencies `h/p` with `p~N` have
Farey spacing `~N^(-2)`, so the classical large sieve over an interval of
only `N` values incurs an `N^2` term.  Even at the conjectural vertical scale,
the resulting pointwise bound is only `N/sqrt(log N)`, worse than the needed
`o(N/log N)`.

Q452 also gives two proof-grade obstructions:

1. The Kummer group and square-factorization conclusions allow artificial
   polynomials with almost every coefficient zero, so those conclusions alone
   cannot constrain coefficient locations.
2. Fixed rational diagonality alone is insufficient: the central binomial
   diagonal has a positive linear-weight interval of primes detected by
   Kummer carries.

The missing theorem must therefore use the Apéry-specific Fourier phases
`U_p(h)` across different residue characteristics.  Q470 asks whether the
exact zero indicator `1-b_r^(p-1)` can turn these phases into a
bounded-complexity trace family or whether its growing power/conductor
rigorously blocks that route.

## 31. The g=1 fixed-sum residual occurs in actual triples

Separating the two q=1 folded branches in the exact scan shows that neither
branch has a uniqueness property.  At `n=11,576`, the direct branch contains

`(p,j)=(8,893,2,683),(9,319,2,257),(11,437,139)`.

For all three,

`g=gcd(p-1,j)=1`,

so with `N=n-1=11,575` they are three distinct solutions of

`a+d=N`, `p=d+1`, `p|b_a`.

Likewise `n=47,066` has three reflected hits, all with `g=1`.  Therefore the
kernel identity, divisor condition, coprimality, and fixed-linear relation do
not force at most one or two candidates.  Any proof of a uniform saving in
this residual must use the Apéry vanishing itself.  Candidate-counting
thresholds cannot remove the term.

## 32. Q448 and Q451: the exact limit of soft q=1 information

Q448 independently confirms that the direct q=1 sum is exactly

`log gcd(b_n,prod_(n/2<p<=n)p)`

and that the strongest presently justified unconditional bound is still

`(1/2+o(1))n`.

It also gives two useful failure certificates.

First, if `m=floor((n-1)/2)`, `d_k=(Delta^k b)_0`, and

`C_n=sum_(k=0)^m binom(n,k)d_k`,

then every direct q=1 bad prime divides `C_n`.  This is a genuine single
integer independent of the unknown bad-prime set.  However the fixed Laurent
polynomial constant-term representation gives

`d_k=CT(Lambda-1)^k >= 4^k`,

because `Lambda-1` has nonnegative coefficients and constant coefficient 4.
Consequently

`C_n >= binom(n,m)4^m`

and

`log C_n >= (log 4)n+O(log n)`.

Thus the most natural exact polynomial interpolation compression is
provably exponential, not a subexponential certificate.

Second, a universal linear finite-window expression derived from the
second-order Apéry recurrence reduces to

`alpha(n)u_n+beta(n)u_(n+1)`.

If it vanishes for every solution with `u_n=0`, the freely specifiable second
state coordinate forces `beta(n)=0`.  Hence such a certificate is only a
multiple of the original value.  A new certificate must be nonlinear,
global, or use genuinely Apéry-specific arithmetic across characteristics.

Q451 performs the corresponding adversarial audit after reflection.  With
only Lucas, reflection, no consecutive zeros, character order, and the
forced continuant factor, an artificial first-block model can make every
top-half prime hit the same `n` while respecting all those restrictions.
Thus these inputs alone cannot improve the leading constant `1/2`.

The two unconditional deletions extracted in Q451 are:

- folded indices `j<n/(log n)^2`, of total weight `O(n/log n)`;
- selected multiplicative characters of order
  `d<=sqrt(n)/(log n)^2`, of total weight `O(n/log^3 n)`, using
  `sum_(d<=D) phi(d)=O(D^2)`.

The remaining obstruction contains primitive-order examples, including the
actual `g=1` triples in Section 31.  The horizontal Fejer proposal in Q451 is
the same exact cross-prime Fourier remainder isolated independently in Q452.
Its outer additive character currently has conductor growing with `p`.

## 33. Folded Newton interpolation: exact divisibility but linear height

Q451 proposes a smaller interpolation degree after folding.  Put

`J=floor((n-1)/3)`,

`c_k=(Delta^k b)_0`,

`F_J(X)=sum_(k=0)^J c_k binom(X,k)`.

For every `0<=j<=J`, `F_J(j)=b_j`.  If a direct folded prime is
`p=n-j`, then

`F_J(n)=F_J(j) = b_j (mod p)`.

If a reflected folded prime is `p=(n+1+j)/2`, then

`F_J(-n-1)=F_J(j) = b_j (mod p)`.

Therefore the whole q=1 radical divides

`gcd(b_n,F_J(n)F_J(-n-1))`

(up to the fixed small primes).  This replaces the moving congruences by an
exact holonomic gcd problem, but it does not win by height.

The reproducible exact-integer calculation in
`projects/Q-series-and-Chan-s-work/q32_newton.py` gives:

| `n` | `log|F_J(n)|/n` | `log|F_J(-n-1)|/n` | sum |
|---:|---:|---:|---:|
| 30 | 1.43485 | 1.52392 | 2.95877 |
| 60 | 1.59477 | 1.69553 | 3.29030 |
| 120 | 1.68659 | 1.79351 | 3.48010 |
| 240 | 1.73836 | 1.84842 | 3.58678 |
| 480 | 1.76714 | 1.87880 | 3.64595 |
| 720 | 1.77754 | 1.88973 | 3.66728 |
| 1000 | 1.78644 | 1.89940 | 3.68584 |

Both evaluations have positive linear logarithmic height; their product is
worse.  Q475 asks whether the three residue-class subsequences admit a
nonzero Ore resultant or low-height Bezout combination with `b_n`, and
whether adjacent evaluations, ratios, or Wronskians cancel the exponential
parts.  Reopen this route only through such a genuine gcd certificate:
interpolation size alone has failed.

Q476 separately audits the toric formulation at the top-half scale.  Since
`p>n/2` and the Laurent polynomial has fixed Newton polytope, only finitely
many multiples of `p-1` occur in the support of `Lambda^(n-1)`.  The question
is whether this finite alias sum can be compressed after substituting
`p=n-j`, rather than merely restating `b_j (mod p)`.

## 34. The Newton gcd is numerically small, but not support-exact

The relevant quantity is the gcd, not the individual interpolation heights.
The exact scan in
`projects/Q-series-and-Chan-s-work/q32_newton_gcd.py` computed

`G_n=gcd(b_n,F_J(n)F_J(-n-1))`

through `n=1200`.  It was nontrivial for 1,103 of the 1,195 tested indices,
mostly because of small fixed factors.  Nevertheless the worst normalized
logarithms on successive dyadic intervals fell as follows:

| interval | max `log(G_n)/n` |
|---:|---:|
| `(10,20]` | 0.46555 |
| `(20,40]` | 0.20174 |
| `(40,80]` | 0.22411 |
| `(80,160]` | 0.09507 |
| `(160,320]` | 0.09565 |
| `(320,640]` | 0.05424 |
| `(640,1200]` | 0.04347 |

This is much stronger numerically than the already sparse q=1 support and
makes the holonomic-gcd formulation worth auditing.  It is not yet a theorem,
and `G_n` is not exactly the desired radical.  For example, at `n=717` it
contains the factors

`5^3*17*59*443*751*821`,

including `751,821>n`.  Thus a proof must either control the full gcd or
isolate only its primes in `(n/2,n]`.  Merely observing that the two
holonomic sequences have different exponential growth constants does not
give a gcd theorem for P-recursive sequences.

Adjacent interpolation cutoffs do not create independence.  For
`K>=j`, a direct candidate `p=n-j` satisfies

`F_K(n)=b_j (mod p)`

as long as `K<p`.  The increment

`F_(K+1)(n)-F_K(n)=c_(K+1)binom(n,K+1)`

is automatically divisible by `p` once `K+1>j`, because the binomial
coefficient has the base-`p` carry.  The reflected evaluation has the same
failure through

`binom(-n-1,k)=(-1)^k binom(n+k,k)`.

Hence taking many neighboring cutoffs only inserts a structural top-prime
factor; it is not a second Apéry congruence.

## 35. Exact finite toric alias formula

The Laurent polynomial used in Q448 is

`Lambda=(1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)`.

Its exponent support lies in `[-1,1]^3` and has 22 monomials.  If

`C_N(e)=[x^(e_1)y^(e_2)z^(e_3)] Lambda^N`,

direct expansion gives the exact one-dimensional sum

`C_N(e)=sum_k binom(N,k)binom(N,k-e_1)
                 binom(2N-k,N-e_2)binom(2N-k,N-e_3)`,

with out-of-range binomial coefficients interpreted as zero.

Let `N=n-1`, `p=n-j`, and `d=p-1=N-j`.  For `0<j<p-1`, torus
orthogonality and `N=j (mod p-1)` give

`sum_(m in Z^3) C_N(md) = b_j (mod p)`.

Q476 sharpens the Newton-polytope step.  For every asymptotically relevant
bad prime, the admissible alias directions are exactly the 22 lattice points

`m in supp(Lambda)`,

not merely a subset of the 27 points in `{-1,0,1}^3`.  The only endpoint
where a `+-2` coordinate could enter has `j=p-1`, but
`b_(p-1)=1 (mod p)` and is never bad.

More decisively, Q476 proves a termwise Frobenius collapse.  Write

`Lambda=sum_(v in S) lambda_v x^v`.

Since `N=(p-1)+j`,

`Lambda^N=Lambda^p Lambda^(j-1) (mod p)`.

Coefficient extraction at `(p-1)v=pv-v` leaves only the Frobenius term
with exponent `pv`, and hence

`C_N((p-1)v)=lambda_v C_(j-1)(-v) (mod p)`.

Summing over the 22 directions is exactly the ordinary convolution identity

`sum_(v in S) lambda_v C_(j-1)(-v)=CT(Lambda^j)=b_j`.

Thus finite aliasing is not just equivalent after summation; each high alias
is the Frobenius lift of one summand in the original constant-term formula.
There is only one scalar vanishing equation.

This finite collapse does not by itself reduce height: the `m=(0,0,0)` term
is `C_N(0)=b_N`, while the full nonnegative alias sum is at most
`Lambda(1,1,1)^N=40^N`.  Its logarithmic height is therefore `Theta(n)`.
The concrete zero `p=17,j=3` has nonzero individual aliases modulo 17, so
the bad event is cancellation and does not force a common alias factor.

The formula also yields a general obstruction
to canceling that central term by changing constant-term representations.
If two fixed Laurent polynomials have the same moments

`CT(Lambda_1^k)=CT(Lambda_2^k)=b_k`,

then their two top-half alias sums are congruent modulo every candidate
prime, regardless of whether `b_j=0`.  Their difference is therefore
divisible by the entire top-half primorial (apart from endpoint exceptions).
Any nonzero such universal cancellation already has logarithmic height at
least `(1/2+o(1))n`, so it cannot be the desired subexponential
bad-prime-specific certificate.

The finite alias route can be reopened only if the additional combination
vanishes specifically when `b_j=0`, rather than producing another universal
Fermat/torus congruence.

## 36. Exact rational interpolation is worse than Newton interpolation

A possible escape from the uniqueness of polynomial interpolation is to find
integer-valued polynomials `U,V` satisfying

`U(j)=b_j V(j)` for every `0<=j<=J`.

For a direct candidate `p=n-j`, reduction at `n=j (mod p)` shows that a bad
prime divides `U(n)` unless it divides `V(n)`; in either case it divides
`U(n)V(n)`.  Thus a rational interpolant with

`log|U(n)V(n)|=o(n)`

would be a valid certificate.

Exact linear algebra was performed in the integer-valued bases
`binom(X,k)`, with `deg U+deg V=J`.  The balanced choice is disastrous:

| `n` | `J` | degrees `(U,V)` | `log|U(n)V(n)|/n` |
|---:|---:|---:|---:|
| 30 | 9 | `(4,5)` | 6.90377 |
| 60 | 19 | `(9,10)` | 15.01463 |
| 90 | 29 | `(14,15)` | 23.62629 |
| 120 | 39 | `(19,20)` | 32.12368 |
| 150 | 49 | `(24,25)` | 40.79897 |

The logarithmic height is quadratic in `n`.  The interpolation equations
contain the exponentially growing values `b_j`; their maximal minors multiply
those heights instead of sharing them.

An exhaustive degree split for `n=30,60,90` found that the smallest product
height always occurs at `deg V=0`, namely the original Newton polynomial.
Thus exact rational interpolation does not beat Section 33; balancing the
degrees makes it vastly worse.  Reopen only if an Apéry-specific
low-displacement-rank determinant identity replaces generic rational
interpolation, not by another choice of degrees or polynomial basis.

The classical irrationality-measure route has a separate exact ceiling.  If
`g_n=gcd(d_n a_n,d_n b_n)`, then the reduced Apéry denominator has exponential
rate

`3+log(lambda)-log(g_n)/n`,

while the rational-approximation error has rate `2log(lambda)`, where
`lambda=17+12sqrt(2)`.  An irrationality-exponent upper bound `mu` yields at
best

`log(g_n)/n <= 3+log(lambda)-2log(lambda)/mu+o(1)`.

Even the conjecturally optimal value `mu=2` gives only the already known
linear ceiling `3+o(1)`.  Therefore sharpening the real irrationality measure
cannot prove `log g_n=o(n)`; the missing input must use the restricted prime
support or a new arithmetic certificate.

## 37. A simpler fixed carrier: the central binomial coefficient

For every prime `n/2<p<=n`, except the odd boundary `n=2p-1`,
Kummer's theorem gives

`p | binom(n,floor(n/2))`.

Indeed, writing `n=p+j` with `0<=j<p-1`, the lower argument lies strictly
between `j` and `p`, so its base-`p` addition has one carry.  At the excluded
boundary `j=p-1`,

`binom(2p-1,p-1)=1 (mod p)`.

However Lucas gives `b_(2p-1)=5 (mod p)` for every `p>=7`, so the missing
carrier prime is never bad.  The only exception is the finite case `p=5`.
Therefore, up to that harmless finite exception, the q=1 radical is exactly
the top-half prime support of

`gcd(b_n,binom(n,floor(n/2)))`.

This carrier is simpler and substantially smaller than the Newton values:

`log binom(n,floor(n/2))=(log 2)n+O(log n)`.

The exact scan in
`projects/Q-series-and-Chan-s-work/q32_central_gcd.py` through `n=3000`
again shows a strongly decaying normalized gcd:

| interval | max `log gcd/n` |
|---:|---:|
| `(10,20]` | 0.25046 |
| `(20,40]` | 0.16866 |
| `(40,80]` | 0.18663 |
| `(80,160]` | 0.07310 |
| `(160,320]` | 0.07288 |
| `(320,640]` | 0.04923 |
| `(640,1280]` | 0.03267 |
| `(1280,2560]` | 0.01442 |
| `(2560,3000]` | 0.01076 |

This is evidence, not a proof.  For the top-half primes the binomial
divisibility is universal, so the carrier supplies no second local
congruence: intersecting with `b_n` is exactly the original q=1 condition.
Its value is that any Apéry-specific gcd identity can now target an elementary
factorial ratio rather than the much larger interpolation sequence.

The logarithmic height cannot be made sublinear by replacing the carrier with
another integer divisible by every top-half prime: the top-half primorial
itself already has logarithm `(1/2+o(1))n`.  A successful use of this carrier
must therefore prove cancellation in the gcd, not merely reduce the carrier's
Archimedean size.

## 38. Exact low-frequency reduction, and why it is still a new theorem

Q470 audited the Fourier zero-indicator route.  If

`Z_p={r mod p:b_r=0 mod p}` and
`U_p(h)=sum_(r in Z_p) exp(-2 pi i h r/p)`,

then the identities

`U_p(h)=-sum_(r mod p) b_r^(p-1) exp(-2 pi i h r/p)`

and

`U_p(h)=p^(-1) sum_(u!=0) sum_(r mod p)
                exp(2 pi i (u b_r-h r)/p)`

are exact for `h!=0`.  They do not produce a bounded-conductor trace
function.  The recurrence realizes `b_r^(p-1)` in the `(p-1)`-st symmetric
power of a two-dimensional transfer system, of rank `p`; the rational
diagonal realization needs `4(p-1)` variables; and the unique polynomial
indicator of `Z_p` has degree at least

`p-|Z_p| >= p-3p^(2/3)-O(1)`.

The 2026 square factorization of the truncated generating polynomial is also
in the wrong algebra: it controls ordinary coefficient convolution, whereas
the zero indicator uses the Hadamard power

`sum_r b_r^(p-1)t^r`.

Thus it does not reduce the growing complexity.

There is nevertheless a clean weaker residual theorem.  Put

`I_n={p prime:n/2<p<=n}`,
`H=floor(n^(2/3)(log n)^2)`,
`V_h(n)=sum_(p in I_n) U_p(h) exp(2 pi i h n/p)`,
and `w_h=1-h/H`.

The positive Fejer kernel gives

`K_1(n) <= H^(-1) sum_(p in I_n)|Z_p|
          +(2/H) Re sum_(1<=h<H) w_h V_h(n)`.

The proved bound `|Z_p|<<p^(2/3)` makes the first term
`O(n/log^3 n)`.  By Cauchy--Schwarz and same-prime Parseval, it is enough to
prove the cross-prime estimate

`|2 Re sum_(h<H) w_h sum_(p<q in I_n)
 U_p(h) conjugate(U_q(h))
 exp(2 pi i h n(1/p-1/q))|
 << n^(8/3)/(log n)^eta`

for any fixed `eta>0`, uniformly in `n`.  The diagonal contribution is
already `O(n^(8/3)/log n)`, so this estimate implies
`K_1(n)=o(n/log n)` and closes q=1.

The derivation and scales check out, but the displayed bound is not supplied
by a standard large sieve.  It is a centered pair-correlation statement for
the moving points `(n-r)/p`, `r in Z_p`, and an adversarial singleton aligned
at one `n` reaches the unsaved `n^(8/3)` scale.  Existing information is
vertical in one characteristic and gives no cancellation between `Z_p` and
`Z_q`.  This route can be reopened only with a genuinely horizontal
cross-prime input; rewriting the zero indicator by Fermat, a diagonal, or the
truncation-square theorem does not provide it.

## 39. Newton/Ore route closed for fixed-order, fixed-degree certificates

Q475 gives a structural audit of Sections 33--34.  On the direct branch,

`p=n-j  =>  F_J(n)=b_j (mod p)`,

and on the reflected branch,

`2p=n+1+j  =>  F_J(-n-1)=b_j (mod p)`.

These hold for every candidate prime in the corresponding geometric branch,
whether or not it is bad.  Thus the two Newton values are alternate
coordinates for the same Lucas/reflection residue, not independent
congruences.

The exact asymptotic rates can be derived from the constant-term Newton
coefficients.  With

`lambda=17+12sqrt(2)` and `gamma=lambda-1`,

the two evaluation heights have limits

`alpha_U=(1/3)log(27 gamma/4)`,
`alpha_V=(1/3)log(256 gamma/27)`.

They exceed the Chebyshev mass of their respective prime branches
(`n/3+o(n)` and `n/6+o(n)`), and their sum exceeds `log(lambda)`.
This explains the numerical positive linear heights rather than leaving them
as a finite-range observation.

There is also an exact algebraic obstruction.  If `a,u,v` denote the Apéry,
direct-Newton, and reflected-Newton coordinates, the branch-union bad locus is

`V(a,u) union V(a,v)`,

whose ideal is

`(a,u) intersect (a,v)=(a,uv)`.

Consequently every universal polynomial certificate obtained only from these
two folded congruences is a multiple of `a` plus a multiple of `uv`; no
linear combination of `u,v` covers both branches.  The product certificate
is algebraically minimal.

Finally, along each residue class modulo 3, the exponential bases of the
Apéry and two Newton slices are multiplicatively independent.  A nonzero
fixed-degree polynomial in the three values with subexponential-height
coefficients has a unique dominant monomial and hence exponential height.
Fixed shifts, ratios, Wronskians, and finite differences do not change this:
either they retain an exponential base, or they lose the pointwise prime
divisibility.

An Ore gcd/resultant only detects common shift-invariant solution modules.
The event here is a common zero at one coordinate; applying an operator
Bezout identity introduces adjacent values not divisible by the same prime.
Moreover, the distinct exponential bases rule out a finite-order rational
intertwiner between the distinguished Apéry and Newton slices.

Therefore the Newton route is closed within natural fixed-order Ore
operations and fixed-degree polynomial combinations.  It can be reopened
only with order/degree growing with `n` together with a new, proved
high-order cancellation mechanism.

## 40. The central-binomial carrier is a whole zero block modulo `p`

The parity-split central binomial sequence

`B_n=binom(n,floor(n/2))`

satisfies

`(m+1)B_(2m+1)=(2m+1)B_(2m)`,
`B_(2m+2)=2B_(2m+1)`.

For an odd prime `p` and every `0<=j<=p-2`, Kummer gives the sharper statement

`v_p(B_(p+j))=1`.

Indeed both denominator factorials in

`B_(p+j)=(p+j)!/
 (floor((p+j)/2)! ceil((p+j)/2)!)`

have arguments below `p`, while the numerator contains exactly one multiple
of `p`.  Moreover,

`B_(p+j)/p =
 ((p+j)!/p)/
 (floor((p+j)/2)! ceil((p+j)/2)!)`

is a `p`-adic unit, with

`B_(p+j)/p =
 -j!/
 (floor((p+j)/2)! ceil((p+j)/2)!) (mod p)`

by Wilson's theorem.

At the final endpoint,

`B_(2p-1)=1 (mod p)`.

Thus modulo `p` the carrier is identically zero on the block

`B_p,B_(p+1),...,B_(2p-2)`,

but not at its final endpoint.

Away from that harmless endpoint, its zero at `n=p+j` is not a second local equation paired with
`b_(p+j)=0`; it is a universal singular factor injected at the beginning of
the block.  In the same block Lucas gives

`b_(p+j)=b_1 b_j=5b_j (mod p)`.

Consequently a pointwise Ore/resultant elimination between the two
recurrences sees the zero central-binomial solution throughout the block and
can recover only the original condition `b_j=0`.  Passing to the normalized
unit `B_(p+j)/p` removes the zero but depends on the unknown candidate `p`;
it supplies no globally defined common-zero certificate.

This does not refute every hypergeometric Bezout identity involving `B_n`,
but it closes the naive common-zero/adjacent-shift argument.  A successful
identity must use the normalized first `p`-adic layer in a way that eliminates
the unknown `p`, or it will merely restate Lucas folding.

## 41. The first Euclidean remainder of the central carrier stays exponential

Let `B_n=binom(n,floor(n/2))` and let `R_(n,d)` be the least absolute
residue of

`n^d b_n (mod B_n)`.

Every divisor of `gcd(b_n,B_n)` divides `R_(n,d)`, so a fixed `d` with
`log|R_(n,d)|=o(n)` would be a valid elementary certificate.  Exact
computation in
`projects/Q-series-and-Chan-s-work/q32_central_remainder.py` gives, for
`d=0`, the rates

| `n` | `log|R_(n,0)|/n` | `log B_n/n` |
|---:|---:|---:|
| 60 | 0.641986 | 0.655195 |
| 120 | 0.661862 | 0.671300 |
| 240 | 0.673770 | 0.680784 |
| 480 | 0.684527 | 0.686245 |
| 720 | 0.686006 | 0.688264 |
| 1000 | 0.688035 | 0.689467 |

The signed remainder approaches the full carrier rate `log 2`, rather than
sublinear height.  Multipliers `n^d` for every fixed `0<=d<=8` behave the
same way; at `n=1000` all measured rates lie between `0.6871` and `0.6886`.

This is an exact finite computation, not an asymptotic impossibility theorem.
It rules out the literal first Euclidean step and provides no evidence for a
low-degree polynomial multiplier.  A useful Bezout construction must exploit
an exact hypergeometric cancellation, not ordinary reduction of `b_n` modulo
the central binomial coefficient.

## 42. Two fixed hypergeometric truncations capture the full q=1 radical

Write the Apéry summand as

`T(n,k)=binom(n,k)^2 binom(n+k,k)^2`

and put

`J=floor((n-1)/3)`,
`L_n=sum_(0<=k<=J) T(n,k)`,
`H_n=sum_(ceil(n/2)<=k<=n) T(n,k)`.

For every prime `n/2<p<=n`, let `r=n-p` and

`j=min(r,p-1-r)`.

Then `0<=j<=J`, and a direct Lucas calculation gives

`L_n=b_j (mod p)`,
`H_n=4b_j (mod p)`.

For the direct branch `r=j`, the only surviving terms in the prefix are
`0<=k<=j`, while the surviving upper block is `k=p,...,p+j`.
For the reflected branch `r=p-1-j`, the terms `0<=k<=j` satisfy the
termwise reflection

`T(p-1-j,k)=T(j,k) (mod p)`,

and again the surviving upper block is `k=p,...,p+j`, with the factor 4
coming from the leading base-`p` digit of `binom(n+k,k)`.

It follows that the q=1 radical is exactly the top-half prime support of

`gcd(L_n,H_n)`.

The reproducible exact audit is
`projects/Q-series-and-Chan-s-work/q32_truncation_gcd.py`.  Through `n=1200`,
the dyadic maxima of `log gcd(L_n,H_n)/n` are

| interval | maximum |
|---:|---:|
| `(10,20]` | 0.46555 |
| `(20,40]` | 0.33540 |
| `(40,80]` | 0.16663 |
| `(80,160]` | 0.11924 |
| `(160,320]` | 0.08280 |
| `(320,640]` | 0.05424 |
| `(640,1200]` | 0.03158 |

At `n=321`,

`gcd(L_n,H_n)=5*179*193*211`,

and the three large factors are exactly the three q=1 hits.  This is strong
numerical evidence for a small gcd, not a theorem.

The full gcd also has unrelated factors.  At `n=717` it is

`5^2*443*751*821`;

only `443` lies in the q=1 interval, while `751,821>n`.  Therefore a theorem
about the full gcd would be stronger than necessary, just as for the Newton
gcd.

The two congruences are still locally dependent:

`H_n-4L_n=0 (mod p)`

for every top-half candidate prime, bad or good.  Hence the full top-half
primorial divides `H_n-4L_n` (Section 46 strengthens this to its square).
The gcd reformulation is therefore

`gcd(L_n,H_n)=gcd(L_n,H_n-4L_n)`,

namely one selective coordinate paired with another universal carrier.  Its
individual heights are also linear:

`log L_n=(log 16+o(1))n`,

while `H_n` contains the dominant portion of the Apéry sum and has rate
`log(17+12sqrt(2))`.

Thus height alone is worse than the central-binomial carrier.  The route is
worth reopening only if creative telescoping or an arithmetic gcd theorem
uses the common hypergeometric origin of `L_n,H_n`; the mod-`p` equations
themselves provide no second condition.

## 43. Q491: exact parity-Ore audit of the central carrier

Q491 independently derived the order-two recurrences for the even and odd
Apéry sections and the order-one hypergeometric recurrences for

`B_(2m)=binom(2m,m)`,
`B_(2m+1)=binom(2m+1,m)`.

The Ore right gcd in `Q(m)[S]` is a unit in each parity.  Substitution of the
binomial hypergeometric ratio into the Apéry operator gives a nonzero
compatibility polynomial of degree 11, with leading coefficient `-80059392`.
This proves that the two recurrences have no common nonzero solution module.
It does not control two different distinguished solutions that merely vanish
at one coordinate.

Every cross-Wronskian factors through the rank-one binomial sequence.  More
generally, reduce a fixed window to formal state coordinates

`(x_0,x_1)` for Apéry and `y_0` for the binomial sequence.

The common-zero locus has ideal

`(x_0,y_0)`.

Thus every fixed-window polynomial certificate is

`x_0 R(x_0,x_1,y_0)+y_0 S(x_0,x_1,y_0)`.

The first part retains the original Apéry coordinate; the second retains the
universal binomial factor.  Removing the latter loses the forced divisibility.
The parity sections have multiplicatively independent exponential bases
`lambda^2` and `4`, so after exact cancellations every nonzero fixed-degree
certificate with subexponential coefficients still has positive exponential
height.

This closes fixed-order Ore subresultants, cross-Wronskians, bounded windows,
and fixed-degree polynomial combinations of the two recurrence modules.  It
also supplied the odd-boundary correction incorporated in Sections 37 and
40: `B_(2p-1)` is a `p`-unit, but `b_(2p-1)=5 (mod p)`, so for `p>=7` that
missing carrier factor is never a q=1 hit.

## 44. Q489 Newton--Pade determinants: exact certificate, bad height

Q489 isolates a degree-growing determinant that is not covered by the
fixed-order no-go results.  Put `M=floor(n/2)` and choose `a,b>=0` with
`a+b=M-1`.  Form the `(M+1)`-row matrix at

`x=0,1,...,M-1,n`

with columns

`binom(x,k), 0<=k<=a`,

followed by

`b_x binom(x,k), 0<=k<=b`.

Let its determinant be `D_(a,b)(n)`.  If a non-boundary top-half bad prime
`p` is given and `r=n-p`, then `r<M`,

`binom(n,k)=binom(r,k) (mod p)`

for every column index, while both `b_n` and `b_r` vanish modulo `p`.
The last row therefore equals row `r` modulo `p`, and

`p | D_(a,b)(n)`.

At `b=0` this is the ordinary degree-`M-1` interpolation-error certificate.
Balanced `a,b` are a genuine degree-growing rational-interpolation/Padé
extension, so they required a separate audit.

Exact Bareiss calculations are in
`projects/Q-series-and-Chan-s-work/q32_pade_determinant.py`.  For every
determinant the computation removes its *entire* `M`-smooth part, which is
more favorable than removing only a predicted systematic factor.  The
results are already strongly negative:

| `n` | degrees `(a,b)` | `log|D|/n` | `log(rough part)/n` |
|---:|---:|---:|---:|
| 20 | `(9,0)` | 3.2241 | 3.1201 |
| 20 | `(4,5)` | 8.6238 | 7.6742 |
| 30 | `(14,0)` | 3.3045 | 3.1353 |
| 30 | `(7,7)` | 11.4884 | 10.4686 |
| 40 | `(19,0)` | 3.3491 | 3.2568 |
| 40 | `(9,10)` | 15.0883 | 13.8355 |

Scanning every degree split for `n=20,30,40` finds the minimum always at
`b=0`.  The balanced rough logarithm is already quadratic in `n`; deleting
all small-prime factors barely changes that scale.  This agrees with the
rational-interpolation failure in Section 36, now for the precise determinant
and after the maximally favorable smooth-factor removal.

This finite computation does not prove a quadratic lower bound on the rough
part.  It makes the proposed balanced small-determinant theorem implausible
and identifies the necessary reopening condition: an exact
Cauchy--Binet/Selberg evaluation must exhibit a cancellation invisible in
these values.  Generic Padé balancing or another degree split should not be
retried.

Q489 also proves fixed-class obstructions consistent with earlier audits:
Gosper antidifferences, bounded contiguous relations, fixed rational gamma
quotients, and fixed-order Pell-scaled corrections cannot change the
incompatible exponential bases enough to leave a nonzero subexponential
integer.  Generic rational approximation with smooth denominators is
circular because its least nonzero cleared remainder is controlled by the
unknown gcd itself.

## 45. Q487 factorial-ratio audit

Q487 confirms the corrected central-binomial support theorem and finds no
published pointwise gcd result for a factorial-ratio sequence against the
Apéry P-recursive sequence.

It checks the second Schmidt--Strehl representation

`b_n=sum_k D_(n,k) f_k`,

where

`D_(n,k)=binom(n,k)binom(n+k,k)`

and `f_k` is the Franel number.  For `n=p+r`, both the direct representation
and the Strehl representation split into the same two Lucas blocks with
multipliers 1 and 2.  Squaring in the direct sum, or multiplying the
Legendre and Franel factors in the Strehl sum, gives the identical total
factor `1+4=5`.  Their difference is the lower-index Strehl identity itself;
there is no second congruence or nonzero resultant.

The symbolic-summation failures are also structural:

- one-term Gosper telescoping would express `b_n` through endpoints of
  exponential base at most 16, contradicting the Apéry base
  `17+12sqrt(2)`;
- fixed edge truncations retain a dominant saddle and hence positive
  exponential height;
- the direct--Strehl edge difference cancels only the leading saddle and
  saves about `0.0091n` in the exponent, while becoming divisible by all
  branch candidates rather than only bad ones;
- fixed-shift bilinear Wronskians have product exponential base and are
  worse than either input.

The one exact symbolic test not covered by these arguments is a
parity-split, central-binomial-adapted two-sieved factorial-basis reduction
using modified Abramov--Petkovsek reduction.  Its success criterion is exact:
the minimal nonsummable residual must have dominant growth constant 1 and
must retain the top-half Apéry residue.  Every standard basis audited so far
has dominant constant greater than 1.  This remains a finite computer-algebra
experiment, not a proved route; the formulas needed to instantiate the
two-sieved basis were not explicit in the returned answer.

## 46. Q501: the truncation difference contains the squared primorial

Q501 independently verifies the full q=1 truncation congruences and proves
the sharper universal statement

`P_n^2 | H_n-4L_n`,

where

`P_n=prod_(n/2<p<=n) p`.

The terms outside the two Lucas-surviving blocks already contain a binomial
coefficient divisible by `p`; the square in `T(n,k)` makes them vanish modulo
`p^2`.  A first-order expansion of the two surviving blocks reduces the
remaining defect to a finite binomial-harmonic identity.  Thus every
top-half prime occurs to depth at least two in the difference, independently
of whether the folded Apéry value vanishes.  The divisibility was also
checked locally by exact arithmetic through `n=120`.

This strengthening still does not select bad primes.  After defining

`D_n=(H_n-4L_n)/P_n`,

one still has `P_n|D_n`.  A bad prime need not divide the further quotient
`(H_n-4L_n)/P_n^2`; that would require a third `p`-adic digit not implied by
`p|L_n`.  Therefore the square factor is a universal obstruction, not a
depth amplification for the desired radical.

The exact local algebra is the line `Y=4X`; certificates derived only from
this relation lie in

`(X,Y-4X)`.

After evaluation, the selective part remains a multiple of `L_n`, while the
other part retains the universal squared primorial.  Since the two
Archimedean bases are `16` and `17+12sqrt(2)`, multiplicative independence
also rules out subexponential height for any nonzero fixed-degree polynomial
combination with subexponential coefficients.

## 47. Full Padé determinantal divisor: tiny numerically, exact localization

The individual Padé determinants in Section 44 are all divisible by every
non-boundary q=1 bad prime.  Instead of choosing one degree split, form the
full `(M+1) x 2M` matrix, `M=floor(n/2)`, with rows

`x=0,1,...,M-1,n`

and columns

`binom(x,k)` and `b_x binom(x,k)`, `0<=k<M`.

A bad top-half prime makes the last row equal row `r=n-p` modulo `p`, so it
divides every maximal minor.  Let `Delta_n` be the gcd of all those minors.

This exponentially large family has a simple exact reduction.  Apply the
unimodular inverse-Pascal transform to the first `M` rows and then use the
identity block of the unweighted binomial columns to eliminate the last row.
The determinantal divisor becomes the gcd of `M` interpolation residuals.
Equivalently, if

`ell_j(n)=(-1)^(M-1-j) binom(n,j) binom(n-j-1,M-1-j)`

is the cardinal Lagrange coefficient for the nodes `0,...,M-1`, then

`Delta_n=gcd_(0<=j<M) ell_j(n)(b_n-b_j)`.

The reason is that evaluation of integer-valued polynomials in the binomial
basis is unimodular on these `M` nodes, while

`q(n)=sum_(j<M) ell_j(n)q(j)`

for every polynomial of degree `<M`.

The exact fast audit is
`projects/Q-series-and-Chan-s-work/q32_pade_minor_gcd.py`.  It avoids
enumerating minors and gives:

| interval | max `log Delta_n/n` through 1200 |
|---:|---:|
| `(10,20]` | 0.38042 |
| `(20,40]` | 0.26361 |
| `(40,80]` | 0.19102 |
| `(80,160]` | 0.09240 |
| `(160,320]` | 0.08674 |
| `(320,640]` | 0.06072 |
| `(640,1200]` | 0.03012 |

At `n=321`,

`Delta_n=2^3*5*179*193*211`,

so its large factors are exactly the three q=1 hits.  At `n=743`,

`Delta_n=2^3*5^2*197*211*631`,

where `631` is the q=1 hit and the other odd factors lie below `n/2`.
This is the strongest small-gcd numerical certificate found in this round.

After dividing out the exact q=1 radical, the remaining nuisance factor is
also polynomial-sized throughout the scan.  Through `n=1200`, the largest
observed value of

`log(Delta_n/R_n)/log n`

is `3.414`, at `n=1093`.  This is evidence only; no uniform polynomial bound
for the small-prime determinantal content has been proved.

The exact formula also exposes why it is not yet a proof.  For a top-half
prime `p`, evaluation at `n` modulo `p` is the cardinal functional at the
single node `r=n-p`:

`ell_r(n)=1 (mod p)`,
`ell_j(n)=0 (mod p)` for `j!=r`.

Therefore

`p|Delta_n <=> p|(b_n-b_r) <=> p|4b_r`

by Lucas, which is exactly the original moving-zero condition.  All the
other residual coordinates vanish because their Lagrange coefficient already
contains the universal candidate prime.  The many determinants do not create
a second local equation; they package the CRT localization perfectly.

The determinant has no prime factors above `n`: if `q>n`, then both
`ell_0(n)` and `ell_1(n)` are `q`-units, so `q|Delta_n` would force

`b_n=b_0=b_1 (mod q)`,

and hence `q|4`, impossible.  The odd midpoint `n=2p-1` is also absent for
`p>=5`, as exact computation confirms and the same two-coordinate argument
shows.  Consequently the prime support of `Delta_n` above `M` is exactly the
q=1 support, with no unrelated primes larger than `n`.  This is cleaner than
the Newton and truncation gcds, which had extraneous factors above `n`.

Thus a theorem `log Delta_n=o(n)` would solve q=1, and the data strongly
suggest it, but proving it is equivalent at top-half primes to controlling
the same aligned moving divisors.  The determinantal description becomes
useful only if a global Smith-normal-form, adjoint-recurrence, or lattice-index
bound controls `Delta_n` without applying the local equivalence prime by
prime.

## 48. Q503: Strehl--Franel lowers a one-certificate exponent to `log 8`

Q503 identified the first genuine height improvement among the positive
Apéry representations.  Let

`F_k=sum_(i=0)^k binom(k,i)^3`

be the Franel number.  The Schmidt--Strehl identity is

`b_n=sum_(k=0)^n binom(n,k)binom(n+k,k)F_k`.

For `J=floor((n-1)/3)`, define

`S_n=sum_(k=0)^J binom(n,k)binom(n+k,k)F_k`.

If `p` is a top-half prime and

`j=min(n-p,2p-1-n)<=J`,

then `S_n=b_j (mod p)` on both folded branches.  On the direct branch this
is immediate from Lucas.  On the reflected branch the two binomial factors
exchange roles and each contributes a sign `(-1)^k`, so their product is
again the lower-index Strehl kernel.  Terms `j<k<=J<p` vanish.  Therefore
every q=1 bad prime divides `S_n`.

At the endpoint `k=n/3+O(1)`, the two-binomial kernel has logarithmic rate
`log 4`, while `F_k` contributes `(log 8)/3=log 2`.  Positivity and
monotonicity near the endpoint give

`log S_n=(log 8)n+O(log n)`.

This improves the earlier direct-prefix bound `log L_n=(log 16)n+O(log n)`,
but it is still linear and cannot prove `R_n=exp(o(n))`.

Q503 also audits a fixed Legendre/binomial translation family containing the
Strehl formula.  Among fixed nonnegative integral translations, the Strehl
center is the height minimum; composition only adds translation parameters.
Whipple merely reverses the original summation index.  Thus repeated positive
hypergeometric transforms do not continue the entropy descent.  Reopening
this route requires a signed WZ or constant-term cancellation that removes
the entire interior saddle *and* preserves the folded residue.

The Notion-to-Markdown extraction omitted Q503's displayed formulas, so only
the claims above that can be reconstructed from the explicit Schmidt--Strehl
identity are accepted.  The identity, both branch congruences, and the gcd
localization are independently checked in
`projects/Q-series-and-Chan-s-work/q32_strehl_gcd.py`.

There is a useful new two-certificate package.  The direct prefix `L_n` and
the Strehl prefix `S_n` both equal `b_j` modulo every top-half candidate
prime, hence

`p | gcd(L_n,S_n) <=> p | b_j`.

Also `S_n-L_n` is divisible by the entire top-half primorial, independently
of badness.  Exact data show the gcd is small (for example at `n=321` it is
`5*179*193*211`, precisely the nuisance factor 5 and the three q=1 hits).
This is strong evidence, not a proof: locally the pair again supplies only
one selective equation, because the difference is universally zero modulo
every candidate prime.  A global subexponential bound for this gcd would
solve q=1 and is a new target distinct from bounding either positive prefix
alone.

## 49. Polynomial content of the full Legendre--Euler family

The formulas omitted by the Markdown exporter were recovered directly from
the Notion equation blocks.  For an integral center `c`, Q503 defines

`Q_n(t)=sum_k L(n,k)t^k=P_n(1+2t)`,

`K_c(n,m)=[y^m]Q_n(c+y)`,

`g_m^(c)=sum_(i=0)^m binom(m,i)(-c)^(m-i)F_i`,

and the prefix

`T_n(c)=sum_(m=0)^J K_c(n,m)g_m^(c)`.

Treating `c` as an indeterminate gives a useful exact expansion.  The
coefficient of `c^0` is the Strehl prefix.  For `d>0`,

`[c^d]T_n(c) = sum_i (-1)^(J-i) L(n,i+d) binom(i+d,i)`

`                 * binom(d-1,J-i) F_i`,

where

`max(0,J+1-d)<=i<=min(J,n-d)`.

This follows by expanding both factors and applying

`sum_(h=0)^a (-1)^h binom(N,h)=(-1)^a binom(N-1,a)`.

More importantly, the actual polynomial congruence at a top-half prime is

`T_n(c) = b_j(1+2c^p) (mod p)`,

not `b_j(1+2c)`.  The latter holds only after evaluating `c` in `F_p`.
Indeed

`Q_(p+r)(t)=(1+2t^p)Q_r(t) (mod p)`,

and the cutoff is below `p`; reflection replaces `Q_r` by `Q_j`.  Therefore,
if `Gamma_n` is the gcd of all monomial coefficients of `T_n(c)`, then

`p|Gamma_n <=> p|b_j`

for every top-half candidate `p`.  The complete q=1 support is again packaged
exactly.

This also kills a tempting false independence argument.  Differentiating in
`c` gives zero modulo `p`, because the derivative of `c^p` vanishes.  All
derivative certificates are universal at candidate primes; they do not
provide a second selective equation.

The exact expansion and congruence are checked by
`projects/Q-series-and-Chan-s-work/q32_legendre_content.py`.  The content is
small numerically, but its available individual coefficients still have
positive linear height.  The constant coefficient has rate `log 8`; the
leading coefficient is

`(-1)^J binom(2n,n)binom(n-1,J)`,

with rate `log 4+H(1/3)`.  Thus merely taking another coefficient does not
give a subexponential certificate.  As with `Delta_n`, only a genuinely
global content/Smith bound could finish the argument.

## 50. Q515 completes and closes the two-sieved factorial-basis test

Q515 supplies the explicit parity bases

`Phi_(epsilon,k)(m)=binom(2m+epsilon,m-k)`,

so `Phi_(0,0)(m)=binom(2m,m)` and
`Phi_(1,0)(m)=binom(2m+1,m)`.  There are unique integral residual
coefficients `e_k,o_k` with

`b_(2m)=sum_(k=0)^m e_k Phi_(0,k)(m)`,

`b_(2m+1)=sum_(k=0)^m o_k Phi_(1,k)(m)`.

The inverse transforms are, for `k>=1`,

`e_k=sum_(m=0)^k (-1)^(k-m) [2k/(k+m)]`

`                   * binom(k+m,k-m)b_(2m)`,

with `e_0=b_0`, and

`o_k=sum_(m=0)^k (-1)^(k-m) [(2k+1)/(k+m+1)]`

`                   * binom(k+m+1,k-m)b_(2m+1)`.

The first values are

`e=(1,71,32711,21263474,16196884679,...)`,

`o=(5,1430,811805,578594525,463454152550,...)`.

These residuals remain exponentially large.  If
`lambda=17+12sqrt(2)`, their dominant root `mu` is characterized by

`mu+mu^(-1)=lambda^2-2`.

Over the rationals it is the dominant root of

`x^4-1150x^3-2301x^2-1150x+1`.

Thus `mu>1` (numerically about 1152), and the growth per original Apéry
index is `sqrt(mu)`, only slightly smaller than `lambda`.  The parity pattern
also leaves infinitely many residual coefficients odd, so no hidden
factorial/lcm content can remove this exponential growth.

For `n=2m+epsilon`, a top-half prime has `r=n-p` and `d=m-r`.  The exact
Lucas transform is the moving tail

`sum_(h=0)^r binom(r,h) r_(d+h)^(epsilon)`

`       = b_n = 5b_r (mod p)`,

where `r_k^(0)=e_k` and `r_k^(1)=o_k`.  Badness therefore does not annihilate
one residual coefficient or a fixed finite combination.  It annihilates a
tail whose start moves with `p`.  Fixed or sublinear carrier subtraction
retains the positive exponential saddle; reaching the saddle requires a
linear-depth subtraction that discards a linear portion of the prime range.

The returned Sage code has been independently translated to exact
pure-Python arithmetic in
`projects/Q-series-and-Chan-s-work/q32_two_sieved_basis.py`.  It verifies the
inverse transform, reconstruction, low values, and every top-half tail
congruence through `n=300`.

This closes fixed-order rational Abramov--Petkovsek reduction and
fixed/sublinear-depth subtraction.  Reopen only with an `n`-dependent signed
growing-order reduction whose moving tails retain the folded residue while
its integer boundary has subexponential height.

## 51. Q521 exact Smith reduction and unconditional global bounds

Q521 independently derives the full integral equivalence

`M_n ~ [I_M 0; 0 c_0 ... c_(M-1)]`,

where `c_j=ell_j(b_n-b_j)`.  Hence the Smith invariants are `M` copies of
one followed by `Delta_n`; there is no hidden recurrence-dependent invariant
factor.

Two coordinates already give a useful global Bezout bound.  Since

`ell_0=+-binom(n-1,M-1)`,

`ell_1/ell_0=n(M-1)/(n-1)`,

putting `L=lcm(|ell_0|,|ell_1|)` and subtracting the scaled residuals gives

`Delta_n | 4L`,

with

`L=binom(n-1,M-1)n(M-1)/gcd(M-1,n-1)`.

Therefore

`log Delta_n <= (log 2)n+O(log n)`.

Q521 also proves the exact valuation statement for every non-boundary prime
`M<p<=n`:

`v_p(Delta_n)=1` if `p|b_n`, and `0` otherwise.

Writing `Delta_n=R_n U_n`, the nuisance quotient `U_n` is `M`-smooth.  Removing
the known top-prime layer from the central binomial factor in `L` yields the
unconditional bound

`log U_n <= (log 2-1/2)n+o(n)`.

This is a real improvement for the smooth quotient, but remains linear.

Finally, for every integer-valued polynomial `g` of degree `<M`,

`T_g(n)=b_n g(n)-sum_(j<M) ell_j b_j g(j)`,

and the set of all such values is exactly `Delta_n Z`.  Abstract Smith or
adjoint theory therefore has no further automatic saving; it must control
the arithmetic content of many residual coordinates simultaneously.

## 52. New obstruction: Q521's proposed single signed adjoint is impossible

Q521 observed that if `max_(j<M)|g(j)|=exp(o(n))` and `g(n)!=0`, then the
`b_n g(n)` term dominates and

`log|T_g(n)|=n log(17+12sqrt(2))+o(n)`.

It therefore proposed searching for `g(n)=0` with a nonzero subexponential
residual.  The exact value lattice shows that this second option is also
impossible.

The evaluation map from integer-valued polynomials of degree `<M` to their
values `z_j=g(j)` on `0,...,M-1` is all of `Z^M`.  The constraint `g(n)=0`
is

`sum_j ell_j z_j=0`.

Because `sum_j ell_j=1`, this kernel is saturated.  On it,

`T_g(n)=-sum_j ell_j b_j z_j`.

The image ideal is the second determinantal divisor of the `2 x M` matrix
with rows `(ell_j)` and `(ell_j b_j)`.  Thus every such residual is a multiple
of

`Theta_n=gcd_(i<j) ell_i ell_j(b_j-b_i)`.

For every non-boundary top-half prime `p`, the Lagrange vector has one
`p`-unit coordinate `ell_(n-p)` and all other coordinates divisible by `p`.
Every `2 x 2` minor is therefore divisible by `p`, independently of whether
the Apéry value vanishes.  Consequently

`prod_(M<p<=n, n!=2p-1) p | Theta_n`,

and the prime number theorem gives

`log Theta_n >= n/2+o(n)`.

Hence no nonzero residual with `g(n)=0` can be `exp(o(n))`, regardless of
degree growth, signs, or coefficient size.  Combined with Q521's dominance
argument for bounded-height `g` with `g(n)!=0`, this closes Q521's proposed
short signed-adjoint target and every bounded-height single-residual strategy.
It does not exclude an exponentially large value vector engineered to cancel
to a small nonzero integer; producing such a vector is essentially the full
global Bezout/gcd problem and requires independent coefficient control.  A
proof of `Delta_n=exp(o(n))` cannot use the proposed `g(n)=0` shortcut.

The determinant formula and universal top-primorial divisibility are checked
in `projects/Q-series-and-Chan-s-work/q32_adjoint_obstruction.py`.

## 53. Two-dimensional lattice-index formulation

The same calculation gives an exact two-dimensional reformulation.  Let

`w_j=ell_j(1,b_j) in Z^2`

and let `L_n` be their integer span.  Its index is

`Theta_n=gcd_(i<j) ell_i ell_j(b_j-b_i)`.

Since `sum_j ell_j=1`, the vector

`(1,C_n)=sum_j w_j`,

where `C_n=sum_j ell_j b_j`, lies in `L_n`.  Smith reduction in dimension
two gives

`L_n=< (1,C_n),(0,Theta_n) >`.

After adjoining `(1,b_n)`, the index becomes `Delta_n`.  Consequently

`Delta_n=gcd(Theta_n,b_n-C_n)`,

and

`[L_n+Z(1,b_n):L_n]=Theta_n/Delta_n`.

For every non-boundary top-half prime, `v_p(Theta_n)=1`: the cardinal
coordinate `ell_r` is a unit, all other coordinates have valuation one, and
one of `b_0=1,b_1=5` differs from `b_r` modulo `p`.  Thus the top-half support
of `Theta_n` is the *entire* candidate primorial, while the top-half support
of `Delta_n` is exactly the bad radical.  The quotient `Theta_n/Delta_n`
therefore contains exactly the good top-half candidates.

This is an exact and potentially useful dual formulation, but a lower bound
for the quotient is meaningful only after controlling its small-prime part.
Otherwise smooth invariant factors can imitate the missing good-prime mass.
The identities and exact top-prime valuations are verified through `n=400`
in `q32_adjoint_obstruction.py`.

## 54. Q533 closes fixed differential/shift operations on the content family

Q533 independently verifies Section 49 and sharpens the universal-factor
description.  If

`T_n(c)=sum_(d=0)^n C_d(n)c^d`

and `P_n` is the non-boundary top-half primorial, then for every `d>0`,

`P_n/gcd(P_n,d) | C_d(n)`.

Indeed, modulo a candidate prime `p`, only degrees `0` and `p` survive:

`C_0=b_j`, `C_p=2b_j`, and `C_d=0` for `d!=0,p`.

Thus every nonconstant coefficient carries almost the full candidate
primorial universally.  Positive derivatives erase the selective direction
because `d(c^p)/dc=0` in characteristic `p`; fixed differential operators
either remain universal or retain only their zero-order multiple of `b_j`.
Finite shifts are likewise rank one:

`T_n(a)=b_j(1+2a) (mod p)`

for integral `a`, so every fixed shift combination has only one selective
scalar.

The simplest second certificate,

`T_n(1)-T_n(0)=sum_(d>0)C_d(n)`,

is selective but has a proved positive linear exponent (Q533 obtains a lower
rate about `2.742n`).  Therefore fixed differential, finite-shift, and Ore
operations do not lower the height.

Q533's proposed “signed global functional” is not itself an advance: proving
a subexponential nonzero functional with the required top-prime localization
is exactly the missing global content theorem.  The only untested extension
left from this family is a genuinely degree-growing algebraic/potential-
theoretic norm construction; Q543 audits that separately.

## 55. Algebraic-center norm probe

If `alpha` is an algebraic integer of degree `d`, coefficient divisibility
gives

`Gamma_n^d | Norm(T_n(alpha))`.

This avoids the denominator penalty of a rational center and replaces the
height by the average over all conjugate centers.  Raw real evaluation is
indeed much smaller for some nonintegral centers: at `n=180`, values near
`c=-1/2` have logarithmic rate about `1.54`, versus about `2.00` at `c=0`.
But multiplying by the denominator contributes `log 2`, making the rational
certificate worse.

`q32_algebraic_centers.py` exactly scans every irreducible monic quadratic

`x^2-sx+t`, `|s|,|t|<=12`,

using the norm formula

`Norm(u+v alpha)=u^2+suv+tv^2`.

The best polynomial throughout `n=60,120,180` is `x^2-x-1`; its averaged
rates are approximately `2.20,2.29,2.31`, all worse than the integral center
`0`.  Complex quadratic centers are also included.  This rules out only the
small quadratic search.  Reopening requires a high-degree algebraic-integer
family whose conjugate distribution beats the relevant weighted capacity;
Q543 is testing whether potential theory permits such a family at all.

A separate 200-digit exploratory evaluation at `n=180` tested
`alpha=zeta_m`, `zeta_m-1`, and `zeta_m+1` for every `3<=m<=100`, averaging
over primitive conjugates.  The best observed rate was about `2.599`, again
worse than `c=0`.  This is not an exact norm proof, but it makes the simplest
cyclotomic and shifted-cyclotomic families unattractive.

## 56. Q539 exact odd-quotient formula for the smooth Smith factor

Write `Delta_n=R_n U_n`, with `R_n` the exact squarefree top-half bad radical.
Q539 completely resolves the valuation mechanism in `U_n`.

For every sufficiently large `n` and every prime

`sqrt(n)<p<=M=floor(n/2)`,

put `q=floor(n/p)` and `r=n mod p`.  Then

`v_p(Delta_n)=1`

if and only if

`q` is odd and `p|b_r`;

otherwise the valuation is zero.  The proof classifies the `p`-unit
coordinates of the Lagrange vector exactly.  In the no-carry case two unit
coordinates force `b_n=1=5 (mod p)`, impossible.  In the carry case the unit
indices are `r,p+r,...`; Lucas makes their Apéry values
`b_0b_r,b_1b_r,...`, so rank loss is exactly `b_r=0`.

For `p<=sqrt(n)`, the divisor `Delta_n|4 lcm(ell_0,ell_1)` gives total
valuation weight `O(sqrt(n))`.  Therefore

`log U_n = sum_(sqrt(n)<p<=M, floor(n/p) odd, p|b_(n mod p)) log p`

`          + O(sqrt(n))`.

Equivalently, it is the sum of the fixed odd quotient slices
`q=3,5,7,...`.  Dropping the Apéry-zero condition gives exactly

`(log 2-1/2)n+o(n)`,

so Q521's constant is the true unresolved odd-slice mass, not a loose lcm
artifact.

Consequently `log U_n=o(n)` is equivalent to the family of fixed-slope
theorems

`sum_(n/(q+1)<p<=n/q, p|b_(n-qp)) log p=o(n)`

for every fixed odd `q>=3`.  This is still a horizontal moving-prime-divisor
problem; local Smith rank, carries, reflection, and zero separation provide
no further condition.  Q546 attacks `q=3` directly.

The exact support and exponent-one statement is independently verified
through `n=1200` in
`projects/Q-series-and-Chan-s-work/q32_delta_smooth.py`.
The only pre-asymptotic failures in `10<=n<22` are the fixed prime `p=5` at
`n=11,21`; the verification covers all 62,441 eligible prime-index pairs for
`22<=n<=1200`.

## 57. Exact `q=3` incidence data

An exact modular recurrence sweep of the first odd slice through `n=20,000`
found at most two simultaneous `q=3` hits.  The first double hit is

`n=743: (p,r,j)=(197,152,44),(211,110,100)`.

Later doubles include `n=1667,5518,7322,8244,10059,12634,19835`.  The
normalized mass of the first double is about `0.01431` and decreases rapidly
in the later examples.  This rules out uniqueness even in the `q=3` slice
but strongly supports the required `o(n)` mass.

The reproducible scanner is
`projects/Q-series-and-Chan-s-work/q32_fixed_q_scan.py`.  The data are
evidence only; as in `q=1`, a slowly growing random maximum cannot be turned
into a pointwise theorem without horizontal arithmetic input.

## 58. Q545: the smooth part of `Theta_n` is exactly the Q539 obstruction

Q545 supplies the full local classification of the lattice carrier.  Put

`U_p(n)={0<=j<M:p does not divide ell_j}`.

Then `p|Theta_n` if and only if the Apéry values `A_j` are constant modulo
`p` on `U_p(n)`.  The unit support itself is described by the simultaneous
Lucas--Kummer conditions

`p does not divide C(n,j)` and
`p does not divide C(n-j-1,M-1-j)`.

For `p>sqrt(n)` and `p>=7`, the resulting classification is:

- in quotient `q=1`, `p|Theta_n` universally, except at the odd midpoint
  `n=2p-1`;
- in every odd quotient `q>=3`, `p|Theta_n` exactly when
  `p|A_(n mod p)`;
- in an ordinary even quotient, `p` does not divide `Theta_n`;
- the even boundary `p=M` is a single universal exception when `M` is
  prime.

Every listed large prime has valuation exactly one.  The fixed primes `2`
and `5` must be absorbed in the small-prime estimate; for example
`(n,p)=(21,5)` violates the large-prime parity classification.  This
exception was found by the local hostile audit, not stated explicitly in
Q545's headline.

There is also an important correction to the informal phrase "full
top-half primorial": when `n=2M+1` and `p=M+1` is prime, all Lagrange
coordinates are units and `A_0=1` differs from `A_1=5`, so this boundary
prime is absent from `Theta_n`.  Earlier scripts already excluded
`n=2p-1`, but future statements must retain the exclusion.

Since `Theta_n` divides the minor from indices `0,1`, Kummer gives

`sum_(p<=sqrt n) v_p(Theta_n) log p=O(sqrt n)`.

After removing the corrected universal top carrier, the primes above
`sqrt n` in the residual factor are therefore exactly the same odd
quotient Apéry-zero slices that occur in `U_n=Delta_n/R_n`, up to the
single factor `M`.  Thus

`log(Theta_n/P_top)=log U_n+O(sqrt n+log n)`.

The apparently polynomial residual seen through `n=1200` is not a new easy
factor.  Proving it polynomial would already give `O_q(1)` hits in every
fixed odd quotient slice.  Q545 therefore closes the hope that the
two-dimensional lattice makes the smooth Smith factor easier.

The complete large-prime classification, exponent-one assertion, midpoint
exception, and the `p=5` cutoff are independently verified through `n=400`
in `q32_adjoint_obstruction.py`.

## 59. A folded double certificate for every fixed odd quotient

Let `q>=1` be fixed and odd, `n=qp+r`, and fold by Apéry reflection:

`j=min(r,p-1-r)`.

Both branches give the exact bound

`(2q+1)j<=n-q`,

so `j<=J_q=floor((n-q)/(2q+1))`.  Define

`D_(n,q)=sum_(k<=J_q) C(n,k)^2 C(n+k,k)^2`

and

`S_(n,q)=sum_(k<=J_q) C(n,k)C(n+k,k)F_k`.

On the branch `r=j`, terms after `j` vanish through `C(n,k)`.  On the
branch `r=p-1-j`, they vanish through the carry in `C(n+k,k)`.  Lucas and
the Strehl--Franel identity consequently give, for every prime in the fixed
slice,

`D_(n,q)=S_(n,q)=A_j (mod p)`.

Thus the entire bad fixed-`q` radical divides `gcd(D_(n,q),S_(n,q))`.
This is a genuine global integer compression, although its individual
heights remain exponential.  For `q=3`, an exact scan through `n=1200`
gives dyadic maxima for `log gcd/n`

`.19851,.09444,.15806,.09840,.06842,.04279,.02685`.

The last maximum occurs at `n=666`, with gcd `58,460,875`.  The trend is
stronger than the single-incidence data, but no resultant or Bezout bound is
yet known.  As in the `q=1` certificate, `D-S` is universally zero modulo
every candidate prime, so the two congruences are not locally independent.

The formulas and local congruences are independently checked in
`q32_fixed_q_certificates.py`.  Q555 asks for a residue-class recurrence,
resultant, or a third transform that can turn the observed small gcd into a
proof.

## 60. Q543 narrows, but does not close, algebraic-center norms

For the q=1 Legendre--Euler polynomial

`T_n(c)=sum_d C_d(n)c^d`, `Gamma_n=gcd_d C_d(n)`,

and an algebraic integer `alpha` of degree `d`,

`Gamma_n^d | Norm(T_n(alpha))`

whenever the norm is nonzero.  Q543 derives the complex saddle equation at
the boundary ratio `beta=1/3`:

`(1+beta)z^2-2(1+2c)z+(1-beta)=0`.

The absolute saddle envelope at `c=0` recovers `log 8`; at `c=-1/2` it is
`(log 136)/3`.  Clearing the denominator of the latter adds `log 2`, so the
rational center does not improve the certified exponent.

Capacity excludes an unbounded-degree algebraic-integer family with every
conjugate confined to `[-1,0]`, whose capacity is `1/4`.  Standard
equidistributed cyclotomic, shifted-cyclotomic, bounded-house Pisot, and
Salem families also retain a positive Mahler-measure floor.  This agrees
with the exact quadratic and exploratory cyclotomic scans.

The unrestricted conclusion is weaker and must not be overstated.  Growing
degree, growing house, or a family tuned to zeros of the moving polynomial
is not covered by fixed-set capacity.  The exact remaining question is the
moving resultant infimum

`delta_n=inf_P log|Res(P,T_n)|/(n deg P)`,

over monic irreducible integer `P` with nonzero resultant.  A positive
`liminf delta_n` would close the algebraic-center strategy negatively; an
explicit family with `delta_n=o(1)` would yield a subexponential norm
certificate.  Neither direction is known.  Reopen this route only with a
moving weighted integer-Chebyshev/resultant estimate or an explicit
growing-degree family, not another bounded quadratic search.

## 61. Full Legendre--Euler content for the `q=3` slice

The Section 49 coefficient formula works with the smaller cutoff

`J=floor((n-3)/7)`.

Call the resulting polynomial `T_(n,3)(c)`.  For every prime
`n/4<p<=n/3`, with folded residue

`j=min(n-3p,p-1-(n-3p))`,

the stronger coefficientwise identity is

`T_(n,3)(c)=A_j Q_3(c^p) (mod p)`,

where

`Q_3(t)=1+12t+30t^2+20t^3`.

There is now a short proof, not only a computation.  The Legendre
polynomials satisfy the coefficientwise Lucas factorization

`Q_(qp+r)(t)=Q_q(t^p)Q_r(t) (mod p)`.

After substituting `t=c+y`, every coefficient of `y^m` with `m<p` acquires
the common factor `Q_q(c^p)`.  On the direct folded branch `r=j`, the
remaining polynomial has degree `j`.  On the reflected branch
`r=p-1-j`, the binomial coefficients above degree `j` vanish and

`Q_(p-1-j)(t)=Q_j(t) (mod p)`.

Finally, the binomial transform defining `g_m(c)` is exactly inverse to the
translation coefficients `[y^m]Q_j(c+y)`, so the complete sum through
degree `j` is `A_j`.  Since `j<=J<p` in the q=3 slice, this proves the
displayed polynomial congruence.

The identity has also been checked for every eligible pair through `n=400`; the
first independent exhaustive check through `n=150` compared every
coefficient with the four predicted degrees `0,p,2p,3p`.  Hence the q=3 bad
radical divides

`Gamma_(n,3)=content_c T_(n,3)(c)`.

This strictly uses more certificates than the direct/Strehl gcd.  Its exact
dyadic maxima through `n=400` are

`.14631,.09444,.13330,.09509,.06293,.03316`.

The last maximum occurs at `n=394`, where

`Gamma_(394,3)=471665=5*17*31*179`.

The generalized Lucas profile is valid more broadly whenever the cutoff is
strictly below `p` and the folded residue is at most the cutoff.  It cannot
silently cover every larger quotient: at `(n,p,q)=(80,11,7)`, the cutoff
equals `p`, a second base-`p` digit enters, and the naive
`A_j Q_q(c^p)` formula fails already in the constant coefficient.  With the
necessary condition `J<p`, the general profile was verified in 813
prime-index pairs through `n=180`.

The implementation is `q32_fixed_q_content.py`.  Q556 is auditing the
general proof and searching for a content/Smith bound.  The exact remaining
obstruction is again global: modulo each candidate, the coefficient vector
has rank-one profile `A_j Q_3(c^p)`, so adding coefficients does not create
local independence.

An attempted exact extension from `400` to `1200` computed each new
coefficient only modulo the gcd accumulated so far.  This avoids constructing
the full integers but was still too slow: fixed small factors often keep the
gcd nontrivial, forcing evaluation of every coefficient, and repeated exact
binomial construction remains quadratic-to-cubic in practice.  Stripping
primes at most `sqrt(n)` did not remove enough unrelated rough factors.
The run was stopped without using partial output.  Reopen the larger scan
with modular binomial recurrences or a C++ multiprecision implementation;
the proof does not depend on extending the numerical range.

## 62. Existing subspace-theorem gcd results do not cover the content sequences

The observed small content suggests comparing a selective coefficient such
as the Strehl prefix with the explicit leading coefficient

`C_n=(-1)^J C(2n,n)C(n-1,J)`.

For most tested `n`, `gcd(C_0,C_n)` already equals the full content; the
exceptions can require many additional coefficients.  This resemblance to
the Bugeaud--Corvaja--Zannier theorem for
`gcd(a^n-1,b^n-1)` prompted a literature audit.

The available unconditional moving-target theorems apply to fixed-degree
polynomials evaluated at moving `S`-unit points, and their recurrence
corollaries apply to algebraic constant-coefficient linear recurrences.
The relevant primary statement is Grieve--Wang, *Greatest common divisors
with moving targets and consequences for linear recurrence sequences*
(arXiv:1902.09109), whose hypotheses explicitly require fixed degree, slow
moving coefficient height, and `S`-unit coordinates.  The related
Corvaja--Zannier/Levin results have the same toric input.

Neither `C_0(n)` nor `C_n(n)` is a constant-coefficient recurrence or an
evaluation at a fixed-dimensional `S`-unit point.  Along residue classes
modulo seven they are P-recursive/diagonal sequences; coefficient extraction
and factorial ratios introduce new primes and do not produce fixed
`S`-units.  No theorem located extends the subexponential gcd conclusion
from C-finite sequences to arbitrary P-recursive diagonals.  Therefore
invoking the subspace theorem here would be an unsupported category jump.

The leading-coefficient shortcut also fails computationally in a concrete
way.  For `3<=n<=160`, `gcd(C_0,C_n)` is strictly larger than the full
content in 46 cases.  At `n=142`, one must include 117 coefficients from the
high-degree end, in addition to `C_0`, before the gcd reaches the true
content.  Thus no fixed collection of leading coefficients is supported by
the data; the small content genuinely uses a degree-growing family.

Reopen this route only with a theorem stated for arithmetic diagonals or
P-recursive sequences, or with an exact reduction of these particular two
sequences to a fixed-dimensional `S`-unit evaluation.

## 63. Q546/Q555: exact factorial collapse of the q=3 double certificate

Q546 does not prove the q=3 estimate.  It independently recovers the two
folded branches, the cutoff `(n-3)/7`, and the direct/Strehl congruences.  Its
useful negative conclusions agree with the local audit: standard block
vectors are rank one over each candidate field, the recurrence crosses the
three singular indices `p,2p,3p`, and the reflection continuant has a forced
factor `p`.  Its claim that no multi-prime q=3 example was available is
false: the exact scan already gives

`n=743: p=197,211`

and

`n=1667: p=431,499`,

with further doubles recorded in Section 57.  Uniqueness cannot be used.

Q555 supplies a stronger exact identity.  Write

`L_k(X)=C(X,k)C(X+k,k)`

as a polynomial in `X`, and

`D_J(X)=sum_(k<=J)L_k(X)^2`,
`S_J(X)=sum_(k<=J)L_k(X)F_k`.

For every integer `0<=m<=J`, both truncated sums become complete and equal
`A_m`.  The symmetry `L_k(-X-1)=L_k(X)` supplies the reflected roots
`-1,...,-J-1`.  Hence

`D_J(X)-S_J(X)=L_(J+1)(X) H_J(X)`

over `Q[X]`, and the denominator audit gives

`(J!)^2(D_n-S_n)=B_n Hcal_J(n)`,

where

`B_n=C(n,J+1)C(n+J+1,J+1)`

and `Hcal_J` is integral.  For `J=floor((n-3)/7)`, every candidate prime
`n/4<p<=n/3` divides `B_n`; its valuation is one except at the possible
reflection midpoint `7p=2n+1`, where it is two.  Therefore, after restricting
prime support to this interval,

`rad gcd(D_n,S_n)=rad gcd(S_n,B_n)`.

This identifies the universal reason that the two residues agree: their
difference contains an explicit factorial carrier.  The target is now the
concrete same-index theorem

`sum_(n/4<p<=n/3, p|gcd(S_n,B_n)) log p=o(n)`.

The normalized quotient does not add a hidden second condition.  At the
first two double hits, exact integer division gives

`Hcal_105(743)=149 (mod 197)`, `=104 (mod 211)`,

and

`Hcal_237(1667)=77 (mod 431)`, `=219 (mod 499)`.

Thus badness does not imply an additional factor of `p` after removing
`B_n`; the proposed p-squared certificate fails already at `n=743`.

Q555 also independently derives the full parameterized
Legendre--Euler family and the coefficientwise profile

`T_(n,3)(c)=A_j Q_3(c^p) (mod p)`,

confirming Section 61.  It does not bound the content.  Q568 and Q569 now
attack respectively the explicit factorial/P-recursive gcd and the Smith
structure of the growing coefficient family.

## 64. The same factorial collapse sharpens the original q=1 formulation

The polynomial factorization in Section 63 is valid for every cutoff.  With

`J=floor((n-1)/3)`,

let `S_n` be the q=1 Strehl prefix and set

`B_n=C(n,J+1)C(n+J+1,J+1)`.

Every prime `n/2<p<=n` divides `B_n`.  Its valuation is one except at the
possible reflection midpoint `3p=2n+1`, where it is two.  The folded Lucas
calculation gives

`S_n=A_j (mod p)`.

Consequently the original q=1 bad radical is exactly the top-half radical
of

`gcd(S_n,B_n)`.

The carrier has exponential rate `log 4`, compared with `log 8` for `S_n`
and `log 16` for the direct positive prefix.  Its rate is comparable to the
earlier central-binomial carrier, but the new identity explains it as the
complete universal factor of the direct--Strehl difference:

`(J!)^2(D_n-S_n)=B_n Hcal_J(n)`.

The normalized quotient again supplies no extra p-adic condition.  At the
three-hit example `n=321`,

`Hcal_106(321)=80,119,41`

modulo `179,193,211`, respectively.  Thus each bad prime occurs only through
the universal factor `B_n`, not through a forced square.

Exact scanning through `n=1200` gives dyadic maxima for the full
`log gcd(S_n,B_n)/n`

`.29263,.19786,.22411,.08374,.10617,.06705,.03553`.

This full gcd contains harmless smooth and off-interval factors, but its
top-half support is exactly the target.  The reduction and all candidate
valuations are verified in `q32_strehl_gcd.py`.  Q572 replaces the failed
Q534 task with a focused attack on this factorial-carrier formulation.

## 65. Q556 confirms the general content law and isolates moving positions

Q556 independently proves the general fixed-quotient formula.  If
`n=qp+r`,

`J=floor((n-q)/(2q+1))`, `j=min(r,p-1-r)`,

then

`T_(n,q)(c)=A_j Q_q(c^p) (mod p)`.

It also records the exact coefficient formula already implemented in
Section 61 and the q=3 Smith profile

`A_j*(1,12,30,20)`

at degrees `0,p,2p,3p`.  All other coordinates are universally zero.  Thus
the content and `gcd(S_n,B_n)` have the same q=3 target support; the extra
coefficients reduce unrelated integer factors but add no local equation.

The only potentially new functional must access the four positions in a
way that depends on `p`.  Fixed derivatives erase them, while the relevant
Hasse derivative order is `p,2p,3p`.  Ordinary resultants also fail: the
effective degree modulo `p` is `3p`, so the selective principal
subresultant has a p-dependent index and its standard determinant height is
`exp(O(n^2))`.

## 66. Growing Euler-derivative functionals: exact modular coset obstruction

There is a nontrivial-looking degree-growing family.  For
`P_n(x) in Z[x]` with `P_n(0)=1`, define

`R_n^+=P_n(theta)T_(n,3)(1)`,
`R_n^-=P_n(theta)T_(n,3)(-1)`,

where `theta=c d/dc`.  Since `P_n(ap)=1 (mod p)`, their local scalars are
respectively

`Q_3(1)=63` and `Q_3(-1)=-1`.

Arbitrary degree therefore preserves selectivity for `p>7`, which escapes
the *fixed-order* wording of the earlier obstruction.

However, every positive moment

`M_k^+=theta^k T(1)`, `M_k^-=theta^k T(-1)`, `k>=1`,

is universally zero modulo each candidate prime: the only nonzero degrees
are `p,2p,3p`, and multiplication by `d^k` kills them.  Hence all such
growing functionals lie in the fixed cosets

`T(1)+G_t^+ Z`, `T(-1)+G_t^- Z`,

where `G_t^±=gcd(M_1^±,...,M_t^±)` contains the full q=3 candidate
primorial.  Degree growth can create Archimedean cancellation only inside
this coarse modular coset.

An exact computation with all moments through degree `n` gives least
nonzero centered-residue rates between `0.713` and `0.815` on representative
`30<=n<=160`, far worse than the candidate-primorial rates and the observed
content.  The moment ideals also contain large unrelated factors.  Thus the
most natural `Z[x]` growing differential family does not numerically approach
subexponential height.

The reproducible audit is `q32_growing_functional.py`.  Q575 searches for a
more structured saddle-cancelling polynomial; Q578 audits whether any such
choice can evade the exact modular coset rather than merely choosing a
different representative of it.

## 67. Literature audit on Apéry prime factors

A targeted primary-source search found no existing theorem controlling the
pointwise top-half radical of `A_n`.

- Malik--Straub, *Divisibility properties of sporadic Apéry-like numbers*
  (arXiv:1508.00297), proves Lucas congruences and studies primes that divide
  no term.  It does not give a pointwise radical bound.
- Delaygue, *Arithmetic properties of Apéry-like numbers*
  (arXiv:1310.4131), gives lower p-adic valuation laws once zero digits are
  present.  In the one-digit top-half situation this recovers the single
  forced factor and does not add cross-prime control.
- The paper *Arithmetic properties of Apéry numbers* proves only
  density-one lower bounds for the number and size of prime factors; its
  quantifiers and direction do not address the required upper bound.
- The published note *A remark on Apéry's numbers* explicitly refutes naive
  extensions of Lucas congruences modulo `p^2` and `p^3`, consistent with the
  direct counterexamples to a forced extra `p` in Sections 63--64.

Thus no located prime-factor theorem closes the factorial-carrier gcd.
Reopen the literature route only with a theorem explicitly bounding a
large-prime radical of one P-recursive/diagonal term for every index, not an
average or density-one largest-prime-factor result.

## 68. Exact polynomial resultant is global but far too tall

For fixed `J`, let

`B_J(X)=(J+1)!^2 L_(J+1)(X)`

be the monic carrier polynomial and

`Sbar_J(X)=J!^2 S_J(X)`.

At the `2J+2` roots

`0,...,J,-1,...,-J-1`

of `B_J`, symmetry gives the values

`J!^2 A_0,...,J!^2 A_J`

twice.  Therefore the resultant is exactly

`Res_X(B_J,Sbar_J)`
`=J!^(4(J+1)) prod_(m=0)^J A_m^2`.

This was symbolically verified through `J=8` in
`q32_factorial_resultant.py`.  Even after ignoring the factorial clearing,
the Apéry product has logarithmic height `Theta(J^2)`; the displayed
integral clearing adds `Theta(J^2 log J)`.  With `J` proportional to `n`,
the ordinary Bezout/resultant certificate is vastly larger than linear,
let alone `o(n)`.

The exact formula also exposes the rainbow obstruction: the resultant pays
for every possible folded root `m<=J`, whereas a fixed `n` uses only the
single root aligned with each different prime.  A useful resultant must
localize those prime-dependent roots without multiplying all Apéry values.
Q574 is auditing whether divided differences or p-orderings can do this.

## 69. Exact branch-separated q=1 data through `n=100,000`

The C++ modular-recurrence scanner `q32_scan.cpp` was compiled with `-O3`
and run exactly through `100,000`.  The q=1 slice has at most three
simultaneous hits in this range.  The dyadic maximum on
`50,000<n<=100,000` is

`mass/n=0.00050856678 at n=63887`,

with three hits.

The maximum three is attained at

`n=321,11576,18444,22101,26164,47066,47859,63887,64555`.

Branch separation refutes two tempting stronger guesses:

- `n=11576` has three direct hits
  `p=8893,9319,11437`;
- `n=321` has three reflected hits
  `p=179,193,211`.

Thus neither branch has a uniqueness or at-most-two theorem.  Mixed triples
also occur.  The direct-branch dyadic maximum on the last range is
`0.00040597238`; the reflected maximum is `0.00033397788`.

The full lower channel reaches six simultaneous hits, so the q=1 sparsity
does not extend verbatim after combining quotient slices.  These data are
strong evidence for `o(n)` but do not imply a bounded hit count; a random
maximum can grow arbitrarily slowly.  Q585 is using the exact branch data
to rule out unsupported uniqueness arguments.

## 70. Euclidean reduction modulo the split factorial carriers fails

The q=1 carrier factors canonically as

`B_n=B_n^- B_n^+`,

where, for `J=floor((n-1)/3)` and `m=J+1`,

`B_n^-=C(n,m)`, `B_n^+=C(n+m,m)`.

Kummer's carry criterion separates the two folded branches: every direct
candidate divides `B_n^-`, while every reflected candidate divides `B_n^+`.
This suggested a much cheaper certificate than a polynomial resultant:
replace the Strehl prefix `S_n` by its least centered residue modulo the
relevant carrier.  Every bad prime in that branch still divides the residue.
If either residue had logarithmic height `o(n)`, that branch would be closed.

Exact computation through `n=1200` decisively rejects the natural
representatives.  At `n=1200`,

`log |center(S_n mod B_n^-)|/n = 0.6325578`,

`log |center(S_n mod B_n^+)|/n = 0.7445801`,

and

`log |center(S_n mod B_n)|/n = 1.3778913`.

The corresponding carrier rates are tending respectively to the positive
entropy constants

`H(1/3)=0.636514...`,

`(4/3)H(1/4)=0.749780...`,

and `log 4=1.386294...`.  Across the last dyadic range the centered residues
remain close to the full carrier scale rather than decreasing.  Thus a
single Euclidean quotient supplies no Archimedean compression; it behaves
like a generic residue.

The reproducible computation is `q32_carrier_remainders.py`.  Reopen this
route only with a structured *family* of Bezout multiples whose attainable
residue lattice is proved to have a subexponential nonzero vector.  Merely
centering `S_n` modulo either binomial factor or their product cannot prove
the target.

## 71. Q568 correction: the carrier has a unique doubled midpoint

Q568 found and the local scripts now verify a valuation correction to
Sections 63--64.  Put `K=J+1`.  Since `K<p` in either fixed slice, Legendre's
formula gives

`v_p(C(n,K))=1_[r<K]`,

`v_p(C(n+K,K))=1_[p-r<=K]`,

where `n=qp+r`.  Thus

`v_p(B_n)=1_[r<K]+1_[p-r<=K]`.

Here the general fixed-odd-q cutoff is

`K=floor((n+q+1)/(2q+1))`.

The inequalities are exact:

`r<K <=> r<=(p-1)/2`,

`p-r<=K <=> r>=(p-1)/2`.

Thus the two indicators overlap only at the reflection-fixed residue
`r=(p-1)/2`, equivalently

`(2q+1)p=2n+1`.

Consequently:

- for `q=1`, the unique possible doubled carrier prime satisfies
  `3p=2n+1`;
- for `q=3`, it satisfies `7p=2n+1`.

All other candidate primes occur in `B_n` exactly once.  The midpoint is at
most one prime for each `n` and contributes only `O(log n)` to the target, so
every radical/support reduction remains valid.  What was wrong was only the
unqualified valuation-one claim.  The exhaustive assertions are now part of
`q32_strehl_gcd.py` and `q32_fixed_q_certificates.py`.

Q568 otherwise does not prove the q=3 bound.  Its useful sharpening is the
exact carrier valuation formula and the observation `p>=2j+1` on both q=3
folded branches.  Since all factorials inside `A_j` are then p-adic units,
the event `p|A_j` is cancellation among units; Kummer carry counting cannot
distinguish bad primes.  Its proposed moving-divisor theorem is another exact
statement of the surviving horizontal obstruction.

## 72. Q569: exact content Smith block, but the target radical survives intact

Q569's matrix claim has been reconstructed from the coefficient formula and
independently verified for both `q=1` and `q=3` through `n=180`.  Let `M_n`
be the integer matrix taking arbitrary input data `x_0,...,x_J` to the
coefficients of the Legendre--Euler truncation, and let `P_(n,q)` be the
squarefree product of all primes in the fixed-q slice.  If `m_(d,i)` is its
entry in degree `d`, then

`P_(n,q) | m_(d,i)-m_(d,0)m_(0,i)`

for every `d>0,i>0`.  Unimodular column and row operations therefore give
the exact integral block form

`M_n ~ [1,0; 0,P_(n,q) W_n]`.

For the Franel input, writing `S_n` for the zero-degree coefficient and

`H_d=sum_(i=1)^J W_(d,i) F_i`,

one gets the exact content identity

`Gamma_(n,q)=gcd(S_n, P_(n,q) gcd_(d>0) H_d)`.

The low `J x J` block of `W_n`, with its columns reversed, is triangular.
Its diagonal entries are

`+- (B_n/P_(n,q)) C(J+1,i)`.

At every non-midpoint candidate prime `p`, these entries are p-units.  Hence
`W_n` is invertible modulo `p`.  Since the Franel tail is nonzero modulo odd
`p` (`F_1=2`), a bad non-midpoint candidate occurs in `Gamma_(n,q)` to
valuation exactly one.  This explains rigorously why the normalized
direct--Strehl quotient is a p-unit at the `n=743` double hit.

The result is structurally useful but does not prove a saving.  Every
nontrivial Smith invariant contains the *entire* candidate primorial before
the special Franel vector is used.  Q569 gives an explicit CRT construction
of a primitive artificial input vector whose content contains all candidates,
so the matrix, its Smith form, and input primitivity alone cannot distinguish
the target.  Ordinary determinant bounds for `W_n` have quadratic logarithmic
height.  The missing input remains cross-prime arithmetic special to the
Franel vector.

The regression audit is `q32_content_smith.py`.  Reopen the content route only
with either a Franel-specific p-adic saturation theorem or a growing signed
dual functional of subexponential nonzero value; abstract Smith reduction
cannot remove the horizontal radical.

## 73. Q575: linear Krawtchouk filters survive formally but fail the first scan

Q575 correctly identifies the exact degree-growing moment lattice.  In
falling-factorial coordinates,

`D_k^eps=sum_d eps^d C(d,k) C_d`,

and degree-`m` integer-valued Euler filters attain precisely the affine
lattice

`D_0^eps + <D_1^eps,...,D_m^eps>_Z`.

At full degree the binomial transform is unimodular, so the shortest nonzero
value is another formulation of the coefficient content.  Every q=3
candidate divides all positive moments through the relevant degree; hence
the affine lattice still contains the universal candidate primorial in its
step size.  Formal Bezout cancellation alone does not prove that its least
nonzero representative is small.

The one named family not excluded by Q575 is the binary Krawtchouk filter

`K_m(x;N)=[z^m](1+z)^(N-x)(1-z)^x`,

with `N` below the least q=3 candidate prime.  Since every candidate then has
`p>N`,

`K_m(ap;N)=K_m(0;N)=C(N,m) (mod p)`,

and the signed value

`R^-(n,N,m)=sum_d (-1)^d K_m(d;N) C_d`

satisfies

`R^-(n,N,m)=-C(N,m) A_j (mod p)`.

The multiplier is a p-unit, so every nonzero such value is a valid selective
certificate.

The plus functional was checked as well:

`R^+(n,N,m)=63 C(N,m) A_j (mod p)`.

Its multiplier is a unit for candidate primes `p>7`.

Q575's stated baseline exponential rate is not correct.  It replaces the
global entropy of `L(n,J)` by the *local endpoint ratio* 48.  At
`beta=J/n=1/7`,

`log L(n,J)/n -> H(1/7)+(8/7)H(1/8)`,

not `log(48)/7`.  Since the binomial transform of the Franel numbers has base
9, the correct signed-baseline rate is

`H(1/7)+(8/7)H(1/8)+log(9)/7`,

a positive constant about `1.154`, consistent with the exact data.

An exhaustive exact scan of every `0<=m<=N<p_min` at

`n=40,60,80,100,120,160,200,240,300,400`

found no improvement at all: in every case the minimum was the unfiltered
choice `m=0` for both signs.  This range is substantially larger than Q575's proposed
`N<=J`; at `n=400` it exhausts `N<=100` while `J=56`.  The unfiltered rate
was `1.11060` for the minus sign and `2.36216` for the plus sign, and every
nonconstant Krawtchouk filter was larger.  This does not prove a positive asymptotic
lower bound, but it removes the only numerical support for the proposed
linear-degree survivor.

The reproducible local-congruence and exhaustive-filter audit is
`q32_krawtchouk_filters.py`.  Reopen this family only with a saddle analysis
showing a negative exponential contribution on a specific linear scaling
`m~eta n,N~nu n`; blind searches at the natural finite sizes show the
opposite behavior.

## 74. Endpoint asymptotic cancellation gains only a polynomial factor

There is another exact family hidden in the factorial carrier.  For fixed
odd `q`, set `K=J+1`,

`S_n=sum_(k<K)L(n,k)F_k`,

`T_n=L(n,K)F_K=B_n F_K`.

Every fixed-q candidate divides `T_n` universally and divides `S_n` exactly
when it is bad.  Hence every nonzero integer combination

`Q(n)S_n-P(n)T_n`

is a valid bad-prime certificate.

The forward ratio of endpoint summands tends to

`rho_q=32q(q+1)`.

Therefore

`S_n/T_n -> 1/(rho_q-1)`;

the constants are `1/63` for `q=1` and `1/383` for `q=3`.  The first natural
integer cancellation is

`E_n=(rho_q-1)S_n-T_n`.

Exact data through `n=1200` shows that `E_n` retains the same positive
exponential rate as `T_n`.  At `n=1200` the rates are respectively

- q=1: `log|E_n|/n=2.06176`, `log T_n/n=2.06738`;
- q=3: `log|E_n|/n=1.12607`, `log T_n/n=1.12978`.

Thus the cancellation is polynomial, not exponential.

The growing-order version faces an exact integrality cost.  Cancelling the
first `t` terms of a `1/n` expansion normally produces a rational
approximant with denominator containing `n^t`.  Multiplying through to obtain
an integer certificate restores precisely the `n^t` factor that the
asymptotic remainder gained.  A successful construction must therefore use
an integral recurrence or finite-difference identity whose coefficient
height grows substantially slower than its cancellation gain; formal
high-order asymptotics alone do not suffice.

The initial exact computation is `q32_endpoint_cancellation.py`; Q594 is
auditing the full residue-class expansions and whether adjacent endpoint
terms can evade the denominator obstruction.

## 75. The complete universal endpoint tail reduces only modulo the primorial

The single endpoint term in Section 74 can be enlarged substantially.  Let
`p_min` be the smallest prime in the fixed-q candidate interval.  Once either
Kummer threshold is crossed at `k=K`, it remains crossed, so every candidate
prime divides

`L(n,k)F_k`

for every

`K<=k<p_min`.

This supplies a linear-dimensional family of universal certificates.  Let
`G_tail` be their gcd and `P_(n,q)` the squarefree candidate primorial.  Since
`P_(n,q)|G_tail`, divide every tail term by the common nuisance factor

`H_tail=G_tail/P_(n,q)`.

The normalized tail terms then have gcd exactly `P_(n,q)`.  Therefore their
full integer span is not a new small lattice:

`<L(n,k)F_k/H_tail : K<=k<p_min>_Z=P_(n,q) Z`.

The best unrestricted affine reduction of `S_n` by the whole normalized tail
is simply its least centered residue modulo `P_(n,q)`.

Exact data shows no unexpected compression.  At `n=1200`:

- q=1: `log P_(n,1)/n=0.490947`, centered-residue rate `0.489069`;
- q=3: `log P_(n,3)/n=0.078090`, centered-residue rate `0.077221`.

Thus the residue is empirically at essentially the full primorial scale.
Adding adjacent endpoint terms removes the large smooth common carrier but
does not distinguish good candidate primes from bad ones.  Proving the
centered residue subexponential would require new arithmetic at least as
strong as a horizontal anti-alignment theorem; lattice dimension alone gives
only the trivial candidate-interval bound.

The regression computation is `q32_tail_lattice.py`.

## 76. Q588: the full coefficient family gives exact branch ideals

Q588 extends the Section 72 Smith block from q=3 to the original q=1 slice
and agrees with the independent reconstruction in `q32_content_smith.py`.
Its genuinely sharper formulation separates the two binomial carriers.  If
`Gamma_n` is the content of the complete q=1 Legendre--Euler coefficient
vector and `m=J+1`, define

`Gamma_n^-=gcd(Gamma_n,C(n,m))`,

`Gamma_n^+=gcd(Gamma_n,C(n+m,m))`.

Then the top-half support of `Gamma_n^-` is exactly the direct bad branch,
and that of `Gamma_n^+` is exactly the reflected bad branch; the unique
possible midpoint belongs to both.  These are the
least positive scalar ideals obtainable by arbitrary integer Bezout
combinations of the full coefficient family together with the corresponding
carrier.  Thus they strictly dominate any fixed collection of evaluations,
derivatives, divided differences, or contiguous cutoffs inside this
coordinate construction.

The enlarged family still has only one selective local direction.  Contiguous
cutoffs differ by terms universally divisible by the branch carrier; fixed
Legendre--Euler translations reduce to a scalar multiple of `A_j`; higher
divided differences are universal.  The low Smith block has a determinant of
quadratic logarithmic height even after removing one candidate primorial per
direction.  Therefore generic determinant, Cramer, Minkowski, or LLL bounds
remain exponential; an exceptional short vector would itself be the new
Franel-specific theorem.

An exact branch scan through `n=400` gives the following last dyadic maxima:

- direct `(320,400]`: rate `0.0345067` at `n=340`;
- reflected `(320,400]`: rate `0.0542410` at `n=321`;
- full content `(320,400]`: the same `0.0542410` at `n=321`.

At `n=321`, the reflected ideal is

`36447085=5*179*193*211`,

so it retains exactly the known three reflected bad primes, plus the small
factor 5.  The rapid numerical decay is encouraging but not a proof; the
large-prime support theorem remains exactly the desired horizontal bound.

The reproducible scan is `q32_branch_content.py`.  Reopen this route with a
Franel-adjoint telescoper or a proof that either branch ideal has
subexponential height, not with another coordinate change inside the same
rank-one lattice.

## 77. The full allowed integer-valued filter lattice stays exponential

Q575 tested named degree-growing filters, but the full allowed class has an
exact one-dimensional affine description.  For sign `eps=+-1`, put

`D_k^eps=sum_(d>=k) eps^d C(d,k) C_d`.

If `m<p_min`, where `p_min` is the least q=3 candidate prime, every
integer-valued polynomial

`P(x)=1+sum_(1<=k<=m) a_k C(x,k)`

satisfies `P(ap)=1 (mod p)` at the four local degrees `0,p,2p,3p`.
Therefore its attainable values are exactly

`D_0^eps + gcd(D_1^eps,...,D_m^eps) Z`.

The least nonzero value can be computed without LLL: it is the least centered
nonzero residue of `D_0^eps` modulo that gcd.  Coefficient height is irrelevant
once this final integer is known; Bezout gives an integer-valued filter
realizing it.  The constraint `m<p_min` is essential.  At full degree the
binomial transform is unimodular and reaches the small coefficient content,
but degrees `m>=p` no longer satisfy `C(ap,m)=0 (mod p)` because the
factorial denominator contains `p`; selectivity is then lost.

The exact optimum over *every* allowed `m`, for both signs, remains
exponential:

- `n=400`, `m_max=100`: minus `0.585610` at `m=97`, plus `0.591297`
  at `m=97`;
- `n=600`, `m_max=150`: minus `0.553965` at `m=149`, plus `0.556884`
  at `m=149`;
- `n=800`, `m_max=210`: minus `0.544535` at `m=199`, plus `0.546000`
  at `m=199`.

This is a real improvement over the unfiltered q=3 rate near `1.15`, but the
data point toward a positive constant near `0.5`, not toward zero.  It also
strictly contains every Krawtchouk filter from Section 73, explaining why a
named-family search cannot do better than this affine optimum.

The reproducible computation is `q32_binomial_filter_lattice.py`; Q599 is
auditing its Smith/Vandermonde structure and whether a positive lower rate can
be proved.

## 78. Q594 confirms fixed-order cancellation; its growing-order closure is conditional

Q594 derives the residue-class Birkhoff expansion of `S_n/T_n` and confirms
that the coefficient after the leading `1/(rho_q-1)` term is nonzero in
every residue class for q=1 and q=3.  This is consistent with the exact
Section 74 data and makes the polynomial-only gain of
`(rho_q-1)S_n-T_n` rigorous in principle.

Its stronger claim that every growing local endpoint jet retains a positive
exponential rate is not accepted as a theorem from the delivered page.  The
displayed constants and inequalities were lost in the Notion export, and the
argument charges exponential coefficient/denominator height as part of the
certificate cost.  For pure divisibility this is too restrictive: once an
integer combination has a small nonzero final value, the sizes of its Bezout
coefficients do not affect the bound.

Section 75 supplies the correct unrestricted lattice audit.  After dividing
the complete universal endpoint tail by its common nuisance factor, its
integer span is exactly `P_(n,q) Z`; all affine values are

`S_n+P_(n,q) Z`.

This proves that extra tail dimension alone gives only centered reduction
modulo the candidate primorial.  It does **not** prove a positive lower bound
for that centered residue.  Thus Q594 closes fixed-order inverse-power
clearing and controlled local Richardson schemes, but it does not close every
unrestricted growing-order integer combination.

## 79. The full q=1 integer-valued filter lattice also stays exponential

The Section 77 affine-lattice computation has a q=1 analogue with a larger
allowed degree.  Since every q=1 candidate satisfies `p>n/2`, the binomial
filter degree may grow to

`m_max=p_min-1~n/2`,

well beyond the coefficient input cutoff `J~n/3`, while retaining
`P(ap)=P(0) (mod p)`.

The exact least nonzero residue over every allowed degree and both signs
still has a stable positive rate:

- `n=300`, `m_max=150`: minus `0.859095`, plus `0.857866`;
- `n=400`, `m_max=210`: minus `0.810745`, plus `0.808752`;
- `n=600`, `m_max=306`: minus `0.820437`, plus `0.824274`.

The optima occur just below isolated drops in the binomial-moment gcd
(`m=149,199,293` in these samples), not at the formal full degree.  Thus even
using the extra interval between `J` and `p_min` does not numerically approach
the small full coefficient content.  Crossing `m=p` remains forbidden because
the integer-valued binomial denominator then contains the candidate prime and
the local unit-scalar congruence fails.

The reproducible computation is `q32_q1_binomial_filter_lattice.py`; Q605 is
auditing whether the apparent positive rate can be proved or escaped.

## 80. Q589 closes recurrence-only close-pair amplification

For a fixed `n`, separate the q=1 hits into the direct and reflected index
sets

`D_n={j: p=n-j is prime, p>2j, p|A_j}`,

`R_n={j: p=(n+1+j)/2 is prime, p>2j, p|A_j}`.

If one branch has `K` hits in `[0,J]` and that interval is partitioned into
blocks of length `H`, the number of same-block pairs is at least

`K^2 H/(2J)-K/2`.

Thus a bound of `H n^(o(1))` close pairs for every `H<=n^(1/3)` would give
`K=n^(2/3+o(1))` and solve the weighted q=1 slice.  This is a genuinely
weaker target than bounding the whole radical directly.

The bare order-two Apéry recurrence does not supply that bound.  Two direct
hits at `j,j+h` live modulo `p=n-j` and `q=p-h`; two reflected hits live
modulo primes differing by `h/2`.  The continuant transfer identity is a
two-zero eliminant only when both endpoint zeros lie in the same field.
Across the two moving primes it instead contains the unknown cross-quotients
`A_j/p` and `A_(j+h)/q`.  Over `F_p x F_q`, the two zero conditions are
independent state lines, so elimination of a general recurrence state gives
the zero ideal.

This is a proof-grade recurrence-level obstruction: by choosing a state line
separately modulo every prescribed prime and applying CRT to the two initial
coordinates, an integer solution of the same cleared recurrence can realize
any finite collection of prescribed moving zeros, even modulo prime squares.
Any successful pair theorem must therefore use the distinguished Apéry
initial direction or its hypergeometric/diagonal arithmetic.

Cross-weighting the two endpoints does give a genuine integer divisible by
`pq`, but its size is governed by the full bulk term `A_(j+h)`, not by the
gap `h`.  The companion Casoratian has the desired low-gap height, but
cross-weighting it to recover `pq` destroys the cancellation; the first
lattice point on the cancellation line contains `pq` universally and loses
all selectivity after division.  This closes ordinary continuants,
Casoratians, Sylvester resultants, transfer matrices, and recurrence-only
prime-power amplification for the close-pair route.

Q607 now asks whether the distinguished binomial/diagonal realization gives
the additional cross-prime constraint that the recurrence module cannot.

## 81. The full coefficient family removes only small carrier nuisance

An exact comparison through `n=400` tested the two q=1 branch ideals

`gcd(Gamma_n,B_n^-)`, `gcd(Gamma_n,B_n^+)`

against the simpler Euclidean ideals

`gcd(S_n,B_n^-)`, `gcd(S_n,B_n^+)`.

They differ in 185 of the 796 branch cases.  The first examples show quotients
`5` and later `41`; for example at `n=165` on the direct branch the full
coefficient ideal is `27935`, while the simple ideal is `1145335=41*27935`.
This confirms that the extra coefficient directions genuinely remove
unrelated carrier factors.  They do not remove any top-half bad prime, by the
exact branch-support theorem in Section 76.

Consequently an adjoint reduction of the Franel input is still worth testing:
it may explain and uniformly bound the small nuisance quotient.  But a proof
that only removes those nuisance primes is not the target theorem; it must
also bound the remaining aligned top-half radical.  Q606 is assigned this
exact adjoint/telescoping audit.

## 82. Local Newton block divisors package close hits exactly, but remain hard

There is a clean interpolation version of the pair-certificate target.  Fix
a block `j=u,...,u+H` and let

`P_m(x)=sum_(k=0)^m Delta^k A_u binom(x,k)`

be the degree-`m` integer-valued Newton interpolant through
`A_u,...,A_(u+m)`.

For the direct branch evaluate at `x=n-u`.  If `j=u+h` is a direct hit, then
`p=n-u-h` and, for every `H<=m<p`,

`P_m(n-u)=P_m(h)=A_(u+h)=0 (mod p)`.

For the reflected branch evaluate at `x=-(n+1+u)`.  A reflected hit satisfies
`2p=n+1+u+h`, hence the evaluation point is congruent to `h` modulo `p`, and
the same conclusion holds for `m<p`.

Therefore every hit prime in the block divides the exact local divisor

`Delta_(n,u,H)^branch=gcd_(H<=m<p_min) P_m(x_branch)`.

This succeeds where the bare two-prime transfer matrix did not: all moving
prime conditions are embedded into divisibility of the same global integers.
The denominators in `binom(x,k)` are p-units precisely because the degree is
stopped below the least candidate prime.

Exact computations in `q32_local_newton_blocks.py` show strikingly clean
localization.  At `n=321`, `H=40`, the three reflected blocks have divisors
`179,193,211`, exactly the three reflected hits.  With shorter blocks the
only observed nuisance is usually the fixed factor 5; for example at
`H=10` the nontrivial reflected divisors are
`5*179,193,5,211`.  At `n=240,400` the only nontrivial sampled reflected
divisors are 5, and at `n=120,600` all sampled divisors are 1.  Direct blocks
in these samples are 1 or 5.

This is not yet a height proof.  The family uses degrees up to a positive
fraction of `n` and Apéry values reaching the bulk.  Its small gcd is another
determinantal divisor phenomenon like `Delta_n`; bounding it by
`exp(O(H polylog n))` would give exactly the missing close-pair theorem, but
generic determinant bounds are much larger.  Stopping after only `O(H)`
degrees retains the full exponential Apéry history.  Thus the construction
is a sharper target and a useful local Smith object, not an unconditional
advance until its smooth/extraneous factors can be bounded independently of
the hit radical.

## 83. Applying the Apéry recurrence to the local interpolant adds no depth

Let `P_m(x)` interpolate the consecutive values
`A_u,...,A_(u+m)`.  Apply the Apéry recurrence operator polynomially:

`R_m(x)=(x+1)^3 P_m(x+1)-P(x)P_m(x)+x^3P_m(x-1)`,

where `P(x)=34x^3+51x^2+27x+5`.

At every interior node `x=u+1,...,u+m-1`, the three interpolated values are
actual consecutive Apéry values, so `R_m(x)=0`.  Hence the exact factorization

`R_m(x)=prod_(h=1)^(m-1)(x-u-h) Q_(u,m)(x)`

holds with `deg Q_(u,m)<=4`.

This looked promising because the large interpolation degree collapses to a
quartic quotient.  At a block hit, however, the universal node product
already contains the prime: evaluating at the direct extrapolation point
gives the factor `n-(u+h)=p`, and the reflected congruence has the analogous
forced factor.  A useful depth amplification would require the quartic
quotient to vanish modulo the bad prime.

The exact audit `q32_newton_recurrence_residual.py` tests the three reflected
hits at `n=321`, using blocks `(u,m)=(0,40),(41,40),(82,24)`.  The quartic
quotient is a p-unit in all three cases; its residues modulo
`p=179,193,211` are respectively `97/103`, `158/121`, and `31/54`, all
nonzero.  Thus the polynomial recurrence residual has exactly the universal
first factor and does not distinguish bad nodes.  This closes the immediate
“low-degree recurrence quotient” rescue of the local Newton divisor.  Reopen
only if a different operator makes badness force an additional factor rather
than merely reproducing the interpolation-node carrier.

## 84. Multi-hit companion cancellation saturates at the target product

The two-point Casoratian obstruction leaves open a genuinely different
possibility: use many hit indices and Siegel's lemma to cancel the rational
Apéry companion while sharing the coefficient cost across the hits.

For hits `(j_i,p_i)`, put `P=prod_i p_i`, choose
`D=lcm(1,...,L)^3` with `L=max j_i`, and define the two integral rows

`c_i=(P/p_i) D a_(j_i)`,

`B_i=b_(j_i)/p_i`.

If `z in Z^K` satisfies `sum c_i z_i=0`, the companion coordinate cancels
exactly, while the Apéry coordinate is

`sum_i (P/p_i)b_(j_i)z_i=P sum_i B_i z_i`.

The possible residual coordinates form the Smith image ideal

`I=B(ker_Z c)`.

Equivalently, its generator is the gcd of the 2-by-2 minors of the matrix
with rows `c,B`, divided by `gcd_i c_i`.  Thus a multi-hit approximation
argument must control not merely the shortest vector in `ker c`, but the
shortest vector outside the codimension-one sublattice `ker c intersect
ker B`.

The exact audit `q32_multihit_companion_lattice.py` gives:

- `n=321`, reflected indices `36,64,100`: `I=Z`, so the remaining integer
  ideal is exactly `P=179*193*211`;
- `n=11576`, direct indices `139,2257,2683`: `I=85Z`, so the remaining
  ideal is `85P`.

This explains both the promise and the trap.  Many-point cancellation really
does preserve all hit primes and removes almost everything else, but it
produces the desired prime product itself, not an additional divisor.
Ordinary homogeneous Siegel lemma can return a very short vector lying in
the larger common kernel, for which the residual coordinate is zero.  The
first vector with nonzero residual is the last relevant successive minimum
and is enormous in the two triples.  For `K=3`, exact minimization along the
rank-one affine solution fiber gives coefficient bit lengths

`[572,716,534]` and `[13690,24469,22302]`,

respectively.  These costs overwhelm the Archimedean Apéry error and agree
with the earlier critical-adelic-ledger obstruction.

This closes the naive claim that increasing the number of hits automatically
divides the companion-relation height by `K`.  A reopening would require a
new upper bound for the *last non-common-kernel successive minimum* of this
special Apéry lattice, not a standard Siegel bound for some nonzero kernel
vector.  Such a bound strong enough to beat the `D`-denominator tax would be
a substantive new place-mixing theorem.

## 85. Q599 completes the safe q=3 affine-lattice audit

For sign `eps` and Taylor moments

`D_k^eps=sum_d eps^d binom(d,k) C_d`,

the degree-`m` affine filter values are exactly

`D_0^eps+g_m^eps Z`, where
`g_m^eps=gcd(D_1^eps,...,D_m^eps)`.

Thus the least nonzero value is the centered residue of `D_0^eps` modulo
`|g_m^eps|`.  Since the gcd chain decreases, these affine lattices enlarge
with `m` and their least nonzero values cannot increase.  The apparent
interior q=3 optimum at `n=400,m=97` was only the first point of a plateau.
An independent exact rerun gives, for the minus sign,

`log(mu_m)/400=0.5856101401068076` for every `97<=m<=100`,

and for the plus sign,

`log(mu_m)/400=0.591296728393659` for every `97<=m<=100`.

The original script records the lexicographically first minimizer and is
therefore correct.

The local evaluation matrix

`(binom(ap,k))_(0<=a<=r,0<=k<=m)`, for `r<=m<p`,

has nonzero Smith invariants `1,p,...,p^r`.  In particular the q=3 profile
has invariants `1,p,p^2,p^3`.  At degree `p` there is an exact rank jump:
`binom(ap,p)=a (mod p)`, so a good candidate can escape and selectivity is
lost.  Globally, the consecutive binomial-basis minor is unimodular, so
there is no generic covolume obstruction and no lattice-only proof of a
positive exponential rate.

Q599 proves neither a positive lower rate nor a subexponential sequence.
The remaining quantity is the Apéry-specific centered Hasse-jet residue.
Its statement that the radical is "exact" must be read after restriction to
the candidate interval; nuisance primes outside that interval can occur.

## 86. Q605 finds a much smaller branch-filtered safe certificate

The unrestricted q=1 safe affine optimum remains exponential through
`n=800`; Q605 reports rates between approximately `0.765` and `0.915` on
`100<=n<=800`.  The important new refinement is to intersect separate safe
affine residues with the two exact branch carriers

`B_n^-=binom(n,K)`, `B_n^+=binom(n+K,K)`,
`K=floor((n-1)/3)+1`.

For a branch, let `p_min` be the least top-half prime dividing its carrier,
take `m=p_min-1`, let `mu_m^eps` be the exact centered affine residue, and
put

`H_branch^eps=gcd(B_n^branch,mu_m^eps)`.

Every branch bad prime divides `H_branch^eps`, while every top-half carrier
prime in it is bad; the remaining factors are small/medium-prime powers.
An independent exact reconstruction of the definition gives the best
per-branch data

- `n=400`: direct `H=3`, reflected `H=3`, product rate `0.0054930614`;
- `n=600`: direct `H=29`, reflected `H=1`, product rate `0.0056121597`;
- `n=800`: direct `H=85`, reflected has rate `0.0102778563`, total product
  rate `0.0158311704`.

These rates are dramatically below the unfiltered affine residue.  They are
finite data, not an asymptotic proof, and minimizing the two branches with
different signs is legitimate only because the two integers are used as
separate certificates.  The reproducible audit is
`q32_branch_filtered_moments.py`.

After reconstructing the data, this is **not** the strongest concrete route.
The full coefficient content `Gamma_n` divides every Taylor moment and hence
every affine value, including `mu_m^eps`.  Therefore

`gcd(Gamma_n,B_n^branch) | H_branch^eps`.

The Section 76 full-content branch ideal always dominates the Q605
certificate.  The domination can be strict:

- `n=400`: both full-content branch ideals are 1, while both Q605 ideals
  are 3;
- `n=600`: both full-content branch ideals are 1, while the Q605 direct
  ideal is 29;
- `n=800`: the full-content ideals are `85` and `1241`, while the Q605
  ideals are `85` and `3723=3*1241`.

Thus Q605 is useful evidence that degree-growing safe filters can nearly
reach the full ideal, but it does not reopen a narrower theorem.  The
stronger target remains

`log gcd(Gamma_n,B_n^-)+log gcd(Gamma_n,B_n^+)=o(n)`.

Q619 and Q620 were already dispatched before this domination was noticed;
their answers must be audited against it and only genuinely new
full-content information retained.

Operational note: invoking `ask-gpt.py --help` does not print help; it
submitted the literal low-value question `--help` as Q616 on dm2.  This was
an operator error.  The poller was stopped locally, but the remote task may
still complete.  Future invocations must read the script header/source for
usage rather than pass conventional help flags.

## 87. Full-content Lucas self-similarity is false beyond the top half

The small full-content branch ideals suggested that the q=1 truncation
polynomial might package the complete lower channel through a digit recursion

`T_n(c) = A_j T_q(c^p) (mod p)`,

for `n=qp+r`, `j=min(r,p-1-r)`, and `p>sqrt(n)`.  If true, this would have
classified every large prime of `Gamma_n` by one lower Apéry zero plus a
small quotient content.

An exact coefficient search disproves it immediately:

- `n=16`, `p=5`, `q=3`, `r=j=1`, and `A_1=5=0 (mod 5)`;
- nevertheless `T_16(c)` has nonzero coefficients modulo 5, beginning
  `C_0=4`, `C_1=3`, `C_5=3`, `C_6=1`.

The failure mechanism is precise.  The q=1 cutoff
`J=floor((n-1)/3)` crosses a p-block as soon as `q>=3`.  In
`Q_n(t)=Q_q(t^p)Q_r(t) (mod p)`, the boundary block is only partially
included, so the high-digit factor cannot be pulled out as a scalar
polynomial.  Even the weaker product guess with the low truncation fails;
for `n=22,p=7`, the coefficient at degree 7 contradicts it.

For `q=2` the cutoff remains below `p` and is complete in the low digit, so
one does have

`T_n(c)=A_r Q_2(c^p) (mod p)`.

Consequently every q=2 Apéry zero can appear as a non-top-half nuisance
prime of the q=1 full-content branch ideal.  Example: from
`179|A_36`, the index `n=2*179+36=394` supplies such a q=2 factor.
Thus bounding the nuisance is already another moving-zero problem, not a
pure small-prime estimate.

Small exact quotient tables also show that q>=3 support depends on the
partial boundary block, not only on `p|A_r`.  For `p=5`, whose zero
positions are `r=1,3`, `p|Gamma_(qp+r)` occurs for `q=1,2,4` but not
`q=3`.  For `p=11,r=5`, it occurs for every `1<=q<11` except `q=3,6`.
For `p=17,r=3`, it fails at `q=3,6,9,12` but holds at the other tested
quotients, while the reflected zero `r=13` holds for every quotient.

This closes any proof that treats old primes in the branch ideal as generic
small nuisance or assumes a one-scalar Lucas profile.  Q624 is assigned the
correct multi-block profile problem.  Reopen the recursive route only after
the partial boundary polynomial is written explicitly and at least two of
its coefficient equations are used.

## 88. The q=2 nuisance profile is exact, and smooth stripping is tautological

There is one quotient beyond the top half for which the full-content profile
does factor cleanly.  Write `n=2p+r`, `0<=r<p`, and retain the q=1 cutoff
`J=floor((n-1)/3)`.  Since

`r<=J<p`,

Lucas gives

`Q_n(t)=Q_2(t^p)Q_r(t) (mod p)`.

For every `m<=J<p`, the `y^p` terms in
`Q_2((c+y)^p)` cannot contribute to the coefficient of `y^m`.  Hence

`K_c(n,m)=Q_2(c^p)K_c(r,m) (mod p)`.

The cutoff contains the complete lower transform because `J>=r`, so summing
over `m` proves

`T_n(c)=A_r Q_2(c^p)
       =A_r(1+6c^p+6c^(2p)) (mod p)`.

For `p>=5`, therefore,

`p|Gamma_n <=> p|A_(n-2p)`

throughout `n/3<p<=n/2`.

The full coefficient profile, including every zero coefficient away from
degrees `0,p,2p`, was independently verified for all 316 q=2 pairs with
`n<=120`.

Kummer also shows that every such q=2 prime divides the direct carrier
`B_n^-=binom(n,K)`, `K=J+1`; the reflected carrier retains the subrange
with `r+K>=p`.  Thus these are genuine moving-zero factors of the
full-content branch ideals.  The check

`n=394=2*179+36`,
`Gamma_394=15215=5*17*179`,

confirms the theorem and explains the factor 179.

One can remove this and every other old prime exactly, but doing so is
tautological.  Let `M=floor(n/2)` and choose `L>log_2(2n)`.  For the direct
carrier define

`B_top^-=B_n^-/gcd(B_n^-,(M!)^L)`.

The reflected carrier `B_n^+=binom(n+K,K)` has a second pollution interval
`n<p<=n+K`.  This is real: at `n=717`, its common factors with `Gamma_n`
include `751` and `821`, both larger than `n`.  Put

`E_n=lcm(1,...,n+K)/lcm(1,...,n)`

and define

`B_top^+=B_n^+/gcd(B_n^+,(M! E_n)^L)`.

The elementary bound `v_p(B)<=log_p(2n)` shows that these purified carriers
have prime support only in `(M,n]`, while every target carrier prime is
unchanged.  For primes above `M`, no higher prime power crosses between the
two lcm cutoffs, so `E_n` removes precisely the above-`n` pollution relevant
to this range.  Consequently

`gcd(Gamma_n,B_top^-)`, `gcd(Gamma_n,B_top^+)`

have exactly the desired direct/reflected supports and no q>=2 or above-`n`
nuisance.  But these gcds are just another exact presentation of the unknown
target radical; the stripping supplies no upper bound for their height.  It
corrects the object but does not advance the asymptotic theorem.

## 89. Q617 proves the transverse companion lattice has exactly zero slack

Q617 confirms the Section 84 Smith formula in arbitrary dimension.  For
integer rows `c,B in Z^K`, with `g_1=gcd_i c_i` and independent rows, put

`Delta_2=gcd_(i<k)(c_i B_k-c_k B_i)`.

Then

`B(ker_Z c)=(Delta_2/g_1) Z`.

Thus the `n=321` and `n=11576` image steps are exactly 1 and 85.  Duplicate
folded indices only add common-kernel directions and can be removed; within
one q=1 branch the relevant folded indices are distinct.

More decisively, if `Lambda_1=ker_Z c` and
`Lambda_2=ker_Z c intersect ker_Z B`, the quotient between their canonical
Bombieri--Vaaler height budgets is exactly the spacing between consecutive
affine layers of `B` on `Lambda_1`.  The exact Apéry error identity gives the
same spacing.  Hence

`one-row Siegel budget
 = common-kernel budget * forced transverse cost`.

There is no missing exponential factor.  Standard Bombieri--Vaaler,
inhomogeneous Siegel, covering-radius, or determinant arguments can return
`K-2` short vectors entirely inside the common recurrence kernel.  A sharp
model with unit common-kernel vectors and arbitrarily large first transverse
vector proves this is not a technical weakness of the estimates.

The Schmidt Subspace Theorem also does not directly apply to the companion
linear form: its decisive coefficient is `zeta(3)`, whose algebraicity is
unknown, and both the forms and dimensions vary with `n`.

The exact reopen condition is now a special *transverse well-roundedness*
theorem for Apéry hit rows: a polynomial-factor bound for the first vector
in `Lambda_1\\Lambda_2`, or equivalently a proof that the geometric-mean
shell contains a relation nonzero in `Lambda_1/Lambda_2`.  Q617 shows such a
bound would give a square-root bound for the number of same-branch hits and
therefore close q=1.  It is false for arbitrary integer recurrence rows, so
it must use arithmetic beyond the order-two recurrence and its Casoratians.

## 90. Q607 isolates the cross-characteristic Mellin theorem, but no proof

Q607 audits the distinguished hypergeometric, toric, K3, and modular
realizations.  For `p>2j`, every summand in the terminating Apéry
`4F3` is a p-adic unit.  Thus `p|A_j` is cancellation among units; Kummer
valuation and Gross--Koblitz slope arguments have no termwise zero to exploit.

For the fixed Laurent polynomial `Lambda`, `A_j mod p` is an exact complete
Kummer character sum on a fixed three-dimensional torus.  In a fixed q=1
row, the selected character can be written uniformly as

`chi_p=omega_p^(n-1)`.

This is bounded-complexity geometry locally at each `p`, but the character
order, trace field, selected prime over `p`, and residue characteristic all
move.  Fixed-family Chebotarev, character-aspect equidistribution at one
prime, and ordinary Sato--Tate therefore have the wrong quantifiers.

The fixed-`n` coefficient-shell identity in Q607 is the same structure
already sharpened in Section 35/Q476.  Q607 initially counts the ambient 27
directions in `{-1,0,1}^3`; the exact support of `Lambda` has only 22, so the
accepted version is the 22-shell formula.  Q476's termwise Frobenius collapse
still applies and shows that the entire shell sum is one scalar vanishing
equation, not 22 independent congruences.  Raw shell heights are
`exp(Theta(n))`.

The modular parametrization expresses `A_j` as a constant term of a
weight-four weakly holomorphic modular form, but its pole order grows with
`j`.  It is not a coefficient of one fixed eigenform.  Clearing the pole
restores linear height.  The central Ahlgren--Ono slice fixes at most one
prime per row and costs only `O(log n)`.

Small character order is again unconditionally negligible: if
`d_p=(p-1)/gcd(p-1,n-1)<=D`, then `p-1=d_p g` with `g|n-1`, giving only
`D n^(o(1))` candidates.  The hard characters may therefore be assumed
almost faithful.  This does not amplify norms because `p` splits completely
in the corresponding cyclotomic trace field.

The weakest new theorem identified by Q607 is an aggregate short-gap
horizontal Mellin non-clustering estimate for same-branch double hits.  With
the Q589 block-packing inequality it would give
`K(n)=n^(2/3+o(1))`.  No cited hypergeometric, Chebotarev, Sato--Tate, or
large-sieve theorem proves it.  Q634 now attacks the exact 22-shell
bivariate holonomic/resultant formulation; any eliminant of
`exp(Omega(n))` height is a failure.

## 91. The evident close-shell pair certificate has full bulk height

The shell formula can be evaluated exactly without expanding the Laurent
power.  If

`C_N(e)=[x^e1 y^e2 z^e3] Lambda^N`,

then

`C_N(e)=sum_k binom(N,k)binom(N,k-e1)
                  binom(2N-k,N-e2)binom(2N-k,N-e3)`.

Put

`E_N(m)=sum_(v in {-1,0,1}^3) C_N(mv)`.

In the top-half range the five directions outside `supp Lambda` vanish
identically, leaving the exact 22-shell sum.  At a q=1 hit `p`,

`p|E_(n-1)(p-1)`.

For two hits `p,q`, the elementary cross-weighted value

`q E_(n-1)(p-1)-p E_(n-1)(q-1)`

is divisible by `pq`.  Exact evaluation on the `n=321` reflected triple
gives shell rates `3.554616`, `3.522692`, and `3.482928` for
`p=179,193,211`.  The pair-certificate rates for gaps `14,18,32` are
respectively `3.571010`, `3.539365`, and `3.571288`.

Thus the first distinguished cross-characteristic certificate is actually
slightly larger than its individual shell values; close spacing produces no
Archimedean gain.  The reproducible audit is
`q32_shell_pair_certificate.py`.  A useful result from Q634 must normalize
the bivariate recurrence deeply enough to replace this bulk rate by
`O(h polylog n)` height.

## 92. Experimental ideal inclusion: the full content divides the Apéry term

An exact scan for every `3<=n<=400` found

`Gamma_n | A_n`.

The same divisibility holds in the separately computed larger samples
`n=600,717,800`.  At `n=717`,

`Gamma_717=5^2*17*443*751*821`,

so the observation includes the above-`n` pollution primes 751 and 821; it
is not merely the known top-half support theorem.  In the tested samples
`Gamma_n` also divides the folded Newton gcd

`gcd(A_n,F_J(n)F_J(-n-1))`.

This is currently an experimental ideal inclusion, not yet a proved
identity.  A local explanation is plausible: coefficientwise vanishing of
the truncated Legendre--Euler polynomial is stronger than the single Lucas
condition `A_n=0`; the converse already fails at `n=16,p=5`.  A global proof
would require showing that `A_n` lies in the integer coefficient ideal, or
using the Franel adjoint recurrence to express the omitted tail in that
ideal.  Q606 is directly relevant.  If proved, the inclusion unifies the
content and Newton gcd packages, but by itself does not give an
`exp(o(n))` upper bound.

The observation is not specific to the q=1 cutoff in the tested small range.
Writing `T_(n,J)(c)` for the same shifted transform truncated at an arbitrary
`0<=J<=n`, exact polynomial construction verifies

`content_c T_(n,J)(c) | A_n`

for every pair `(n,J)` with `n<=35`.  The analogous statement is false for
arbitrary input vectors in the two inverse binomial bases, so any proof must
use the special Legendre kernel and Franel input rather than formal
translation invariance alone.

## 93. Q619 removes centered-residue overshoot but returns to moving zeros

For sign `eps` and branch carrier `B`, Q619 replaces the centered certificate
by the strictly smaller core

`H_core(n,m,eps)=gcd(B,D_0^eps,D_1^eps,...,D_m^eps)
                =gcd(B,D_0^eps,g_m^eps)`.

Every common divisor of `B,D_0,g_m` divides the centered representative, so
`H_core|H_centered`.  The core has the same bad top-prime support and is
monotone decreasing in `m`; the centered certificate is not monotone and can
acquire accidental factors from the nearest-residue quotient.  Exact
examples include:

- `n=400`, reflected minus: centered 151, core 1;
- `n=600`, direct centered values 101 and 29, both cores 1;
- `n=800`, centered values contain extra factors `245` or `7`.

The full coefficient branch ideal still dominates this refinement:

`gcd(Gamma_n,B) | H_core | H_centered`.

Thus the logical hierarchy found locally is confirmed; Q619's core is the
right safe-moment object but not stronger than the full content.

Q619 also gives an exact prime-power interpretation:

`ell^a|H_core`

iff the signed truncation polynomial has a zero of Taylor multiplicity
greater than `m` at the chosen center modulo `ell^a`.  Small prime powers
up to `Y=n^(2/3)` have total carrier weight `O(n^(2/3))`; above `Y` all
carrier valuations are one.

For a large carrier prime `ell`, write `n=q ell+r` and decompose both the
cutoff and the moment order into base-`ell` blocks.  Lucas and the Franel
Lucas congruence factor the truncation into:

1. the genuine lower-digit factor `A_r`;
2. a smaller quotient Taylor polynomial.

The quotient-residual alternative for fixed quotient data divides a nonzero
integer of height `exp(O(q))`.  Summed over `q<n^(1/3)`, these residual
factors cost `O(n^(2/3))`.  What is **not** removed is the genuine
lower-quotient condition `ell|A_(n mod ell)`.  The factor 73 at `n=800` is
exactly such an old-slice zero.

Therefore the core admits a clean

`moving Apéry-zero support + O(n^(2/3))`

description, but after separating the top-half q=1 radical it still contains
lower-quotient moving zeros.  It is an equivalent-plus-contamination
reformulation, not a proof and not a theorem logically weaker than the
horizontal target.  Q638 now tests the stronger experimental inclusion
`Gamma_n|A_n` and the exact ideal relations.

## 94. Stronger tail-ideal conjecture behind `Gamma_n|A_n`

For an arbitrary cutoff write

`T_(n,J)(c)=sum_d C_d c^d`,

`C_0=sum_(k<=J)L(n,k)F_k`,

so

`A_n-C_0=sum_(k>J)L(n,k)F_k`.

Let

`H_(n,J)=gcd_(d>0) C_d`.

Exact construction of every polynomial `T_(n,J)` verifies the stronger
divisibility

`H_(n,J) | A_n-C_0`

for all `(n,J)` with `n<=60`.  This immediately implies

`content(T_(n,J))=gcd(C_0,H_(n,J)) | A_n`.

The statement is not termwise: already at `(n,J)=(6,1)`,
`H=210` does not divide the omitted `k=5` summand (remainder 84).  The
divisibility is created only after summing the whole Franel tail, pointing
directly to a formal-adjoint/telescoping identity rather than Kummer
valuation of individual terms.

The exact Lucas block decomposition for `J=A p+B` is

`T_(n,J)(c)
 = A_r T_(q,A-1)(c^p)
   + K_(c^p)(q,A) g_A^(c^p) T_(r,B)(c)  (mod p)`,

with `n=qp+r`; empty prefixes are interpreted as zero.  This formula explains
the `n=16,p=5` boundary-block counterexample and may support an induction for
the radical version, but it does not yet prove the integer prime-power ideal
inclusion.  The missing global step is an explicit Franel-adjoint expression
of the full omitted tail in the nonconstant coefficient ideal.

The broad Q606 Franel-adjoint task failed at the bridge without an answer.
It was replaced on dm3 by Q643, which asks only for this exact tail-ideal
identity and its saturation denominators.

## 95. The order-two rational adjoint is exactly unsaturated at the target

Because the Franel recurrence has order two, regard its two initial values
as independent formal variables.  For fixed `(n,J)`, the omitted tail

`R_(n,J)=sum_(k>J) binom(n,k)binom(n+k,k)F_k`

and every nonconstant coefficient `C_d` are rational linear forms on this
two-dimensional solution space.  Whenever `C_1,C_2` are independent there
are therefore unique rational numbers `x,y` such that

`R_(n,J)=x C_1+y C_2`

for every formal solution of the Franel recurrence.

Exact `Fraction` arithmetic for the q=1 cutoff gives:

| `n` | `Gamma_n` | common denominator of `x,y` | gcd |
|---:|---:|---:|---:|
| 20 | 17 | 2898800997129995958000 | 17 |
| 30 | 85 | 67654179405461703371412117122100000 | 85 |
| 40 | 5 | 17652791220399897506969952735758163573360000 | 5 |
| 60 | 11 | 290143197036476511651301615440641748780830243137544853110204448792000000 | 11 |
| 80 | 1 | 4690739093463869359851560547541227687761967296376413360960313710436794478218487310663157559096000 | 1 |
| 100 | 1 | 7164219947249817222359995293087967137516077298241604527238549848953517594164223285321142115673542340048930834078550812500 | 1 |

Thus, in every tested case with nontrivial content, the two-boundary
rational adjoint denominator contains the full observed target content.
Clearing denominators loses exactly the desired divisibility rather than
proving it.  This closes the naive order-two reduction: the required theorem
is an **integral saturation** statement involving the whole coefficient
lattice, not merely a rational recurrence identity using two coefficients.

The stronger tail divisibility is also not termwise.  For example, at
`(n,J)=(6,1)`, the nonconstant coefficient gcd is `210`, whereas the gcd of
the individual omitted summands is only `84`.  Further samples have the same
failure.  Any proof must telescope or saturate the *summed* tail; it cannot
show that each omitted hypergeometric term lies in the coefficient ideal.

## 96. Q632: no existing holonomic gcd theorem reaches the target

Q632 audited the known Subspace-Theorem gcd literature against the exact
Apéry, Newton, Strehl-prefix, and branch-carrier sequences.  The usable
theorems require an exact finite power sum/Binet representation, or points in
a fixed finite-rank multiplicative group/fixed `S`-unit set.  Passing to the
three residue classes modulo 3 leaves genuine P-recursive/rational-diagonal
sequences with negative powers in their saddle asymptotics; it does not
produce constant-coefficient recurrences.

The favorable multiplicative independence of the dominant algebraic bases
therefore checks only one hypothesis.  It does not supply the missing exact
finite-torus representation.  Standard G-function denominator theorems are
also irrelevant because the sequences here are already integral and the
problem is a numerator gcd.

For the coefficient content there is an additional mismatch: `Gamma_n` is
the gcd of a coefficient family whose dimension grows linearly with `n`, not
the gcd of two fixed recurrence sequences.  A full `exp(o(n))` gcd theorem
would be sufficient despite support pollution, but no published theorem
applies.

The weakest reusable missing theorem is consequently a carrier-saturated
holonomic gcd theorem for one fixed arithmetic diagonal and one fixed
factorial-ratio sequence with multiplicatively independent critical values.
Proving such a theorem would be substantially new; merely computing Ore
operators or their right gcds does not imply a value-gcd bound.

## 97. Q634 closes the generic close-shell resultant

Q634 corrects the shell carrier used in the first formulation.  The q=1 hit
carrier is the **unweighted** sum over the 22 nonzero Laurent directions.
Weighting again by the Laurent coefficients is false: already at the hit
`(p,r)=(17,3)` the weighted low shell is `3977=16 (mod 17)`, while the
unweighted shell vanishes.

For the corrected shell, every direction is a one-fold proper
hypergeometric sum and has an explicit bivariate Horn/Ore system.  At two
hits `p=m+1` and `q=m-h+1`, however, the exact universal ideal is

`(p,X) intersect (q,Y)=(pq,pY,qX,XY)`.

After the universal factor `pq` is removed, two independent quotient
coordinates remain.  Ore recurrences relate neighboring shell values but do
not increase the codimension of either hit condition.  Positive or
fixed-sign nonuniversal certificates retain full Apéry-scale exponential
height.  Thus generic scalar resultants and cross-weighted determinants
cannot yield the desired close-pair saving.

Reopen this route only with a genuinely arithmetic relation between the two
quotient residues, a single global fixed-gap resultant of sublinear height,
or the horizontal fixed-gap pair-count theorem itself.

## 98. A much stronger experimental Franel tail-lattice lemma

The tail-ideal experiment suggests an integral statement independent of the
outer Legendre weights.  Fix `J` and define, for `k>J`,

`a_(k,d)=(-1)^(J-i) binom(k,i)binom(k-i-1,J-i)F_i`,

where `i=k-d` lies in `[0,J]`, and put the entry equal to zero otherwise.
Exact Smith computations first showed that

`(F_(J+1),...,F_N) in A_(J,N) Z^N`.

There is a particularly sharp candidate solution.  Set
`z_1=...=z_J=0` and recursively require

`F_k=sum_(i=0)^J (-1)^(J-i) binom(k,i)`
`                 * binom(k-i-1,J-i) F_i z_(k-i)`.

The coefficient of `z_k` is
`(-1)^J binom(k-1,J)`, so this rational recursion has no formal reason to
remain integral.  Nevertheless exact arithmetic verifies every division and
identity for

`0<=J<k<=250`.

Examples are

- `J=1`: `z_2,z_3,...=-10,-58,-270,-1238,-6008,...`;
- `J=2`: `z_3,z_4,...=56,414,1512,2566,-2196,...`.

The audit is implemented in
`projects/Q-series-and-Chan-s-work/q32_franel_tail_lattice.py`.

If integrality holds for all `J,k`, then for **every** integer outer weight
vector `(q_k)_(k>J)`,

`sum_(k>J)q_k F_k`

is an integer linear combination of the nonconstant coefficients

`sum_(k>J)q_k a_(k,d)`.

This immediately proves the arbitrary-cutoff tail-ideal theorem and hence
`content(T_(n,J))|A_n`.  It also explains why the two-coordinate rational
adjoint in Section 95 sees target primes in its denominator: the saturation
uses the entire growing family of degree coordinates.

The exponential-gcd target would still require a height or prime-support
bound after this structural lemma, so this is not yet the final
`exp(o(n))` estimate.  Q645, Q648, and Q649 now seek respectively a
constructive proof, a generating-function/combinatorial formula, and a
prime-power Smith proof of the integrality.

## 99. EGF transvectant form and experiments narrow the integrality source

Let

`\widehat F(x)=sum_(k>=0)F_k x^k/k!`,

`\widehat Z_J(x)=sum_(d>=1)z_d x^d/d!`.

Direct coefficient extraction rewrites the entire left-inverse recurrence as

`\widehat F(x)-sum_(k=0)^J F_k x^k/k!`

` = [u^J] \widehat F(xu) \widehat Z_J(x(1-u))/(1-u)`.

Indeed the coefficient of `u^J x^(i+d)/(i!d!)` on the right is

`(-1)^(J-i) binom(d-1,J-i) F_i z_d`,

and multiplication by `(i+d)!` supplies `binom(i+d,i)`.  Thus the candidate
integer inverse is an exact EGF transvectant/associated-convolution problem,
not an artifact of the Legendre outer weights.

Further exact tests sharpen what must be proved:

- the Franel solution is integral for every tested fixed
  `J in {1,2,3,4,5,7,10,16,25,40,64,100}` through `k=800`;
- replacing `F_k` by the generalized Franel sum
  `sum_a binom(k,a)^s` also passes for every `0<=s<=8` through `k=100`;
- arbitrary integer input sequences fail immediately;
- generic factorial ratios fail: `binom(3k,k)` fails already at
  `(J,k)=(1,3)`, and the Catalan sequence fails at the same point;
- a single exponential-polynomial sequence `P(k)c^k` passes in the tested
  cases, while nontrivial integer linear combinations of distinct
  exponentials generally fail.

The evidence therefore points to a special associated-binomial inversion for
the generalized Franel family, plausibly accessible from

`sum_a binom(k,a)^s`

or its product-of-binomial constant-term model.  It is not a generic
integrality property of diagonals, moments, P-recursive sequences, or
factorial ratios.

An attempted full recomputation of all q=1 branch contents through `n=400`
was stopped after several minutes: the current routine repeatedly reduces
huge coefficient sums modulo a still-large running gcd.  This was a
performance failure only, not mathematical evidence.  If the tail-lattice
lemma is proved, initializing the gcd with `gcd(C_0,A_n)` will be a valid
substantial optimization; using that initialization before proof would
silently assume the conjectured divisibility and is therefore not used in
the audit.

## 100. The tail-lattice lemma is now proved for cutoffs 1 and 2

Write `F(x)=sum_(k>=0)F_k x^k`.  Its exact differential equation is

`2(1+4x)F+(-1+14x+24x^2)F'`
` +x(1+x)(-1+8x)F''=0`.

For fixed `J`, the inverse recurrence has the ordinary generating-function
form

`L_J Z_J=F-sum_(k=0)^J F_k x^k`,

where, with `theta=x d/dx`,

`L_J=sum_(i=0)^J (-1)^(J-i)F_i x^i`
`    * binom(theta+i,i)binom(theta-1,J-i)`.

For `J=1`, exact reduction modulo the Franel differential equation gives

`Z_1=(A_1F+B_1F'+C_1)/(1-2x)^2`,

with

`A_1=(1+x)(1-8x)`,

`B_1=x(1+x)(1-8x)`,

`C_1=-1+3x-4x^2`.

All numerator series are integral and the denominator has constant term one,
so this proves every `z_k^(1)` is an integer.

For `J=2`, put `P_2=1-4x+10x^2`.  The exact identity is

`Z_2=(A_2F+B_2F'+C_2)/P_2^2`,

where

`A_2=(1+x)(8x-1)(4x^2+6x-1)`,

`B_2=x(1+x)(8x-1)(8x^2+14x-3)/2`,

`C_2=-100x^4+26x^3-9x^2+8x-1`.

Every positive-index Franel number is even (pair the summands; the middle
summand, when present, is also even), so `F'(x)/2` is integral.  Again
`P_2(0)=1`, proving integrality for all `z_k^(2)`.

The symbolic certificate
`projects/Q-series-and-Chan-s-work/q32_franel_tail_ogf.py`
checks both operator identities exactly, not by series guessing.

Series reconstruction for `J=3,4,5` gives the same shape

`Z_J=(A_JF+B_JF'+C_J)/P_J^2`,

`P_J(x)=sum_(i=0)^J(-1)^i binom(J,i)F_i x^i`,

with `P_J(0)=1`.  The observed scalar denominators of `A_J` are

`J!/2^floor(J/2)`

for `J=1,...,5`, and those of `B_J` are twice as large.  The numerator
combination is integral in every computed coefficient, but a uniform
divisibility proof for these growing factorial denominators is still
missing.  The next target is a general Darboux/contiguous formula for
`A_J,B_J,C_J` plus its coefficientwise saturation theorem.

## 101. The moving Jensen polynomial does not give a small resultant

The new denominator

`P_J(x)=sum_(i=0)^J(-1)^i binom(J,i)F_i x^i`

is a natural growing-degree algebraic-center candidate not covered by the
earlier bounded-degree scans.  Since `P_J(0)=1`, it is primitive, and for the
q=1 coefficient polynomial `T_n`,

`Gamma_n^J | Res_x(P_J,T_n)`.

If the resultant had logarithmic height `O(n log n)`, this would immediately
give `log Gamma_n=O(log n)`.  Exact resultants show the opposite.  For
`7<=n<=35`, `J=floor((n-1)/3)`, every resultant is nonzero, but

`log |Res(P_J,T_n)|/(nJ)`

increases from approximately `2.82` to `3.81`.  Representative values are:

| `n` | `J` | digits of resultant | normalized rate |
|---:|---:|---:|---:|
| 13 | 4 | 75 | 3.316606 |
| 20 | 6 | 184 | 3.518273 |
| 27 | 8 | 342 | 3.638968 |
| 35 | 11 | 638 | 3.814668 |

Thus the resultant has the generic `exp(Theta(nJ))` height and, after taking
the `J`-th root, gives only another positive-linear bound for
`log Gamma_n`.  The left-inverse denominator is structurally meaningful but
its roots are not low-height evaluation centers for `T_n`.

Reopen this specific algebraic-center route only if an exact factorization
removes a universal `exp(Theta(nJ))` factor while retaining the full
`Gamma_n^J` divisibility.  The raw resultant is closed.

## 102. The tail-left-inverse system is exactly an ordinary convolution

The banded recursion in Section 98 has a much simpler exact form.  The
binomial coefficient in its `(k,i)` entry satisfies

`binom(k,i)binom(k-i-1,J-i)`

` = k binom(k-1,J) binom(J,i)/(k-i)`.

Put

`P_J(x)=sum_(i=0)^J (-1)^i binom(J,i)F_i x^i`

and

`W_J(x)=sum_(d>J) z_d x^d/d`.

Then the coefficient of `x^k` in `P_JW_J` is

`sum_(i=0)^J (-1)^i binom(J,i)F_i z_(k-i)/(k-i)`.

Multiplying this identity by
`(-1)^J k binom(k-1,J)` shows that the entire tail-left-inverse
recursion is equivalent to the single formal-series identity

`P_J(x)W_J(x)`

` =(-1)^J sum_(k>J) F_k x^k/[k binom(k-1,J)]`.                 `(102.1)`

Because `P_J(0)=1`, this determines `W_J` uniquely over the rationals.  If
`Z_J=sum_(d>J)z_d x^d`, then

`Z_J=theta W_J`

` =(-1)^J theta((1/P_J)`
`   sum_(k>J)F_k x^k/[k binom(k-1,J)]).`                      `(102.2)`

This proves, for every cutoff, why the experimental closed form for `Z_J`
has denominator `P_J^2`: one factor is already present in `(102.1)`, and
the Euler derivative introduces the second.  It also removes the need to
guess a separate Darboux formula for each `J`.

The remaining integrality statement is now precise:

`k [x^k] ( (1/P_J)`
`  sum_(m>J)F_m x^m/[m binom(m-1,J)] ) in Z`                  `(102.3)`

for every `k>J`.  This is equivalent to the tail-lattice lemma, rather than
merely sufficient for it.  Generic integer sequences fail `(102.3)`, so
the unit constant term of `P_J` alone does not prove it; cancellation must
come from the Franel numerators.

The finite exact audit
`projects/Q-series-and-Chan-s-work/q32_franel_tail_convolution.py`
checks all signs, shifts, and rational divisions for `0<=J<k<=160`.

The literature check of Strehl's *Recurrences and Legendre Transform*
does not contain `(102.1)` or its required saturation.  Its inverse
Legendre matrix and Franel--Apéry operator conjugacy concern the standard
kernel `binom(n,k)binom(n+k,k)`, not the moving-cutoff kernel above.
Likewise the Jin--Dickinson--Zudilin inverse-Legendre divisibility theorem
controls a different very-well-poised transform.  These sources motivate
the binomial manipulation but do not close `(102.3)`.

## 103. Parameter differentiation isolates a stronger Franel integrality

Introduce the formal-parameter series

`Phi(alpha,x)=sum_(m>=0)(-1)^m binom(alpha,m)F_m x^m`

and put

`D_J(x)=partial_alpha Phi(alpha,x)|_(alpha=J)`.

For `m>J`,

`partial_alpha binom(alpha,m)|_(alpha=J)`

` =(-1)^(m-J-1)/[m binom(m-1,J)]`.

Consequently

`(D_J)_(>J)=(-1)^(J+1)`
`  sum_(m>J)F_m x^m/[m binom(m-1,J)]`.

Combining this with `(102.1)` gives the exact identity

`W_J=-D_J/P_J + D_J^(<=J)/P_J`,                         `(103.1)`

where `D_J^(<=J)` is the degree-`J` truncation.  In other words, the
tail-left-inverse problem is the difference of two parameter-logarithmic
derivatives:

`Z_J=-theta(D_J/P_J)+theta(D_J^(<=J)/P_J)`.              `(103.2)`

Exact computation reveals a stronger phenomenon:

`theta(D_J/P_J) in Z[[x]]`

and

`theta(D_J^(<=J)/P_J) in Z[[x]]`

separately for every tested `0<=J<k<=180`.  The first few coefficients of
the full logarithmic derivative for `J=1` are

`-2, 2, 34, 206, 1078, 5624, 30616,...`.

This is not a formal consequence of `P_J(0)=1`: random integer sequences
fail the full assertion, for example already at `(J,k)=(1,5)`.  Nor is it
a generic Laurent-period property: the constant terms of powers of
`x^-2+x^-1+2+x+x^2` fail the full assertion at `(J,k)=(2,5)`.

For Franel numbers there is a precise two-variable hypergeometric form:

`Phi(alpha,x)=sum_(r,s>=0)`
` (-alpha)_(r+s)(r+s)!^2 x^(r+s)/(r!^3 s!^3)`.           `(103.3)`

Thus the missing theorem can be attacked as an integrality theorem for a
parameter derivative of this particular `A`-hypergeometric series.  It is
reminiscent of mirror-map integrality but is not covered verbatim by the
standard statement `exp(G/F) in Z[[x]]`: direct computation shows
`exp(D_J/P_J)` develops unbounded prime denominators for `J>=1`, even
though its Euler logarithmic derivative is integral.

The audit
`projects/Q-series-and-Chan-s-work/q32_franel_parameter_derivative.py`
checks both separate integrality assertions.  A proof of just the first
one together with a proof of the finite-numerator assertion would close
the Franel tail-lattice lemma for every cutoff.

## 104. Correction: the reduced scalar denominator is an lcm, not a factorial

Section 100 reported the unreduced pattern
`J!/2^floor(J/2)` from the first five reconstructed formulas.  Extending
the exact reconstruction and reducing every polynomial coefficient shows
that this extrapolation was false.

In the representation

`theta(D_J/P_J)=(A_J F+B_J F'+C_J)/P_J^2`,

the exact reduced scalar denominators for `J=1,...,10` are

`den(A_J)=1,1,3,6,30,30,210,420,1260,1260`,

`den(B_J)=1,2,6,12,60,60,420,840,2520,2520`,

and `den(C_J)=1` throughout.  Thus, for `J>=2`,

`den(A_J)=lcm(1,...,J)/2`,

`den(B_J)=lcm(1,...,J)`.                                  `(104.1)`

This is exactly the denominator scale predicted by the harmonic
coefficients

`partial_alpha binom(alpha,m)|_(alpha=J)`

` =binom(J,m)(H_J-H_(J-m))` for `m<=J`,

not a factorial-scale saturation problem.  The earlier factorial claim
must not be used.

The parameter family also satisfies the exact third-order equation

`[theta^3`
` -x(theta-alpha)(7theta^2+7theta+2)`
` -8x^2(theta+1)(theta+1-alpha)(theta-alpha)]Phi=0`.       `(104.2)`

This follows directly by multiplying the Franel recurrence by
`(-1)^m binom(alpha,m)`.  At `alpha=J`, `P_J=Phi(J,x)` is its terminating
polynomial solution and `D_J` satisfies the parameter-differentiated
inhomogeneous equation.  Consequently

`Q_J=theta(D_J/P_J)`

is the logarithmic Frobenius quotient attached to the terminating
solution.  Its observed integrality is now a prime-power/lcm congruence
problem.  Standard mirror-map theorems do not apply verbatim because
`exp(D_J/P_J)` is not integral for `J>=1`; only its Euler derivative is.

There is also a universal contiguity relation

`(alpha+1-theta)Phi(alpha+1,x)=(alpha+1)Phi(alpha,x)`.      `(104.3)`

Differentiating `(104.3)` relates consecutive `D_J` and `P_J`.  This is a
possible induction mechanism, but the resulting coefficient equation has
the singular divisor `k-(J+1)` at degree `k`; integrality does not follow
formally without an additional Franel congruence.

## 105. Q664 closes positive multicutoff compression

Q664 gives a clean obstruction to the proposed dyadic prefix/block route.
For a q=1 hit with folded index `j`, Lucas gives

`S_K=0 mod p` for every `K>=j`.

For a block `S_b-S_(a-1)`, this forces divisibility only when `j<a`; a
non-prefix block does **not** isolate the middle band `a<=j<=b`.  More
generally, a prefix linear combination which vanishes universally for all
`j<=h` can use only cutoffs `K>=h`.  Thus the exact universal module is
generated by

`S_h,S_(h+1),...,S_J`

and, if the full sum is included, by the disjoint tail `A_n-S_J`.

The summand at `k=beta n` has exponential rate

`Phi(beta)=(1+beta)log(1+beta)-(1-beta)log(1-beta)`
`          -2 beta log beta+beta log 8`,

with

`Phi'(beta)=log(8(1-beta^2)/beta^2)>=log 64`

on `0<beta<=1/3`, and `Phi(1/3)=log 8`.  The height cost of every positive
block is therefore larger than the Chebyshev mass of the candidate-prime
interval it covers.  Deleting the block and charging those primes directly
never worsens the bound.  This closes prefixes, positive blocks,
overlapping/dyadic products, and nonnegative geometric-mean weights: their
optimum is the trivial `exp((1/2+o(1))n)` candidate-interval bound.

Arbitrary signed combinations evade this entropy comparison only by
computing the gcd/shortest vector of the universal module.  That is exactly
the original branch-purified holonomic gcd problem, not a new height
certificate.  Reopen multicutoffs only with a genuinely global signed
functional whose small height is proved independently of local Lucas
threshold geometry.

## 106. Congruence-preserving polynomials are the exact integrality target

For fixed coefficient degree `n`, define

`H_n(alpha)=[x^n] theta log Phi(alpha,x)`.

Then

`[x^n] theta(D_J/P_J)=H_n'(J)`.

Writing the integer-valued polynomial in its Newton basis,

`H_n(alpha)=sum_(k=0)^n c_(n,k) binom(alpha,k)`,

exact computation through `n=16` shows

`lcm(1,...,k) | c_(n,k)` for every `k`.                   `(106.1)`

Condition `(106.1)` is the standard Newton-basis criterion for
`H_n:Z->Z` to preserve all congruences:

`a-b | H_n(a)-H_n(b)`.

It is also the natural sufficient condition for `H_n'` to remain
integer-valued.  Thus the infinite tail-lattice divisibility is reduced to
one uniform statement:

**Franel congruence-preserving logarithmic derivative.**
Every coefficient polynomial of `theta log Phi(alpha,x)` preserves
congruences on the integer parameter.

The first Newton coefficient rows are

`n=3: 0,-8,72,-36`,

`n=4: 0,-16,216,-720,-72`,

`n=5: 0,-32,240,-4680,3600,1800`,

and each entry in column `k` has the predicted `lcm(1,...,k)` factor.

A tempting stronger explanation is false.  Although

`theta log(P_a/P_b) in (a-b) Z[[x]]`

is exactly the desired congruence, it does not follow because
`P_a/P_b` is an `(a-b)`-th power in `Z[[x]]`: already the formal square
root of `P_2=1-4x+10x^2` has coefficient `15/2` at `x^4`.
Any proof must work at the logarithmic-derivative level.

Q671--Q673 were dispatched with this exact reformulation, respectively
targeting the Newton/lcm proof, the Wronskian prime-power sum, and the
congruence-preserving polynomial theory.

The exact triangular audit
`projects/Q-series-and-Chan-s-work/q32_franel_congruence_preserving.py`
checks `(106.1)` for every `0<=k<=n<=140`, rather than only evaluating
derivatives at a finite set of integer parameters.

## 107. Frobenius-twisted Dwork congruence proves the full tail lattice

The missing prime-power mechanism is present in the coefficient-ring
version of the Samol--van Straten/Mellit--Vlasenko congruence.

Let

`Lambda(u,v)=(1+u)(1+v)(1+1/(uv))`,

so that `F_m=CT Lambda^m`, and put

`P_J(x)=CT_(u,v)(1-x Lambda(u,v))^J`

`      =sum_(m=0)^J (-1)^m binom(J,m)F_m x^m`.

The Newton polygon of `1-x Lambda`, considered only in the Laurent
variables `u,v`, has the origin as its unique interior lattice point.
Work over the coefficient ring

`R=Z_p[x]`, with Frobenius lift `phi(x)=x^p`.

The ghost-term proof of the D3 congruence works over any p-torsion-free
ring equipped with a Frobenius lift.  In this setting it gives, for
`r>=1` and `J,m>=0`,

`P_(J+m p^r)(x) P_floor(J/p)(x^p)`

` == P_J(x) P_(floor(J/p)+m p^(r-1))(x^p) mod p^r`.       `(107.1)`

Here is the coefficient-ring extension explicitly.  For a Laurent
polynomial `A` over `R`, replace the ordinary ghost terms by

`R_s(A)=A^(p^s)-phi(A)(u^p)^(p^(s-1))`

at the first level and iterate the same telescoping construction with
`phi`.  The two properties used in the original proof remain unchanged:

`R_s(A) in p^s R[u^+-1,v^+-1]`

and

`Newt_(u,v)(R_s(A)) subset p^s Newt_(u,v)(A)`.

The indecomposable-digit decomposition therefore gives twisted blocks

`c_(n1) phi^ell(n1)(c_(n2)) ...`

with `c_n in p^(ell(n)-1)R`.  In the proof of D3, the same good-partition
bijection preserves the absolute base-p positions of all blocks, hence
also their Frobenius powers.  Every unmatched partition is still killed
by the same sum of `ell(block)-1`, which is at least `r`.  This proves
`(107.1)` without specializing `x`.  It is important that `x` is in the
coefficient ring and that its Frobenius is retained; applying the scalar
version pointwise and then claiming coefficientwise divisibility would
not be valid.

Define

`H(J,x)=theta log P_J(x)`.

All four factors in `(107.1)` have constant term one.  Logarithmically
differentiating the congruence is therefore legitimate modulo `p^r` and
gives

`H(J+m p^r,x)-H(J,x)`

` == p[H(floor(J/p)+m p^(r-1),x^p)`
`       -H(floor(J/p),x^p)] mod p^r`.                    `(107.2)`

The extra factor `p` is exactly what was missing from the untwisted
attempt.  Induction on `r`, starting with the tautological modulus
`p^0`, yields

`H(J+m p^r,x)-H(J,x) in p^r Z_p[[x]]`.                  `(107.3)`

Consequently, for all nonnegative integers `a,b`,

`H(a,x)-H(b,x) in (a-b) Z[[x]]`.                        `(107.4)`

Thus every coefficient polynomial

`H_n(alpha)=[x^n]theta log Phi(alpha,x)`

is congruence-preserving.  This proves the experimental criterion
`(106.1)`.

There is a short p-adic proof that the full parameter derivative is
integral.  For every prime `p`,

`[H_n(J+p^r)-H_n(J)]/p^r in Z_p`.

As `r` tends to infinity, this difference quotient tends p-adically to
`H_n'(J)`.  Hence `H_n'(J)` belongs to every `Z_p`, and being rational it
is an integer.  Therefore

`theta(D_J/P_J) in Z[[x]]`.                              `(107.5)`

The finite-numerator half follows from the same congruence; no second
Dwork theorem is needed.  Define the degree-`J` interpolation

`T_J(alpha,x)=sum_(m=0)^J`
` (-1)^m binom(alpha,m)F_m x^m`.

For every integer `0<=k<=J`, `T_J(k,x)=P_k(x)`.  Put

`G_J(alpha,x)=theta(T_J(alpha,x)/P_J(x))`.

Then `G_J(J,x)=0` and

`partial_alpha G_J(J,x)=theta(D_J^(<=J)/P_J)`.

The endpoint derivative formula for a polynomial of degree at most `J`
is

`partial_alpha G_J(J,x)`

` =sum_(r=1)^J (-1)^r binom(J,r) G_J(J-r,x)/r`.          `(107.6)`

By `(107.4)`,

`theta log(P_(J-r)/P_J) in r Z[[x]]`.

Since `P_(J-r)/P_J` itself is an integral unit series,

`G_J(J-r,x)=theta(P_(J-r)/P_J) in r Z[[x]]`.

Every summand in `(107.6)` is integral, proving

`theta(D_J^(<=J)/P_J) in Z[[x]]`.                        `(107.7)`

Finally `(103.2)`, `(107.5)`, and `(107.7)` give

`Z_J in Z[[x]]`

for every cutoff `J`.  Equivalently, all coefficients in the exact
Franel tail-left-inverse recurrence are integers.  This closes the
tail-lattice conjecture unconditionally.

The audit
`projects/Q-series-and-Chan-s-work/q32_franel_twisted_dwork.py`
checks `(107.1)` and `(107.3)` coefficientwise for primes `2,3,5,7`,
all prime powers and parameters up to `90`, and degrees up to `90`.
The proof has also been extracted as
`projects/Q-series-and-Chan-s-work/Q32_FRANEL_TAIL_DWORK_PROOF.md`.

Primary sources checked:

- K. Samol and D. van Straten, *Dwork Congruences and Reflexive
  Polytopes*, Theorem 3.3 and its ghost-term proof.
- A. Mellit and M. Vlasenko, *Dwork's Congruences for the Constant
  Terms of Powers of a Laurent Polynomial*, especially the
  indecomposable-digit decomposition.  Their printed theorem uses
  coefficients in `Z_p`; `(107.1)` is the Frobenius-ring extension of
  the proof, not a verbatim quoted statement.

## 108. Q671 independently confirms the relative Dwork mechanism

Q671 independently found the same coefficient-ring theorem `(107.1)`,
the logarithmic induction `(107.2)`, and the congruence-preserving
conclusion `(107.4)`.  It also identified the Newton-basis theorem of
Cégielski--Grigorieff--Guessarian as an alternative route from integral
difference ratios to

`lcm(1,...,k) | Delta^k H_n(0)`.

This independently confirms the full-derivative half `(107.5)`.

One overclaim in Q671 must not be copied: it says that
`theta(D_J/P_J)` alone proves the tail-left-inverse integrality.  The exact
identity `(103.2)` also contains

`theta(D_J^(<=J)/P_J)`.

Q671 does not prove that finite-numerator term.  The missing half is supplied
by the endpoint interpolation argument `(107.6)`--`(107.7)`, which again
uses the stronger difference-ratio theorem.  Thus Q671 is a valuable
independent audit of the load-bearing Dwork step, but §107 remains the
complete tail proof.

Q684 was dispatched to test whether multiplying the contents over all
cutoffs and using row-subtracted coefficient determinants amplifies each
bad prime enough to turn the final height bound into `o(n^2)`.

## 109. Q645 gives a simpler canonical lattice vector (but not the zero-prefix one)

Q645 independently proves the tail **lattice membership** using only the
full logarithmic derivative.  If

`Q_d(alpha)=d[x^d]log Phi(alpha,x)`,

then the canonical choice

`z_d^can=-Q_d'(J)`

is integral by `(107.4)` and satisfies every matrix row

`F_k=sum_d A_(k,d) z_d^can`, `k>J`.

For example,

- `J=1`: `z^can=2,-2,-34,-206,...`;
- `J=2`: `z^can=2,-14,-88,-74,...`.

This is a different solution from the sharp triangular solution of §98:
it generally has nonzero entries for `d<=J`.  Thus Q645 does **not** prove
the stated normalization `z_1=...=z_J=0`; its phrase "the canonical
sequence" must not be confused with uniqueness of the underdetermined
column system.  It nevertheless suffices immediately for the tail-ideal
theorem and `content(T_(n,J))|A_n`.

The endpoint-interpolation half `(107.6)`--`(107.7)` proves the strictly
stronger statement that the special solution with the first `J` coordinates
zero is also integral.  This resolves the apparent conflict between Q645's
first values and the earlier values

`J=1: 0,-10,-58,-270,...`.

Q645 supplies no useful height saving: the canonical coordinates are
logarithmic-derivative coefficients of the fixed denominator `P_J` and
still grow exponentially with the coefficient index.  Their integrality
solves saturation, not the remaining Archimedean-size problem.

## 110. All-cutoff determinant amplification loses its multiplicity universally

The first exact audit of Q684's proposed all-cutoff determinant finds a
structural obstruction before any height estimate.

Fix the maximal q=1 cutoff `H=floor((n-1)/3)` and form a square matrix from
the coefficient/evaluation rows

`T_(n,0),T_(n,1),...,T_(n,H)`.

Its determinant is divisible by

`prod_(J=0)^H Gamma_(n,J)`.

Row subtraction replaces these rows, without changing the determinant, by

`h_J(c)=T_(n,J)(c)-T_(n,J-1)(c)`

`      =K_c(n,J)g_J^c`.

Now fix a q=1 candidate prime `p` with folded index `j`.  For every cutoff
`J>=j`, the Lucas block is already complete and

`T_(n,J)(c)=A_j(1+2c^p) mod p`.

Consequently

`h_J(c)=0 mod p` for every `J>j`,

**whether or not** `p|A_j`.  Thus every such determinant contains the
universal factor

`p^(H-j)`

before the selective bad-prime condition is imposed.  If `p|A_j`, the
product of contents asks for exponent `H-j+1`; after the universal factor is
removed, only one selective copy of `p` remains.

Therefore multiplying over cutoffs does not amplify the useful condition.
It amplifies a universal completed-block zero by exactly the same
multiplicity.  This is the determinant analogue of Q664's positive-prefix
obstruction.

Exact consecutive coefficient minors support the diagnosis.  For rows
`h_J`, `0<=J<(n-1)/3+1`, the best consecutive minor through `n=40` always
uses the lowest coefficient columns, but

`log|det|/(n * number_of_rows)`

stays between about `1.0` and `1.44` and increases along the fixed residue
classes.  The top-coefficient minor has rate about `1.7`--`2.0`.
These are full quadratic-height determinants dominated by the universal
factors above.

Reopen the all-cutoff route only if a quotient construction removes the
universal `p^(H-j)` for every candidate while provably retaining the final
selective copy.  Merely choosing a different coefficient minor, evaluation
basis, or unimodular row basis cannot change this local Smith multiplicity.

The exact audit is
`projects/Q-series-and-Chan-s-work/q32_all_cutoff_determinant.py`.

## 111. Rational kernel of the tail matrix is explicit

The convolution `(102.1)` also identifies the entire rational kernel of the
tail matrix.  Write a prospective homogeneous column vector as

`U(x)=sum_(d>=1)u_d x^d/d`.

Then `A u=0` in every row `k>J` if and only if

`P_J(x)U(x)=R(x)`

for a polynomial `R` of degree at most `J` with `R(0)=0`.  Hence

`U(x)=R(x)/P_J(x)`.                                      `(111.1)`

The rational kernel has dimension exactly `J`.  A simple integral
sublattice is obtained from

`R=x^r`, `1<=r<=J`,

because `P_J(0)=1` implies `1/P_J in Z[[x]]`, and therefore

`u_d=d[x^d]x^r/P_J in Z`.

This explains conceptually why the canonical Dwork solution and the
zero-prefix solution can differ: their difference is one of the logarithmic
kernel vectors `(111.1)`.

The full integral kernel may be larger than the displayed sublattice,
because `u_d in Z` requires only

`d[x^d]R/P_J in Z`,

not `[x^d]R/P_J in Z`.  Determining that finite-index enlargement and its
Smith invariants is Q690.  Even the rational reduction leaves `N-J`
independent coefficient coordinates, so the kernel formula alone does not
yet give a bounded-dimensional or subexponential height certificate.

## 112. Q672 independently proves the logarithmic-Wronskian integrality

Q672 gives a second prime-power proof of the canonical Franel tail vector.
In the notation of Sections 103--107, it proves

`theta(D_J/P_J) in Z[[x]]`

by applying the Frobenius-twisted Dwork congruence to the same
reflexive-hexagon Laurent polynomial over `Z_p[[x]]`, logarithmically
differentiating, and taking the p-adic limit of integral difference
quotients.  It also packages the coefficient statement as integrality of
the Wronskian numerator

`P_J theta D_J-D_J theta P_J`.

This independently confirms the prime-power cancellation in every apparent
denominator

`m binom(m-1,J)`.

As in Q645, Q672 constructs the canonical lattice solution and therefore
proves tail *membership*.  It does not impose the sharp zero-prefix
normalization.  The latter still needs the finite-numerator endpoint
interpolation in `(107.6)`--`(107.7)`.  It also gives no Archimedean height
saving, so it confirms the auxiliary theorem but does not change the final
moving-zero barrier.

Operational note: a conventional `ask-gpt.py --help` invocation was
accidentally repeated and submitted the literal question `--help` as Q691
on dm6.  The script has no help mode; its source header is the usage
documentation.  Q691 must be allowed to clear and dm6 then refilled with the
intended content-to-height question.  This repeats the already recorded Q616
operator failure and yields no mathematical information.

Q692 was sent to dm4 to test a precise multiscale-cutoff quotient:
remove the universal completed-block multiplicities integrally, retain one
selective bad-prime copy, and prove a sublinear normalized determinant
height.  A valid answer must address both the local Smith quotient and the
Archimedean height; tail integrality alone is insufficient.

## 113. Primitive row saturation makes the all-cutoff divisor small but nonselective

Section 110 measured only individual raw minors.  A stronger normalization is
available because

`h_J(c)=K_(n,J)(c) g_J(c)`

and `g_J` is primitive: its top coefficient is `+-F_0=+-1`.  Gauss's lemma
therefore gives

`content(h_J)=content(K_(n,J))`,

an entirely universal row factor.  Divide every cutoff-difference row by
this content and take the gcd of all maximal minors, equivalently the product
of the nonzero Smith invariants of the primitive row matrix.

This primitive determinantal divisor is dramatically smaller than the raw
minors.  For example its values for `n=20,30,39,42` are respectively

`1105, 8500, 53476453125, 4599472500000`.

Every aligned q=1 bad prime survives through `n=42`.  However exact
factorization gives decisive contamination:

- `n=23`: both `17` and `19` divide the primitive divisor, although neither
  is an aligned bad prime for this `n`;
- `n=38`: the above-`n` prime `47` divides it;
- `n=42`: the good candidate `29` and above-`n` prime `53` divide it.

Thus universal row-content saturation does not retain only the final
selective copy.  It detects additional Franel/structural rank losses at
unrelated positions.  Intersecting with a q=1 carrier removes the above-`n`
factors but not good candidates such as `17,19` at `n=23`.

This does not rule out a more delicate multiscale quotient or an intersection
with a second independently controlled lattice, but it closes the direct
claim that the primitive all-cutoff Smith divisor itself has the target
top-half support.

The exact audit is
`projects/Q-series-and-Chan-s-work/q32_primitive_cutoff_smith.py`; it checks
the target divisibility and the stated counterexamples through `n=42`.

## 114. Q684 proves the exact all-cutoff stabilization obstruction

Q684 independently derives

`h_J(c)=K_c(n,J)g_J(c)`

and proves that, for a q=1 candidate prime with folded index `j_p`, every
stabilized difference row `h_J`, `J>j_p`, vanishes modulo `p` without any
Apéry-zero hypothesis.  For an arbitrary subset of cutoffs, the universal
valuation is exactly the number of selected stabilized differences above
`j_p`; badness supplies at most one further copy.  Sparse, dyadic, and
overlapping cutoff sets therefore do not restore a weighted multiplicity
gain.

It also evaluates the universal q=1 factor asymptotically:

`log U_(1,n)=n^2/36+o(n^2)`.

There is a disjoint q=2 stabilization factor with

`log U_(2,n)=n^2/300+o(n^2)`.

Thus dividing only the q=1 factor cannot produce a small global determinant.
Q684 further gives an exact product formula for the top-coefficient minor and
shows a strictly positive raw normalized quadratic-height rate.

The Smith determinantal divisor is the intrinsic optimum over all unimodular
coefficient/evaluation functional bases, but it contains the same universal
factors.  After all explicit stabilization factors are divided out, the
sufficient statement is again a new top-localized Smith theorem of
sublinear logarithmic height.  This is not supplied by tail integrality.

This confirms Section 110 and sharpens its reopen condition: a future
determinant construction must remove *all* Lucas digit-stabilization factors
integrally and then prove an `o(n)` bound for the residual Smith divisor.

Q699 was sent to dm7 for the exact dual Smith transform of the Franel tail
matrix and the transformed Legendre outer row.  Q700 was sent to dm1 for the
remaining growing-degree algebraic-center problem: construct a moving monic
integer polynomial with sublinear normalized resultant, or prove the
corresponding weighted-capacity obstruction.

Operational correction to Section 112: Q691 was not allowed to consume dm6.
After an accidental read of the bridge's task-claim endpoint moved it from
pending to processing, the task was explicitly marked
`operator-cancelled-q691` through its own completion endpoint, without
clearing or altering any other queued task.  The later failure notification
is stale.  dm6 was refilled with Q695, the Apéry-truncation square-factor
attack.

## 115. Q692 closes fixed-order multiscale determinants; cutoff-content intersection survives

Q692 gives the exact local Smith ledger for arbitrary cutoff determinants.
For cutoffs

`J_0 < ... < J_r`,

rows beyond the folded index `j_p` are all equal modulo a q=1 candidate
prime `p`.  If `q` selected rows lie beyond `j_p`, every determinant has
the universal factor `p^(q-1)`.  The extra hypothesis `p|A_(j_p)` raises
the forced valuation only to `q`.  Thus after exact universal saturation
there is only **one** selective copy of `p`, independent of the number of
rows.  A model matrix shows this is sharp.

For the canonical constant/linear two-coefficient minor, Q692 computes the
exact exponential rate.  If the top cutoff is `B=beta n+O(1)`, then after
universal saturation the quotient still has rate

`Phi(beta)=2(1+beta)log(1+beta)-4 beta log beta
           -2(1-beta)log(1-beta)`.

At `beta=1/3`, `Phi(beta)=log 8`.  Repeating over many lower cutoffs
therefore repeats one target copy at cost at least `exp((log 8+o(1))n)`
per copy.  Combining a partial cutoff with the untreated candidate tail
recovers only the trivial total constant `1/2`; branchwise it recovers
`1/3` for the direct branch and `1/6` for the reflected branch.  The same
positive-rate obstruction holds for every fixed set of coefficient degrees.

This closes raw/fixed-order multiscale determinants.  It does not close
degree-growing signed functionals.  The precise reopen condition is a
subexponential geometric-mean bound for branch-purified saturated
determinants, not merely an `o(nr)` bound on one `r`-row determinant.

A different use of all cutoffs avoids determinants entirely.  Let

`Gamma_(n,J)=content(T_(n,J))`,

`H=floor((n-1)/3)`, `M=floor(n/2)`,

and

`I_n=gcd_(H<=J<=M) Gamma_(n,J)`.                    `(115.1)`

Every q=1 bad prime divides `I_n`: its folded index is at most `H`, and
all cutoffs in `(115.1)` remain strictly below the prime.

There is an exact efficient reduction.  Consecutive transforms satisfy

`T_(n,J)-T_(n,J-1)=K_(n,J) g_J`.

The polynomial `g_J` is primitive because its top coefficient is
`+-F_0=+-1`.  Gauss's lemma gives

`content(T_(n,J)-T_(n,J-1))=content(K_(n,J))=:kappa_(n,J)`,

where

`kappa_(n,J)=gcd_(J<=k<=n)
 binom(n,k)binom(n+k,k)binom(k,J)`.                 `(115.2)`

For integer coefficient vectors, the ideal generated jointly by two
vectors `u,v` equals the ideal generated by `u,v-u`.  Iterating this
identity proves

`I_n=gcd(Gamma_(n,H), kappa_(n,H+1),...,kappa_(n,M))`. `(115.3)`

This is not a nesting claim about the individual contents.

Exact exhaustive computation for every `n<=400` found no failure of target
divisibility.  The maximal `log(I_n)/n` on dyadic windows was

```
(10,20]    0.217940
(20,40]    0.148090
(40,80]    0.137889
(80,160]   0.083736
(160,320]  0.058712
(320,400]  0.049227
```

The corresponding maximal nuisance rates after dividing by the q=1 target
were

`0.217940, 0.076640, 0.083220, 0.083736, 0.042420, 0.016510`.

Representative exact values are

```
n=20   I_n=17                    target=17
n=30   I_n=85                    target=17
n=142  I_n=5*163*179             target=1
n=200  I_n=5*139*181             target=139*181
n=321  I_n=179*193*211           target=I_n
```

The apparent purification is strong, but it is not yet a proof:
`log I_n=o(n)` already contains the hard q=1 moving-zero statement.
The value of `(115.3)` is structural: it converts the target into the
intersection of one Franel-tail content and explicit binomial tail
contents, and it automatically removes the q=2 stabilization range.

Modulo a prime `ell`, write `n=sum n_i ell^i` and

`s_i=min(n_i,ell-1-n_i)`.

Lucas gives the nonzero exponent box of the Laurent coefficient sequence

`0<=k_i<=s_i`.

Consequently the `J`-th Hasse derivative is nonzero modulo `ell` exactly
when `J_i<=s_i` digitwise.  Equivalently,

`ell|kappa_(n,J)` iff some base-`ell` digit `J_i>s_i`. `(115.4)`

Thus a prime survives the intersection `(115.3)` precisely when its digit
support box avoids the whole interval `[H+1,M]`, in addition to dividing
the initial content.  For q=1 candidates this avoidance is automatic; for
q=2 candidates the point `J=p` lies in the safe interval and removes them.

The executable audit is
`projects/Q-series-and-Chan-s-work/q32_multicutoff_content_intersection.py`.
It compiles, verifies the stated samples, and checks target divisibility.

## 116. The central-third intersection exactly purifies the direct large-prime support

The common safe interval can be enlarged if only the direct q=1 branch is
retained.  Put

`H=floor((n-1)/3)`,

`U=n-H-1=floor(2n/3)`,

and

`D_n=gcd_(H<=J<=U) Gamma_(n,J)
    =gcd(Gamma_(n,H),kappa_(n,H+1),...,kappa_(n,U))`. `(116.1)`

A direct candidate has `p=n-j`, `j<=H`.  Hence every cutoff in `(116.1)`
satisfies

`j<=J<=n-H-1<=n-j-1=p-1`,

so every direct bad prime divides `D_n`.

There is an exact large-prime support theorem for the difference part of
`(116.1)`.

**Central-third Lucas lemma.**  Let `p<=n` be prime and `p^2>2n`.  Then
`p` divides every

`kappa_(n,J)`, `ceil(n/3)<=J<=floor(2n/3)`,

if and only if

`n=p+r` with `0<=r<=H`.                              `(116.2)`

To prove it, write `n=qp+r`.  The size hypothesis gives
`q< p/2`, hence the high base-`p` support digit is

`min(q,p-1-q)=q`.

Writing `j=min(r,p-1-r)`, Lucas gives the complete nonzero coefficient
support

`union_(0<=b<=q) [bp,bp+j]`.                         `(116.3)`

If `q>=2`, the real interval

`[n/(3p),2n/(3p)]`

contains an integer `b` with `1<=b<=q`: take
`b=ceil(n/(3p))`; for `x=n/p>=2`, `ceil(x/3)<=2x/3`.
Thus `bp` belongs both to `(116.3)` and to the central third, so the
corresponding `kappa` is nonzero modulo `p`.

If `q=1`, the reflected case `r>(p-1)/2` has `p` itself in the central
third and is likewise removed.  In the direct case, the support blocks are
`[0,r]` and `[p,p+r]`; they miss the central third exactly when `r<=H`.
This proves `(116.2)`, including the floors.

Consequently, after intersecting `D_n` with the direct binomial carrier,
every prime larger than `sqrt(2n)` is an actual direct q=1 candidate; the
carrier removes primes above `n`.  Primes at most `sqrt(2n)` have only
`exp(o(n))` total possible carrier weight by Kummer: each valuation is at
most the number of base-`p` carries, hence contributes `O(log n)`, and
there are `O(sqrt(n))` possible primes.

This is an exact support purification, but it still does **not** bound the
large direct factors.  In fact

`gcd(D_n,binom(n,H+1))`

has the direct target as its large-prime support, so asserting that this gcd
is subexponential is another equivalent form of the remaining direct
moving-zero problem.  The central-third lemma removes every q>=2 and
reflected large-prime nuisance but supplies no Archimedean bound for the
surviving direct factors.

Exact examples are

```
n=56   D_n=61
n=142  D_n=5*163*179
n=200  D_n=139*181       (both are direct bad primes)
n=321  D_n=1             (the three reflected hits are all removed)
n=394  D_n=5*17
n=400  D_n=1
```

The first scan used the slightly smaller upper endpoint
`floor((2n-2)/3)` and reported `D_56=37*61`.  This was not the maximal
common safe interval: at `n=56`, the correct endpoint is `37`, whose final
`kappa` removes the factor `37`.  The script and the values above use the
correct identity `U=n-H-1=floor(2n/3)`.

The updated executable audit proves the finite support equivalence through
`n=2000` in addition to the exact content samples.  Q709 asks for a hostile
independent audit and for any non-tautological height consequence.

## 117. Q701 classifies the common intersection; Q700 narrows resultants

Q701 independently proves the Euclidean reduction `(115.3)` and the exact
base-`p` Hasse-support box `(115.4)`.  For the original interval
`[H+1,M]`, its complete large-prime classification has exactly three
surviving slopes:

Terminology correction: the integer `gcd_J Gamma_(n,J)` generates the sum
of the corresponding principal ideals; their ideal-theoretic intersection
is generated by the lcm.  Hereafter “intersection” means only the
common-content gcd operation, not an intersection of ideals.

1. `q=0`: primes above `n` whose folded degree is at most `H`;
2. `q=1`: the desired two top-half branches;
3. one reflected `q=3` band, of asymptotic prime-interval width `n/60`.

The q=2, q=4, q=5, and every q>=6 large-prime band meets the digital
support and is removed by one increment.  Q701 further uses the recursive
cutoff profile to show that divisibility of the initial content on the
three surviving plateaux is respectively an Apéry-zero condition at the
associated folded index.  The large carrier is squarefree; small-prime
valuations contribute only `O(sqrt(n) log n)`.

Thus the raw common-content gcd asks for a three-slope moving-zero theorem,
strictly stronger than q=1.  Its q=0 branch is permanent under later
cutoffs: once the reduced polynomial has degree below the starting cutoff,
every subsequent increment vanishes.  Extending an undifferentiated tail
of cutoff contents cannot remove it.

Section 116 avoids the reflected q=3 branch by extending to the entire
central third and retaining only direct q=1 primes.  It still needs a
non-cutoff selective height estimate; the exact support theorem by itself
does not provide one.

Q700 analyzes the growing-degree monic-resultant route.  If `P_n` is monic
of degree `d_n`, then

`Gamma_n^d_n | Res(P_n,T_n)`.

The route would solve q=1 if a nonzero resultant of size
`exp(o(n d_n))` could be constructed.  The canonical algebraic-integer
construction on a polynomial lemniscate reaches only the positive rate

`log 3+(4/3)log 2 = 2.0228...`.

This is the exact threshold for the strategy that places every conjugate
inside one sublevel set of `|T_n|`.  Classical capacity does not exclude a
moving exceptional family below that threshold because both the polynomial
and compact set vary with `n`; finite splitting conditions also cannot
force a unit resultant against all primes.  The remaining question is an
arithmetic weighted-equilibrium problem over complete conjugate measures,
including an escape-of-house variable.  Q700 neither constructs the
required subexponential resultants nor proves a uniform positive lower
bound.  It closes only the canonical single-lemniscate implementation.

## 118. Apéry truncation squares do not eliminate congruence classes; naive mod-p-square reflection fails

The source of Caruso--Fürnsinn--Vargas-Montoya--Zudilin was checked directly.
For every odd prime `p`, their theorem states

`A_p(t)=B_p(t)^2`

for `p=1,5,7,11 mod 24`, and

`A_p(t)=(t^2-34t+1)B_p(t)^2`

for `p=13,17,19,23 mod 24`.

The proof obtains this from the rational substitution relating the Apéry
and Franel series and from the quadratic extension with discriminant
`t^2-34t+1`.  It does not provide a second low-coefficient vanishing
condition when one coefficient `A_j` is zero.  Since the constant term of
the square factor is a unit, its coefficients through degree
`j<p/2` are recursively determined by the same initial Apéry coefficients;
a zero coefficient is still one scalar convolution equation.

An exact scan of every zero pair for `p<20000` found zeros in all eight
admissible odd residue classes modulo 24, with nearly balanced counts:

```
class:  1   5   7  11  13  17  19  23
zeros: 308 282 268 267 322 289 296 258
```

Direct and reflected counts are equal up to central endpoints, as expected
from the universal palindromy

`A_(p-1-j)=A_j mod p`.

The most naive lift of this symmetry is false.  Computing both terms
modulo `p^2`, the selective condition `p|A_j` generally does **not** imply

`p^2 | A_(p-1-j)-A_j`.

Already `(p,j)=(5,1)` gives the normalized difference `3 mod 5`;
further counterexamples include `(17,3)` giving `16 mod 17`,
`(19,8)` giving `17 mod 19`, and `(73,2)` giving `9 mod 73`.
The accidental zero at `(11,5)` is not a uniform theorem.

This rules out using the uncorrected reflection difference as a second
p-adic equation.  A harmonic correction could still exist; Q715 is tasked
with deriving and auditing the exact mod-`p^2` reflection formula.  Q716
tests whether the separate high-cutoff plateau of reflected q=1 primes can
yield a branch-separating integral family.

The harmonic correction can in fact be derived directly.  For
`0<=j<=(p-1)/2`, put

`W_j=sum_(k=0)^j binom(j,k)^2 binom(j+k,k)^2
                (H_(j+k)-H_(j-k))` in `Z_(p)`.       `(118.1)`

Expanding

`binom(p-1-j,k)`

and

`binom(p-1-j+k,k)`

to first order in `p`, while observing that every `k>j` term contains a
squared factor `p`, proves the exact congruence

`A_(p-1-j) = A_j-2p W_j (mod p^2)`.                  `(118.2)`

Hence under `p|A_j`,

`(A_(p-1-j)-A_j)/p = -2W_j (mod p)`.                 `(118.3)`

The right side is generally nonzero; `(118.2)` explains every numerical
counterexample above.  It supplies a corrected reflection identity but no
second vanishing condition: forcing `(118.3)` to vanish would add the new
independent hypothesis `p|W_j`, which is absent in the target.

The executable proof audit
`projects/Q-series-and-Chan-s-work/q32_apery_reflection_p2.py`
checks `(118.2)` for every prime `p<=500` and finds 47 selective
counterexamples.  It compiles and passes `git diff --check`.

## 119. Literature recheck: available p-adic and prime-factor theorems have the wrong direction

A fresh search for an existing gcd theorem between Apéry values and
factorial ratios found no result beyond the obstruction already recorded in
Q632.

Delaygue's *Arithmetic properties of Apéry-like numbers*
(`arXiv:1310.4131`) proves lower p-adic valuation bounds from zero digits in
the base-`p` expansion, together with Lucas criteria for factorial-ratio
multisums.  In the present q=1 case this propagates a known zero
`p|A_j` to indices having digit `j`; it does not upper-bound the set of
primes for which the first zero occurs.

Luca--Shparlinski, *Arithmetic properties of Apéry numbers*,
J. London Math. Soc. 78 (2008), proves lower bounds for the number and size
of prime factors of `A_n` on a density-one set of indices.  These are
vertical/typical-index results and do not bound, for every fixed `n`, the
diagonal primes `p=n-j` or `(n+1+j)/2` dividing `A_j`.

The searches for “G-function gcd factorial ratio”, “holonomic gcd binomial
coefficient”, and “Apéry zero set modulo primes” returned congruence and
supercongruence papers, but no pointwise subexponential gcd theorem of the
required type.  Therefore no literature theorem can currently be inserted
at the last step without proving a new horizontal nonalignment result.

Q718 tests whether Delaygue's zero-digit multiplicity can be compressed into
a height-efficient global certificate.  Literal digit packing is suspect:
placing the same nonzero digit in `m` base-`p` positions forces the index to
be at least `p^(m-1)`, so direct evaluation of the exponentially growing
Apéry sequence loses far more height than the gained p-adic multiplicity.

## 120. Exact q=1 all-cutoff profile: two plateaux and one gap

The arbitrary-cutoff transform has substantially more structure than the
single low-cutoff congruence.  Write

`n=p+r<2p`, `0<=r<p`, `j=min(r,p-1-r)`,

and put `x=c^p`.  For every odd prime `p`, Lucas factorization gives

`T_(n,J)(c)=(1+2x)T_(r,J)(c) mod p` for `J<p`,        `(120.1)`

while for `J=p+B`,

`T_(n,p+B)(c)
 =A_r(1+2x)+2(2-x)T_(r,B)(c) mod p`.                `(120.2)`

Here the factor `2-x` comes from the shifted Franel moment:

`S_(p+s)(c)=(2-c^p)S_s(c) mod p`.

Indeed, in its Lucas expansion the high digit of the summation index is
zero or one; these two pieces contribute respectively `(-c)^p S_s(c)` and
`F_1 S_s(c)=2S_s(c)`.

Modulo `p`, the effective degree of the Legendre kernel `Q_r` is exactly
`j`.  In the reflected case `r=p-1-j`, all coefficients above degree `j`
vanish because `binom(r+k,k)` carries.  Thus

`T_(r,J)=A_r=A_j mod p` once `J>=j`.                 `(120.3)`

For `J<j`, its degree-`j` coefficient is

`(-1)^J [t^j]Q_r(t) binom(j-1,J)`,                  `(120.4)`

which is nonzero modulo `p`.  Equations `(120.1)`--`(120.4)` prove the
exact selective profile:

**Two-plateau theorem.**  If `p|A_j`, then

`p|Gamma_(n,J)`

if and only if

`J in [j,p-1] union [p+j,n]`.                       `(120.5)`

The middle gap `[p,p+j-1]` has exactly `j` cutoffs.  In the direct branch
`r=j`, the high plateau is only the endpoint `J=n=p+j`.  In the reflected
branch `r=p-1-j`, the high plateau has `p-2j` cutoffs:

`p+j<=J<=2p-1-j=n`.

The executable audit `q32_q1_all_cutoff_profile.py` verifies the complete
polynomial identities, not only contents, at all 5,110 triples for odd
primes through `31`; it also checks 237 cutoff instances conditioned on an
actual Apéry zero.

This explains exactly why extending the low-cutoff common-content gcd into
the central third removes reflected primes: their gap begins at `J=p`,
whereas direct primes have `p>floor(2n/3)`.  It also identifies a possible
new source of branch separation, namely the different high-plateau lengths.
It is not yet a height bound.  The endpoint `Gamma_(n,n)=A_n` is
exponential, and multiplying all plateau contents gives only linear
selective multiplicity per prime at a generally quadratic total height.
Any successful use of `(120.5)` must remove the universal cutoff factors or
find a signed/Smith quotient whose Archimedean height is subquadratic.

Q727 asks ChatGPT for a hostile height analysis of the two-plateau
signature.  Q690 failed at the bridge level without a git or Notion drop;
its detached watcher remains alive, so Q727 is a distinct replacement
rather than a resend.

## 121. Q715 confirms that the reflection derivative is independent

Q715 independently derives `(118.2)` in the notation

`D_j=2W_j`,

so that

`A_(p-1-j)=A_j-pD_j mod p^2`.                        `(121.1)`

It identifies `D_j` as the derivative in the **index parameter** of the
truncated hypergeometric continuation.  This is not the derivative in the
generating variable appearing in the CFVMZ square factorization.  Thus a
coefficient zero of the latter polynomial gives no reason for `D_j` to
vanish.

The strongest counterexample is `(p,j)=(17,3)`.  Here

`A_3=1445=5*17^2`,

but `D_3=1 mod 17`, so the reflected normalized value is nonzero.  Because
the lower member already vanishes modulo `17^2`, no choice of a scalar
multiplier `lambda` can make

`A_(p-1-j)-lambda A_j`

vanish modulo `p^2`.  This closes every scalar-corrected reflection
combination, not only the uncorrected choice `lambda=1`.

The fast differentiated-recurrence scan from Q715 was rerun locally and
incorporated into `q32_apery_reflection_p2.py`.  Through every odd prime
`p<=5000`, excluding the central reflection fixed point, it finds:

```
303 lower-half zero pairs at 247 primes,
no zero with D_j=0,
no p^2 reflection coincidence,
only (p,j)=(17,3) with p^2|A_j,
and no reflected member divisible by p^2.
```

The termwise harmonic formula is still checked separately through `p<=500`.
The enlarged script compiles and completes in about 18 seconds.

The observed simplicity of all 303 roots is evidence, not an all-prime
squarefreeness theorem.  Even proving it for all primes would only certify
that reflection supplies no second zero; it would not give the desired
horizontal prime-weight bound.

After Q715 freed dm1 and Q695 failed without a drop on dm6, the two tabs
were refilled with distinct tasks:

- Q729: recurrence-aware direct-branch interpolation or a sharp no-go
  theorem for that class;
- Q730: a Cartier/Picard--Fuchs/finite-field attack that must explicitly
  bridge per-prime information to the pointwise additive diagonal.

## 122. Q709 gives the exact central-third carrier and its positive height

Q709 independently verifies the floors in §116 and strengthens the support
statement to an exact exponent classification.  Put

`J0=H+1=ceil(n/3)`, `U=floor(2n/3)`,

and

`K_n=gcd_(J0<=J<=U) kappa_(n,J)`.

For every prime `ell>sqrt(2n)`,

`ell|K_n` if and only if `U<ell<=n+J0`,                  `(122.1)`

and every such prime occurs to exponent exactly one.  Thus

`K_n=S_n product_(U<ell<=n+J0) ell`,                    `(122.2)`

where `S_n` is supported on primes at most `sqrt(2n)` and
`log S_n=o(n)`.  The interval splits into:

- the within-`n` direct carrier `U<ell<=n`;
- an unavoidable above-`n` q=0 carrier `n<ell<=n+J0`.

The proof is a two-digit Lucas/Kummer box calculation.  If `ell<=U`, the
middle third contains a multiple of `ell`, producing a unit coefficient in
one `kappa_(n,J)`.  If `U<ell<=n`, every middle-third cutoff lies strictly
between the two nonzero support blocks.  If `ell>n`, the smallest cutoff
`J0` gives exactly the band through `n+J0`.  A boundary coefficient then
shows valuation one.

This corrects any hope that the common increment gcd might itself have
small height.  By the prime number theorem, `(122.2)` has a strictly
positive linear logarithmic rate.  Intersecting with the initial content
`Gamma_(n,H)` prunes the within-`n` carrier exactly by the original
condition `ell|A_(n-ell)` and prunes the above-`n` carrier by another
Apéry-zero condition.  The direct binomial carrier

`B_n^-=binom(n,J0)`

removes the above-`n` band but contains every within-`n` carrier prime
regardless of Apéry badness.  Therefore the purified gcd is an exact
ordinary-gcd reformulation of the direct target, with only `o(n)` small
prime nuisance, but has no non-tautological height saving.

This closes the raw central-third increment/carrier route.  Reopen it only
if a new integral row prunes within-`n` carrier primes independently of
`Gamma_(n,H)`.

## 123. Q716: arbitrary-cutoff blocks and exact reflected-branch jumps

The full Notion page was fetched because its local Markdown conversion
dropped displayed formulas.  Its key identity is valid for every quotient.
Write

`n=qp+r`, `J=Ap+B`, `0<=r,B<p`.                       `(123.1)`

Then polynomial Lucas and Franel Lucas give

`K_(qp+r,ap+b)(c)
 =K_(q,a)(c^p)K_(r,b)(c) mod p`,

`g_(ap+b)(c)=g_a(c^p)g_b(c) mod p`,

and hence

`T_(qp+r,Ap+B)(c)
 =A_r T_(q,A-1)(c^p)
  +K_(q,A)(c^p)g_A(c^p)T_(r,B)(c) mod p`.            `(123.2)`

The convention is `T_(q,-1)=0`.  This is a direct block factorization, not
an experimental cancellation.  The new executable audit
`q32_allq_cutoff_block.py` verifies 2,984 complete polynomial identities
for `p=5,7,11,13` and every `q<=3`, residue, and cutoff.

For a reflected low digit `r=p-1-j`, the low family enters at cutoff `j`,
while the final base-`p` plateau begins at

`S_q(j)=qp+j`.                                        `(123.3)`

Q716 turns this into a branch separator.  Define

`H=floor((n-1)/3)`, `Hsharp=floor((n-3)/3)`,

`L_h^flat=rad gcd_(h<=J<=H) Gamma_(n,J)`,

`lambda_h=L_h^flat/L_(h-1)^flat`,

and

`sigma(h)=ceil((n+1+3h)/2)`,

`U_h^flat=rad gcd_(sigma(h)<=J<=n-1) Gamma_(n,J)`,

`upsilon_h=U_h^flat/U_(h-1)^flat`.                    `(123.4)`

These are genuine integer jump quotients because both gcd families are
nested.  A reflected quotient-`q` bad prime enters the low jumps at `h=j`
and the high jumps at

`rho_q(j,p)
 =floor((2S_q(j)-n)/3)
 =j+floor(((q-1)p+1)/3)`.                             `(123.5)`

Under the extra hypothesis `p>H`, Q716 states that `rho_q<j` for q=0,
`rho_q=j` for q=1, and `rho_q>j` for q at least 2.  Above-`n` q=0 factors
are already present in every high-tail gcd and make no positive jump.
Direct q=1 factors have no high plateau before `J=n`.  It then claims that

`product_h gcd(lambda_h,upsilon_h)`                    `(123.6)`

has exactly the reflected q=1 large-prime support.  Section 125 below
records a hostile counterexample to this extrapolation beyond `p>H` and a
carrier-based repair.

This is a real support purification, but not yet a height theorem.  The
trivial prime-interval bound is still

`log (123.6) <= (1/6+o(1))n`,

and each surviving bad prime appears only once.  Q716 promised a Section 9
with a weighted diagonal-jump estimate, but the delivered page ends at
Section 8.  Q733 explicitly requests the missing height analysis; Q734
separately studies the exact loss at cutoff `J=p`.

## 124. Larger diagonal scan: strong sparsity, still no deterministic bridge

The existing C++ diagonal scanner was compiled with optimization and rerun
through `n=50000`.  It computed 5,172 Apéry zero pairs over all primes in
range.  In the q=1 slice:

- the maximum number of hits at one `n` is still only three;
- the maximum direct count is three, attained at `n=11576,26164`;
- the maximum reflected count is three, attained at `n=321,47066`;
- on the last dyadic window `(25000,50000]`, the worst q=1 logarithmic
  mass ratio is `0.0011467269`, at `n=26164`;
- the seven three-hit outer indices through 50,000 are
  `321,11576,18444,22101,26164,47066,47859`.

The pure triples include both extremes: `n=11576` has three direct hits,
while `n=321` and `n=47066` have three reflected hits.  Thus neither branch
can be dismissed as an endpoint artifact.

This strongly suggests that even a very weak uniform incidence bound, for
example

`# {p in (n/2,n] : p|A_(n mod p)} = n^(o(1))`,        `(124.1)`

would finish q=1 immediately, since each prime contributes `O(log n)`.
But the scan supplies no proof of `(124.1)`.  Per-prime zero-fiber bounds,
including the unconditional `O(p^(2/3))` continuant bound, do not control
concentration on one additive diagonal.  This route should be reopened only
with a genuine cross-prime incidence theorem, not a larger computation.

## 125. Hostile correction to Q716: a q=1 mod 3 alias

Q716's intrinsic purification claim is false as stated.  Its Section 8
starts with `p>H`; for a reflected q=3 prime one generally has `H>p`, and
for q=4 the hypothesis certainly fails.  The conclusion silently removes
this assumption.

The block identity `(123.2)` gives the exact correction.  Suppose

`n=qp+p-1-j`, `p|A_j`,

with `q` fixed and small relative to `p`.  In every base-`p` cutoff block,
the zero plateau begins at low digit `j`.  Let

`a=floor(H/p)`.

When the low tail is nonempty at this prime, its jump is at

`h_low=ap+j`.                                          `(125.1)`

The high jump remains

`h_high=j+floor(((q-1)p+1)/3)`.                       `(125.2)`

For `q=3s+1`,

`a=s`

and the right side of `(125.2)` is exactly `sp+j`.  Thus the two jump
indices agree for **every reflected quotient**

`q=1,4,7,10,...`,                                     `(125.3)`

not only q=1.  For the genuine zero `17|A_3`, exact polynomial computation
gives:

```
q=1, n=30:  low jump=3,  high jump=3
q=4, n=81:  low jump=20, high jump=20
q=7, n=132: low jump=37, high jump=37
```

The q=2,3,5,6 cases have distinct jump indices.  The executable audit
`q32_reflected_jump_alias.py` checks every cutoff polynomial for q=1 through
q=7 and reproduces this periodic alias.

There is a clean support repair, but it is no longer intrinsic to the jump
system.  Intersect `(123.6)` with the lcm-interval carrier

`C_n^R
 =lcm(1,...,floor(2n/3))/lcm(1,...,floor(n/2))`.       `(125.4)`

For primes above `sqrt(n)`, `(125.4)` has exactly the support

`n/2<p<=2n/3`.

It therefore keeps reflected q=1 and removes the q=4,7,... aliases; prime
powers at most `sqrt(n)` contribute only `o(n)`.  This repairs the exact
large-prime support theorem, but its PNT height is still

`log C_n^R=(1/6+o(1))n`.                              `(125.5)`

So Q716 supplies a useful quotient-congruence filter only after the
interval carrier is added, and the repaired construction still gives no
sublinear height.  Q733 was already running from the uncorrected premise;
its answer must be audited against `(125.3)` before use.

Operationally, Q699 failed without a git or Notion drop.  Its watcher was
left alive, and dm7 was refilled with the distinct Q737 irrationality/
middle-prime localization audit.

## 126. A cleaner all-cutoff quotient filter: count loss transitions

The failure in §125 suggests using all repeated block transitions rather
than matching one low and one high jump.  This gives an exact new quotient
statistic.

Put

`R_(n,J)=rad Gamma_(n,J)`

and define the squarefree loss quotient

`Lambda_(n,J)
 =R_(n,J-1)/gcd(R_(n,J-1),R_(n,J))`, `1<=J<=n`.       `(126.1)`

It is an integer.  Let

`P_n=product_(J=1)^n Lambda_(n,J)`.                   `(126.2)`

For a prime `p`, `v_p(P_n)` is exactly the number of zero-to-nonzero
transitions in the all-cutoff content profile.

Write `n=qp+r`, `0<=q,r<p`, and

`j=min(r,p-1-r)`.

Suppose `p|A_j`.  The effective degree of `Q_q(t)` modulo `p` is

`e=min(q,p-1-q)`.                                     `(126.3)`

Indeed, the coefficients through `e` are nonzero and every later
coefficient has a carry.  In block `A`, equation `(123.2)` reduces to

`T_(n,Ap+B)
 =K_(q,A)(c^p)g_A(c^p)T_(r,B)(c) mod p`,             `(126.4)`

because `A_r=A_j=0 mod p`.  The factor `g_A` is always nonzero because it
is primitive.  The factor `K_(q,A)` is nonzero exactly when `A<=e`: its
top coefficient comes from the effective top degree of `Q_q`.  Finally,
`T_(r,B)` is nonzero for `B<j` and zero for `B>=j`.

Consequently, for `A=0,...,e`, each base-`p` block has a terminal zero
plateau, followed at the next block boundary by a loss; once `A>e`, the
whole remaining transform is identically zero modulo `p`.  Therefore

`boxed(v_p(P_n)=min(q,p-1-q))`.                       `(126.5)`

This counts selective losses without any q=1-mod-3 alias.  In particular,
a bad large prime occurs exactly once in `P_n` if and only if

`q=1` or `q=p-2`.                                     `(126.6)`

The second case is harmless pointwise.  It forces

`p^2-2p<=n<=p^2-p-1`,

and these intervals are disjoint as `p` varies, so at most one such prime
exists for a fixed `n`; its contribution is `O(log n)`.

Equivalently, define the radical of the exactly-once part by

`Once(P_n)
 =rad(P_n)/rad(gcd(P_n,P_n/rad(P_n)))`.               `(126.7)`

Up to small primes and the single near-square-root case in `(126.6)`, the
large-prime support of `(126.7)` is exactly the full q=1 bad radical, with
direct and reflected branches simultaneously.

Intersecting `(126.7)` with

`lcm(1,...,n)/lcm(1,...,floor(n/2))`

removes even that near-square-root exception: above `sqrt(n)` this carrier
has exactly the prime support `(n/2,n]`.

The updated `q32_reflected_jump_alias.py` checks the actual cutoff
polynomials for the zero `17|A_3`.  For q=1 through q=7 it finds losses
exactly at

`p,2p,...,qp`,

so their counts are `1,2,...,7`, even though the old same-index construction
aliases q=1,4,7.

This is the cleanest support purification obtained so far, but it still
does not prove a height bound.  The operation `(126.7)` merely extracts the
squarefree-once part of a potentially huge product of content-transition
quotients.  No current estimate shows

`log Once(P_n)=o(n)`.                                 `(126.8)`

The new narrow route is therefore to exploit the *global transition
lattice* to bound the exactly-once radical, rather than bounding individual
contents or matching two selected jumps.

## 127. Q733 and Q734 audit: useful transition algebra, no height theorem

Q734 independently derives the exact loss/gain quotients for adjacent cutoff
contents and verifies the q=1 transition profile:

- a bad q=1 prime is lost at cutoff `J=p`;
- a direct prime is regained only at `J=n`;
- a reflected prime is regained at `J=p+j`;
- the product of prime-cutoff loss factors has exactly the q=1 bad support.

It also proves that an lcm of diagonal losses over a terminal cutoff interval
recovers the full large moving-zero set, up to a square-root band of total
weight `O(sqrt(n))`.  This is an exact re-encoding, not a height compression.
Its Smith, resultant, factorial-carrier, and telescoping audits all have the
same conclusion: after artificial repeated losses are squarefree-saturated,
the remaining radical is the original moving-zero radical.

Q734's raw-product multiplicity is the number of relevant cutoff multiples in
its chosen terminal interval.  Formula `(126.5)` is a stronger all-cutoff
refinement: over *every* cutoff, the loss multiplicity is exactly
`min(q,p-1-q)`.  The exactly-once operation `(126.7)` therefore removes all
ordinary `q>=2` slices and leaves q=1 plus the isolated `q=p-2` edge.

Q733 analyzes the same-index low/high jump system of Q716.  Its support claim
is not valid as written: it repeats Q716's assertion that equal jump indices
characterize reflected q=1, without checking the q=4 and q=7 counterexamples
in §125.  In particular, the exact examples

```
17 | A_3,  n=81  (q=4),   low jump=high jump=20,
17 | A_3,  n=132 (q=7),   low jump=high jump=37
```

contradict the claimed intrinsic purification.  After intersecting with the
interval carrier `(125.4)`, its height/no-go discussion remains informative:
the carrier already costs `(1/6+o(1))n`, first Smith divisors retain only one
selective copy, mixed low/high minors retain at most two selective copies
after universal-rank saturation, and telescoping controls total marginal
jump mass rather than diagonal overlap.  Thus Q733 contributes a barrier
analysis, not a new unconditional bound.

## 128. Reflection-to-distinct-indices does not combine with density one

For any q=1 hit `(n,p)`, reflection sends it to the second outer index

`m_p=3p-1-n`.

The same prime `p` is a q=1 hit at `m_p`, and distinct primes give distinct
`m_p`.  This initially suggests using the unconditional density theorem on
the collection of partner indices.

The quantifiers do not match.  The density theorem controls indices `m` for
which

`K(m)=#{sqrt(m)<ell<=m : ell is bad at m}`

is large enough to give linear logarithmic mass.  Each distinct partner
`m_p` is only known to have the single hit `p`; hence all of them may be
nonexceptional.  The proven first moment

`sum_(m in a dyadic shell) K(m) << X^(5/3)/log X`

does not bound the support `{m:K(m)>=1}` below the trivial `O(X)` scale.
In particular, it gives no improvement over the `O(n/log n)` number of
candidate primes attached to one fixed `n`.

Graphically, a large q=1 multiplicity at `n` is a star whose leaves are the
distinct `m_p`.  Density one forbids many high-degree vertices but does not
forbid one high-degree center with degree-one leaves.  Turning this star into
a contradiction needs a new collision/expansion statement; the existing
zero-fiber bound and exceptional-index estimate do not provide it.

## 129. Q729: exact recurrence-aware Padé boundary

Q729 audits direct q=1 interpolation after incorporating the Apéry
recurrence rather than treating the values as arbitrary data.

The fixed-window obstruction is exact.  Every bounded-shift rational
combination reduces to the two-dimensional recurrence state

`U(n)A_n+V(n)A_(n+1)`.

At a direct bad prime, the first coordinate vanishes and the second does not
(two consecutive Apéry values cannot vanish).  Thus divisibility of the new
combination is exactly divisibility of the scalar `V(n)`.  Fixed-order Ore
adjoints, recurrence residuals, and Christoffel--Darboux combinations do not
create a second selective coordinate.

The mod-`p^2` Hermite version remains rank one: the derivative term refines
the same scalar, and the Q715 computation found no simultaneous value and
index-derivative zero through `p<=5000`.

For degree-growing Newton--Padé interpolation, Q729 proves:

- the ordinary Newton interpolant is exponentially large, with no hidden
  cancellation because its binomial-basis coefficients are positive;
- any nonzero denominator with nonnegative binomial coefficients is
  impossible in an underdegree interpolation, because Apéry numbers and
  products of binomial-positive sequences remain binomial-positive;
- exact signed square-Padé searches through `n<=60` (selected splits through
  `n<=120`) found no numerator smaller than the denominator-free Newton
  choice;
- at full degree, arbitrary integer denominator values give a complete
  numerator ideal generated by explicit Lagrange-weighted Apéry values.
  After intersection with the direct lcm carrier, its large-prime support is
  exactly the original direct q=1 radical, with only `O(sqrt(n))` small-prime
  contamination.

Thus full-degree Padé is tautologically complete rather than compressive.
The one interpolation class not ruled out theoremically is a primitive
**signed**, mildly under/overdetermined pair with total degree
`H+O(1)`, `H=floor((n-1)/3)`, whose nonzero value at `n` has
subexponential height.  Q754 asks for either an explicit asymptotic
construction of this kind or a uniform exponential lower bound for every
signed split.

## 130. The signed Padé denominator-degree-one case is closed

The first genuinely signed case left by Q729 can be handled exactly.  Put

`c_k=Delta^k A_0`.

Binomial transformation of the Apéry differential equation gives the
four-step recurrence

```
(k+4)^3 c_(k+4)
 = (2k+7)(15k^2+105k+184)c_(k+3)
 + (k+3)(95k^2+570k+864)c_(k+2)
 + 48(k+2)(k+3)(2k+5)c_(k+1)
 + 32(k+1)(k+2)(k+3)c_k.
```

Its characteristic polynomial at infinity is

`(x+1)^2 (x^2-32x-32)`.

Thus the unique dominant root is `gamma=16+12sqrt(2)`, with indicial
exponent `-3/2`, and the standard full Poincare expansion gives

`c_k=kappa gamma^k k^(-3/2)(1+alpha/k+O(k^-2))`.

Let `H=floor((n-1)/3)` and take a degree-one denominator
`Q(x)=q_0+q_1 x`.  If `d_k` are the Newton coefficients of `A_x Q(x)`,
the condition that the interpolating numerator have degree at most `H-1`
is `d_H=0`.  The primitive solution is

`q_0=H(c_H+c_(H-1))/g`, `q_1=-c_H/g`,

where `g=gcd(H(c_H+c_(H-1)),c_H)`.  With

`rho_k=c_(k-1)/c_k`, `f(k)=k(1+rho_k)`,

one obtains the exact identity

`d_k=(c_H c_k/g)(f(H)-f(k))`.

The asymptotic expansion gives

`f(k)=(1+1/gamma)k+3/(2gamma)+O(1/k)`,

so `f` is eventually strictly increasing.  Hence, for every sufficiently
large `H`, all numerator Newton coefficients `d_k`, `k<H`, have the same
positive sign.  Taking `k=floor(H/2)` gives

`|P(n)| >= binom(n,floor(H/2))`

and therefore

`liminf log|P(n)|/n >= h(1/6)=0.450561...`.

So a signed Padé denominator of degree one is exponentially too large and
cannot close the direct q=1 branch.  The proof is recorded in
`Q32_SIGNED_PADE_DEGREE_ONE_NO_GO.md`.  The remaining Padé gap begins at
denominator degree two and, more importantly, at degree growing with `H`.

## 131. Q747/Q748: exact countermodels for the two newest structural routes

Q747 proves that the failure in §128 is logical, not merely quantitative.
It constructs abstract reflection-symmetric zero sets with:

- at most two zeros per prime;
- no consecutive zeros;
- the exact forced reflection-continuant factor;
- the all-cutoff loss-count and Once support laws;
- a dyadic incidence first moment as small as `O(X/log X)`;

but with

`K_1(N_t) asymp N_t/log N_t`

for infinitely many centers `N_t`.  Thus reflection, vertical zero-fiber
sparsity, density one, continuant rigidity, and the Once invariant together
still permit a maximal star.  A successful graph route needs a genuinely
off-diagonal second-neighborhood theorem: the degree of a heavy center must
propagate to incidences involving *other* primes at its leaves.  Reflection
itself supplies only the edge returning to the center.

Q748 independently validates the exactly-once theorem of §126 and sharpens
its support bookkeeping.  Above `sqrt(n)`, Once contains exactly:

- every q=1 bad prime; and
- at most one `q=p-2` bad prime.

The `q=p-2` intervals are pairwise disjoint.  Intersecting with the top-half
lcm carrier removes that prime, while all remaining small-prime
contamination has radical weight `O(sqrt(n))`.  Using the literal top-half
prime product gives an exact identity.

Q748 also derives exact pairwise-gcd/Möbius and Smith formulas for Once.  In
the diagonal loss matrix, Once is the ratio of the last two Smith invariant
factors.  This does not bound it: an arbitrary squarefree integer can occur
as precisely that ratio even for a diagonal squarefree matrix.  Homogeneous
resultants and discriminants either retain a predictable content power at
proportional Archimedean cost, or divide out the content and lose the
selective prime.  A formal radical-content model can make every top-half
candidate occur exactly once while satisfying all transition axioms.

Therefore the loss calculus gives the cleanest exact q=1 extractor but no
height compression without additional arithmetic of the actual Apéry
coefficients.

## 132. A proposed total-positivity closure of signed Padé pairs (refuted)

For the Q729 Padé system define

`M_(k,l)=binom(k,l) Delta^(k-l) A_l`.

In terms of `c_m=Delta^m A_0`,

`M_(k,l)=binom(k,l) sum_(t=0)^l binom(l,t)c_(k-l+t)`.

There is a new exact factorization.  If `B_(i,j)=binom(i,j)` is the Pascal
matrix and each row is reversed by putting `r=k-l`, then

`M_(k,k-r)=sum_(m=r)^k binom(k,m)c_m binom(m,r)`

and hence

`M_(k,k-r)=(B diag(c) B)_(k,r)`.                     `(132.1)`

Thus the row-reversed triangle is totally nonnegative.  The unresolved
question is whether undoing this *row-dependent* reversal preserves the
flag minors used by Padé, or even makes the original matrix totally
nonnegative.

The initial evidence looked strong:

- the reproducible script checks every minor only through order 3 with rows
  below 14 and columns below 8; the earlier order-4 claim was not encoded in
  the script and must not be cited as verified;
- 2,000 deterministic random minors of orders 2 through 9, with rows below
  55 and columns below 18, contain no negative determinant;
- for every denominator degree `1<=b<=6` and `H<=40`, the primitive kernel
  alternates in sign and every numerator Newton coefficient has one common
  sign.

The executable audit is `q32_pade_total_positivity.py`.

If the required flag-minor positivity holds for all `H,b`, choose the common
sign so every numerator Newton coefficient is positive.  Then

`P(n)>=P(H)=A_H Q(H)>=A_H`,

because `Q(H)=P(H)/A_H` is a positive integer.  Since
`log A_H asymp H asymp n`, **every** signed square Padé pair would be
exponentially large, including the dangerous regime `deg P=o(n)`.  This
would completely close the last Padé class left by Q729.

Q762 asks for an exact total-positivity or flag-minor proof (or the first
counterexample); Q763 separately audits the proportional-degree saddle and
primitive-normalization problem.

### Sol-max audit correction

The proposed total positivity is false.  The first structured
counterexample occurs at `H=6`, denominator degree `b=5`, and numerator
degree `a=1`.  With rows `1,...,6` and columns `0,...,5`,

`det M[1,...,6;0,...,5] = -2421161406987687811000`.

The primitive Padé kernel is

```
[1267742334817618036530,
 -1016615029261082117035,
  782522391112049501080,
 -564612223198400602008,
  362069276292219883472,
 -174113775109881764984]
```

and the numerator has Newton coefficients

`[1267742334817618036530, -12105807034938439055]`.

Thus the required flag minor is negative and the numerator coefficients
have opposite signs.  The random audit missed a rare consecutive minor, and
the original fixed-degree Padé scan excluded the boundary `H=b+1`.  The
row-reversed factorization (132.1) remains exact, but it gives no sign
control after the row-dependent reversal.  This removes the claimed
all-degree no-go theorem.  The sign failure occurs precisely in the
low-numerator-degree regime, so that regime must be audited by height rather
than positivity.

## 133. Sol-max audit of the 23 July proof claims

All conclusions produced under the earlier main model are provisionally
untrusted until independently rederived.  The first audit pass gives:

1. **Refuted:** the original matrix `M_(k,l)` is totally nonnegative, and
   therefore the proposed one-line closure of every signed Padé split.
   The exact order-six counterexample is recorded above and in
   `q32_pade_total_positivity.py`.
2. **Corrected and repaired:** the degree-one-denominator asymptotic
   argument wrote
   `f(k+1)-f(k)=1+gamma^(-1)+O(k^(-2))` from an
   `O(k^(-1))` expansion of `f`.  Without a fuller expansion the justified
   error is only `O(k^(-1))`; this still tends to zero and is enough for
   eventual positivity.  The missing asymptotic justification is supplied
   directly by Edgar's positive moment density:
   `c_k=int_0^(17+12sqrt2)(x-1)^k phi(x) dx`, and the right-endpoint
   Frobenius expansion followed by Watson's lemma gives the required full
   expansion with dominant base `16+12sqrt2` and exponent `-3/2`.
   Independently, the transformed recurrence is exact through `k=145` and
   `f(k)` is strictly increasing through `k=149`.  Thus the
   degree-one-denominator no-go is restored; the finite checks are
   corroboration, not part of its asymptotic proof.
3. **Independently rechecked algebraically:** the arbitrary-cutoff Lucas
   identity passed 2,984 exact polynomial cases; the q=1 cutoff profile
   passed 5,110 cases; and the q=1 mod 3 alias has the exact examples
   `(q,n,h)=(4,81,20),(7,132,37)`.
4. **Strongly verified but citation boundary still open:** the
   Frobenius-twisted Franel Dwork congruence passed 22,448 tests including
   arbitrary `m`, and the tail inverse is integral for
   `0<=J<=30`, `d<=200`.  The document's coefficient-ring D3 paragraph is
   only a proof sketch, not yet a self-contained proof or a verbatim
   application of the cited scalar theorem.  Until that extension is
   written carefully, `Q32_FRANEL_TAIL_DWORK_PROOF.md` is classified as a
   highly supported proof draft rather than a closed proof.

This audit distinction is essential: exact identities and explicit
counterexamples survive; numerical searches are evidence only; literature
theorems survive only after their hypotheses match the coefficient ring and
Frobenius used here.

## 134. Targeted Padé minimality after total positivity fails

The failure of total positivity does not yet produce a useful small
certificate.  A fresh exact scan used every split `a+b=H`, the target
`n=3H+1`, and primitive integral normalization.  For every
`4<=H<=30`, the smallest absolute numerator value among all splits is still
the denominator-free Newton interpolant (`a=H,b=0`).  Thus the original Q729
phenomenon survives well beyond its first scan even though the proposed
explanation by common coefficient signs is false.

Two cautions:

1. This is exact finite evidence, not an asymptotic theorem.
2. The ordinary Newton value does **not** divide the other primitive Padé
   numerator values: direct checks for every `4<=H<=20` fail divisibility
   for every nontrivial split.  Hence the observed minimality cannot be
   proved by a hidden elementary divisibility relation.

The new precise Padé target is therefore either:

- prove the evaluation inequality
  `|P_(H,a)(3H+1)| >= N_H(3H+1)` for every primitive split and all
  sufficiently large `H`; or
- prove only a uniform exponential lower bound for the left side.

The low-numerator-degree regime must use a primitive-height/Smith argument,
not sign regularity.  Q776 asks ChatGPT for this targeted theorem after
including the exact order-six counterexample.

The exhaustive minimality scan has now been extended through `H=40`; the
ordinary Newton split remains the unique minimum in every case
`4<=H<=40`.  This remains evidence only.

## 135. Q764: ray amplification is exact, but exponent and height scale together

Q764 tested the cross-index ray

`N_s=n+s p`

for a fixed q=1 bad prime `p`.  The low folded digit stays fixed, while the
proved loss-count formula gives multiplicity

`min(s+1,p-2-s)`.

Thus one really can manufacture arbitrarily many local copies of the same
bad prime within the uniform linear range.  Q764 also constructs:

- a single lcm-based global selector for one ray;
- the product over many rays;
- a coefficient determinant and its full Smith/Fitting determinantal
  divisor.

All three have the same obstruction.  A target exponent `t` costs a
universal carrier of logarithmic height `Theta(t n)`.  Multiplying rays
raises both exponent and height by the same quadratic factor.  The natural
coefficient matrix is lower triangular, so its determinant literally
collapses to a product of scalar carriers; there is no determinant
cancellation hidden in that construction.  Resultants and discriminants
have the same homogeneous-degree accounting.

The structural reason is that every propagated row is the same local scalar
`A_j mod p` times a known block polynomial.  The ray repeats one condition;
it does not create independent conditions on `A_j`.

This route is therefore recorded as an **exact amplification/no compression**
result.  Reopen only with an order-`t` cross-cutoff identity whose local
intersection multiplicity is `t` but whose single nonzero boundary has
height `o(t n)`.  Q777 asks specifically whether finite differences,
discrete Wronskians, confluent determinants, or factorial/Vandermonde
normalization can realize this escape already for `t=2`.

## 136. Q754: exact Smith quotient and an adjacent-Padé second filter

Q754 rewrites every signed square Newton--Padé split `a+b=H` as an
explicit integral finite-difference kernel problem.  If `M^(a,b)` is the
`b x (b+1)` denominator matrix, its primitive kernel vector is the
cofactor vector divided by the last determinantal divisor
`delta_(a,b)`.  Appending the evaluation row at `n` gives the exact
identity

`|P_(a,b)(n)| = |det Mhat^(a,b)(n)| / delta_(a,b)`.

Thus an Archimedean saddle estimate of the augmented determinant is not
enough: a proof must simultaneously control the Smith saturation factor.
For fixed excess `s`, the exact replacement is the image ideal of
evaluation on the saturated kernel lattice; it is the gcd of the
`s+1` values obtained from any saturated integral kernel basis.  This is
a useful exact reduction, not a bound for that ideal.

There is also a new adjacent-table identity.  For primitive square pairs
`(P_b,Q_b)` and `(P_(b+1),Q_(b+1))`,

`Omega_b(x)=P_b(x)Q_(b+1)(x)-P_(b+1)(x)Q_b(x)`

vanishes at every interpolation node `0,...,H` and has degree at most
`H+1`.  Hence

`Omega_b(x)=kappa_b binom(x,H+1)`

for an explicit integer `kappa_b` (equivalently, the coefficient of the
monic falling factorial is `kappa_b/(H+1)!`).  At a direct prime
`p=n-j>2H`, the
value identity together with its derivative at `x=j` shows that any
*good* direct prime common to the two adjacent numerator certificates
must divide `kappa_b`.  A fixed adjacent chain therefore separates the
target bad primes from good-prime pollution up to the prime support of
the corresponding finite list of `kappa` factors.

This does **not** yet prove that the branch-purified gcd has
subexponential height: neither `kappa_b` nor the saturated evaluation
ideal is controlled uniformly.  It is also fully consistent with the
order-six negative minor in section 133; Q754 explicitly explains why
the finite-difference matrix need not inherit Hankel total positivity.
The next precise task is to compute or bound `kappa_b` after primitive
Smith normalization and test whether a fixed adjacent chain has
subexponential purified gcd.

## 137. Exact adjacent-`kappa` audit: striking small gcds, but the edge is tautological

The executable audit `q32_adjacent_pade_kappa.py` constructs every
primitive square pair, verifies

`P_b Q_(b+1)-P_(b+1)Q_b
 = kappa_b binom(x,H+1)`

at exact integer points, and also verifies the leading-coefficient formula

`kappa_b=binom(H+1,b+1) p_(b,H-b) q_(b+1,b+1)`,       `(137.0)`

where the last two factors are the top Newton coefficients of `P_b` and
`Q_(b+1)`.  It measures adjacent and three-term numerator gcds at
`n=3H+1` and passed every split for `2<=H<=30`.

There is an equivalent determinantal form.  Let

`E_b=det M[H-b,...,H; 0,...,b]`

and let `delta_b` be the maximal-minor gcd of the `b x (b+1)` kernel
matrix for type `(H-b,b)`.  Cofactor expansion gives

`p_(b,H-b)=+/- E_b/delta_b`,
`q_(b+1,b+1)=+/- E_b/delta_(b+1)`,

and hence

`|kappa_b|
 = binom(H+1,b+1) E_b^2/(delta_b delta_(b+1))`.       `(137.0a)`

This exact square/Smith formula may help audit prime support, but by itself
it supplies no small-height bound.

The numerical signal is strong but not yet a theorem.  The minimum
adjacent gcd has `log(gcd)/n<=0.219` throughout this range and is often
`1`, `5`, or another tiny integer; a three-term chain is no larger in the
reported minima.  The minimizing pair is usually the edge `b=0,1`.

That edge can be simplified exactly.  If

`c_H=Delta^H A_0`,
`g=gcd(H(c_H+c_(H-1)),c_H)`,

then

`Q_1(x)=(H(c_H+c_(H-1))-c_H x)/g`

and comparison of top Newton coefficients gives

`kappa_0=-(H+1)c_H^2/g`.

Consequently

`gcd(P_0(n),P_1(n))
 = gcd(P_0(n), kappa_0 binom(n,H+1))`.                `(137.1)`

This explains why the adjacent gcd is usually small, but it does not
compress the direct bad primes.  For a candidate `p=n-j` with
`0<=j<=H<p`, Lucas already makes `binom(n,H+1)=0 mod p`.  Hence every bad
prime can enter the right side of `(137.1)` through the same automatic
binomial carrier, without dividing `kappa_0`.  Proving that the relevant
part of `(137.1)` is subexponential is therefore another formulation of
the moving-zero problem, not a consequence of the adjacent identity.

The useful surviving content of Q754 is narrower: `kappa_b` filters
*good* common candidate primes.  To become a proof, a fixed chain must
also supply a non-tautological bound on the part supported by the
binomial carrier.  Small computed gcds alone do not provide that bound.

## 138. Q763: proportional Padé reduces to an `H`-scale signed asymptotic

Q763 does not prove a positive lower bound or construct a zero-rate split.
Its main exact clarification is nevertheless useful.  Let `D_l` be the
signed maximal cofactors of the square denominator matrix, let
`g=gcd_l D_l`, and assume `D_b!=0`.  Normalize the rational kernel by
`qhat_b=1`.  Then

`q_primitive=(D_b/g) qhat`,
`P_primitive=(D_b/g) Phat`.

Since `D_b/g` is a nonzero integer, a positive exponential lower bound for
the projectively normalized value `|Phat(n)|` automatically survives
primitive Smith normalization.  Smith saturation remains essential for
constructing a *small* certificate, but it is not an obstruction to a
projective no-go theorem.  This elementary scaling statement was checked
directly.

Using Edgar's moment density, Q763 obtains a uniform endpoint asymptotic
for the matrix entries over the full triangular range.  Replacing the
endpoint power-law factor by its leading weighted-Pascal symbol gives an
exactly solvable Krawtchouk-type model whose projective evaluation rate is
strictly positive on every compact proportional range
`b/H in K subset (0,1)`.

This is a model calculation, not a transfer theorem for the Apéry matrix.
For the actual cofactors, the common `H^2` free energy cancels in their
ratio.  The endpoint `k^(-3/2)` factor, alternating subset phase, discrete
entropy, and possible Stokes cancellation all affect the next `H` scale,
which is exactly the scale of `log|P(n)|`.  Entrywise asymptotics and raw
determinant bounds therefore stop one order too early.

The precise missing result is a uniform strong asymptotic, through the
`H` term and with nonvanishing, for the signed Krawtchouk deformation.
Until such a theorem is proved, the proportional Padé route remains open;
the solvable positive-rate model is evidence favoring a no-go, not an
unconditional conclusion.

## 139. Q777: cross-ray Hasse amplification exists, but the canonical determinant is too large

Q777 follows the ray amplification of section 135 into the first genuine
exterior-power construction.  For order `t`, it forms the Hasse-jet
determinant of the `t` ray rows at derivative orders

`0,J,2J,...,(t-1)J`.

At `J=p>2t`, the all-cutoff block profile makes every entry contain the
same factor `A_j mod p`.  Thus a q=1 bad prime contributes `p^t`, while
the constant Hasse determinant is a `p`-unit for a good prime.  This is
real selective amplification rather than a factorial artifact.

The full polynomial determinant nevertheless has a nonzero leading
coefficient.  Its top term is a product of central-binomial/cutoff factors
times a positive Pascal minor, so its logarithmic height is at least

`(log 4)(t n + J t(t-1)/2) - O(t log(n+tJ))`.

The standard normalizations do not change this conclusion:

- ray finite differences are a lower-unitriangular row operation;
- the ordinary ray Vandermonde and divided differences remove only
  `t`-dependent factors;
- ordinary derivatives add universal factorial powers, and Hasse
  normalization removes exactly those powers;
- the outer-index factor `J^(t(t-1)/2)` is not an integral content divisor,
  since the constant term is nonzero modulo a good prime `J=p`.

Linear WZ telescopes are also locally rank one: after saturation they retain
at most one selective factor unless a genuinely new congruence modulo
`p^2` controls the quotient rows.  Passing to the `t`-th exterior power
raises selective order and determinant degree together.

Hence the full cross-ray Wronskian route is recorded as
**amplification without compression**.  The narrow survivor is a selected
constant Hasse minor, beginning at order two, where signed cancellation
could conceivably reduce the aggregate height.  Q777 does not prove such a
bound.  Reopen only after an independent exact audit finds a selected minor
with `p^t` selectivity and total logarithmic height `o(tn)`, or after a new
mod-`p^t` supercongruence supplies extra local order without exterior
degree.

## 140. Q778: second-neighborhood expansion needs cross-characteristic coherence

Q778 classifies every off-diagonal q=1 incidence at a reflection leaf.  If
`p` is a star prime at `n`, its leaf is

`m_p=3p-1-n`.

If a different prime `ell`, with folded zero `k`, is bad at `m_p`, then
exactly one of the two affine equations holds:

`3p=n+1+ell+k`,

`3p=n+2ell-k`.

For a fixed pair `(ell,k)`, these equations hit at most two leaves.  The
original folded zero that creates the spoke through `n` cannot create a
new off-diagonal leaf incidence; another zero of `ell` is required.
Moreover,

`m_p-m_q=3(p-q)`,

so the leaf Vandermonde is coprime to every target top-half prime.  The
obvious product-of-differences and discriminant certificates therefore
contain none of the desired radical, rather than merely having excessive
height.

Q778 also gives a formal disjoint-star incidence model satisfying all
prime-local zero-set consequences currently in use: reflection, no
consecutive zeros, continuant/Casoratian identities, Lucas blocks, and
vertical sparsity.  It permits degree `asymp n/log n` with degree-one
leaves.  Each individual spoke is realizable in a recurrence module when
the initial state may depend on the prime.  What the model deliberately
does not realize is the same global Apéry initial state `(1,5)` in every
residue characteristic.

Thus no superlinear expansion follows from the present prime-local
toolkit.  A sufficient reopening theorem is an off-diagonal
second-neighborhood/correlation lower bound with exponent
`delta>2/3`; the shell incidence upper bound would then force the star
degree to be sublinear enough for the moving-zero estimate.  Any proof of
such a theorem must exploit cross-characteristic arithmetic coherence of
the fixed initial state, not reflection geometry alone.

## 141. Q786: fixed adjacent Padé chains saturate on the binomial carrier

Q786 confirms and generalizes the exact adjacent formula of sections
136--137.  Let

`tau_b=det M[H-b,...,H;0,...,b]`

be the solid Newton minor and `delta_b` the last determinantal divisor of
the denominator kernel matrix.  In the normal full-rank case,

`|kappa_b|
 = binom(H+1,b+1) tau_b^2/(delta_b delta_(b+1))`.

This is the independently checked Smith-normalized square formula already
implemented in `q32_adjacent_pade_kappa.py`.  Dodgson/Toda identities for
the raw solid minors introduce shifted minors, while consecutive
`kappa` ratios retain a two-step Smith ratio; there is no closed small
recurrence for the `kappa_b`.

For any fixed chain,

`G_(b,r)(n)=gcd(P_b(n),...,P_(b+r)(n))`

divides

`binom(n,H+1) gcd(kappa_b,...,kappa_(b+r-1))`.

The `kappa` gcd supports the extra common zeros coming from good candidate
primes.  After those are removed, the residual prime support is exactly the
original direct bad-prime set.  Every direct candidate divides
`binom(n,H+1)` exactly once, and a local `p`-adic model realizes valuation
one in every numerator while every `kappa_b` remains a unit.  Exact
saturation by the binomial factor therefore removes the entire forced bad
prime; products and Plücker/Dodgson identities only repeat this balance.
The carrier itself has positive linear entropy:

`log binom(n,H+1) asymp n h(1/3)`.

Ordinary resultants reintroduce products of the interpolation values and
have quadratic raw height.  Thus the universal adjacent-identity/Smith
calculus cannot explain the small computed adjacent gcds.

The narrow empirical route is still logically open: prove directly, for
one fixed interior adjacent pair, that its large-prime gcd has logarithm
`o(H)`.  Such a result would be a new arithmetic theorem at least as strong
as the direct q=1 moving-zero estimate; it is not implied by the adjacent
identity.  Reopen this route only with input about the actual Apéry
holonomy/moment sequence beyond interpolation, Smith normalization,
Plücker relations, and the universal binomial carrier.

## 142. Canonical repository consolidation and next hostile batch

The active campaign has been consolidated into the canonical GitHub
repository `xiangyazi24/Ramanujan_Challenge`.  Commit `710a25d` contains
the historical questions and answers, proof drafts, failed-route log, and
reproducible scripts.  The Mini-synchronized Lean commits were preserved
without modification; Lean is not the current research priority.  Commit
`6aea822` adds the exact new prompt texts.  The stale Q-series checkout and
the non-Git snapshot remain only as recovery copies.

The next ChatGPT batch is deliberately split across seven nonoverlapping
paper-proof targets:

- Q800: exploit the common Apéry initial state `(1,5)` across residue
  characteristics to prove cross-prime expansion;
- Q801: compute the selected order-two constant Hasse minor and determine
  whether signed cancellation beats the ray height barrier;
- Q802: attack one fixed interior adjacent Padé gcd using actual Apéry
  arithmetic beyond universal interpolation identities;
- Q803: prove or refute the missing `H`-scale strong asymptotic for the
  signed proportional Krawtchouk deformation;
- Q804: derive the mod-`p^2` ray-quotient formula needed to obtain a second
  selective factor without exterior degree;
- Q806: search for, or sharply rule out, a pointwise moving-linear
  prime-divisor theorem for
  `sum_(j<=H, n-j prime, n-j|A_j) log(n-j)`.
- Q817: audit from first principles whether the fixed-outer-index quantity
  `A_n mod p` is a Frobenius trace or Hasse invariant of one fixed
  `n`-dependent fiber, and whether any available uniform Frobenius sieve
  survives the growing height and the moving trace multiplier.

Every answer remains untrusted until its load-bearing identity is
independently derived or exactly checked.  Failed attempts must be recorded
with the precise saturation mechanism and the condition under which they
could be reopened.

## 143. Q800: the fixed initial state gives an exact star lattice, not expansion

Q800 supplies a clean integral formulation of the cross-characteristic
coherence that Q778 deliberately omitted.  For arbitrary initial data
`u_0=x,u_1=y`, put

`X_m=(m!)^3 u_m`.

Then the Apéry recurrence becomes

`X_(m+1)=P(m)X_m-m^6X_(m-1)`,

where `P(m)=34m^3+51m^2+27m+5`.  If
`X_m=r_m x+s_m y`, direct transfer-matrix induction gives

`det((r_m,s_m),(r_(m+1),s_(m+1)))=(m!)^6`.          `(143.1)`

Hence for a candidate prime `p>m` the row `(r_m,s_m)` is nonzero modulo
`p`.  At an actual folded Apéry zero `p|A_m`, that row annihilates
`(1,5)`, so it is necessarily a nonzero multiple of `(-5,1)`.  Therefore
an arbitrary integral initial state vanishes at that same local zero if
and only if

`y=5x (mod p)`.                                      `(143.2)`

Let `R_1(n)` be the product of the q=1 star primes for the fixed outer
index `n`.  Intersecting `(143.2)` over those distinct primes gives the
exact simultaneous lattice

`Lambda_n=Z(1,5)+Z(0,R_1(n))`,

with Smith invariants `(1,R_1(n))`.  In particular, every independent
comparison state synchronized with all star zeros has determinant with
`(1,5)` divisible by `R_1(n)`, and its Euclidean norm is at least
`R_1(n)/sqrt(26)`.  More generally, any finite family whose synchronized
prime subsets cover the star carries the full radical in the product of
its comparison determinants.

This is exact, but it is a saturation theorem rather than a height saving.
The fixed Apéry state supplies the short first lattice direction; the
unknown radical is precisely the second successive minimum.  The rational
companion cannot supply a second direction because the Apéry Casoratian is
a unit at every low candidate prime.

Q800 also refines the Q778 countermodel by a constant integral conjugacy:
one may CRT-lift prime-local gauges so that a single coordinate vector
realizes any prescribed finite collection of spokes.  This correctly
rules out arguments invariant under constant gauge and depending only on
monodromy or the existence of a fixed state.  It does **not** model the
canonical Apéry scalar sequence at bounded height: the conjugating matrix
and transformed observable encode the prescribed radical.  The remaining
load-bearing feature is therefore the canonical low-height integral gauge,
not fixed-state coherence by itself.

The route can be reopened in either of two genuinely new ways:

1. construct a canonical independent synchronized section, or a finite
   family of such sections, whose determinant product has logarithm
   `o(n)`; or
2. prove an actual cross-characteristic off-diagonal incidence lower bound
   for the canonical Apéry coefficients, for example the Q778/Q800
   condition
   `E_off(n)>=K_1(n)^(1+delta)n^(-o(1))` for some `delta>2/3`.

Q822 asks whether the canonical Picard--Fuchs/modular integral structure
has an endomorphism or quantitative row constraint that survives this
gauge audit.  Until such input appears, Q800 closes the fixed-state,
finite-comparison, product-formula, and gauge-invariant monodromy routes
without advancing the unconditional bound.

## 144. Q804: the first Frobenius jet is a second state, not a second zero

Q804 derives the first genuinely outer-index mod-`p^2` formula for the
cutoff row at `J=p-1`.  Put

`S_q(p,r)=[c^0]T_(qp+r,p-1)(c)`

and

`Acal_p(X)=sum_(m=0)^(p-1) L(X,m)F_m`.

Since `S_q(p,r)=Acal_p(r+qp)` and `Acal_p(r)=A_r`, ordinary Taylor
expansion over `Z_(p)` gives the exact two-state congruence

`S_q(p,r)=A_r+qp Dcal_(p,r) (mod p^2)`,             `(144.1)`

where `Dcal_(p,r)=Acal_p'(r) mod p`.  On a bad digit, writing
`u_r=A_r/p mod p`,

`S_q(p,r)/p=u_r+q Dcal_(p,r) (mod p)`.              `(144.2)`

Thus two ray rows recover the two first-jet coordinates `(u_r,Dcal)`.
Badness kills the zeroth reduction `A_r mod p`; it imposes no vanishing
condition on either quotient coordinate.  The first ray difference has
one universal factor `p` and quotient `Dcal`; the second difference has a
universal factor `p^2`.  Saturating those universal powers leaves no
selective amplification.

The full coefficientwise calculation has the same structure.  For
`k=ap+b`, Q804 defines the integral first Lucas defect

`lambda=[L(qp+r,ap+b)-L(q,a)L(r,b)]/p`

and the integral binomial-translation defect

`mu=[B_p(ap+b;c)-c^(ap)F_b]/p`.

Substitution gives

`T_(qp+r,p-1)=A_r Q_q(c^p)+pE_(q,r,p)(c) (mod p^2)`, `(144.3)`

with an explicit finite double sum for `E`.  In the no-carry range the
Lucas defect is

`lambda=L(q,a)L(r,b) Phi (mod p)`,

`Phi=(q+a)H_(r+b)-(q-a)H_(r-b)-2aH_b`.

The translation defect is the parameter derivative of
`binom(ap+b,m)` and includes the simple-zero terms with `m>b`.  Those
terms are important: truncating the derivative to `m<=r` would give a
false formula for `Dcal`.

Independent exact checks completed here:

- 22,731 no-carry tuples for `p=3,5,7,11,13,17,19,23` satisfy the
  displayed harmonic Lucas-defect formula;
- 13,205 binomial-translation tuples over the same primes satisfy Q804's
  two-branch derivative formula;
- the constant-row values at `(p,r)=(5,1)` are exactly
  `1,188,265` and `159,262,225`, and their divided first difference is
  `2 mod 5`, not zero;
- 1,023 tuples with `p<=19`, `0<=q,r<p` satisfy the published
  Gessel--Lucas congruence modulo `p^2` after rational denominators are
  interpreted in `Z_(p)`;
- for the stronger counterexample `A_3=5*17^2`, one has
  `A_20=85 mod 17^2` and `A_13=-17 mod 17^2`, so propagation and
  reflection both fall back to valuation one.

Straub's Theorem 5.1 (published as Theorem 1.3) gives

`A_(qp+r)=A_r A_q+pq A'_r A_q (mod p^2)`             `(144.4)`

for `0<=r<p`; the author's corrigendum explicitly repairs the published
misprint in this digit range.  Formula `(144.4)` has the same two-state
obstruction.  The source and corrigendum were checked against the author's
preprint page and arXiv preprint, rather than accepted only from Q804.

The order-two constant Hasse minor still has selective `p^2`, but only
because it is an exterior product of two rows already divisible by `p`.
Q804 supplies neither a third factor nor an improved Archimedean bound.
Therefore the proposed new mod-`p^2` linear escape is closed.

Q824 is testing the next and final natural version of this idea: compute
the second jet modulo `p^3` and either find a jet invariant whose selective
order beats homogeneous degree, or prove an all-order jet-algebra
saturation theorem.  Reopen the linear-ray route only if a higher jet
creates selective order not paid for by the same universal finite
difference or exterior degree.

## 145. Q803: exact checkerboard reduction verified; fixed-strip proof still under audit

Q803 rewrites the signed Newton--Padé cofactors as finite checkerboard
orthogonal-polynomial ensembles.  For

`w_H(s)=A_s/[s!(H-s)!]`

and

`Z_m^-(w)=sum_(|S|=m)(-1)^(sum S) Vand(S)^2 product_(s in S)w(s)`,

the proposed projective formula is

`|Phat(n)|=(n)_(H+1)/b! *
 |Z_(b+1)^-(w_H/(n-s))/Z_b^-(w_H)|`.                `(145.1)`

This is not merely a model.  Q803's own exact program verifies the
denominator cofactor through `H=8`.  The independent script
`research/scripts/q32_checkerboard_audit.py` constructs the original
Padé kernel directly from

`M_(k,l)=binom(k,l) Delta^(k-l)A_l`

and verifies the full numerator/denominator ratio `(145.1)` in all 65
cases `1<=H<=10`, `0<=b<=H`.  The signed-to-positive ratio is also the
standard finite determinantal-process identity

`Z_m^-/Z_m^+=det(I-2 K_(m,w) chi_odd)`.

Thus the exact proportional problem really contains an order-`H`
checkerboard full-counting-statistics determinant.  Positive equilibrium
theory controls the unsigned free energy but does not control this
parity determinant, which can vanish precisely when the odd-site
compression has eigenvalue `1/2`.  This is a sharper exact statement of
the signed `H`-scale obstruction from Q763.

Q803 further claims that every fixed denominator degree `b` has the
asymptotic

`Phat_(H,b)(3H+1)
 ~(-1)^b (lambda/gamma)^b c_(H-b) binom(3H+1,H-b)
   (1-1/(2gamma))^(-b-1)`,                           `(145.2)`

where `lambda=17+12sqrt(2)`, `gamma=lambda-1`, and
`c_k=Delta^k A_0`.  Exact computations for `b<=8` and
`H=40,80,120,180` strongly support `(145.2)`: at `H=180` the absolute
exact/predicted ratio lies between `0.991342844761` and
`1.000041243513`.

The mechanism is transparent.  For fixed `b`, the normalized cofactor
kernel approaches the weighted-Pascal polynomial with consecutive roots
`H-b+1,...,H`.  At `k=H-b-r` its normalized value tends to

`(-1)^b (lambda/gamma)^b binom(r+b,b)`.

The endpoint ratios contribute `(2gamma)^(-r)`, and summing
`binom(r+b,b)(2gamma)^(-r)` gives the last factor in `(145.2)`.

Nevertheless the proof as returned is not yet theorem-grade.  Its
displayed equations `(5.4)`--`(5.8)` were absent from the delivered
answer, and it only asserts, without a written bound, the load-bearing
facts that:

1. sufficiently many terms of the Poincare expansion survive the
   consecutive-row determinant cancellation uniformly for
   `r<=C_B log H`;
2. the rational cofactor kernel has a polynomial-in-`H` bound uniform in
   the remaining `r`; and
3. this bound and a uniform ratio bound for `c_k` make the tail
   negligible.

These facts are plausible and finite-dimensional for fixed `b`, but the
claim “fixed `b` is closed” is withheld until these three estimates are
written explicitly.  Q825 failed at delivery; Q828 is performing that
hostile audit.

The reciprocal complement formula supplies a second exact boundary.  For
fixed numerator degree `a=1`, put

`B_r(H)=sum_(s=0)^H (-1)^s binom(H,s)s^r/A_s`.

Then `(145.1)` reduces exactly to

`|Phat_(H,H-1)(n)|
 =|H(nB_0-B_1)/(B_0B_2-B_1^2)|`.                   `(145.3)`

The independent script verifies `(145.3)` for `2<=H<=10`.  This isolates
the smallest unresolved reciprocal boundary as three alternating binomial
transforms of `1/A_s`.  A finite Poincare expansion of `1/A_s` does not
control them: its unsigned error has base `1+lambda^(-1)`, while the
expected checkerboard answer has base `1-lambda^(-1)`.  Reopen this
boundary only with analytic continuation/strong asymptotics for the
reciprocal-Apery generating function, an exact recurrence for the three
transforms, or a direct nonvanishing/Turan estimate for the determinant
in `(145.3)`.  Q826 did not remain attached to a live tab; Q830 and Q831
replace it with the analytic-interpolation and primitive-lattice attacks.

## 146. Q801: the selected order-two Hasse minor is exactly the target squared

Q801 completes the audit of the first exterior-power construction.  With

`C_d(N,J)=[c^d]T_(N,J-1)(c)`,

put

`D_2(n,J)=C_0(n,J)C_J(n+J,J)-C_J(n,J)C_0(n+J,J)`. `(146.1)`

For a top-half prime `J=p>4`, write `n=p+r` and
`j=min(r,p-1-r)`.  The two cutoff rows reduce coefficientwise to
`A_j Q_1(c^p)` and `A_j Q_2(c^p)`.  Taking the columns of degrees `0`
and `p` gives the exact unit test

`D_2(n,p)=4A_j^2 (mod p)`.                           `(146.2)`

If `p|A_j`, every entry of both rows is divisible by `p` over the
integers, so `p^2|D_2(n,p)`.  Conversely `(146.2)` makes `D_2(n,p)` a
`p`-unit whenever `p` is good.  Therefore

`product_(n/2<p<=n) gcd(p^2,D_2(n,p))=R_1(n)^2`      `(146.3)`

after omitting only the fixed primes at most `3`.  This proves that the
prime-only “normalized aggregate” is not merely comparable to the
unknown q=1 radical: it is exactly its square.  Any `o(2n)` bound for
the left side is precisely the original `o(n)` moving-zero theorem
multiplied by two.

Q801 also expands `(146.1)` into four single sums and a double sum.  Its
Legendre determinant kernel is strictly positive in the top-half range,
but the Hasse-column weight contains a moving checkerboard sign.  Hence
the full transform has no sign supplied by total positivity or Andreief.
There is no universal `J`-content to remove, since `(146.2)` is nonzero
at every good prime.  Desnanot--Jacobi, Christoffel--Darboux, WZ
telescoping, Smith saturation, and the lcm over composite cutoffs only
change representations or add contamination; they cannot reduce the
selective-degree/height ratio.

This agrees with the independent local calculation made before Q801
returned.  The selected minor is a valid order-two amplifier but a
tautological compression.  Reopen this route only if one obtains either:

1. a genuinely uniform diagonal recurrence/asymptotic for the alternating
   double sum that proves the original radical bound outright; or
2. a mod-`p^2` quotient-row congruence producing selective order greater
   than exterior degree.

The first option is not a weaker subproblem, and Q804 already shows that
the first jet does not provide the second.

## 147. The reciprocal `a=1` boundary splits into a tiny projective value and a huge primitive lattice

The three moments in `(145.3)` simplify further.  Define

`U_H=sum_(s=0)^H (-1)^s binom(H,s)/A_s`.

Elementary binomial identities give

`B_0=U_H`,

`B_1=H(U_H-U_(H-1))`,

`B_2=H^2U_H-H(2H-1)U_(H-1)+H(H-1)U_(H-2)`.

Consequently, for `n=3H+1`,

`Phat_(H,H-1)(n)
 =+-[(2H+1)U_H+H U_(H-1)]`

`   /[(H-1)U_HU_(H-2)+U_HU_(H-1)-H U_(H-1)^2]`.    `(147.1)`

The rational numerator line itself is, up to scale,

`P_H(x)=H(U_H-U_(H-1))-xU_H`,                       `(147.2)`

with root

`x_H=H(1-U_(H-1)/U_H)`.                             `(147.3)`

The script now clears all coefficients of `(147.2)` and of the Newton
interpolant of `P_H(s)/A_s`, divides their total gcd, and verifies the
result against the original primitive Padé kernel for every `2<=H<=10`.
Thus it computes the actual primitive pair, not only the projective
quotient.

The distinction is dramatic.  The exact projective rates remain small
and irregular (about `0.01647` at `H=300`), but the primitive numerator
has:

| `H` | bits of `|P_primitive(3H+1)|` | `log|P|/(3H+1)` |
|---:|---:|---:|
| 40 | 3,686 | 21.10956 |
| 80 | 15,158 | 43.59374 |
| 100 | 23,978 | 55.21585 |
| 120 | 34,868 | 66.94797 |
| 140 | 47,730 | 78.58300 |

This is numerical evidence of quadratic logarithmic height, not a
theorem.  It does show that a zero projective rate would not by itself
produce a useful integral certificate.

There is a simple sufficient sign lemma.  If `U_(H-1)U_H>0`, then
`x_H<H`.  Since the primitive integer-valued denominator makes

`A_s | P_primitive(s)` for every `0<=s<=H`,

and a nonzero linear polynomial cannot vanish at both `H-1` and `H`,
one gets

`|P_primitive(3H+1)|>=A_(H-1)`.                     `(147.4)`

Thus eventual consecutive same sign would already close `a=1`
exponentially, without any lower bound for `(147.1)`.  Exact computation
finds `U_H>0` for every `84<=H<=2000`, but this is not yet a proof and
the earlier signs occur in long blocks.

A new analytic route is available from the primary source Golyshev--Zagier,
*Interpolated Apéry Numbers*.  It proves that

`A(z)=sum_(k>=0) binom(z,k)^2 binom(z+k,k)^2`

is entire, agrees with `A_s` at nonnegative integers, and satisfies
`A(z)=A(-z-1)` together with an explicit inhomogeneous Apéry recurrence.
Therefore `1/A(z)`, rather than the ordinary reciprocal generating
function, is the natural input to a Nörlund--Rice integral for `U_H`.

A high-precision, non-rigorous search found the simple zeros

`z=0.14598677311801086651 +- 0.58242614024485315072 i`

with residual about `1.4e-26`, plus the conjugate and
`z -> -z-1` partners.  Poles of `1/A(z)` at such zeros would contribute
oscillatory powers `H^z` to the Rice expansion.  This explains why a
finite Poincare expansion of `1/A_s` predicts the wrong exponential
scale and suggests that the long sign blocks are logarithmic
oscillations.

This is promising but conditional on load-bearing analytic facts:

1. locate all zeros of `A(z)` in the rightmost strip and certify their
   simplicity;
2. prove the vertical growth needed to shift the Rice contour;
3. bound the remainder and the sparse integer `H` near zeros of the
   leading cosine; and
4. convert the resulting root control through `(147.3)` into `(147.4)`,
   including sign-change indices.

Q830 attacks this Rice/zero problem.  Q831 attacks the independent
primitive denominator lattice and lcm/gcd route.

## 148. Q802/Q806: every adjacent gcd contains one fixed Apéry--binomial core

For `n=3H+1`, put

`B_H=binom(3H+1,H+1)`,

`T_H=rad_(p>2H) gcd(A_(3H+1),B_H)`.                 `(148.1)`

The large-prime support of `B_H` consists exactly of primes
`p=3H+1-j`, `0<=j<=H`, each with valuation one.  Lucas gives

`A_(3H+1)=5A_j (mod p)`.

Hence

`T_H=product_(0<=j<=H, p=3H+1-j prime, p|A_j) p`,   `(148.2)`

so `log T_H` is exactly the direct q=1 sum.  In the general notation
`H=floor((n-1)/3)`, Q806 records the equivalent identity

`S(n)=log rad_(p>2n/3) gcd(A_n,binom(n,H+1))`.       `(148.3)`

This is a useful fixed two-integer formulation, but the carrier still has
positive linear entropy and no audited existing G-function, Subspace
Theorem, recurrence-gcd, affine-sieve, primitive-divisor, or
greatest-prime-factor theorem gives the required pointwise radical bound.
Q806's block-determinant compression lemma is a correct sufficient
condition, but no construction satisfying it is supplied.

Q802 proves a sharper exact obstruction for every adjacent Padé pair.
Let

`R_H(b)=rad_(p>2H) gcd(P_b(n),P_(b+1)(n))`.

Then

`R_H(b)=T_H N_H(b)`, `gcd(T_H,N_H(b))=1`,           `(148.4)`

and every prime of `N_H(b)` divides the adjacent constant `kappa_b`.
The only delicate case is a common prime that also divides `B_H` but not
`A_n`.  It corresponds to a node `j` where both denominator polynomials
vanish modulo `p`.  Differentiating

`P_bQ_(b+1)-P_(b+1)Q_b=kappa_b binom(x,H+1)`

at `x=j` shows `p|kappa_b`, because

`[d/dx binom(x,H+1)]_(x=j)
 =(-1)^(H-j)j!(H-j)!/(H+1)!`

is a `p`-adic unit.  Thus `(148.4)` is fully exact, not a heuristic
carrier statement.

It follows that a subexponential adjacent gcd at any fixed slope is
equivalent to two separate statements: the original core bound
`log T_H=o(H)` and a Padé-specific pollution bound `log N_H(b)=o(H)`.
The first is slope-independent and survives every adjacent chain.
Therefore the adjacent-gcd route cannot bypass Problem 3.2; it can only
add a second problem.  Q832 now attacks the fixed core `(148.1)` directly
through hypergeometric remainders and Bézout certificates.

## 149. Q822: the canonical endomorphism route is scalar; the exact cubic is the remaining datum

Q822 audits whether the canonical low-height Apéry gauge has an
endomorphism, contiguous operator, or modular correspondence that could
create a second synchronized integral section.  The answer is negative.
For

`T(t)=[[P(t),-t^6],[1,0]]`,

`P(t)=34t^3+51t^2+27t+5`,

every rational horizontal endomorphism satisfying

`E(t+1)T(t)=T(t)E(t)`

is scalar.  The proof reduces a non-scalar endomorphism to a rational
solution of

`z(t+1)=P(t)-t^6/z(t)`.

Writing `z=A/B` forces

`z=c t^3 B(t+1)/B(t)`.

The coefficients at infinity then give

`c+c^(-1)=34`,

`deg B=(51-3c)/(c-c^(-1))=-3/2`,

an impossibility.  Hence every forward `h`-shift intertwiner is a scalar
multiple of the canonical transfer product.  The determinant pairing is
the unique adjoint/self-duality map up to scalar and does not create an
independent solution.

The modular picture has the same outcome.  The rank-two
`Gamma_1(6)` local system has scalar commutant, as does its symmetric
square.  Hecke and Atkin--Lehner correspondences move the base point or
level structure; any same-base descent is scalar.  Thus the
Picard--Fuchs origin supplies no hidden second holomorphic integral
section.

Q822 also quantifies the cost of abandoning the canonical gauge.  If an
integral matrix `G` sends `(1,5)` to an independent state synchronized
with a star subset of radical `R_S`, then

`||G||_infinity>=R_S/36`,

and already

`max(||G^(-1)T(1)G||,||G^(-1)T(2)G||)>=R_S/72`.

Thus Q800's CRT conjugacy necessarily writes the unknown radical into the
initial coefficient height.

The exact surviving same-gauge detector is

`det(row_a,row_b)=(a!)^6 N_(b-a)(a)`,               `(149.1)`

where `N_h` is the gap continuant.  If `p>b` and `p|A_a`, then

`p|det(row_a,row_b) <=> p|A_b`.                     `(149.2)`

This converts a second zero under the same prime into a canonical integer,
but supplies no frequency lower bound.  Q822 formulates one sufficient
cross-incidence estimate: if the star has size `K` and its ordered
off-diagonal incidence count is at least

`K^(1+delta)n^(-o(1))`

for some `delta>2/3`, then the known vertical zero bound forces
`K=o(n/log n)`.

Finally, a hostile recurrence can preserve the fixed initial state,
coefficient size `O(n^3)`, off-diagonal term `-m^6`, adjacent determinant
law, and arbitrary prescribed reflected zero pairs.  What it changes is
exactly the diagonal cubic `P(m)`.  Therefore bounded height and the
Casoratian law do not imply the needed correlation; a successful proof
must use the exact cubic or its modular/Dwork meaning.

Q822 does not improve the unconditional radical bound, but it closes the
canonical symmetry/endomorphism/constant-conjugacy route cleanly.  Q836
now attacks the explicit gap continuants in `(149.1)` for strong
divisibility, resultants, or a genuine exact-cubic correlation theorem.

## 150. Q824: higher Frobenius jets saturate at homogeneous degree

Q824 independently rederived the first ray defect and continued it through
`p^3`.  For the constant coefficient

`S_q(p,r)=[c^0]T_(qp+r,p-1)(c)`,

ordinary Taylor expansion of the integral polynomial

`Acal_p(X)=sum_(m=0)^(p-1)L(X,m)F_m`

gives

`S_q=A_r+qp Acal_p'(r)+q^2p^2 Acal_p''(r)/2 (mod p^3)`. `(150.1)`

The first and second forward differences gain `p` and `p^2`
universally, whether or not `p|A_r`.  At the first bad fiber `(p,r)=(5,1)`
the normalized first three scalar jets are all units.  The natural
Schwarzian numerator has valuation exactly its universal value `4`.

The full coefficient vector is even less degenerate.  After dividing the
three bad rows at outer indices `6,11,16` by `5`, their columns of degrees
`0,5,10` reduce to

`[[3,4,0],[0,4,3],[2,2,0]] (mod 5)`,

whose determinant is `1`.  The corresponding order-three Hasse minor has
valuation exactly `3`, and the order-two subminor has valuation exactly
`2`.  The final checker in Q824 was extracted and run locally.  It
verified the scalar `p^3` Taylor law for every odd prime through `43` and
every residue, as well as all quoted `p=5` valuations and minors.

The resulting formal saturation statement is precise.  After universal
finite-difference powers are removed, one equation `A_r=0 (mod p)` can
force at most one additional factor per homogeneous row degree.  A larger
order is equivalent to a new algebraic relation among the normalized
singular jets.  This is a theorem in the universal jet algebra; ruling out
every specially designed invariant on the actual Apéry family would still
require the horizontal singular-jet noncontainment statement formulated in
Q824.  The explicit `p=5` and `p=17` fibers already rule out the natural
derivative, reflection, low-order determinant, Hessian, and Schwarzian
relations.

Thus higher same-prime rays cannot improve the selective-order/height
ratio by formal differentiation alone.  Q842 instead attacks the
different-modulus short-gap pair bound isolated by Q824.

## 151. Q828: the fixed-denominator strip asymptotic, with one repair

After refetching the Notion equation blocks, Q828 supplies the missing
uniform proof of the Q803 fixed-`b` asymptotic:

`Phat_(H,b)(3H+1)`

` ~(-1)^b (lambda/gamma)^b c_(H-b) binom(3H+1,H-b)`

`   *(1-1/(2gamma))^(-b-1)`                         `(151.1)`

uniformly on every fixed bounded set of denominator degrees.

The primary input is genuine.  Edgar, *The Apéry Numbers as a Stieltjes
Moment Sequence*, Proposition 11, represents `A_m` as the moments of a
positive density on `[0,lambda]`; Proposition 25 and the normalized
Frobenius solution give a convergent square-root expansion at the right
endpoint.  Shifting the moment variable by one yields

`c_k=int_(-1)^gamma y^k phi(y+1)dy`

and hence the full differentiated endpoint expansion

`c_k=kappa gamma^k k^(-3/2)(1+alpha_1/k+...)`.       `(151.2)`

Q828 then converts each normalized multiplication column into a
fixed-degree polynomial plus a symbol remainder.  Divided differences on
the consecutive high rows expose the nonzero Vandermonde term uniformly
over an `O(log H)` window.  Cramer's rule supplies a global polynomial
cofactor bound, while the exact binomial ratio supplies a factor at most
`2^(-r)`.  Choosing the logarithmic cutoff to pay for both polynomial
losses makes the remaining tail `O(H^(-2))`.

There is one false displayed estimate in the returned proof:

`|c_k-C_+(k)| <= (1-epsilon)^k`.                    `(151.3)`

The omitted negative interval reaches `y=-1`, and Edgar's density has a
logarithmic singularity there, so its absolute contribution need not
decay exponentially.  The proof only needs a relative estimate, and that
repair is immediate:

`|c_k-C_+(k)|<=int_0^1 phi(x)dx<=1`,

whereas `C_+(k)` has size `gamma^k k^(-3/2)`.  Hence the relative error,
and consequently the error in every fixed shifted ratio, is
`O(gamma^(-k) k^(O(1)))`, which is exponentially small.  The determinant
and tail argument is unchanged.  Positivity of the finitely many early
`c_k` should also be cited or included as an exact finite verification
rather than silently inferred from the eventual asymptotic.

An independent exact computation with the existing checkerboard script
confirmed the predicted ratios for every `b<=8` at `H=200,250,300`.  At
`H=300` they range from `1.000024639128` for `b=0` to
`0.994738599821` for `b=8`.

This closes the fixed-denominator boundary, not the proportional
checkerboard regime and not the reciprocal fixed-numerator boundary.
Q845 is performing a final hostile audit with every formula forced into
plain-text fences.

## 152. The reciprocal boundary: the first complex zero is not a dominant term

For

`U_n=sum_(k=0)^n (-1)^k binom(n,k)/A_k`,

there is an exact ordinary generating-function identity

`sum_(n>=0) U_n z^n`

` =(1-z)^(-1) R(-z/(1-z))`,

`R(w)=sum_(k>=0) w^k/A_k`.                          `(152.1)`

Thus the ordinary reciprocal generating function remains a natural
object even though the Golyshev--Zagier interpolation also gives a formal
Nörlund--Rice representation.

For a simple zero `alpha` of that interpolation, the residue contribution
to the finite difference, with the standard Rice orientation, is

`Gamma(n+1)Gamma(-alpha)`

` /(Gamma(n+1-alpha) A'(alpha))`.                   `(152.2)`

The sign was checked against the exact model `f(z)=1/(z+a)`.
Using the numerical zero

`alpha=0.1459867731180108665+0.5824261402448531507 i`,

the conjugate-pair contribution was compared with exact or
high-precision values through `n=2000`.  It does not approximate the
answer in this range.  For example, at `n=100` the true value is
`0.2485554104` while this pair contributes `-1.169866351`; at `n=2000`
the values are `0.06275541347` and `0.1772231020`.  Sampled true values
remain positive through `n=10000` but oscillate slowly rather than
converging to the initially suspected `1/16`.

This does not refute eventual residue dominance.  It proves that the
single-zero heuristic is not an audited argument: one must locate and
sum the other poles and prove that the Rice contour remainder is smaller.
The reproducible calculation is
`research/scripts/q32_reciprocal_boundary_audit.py`.  Q830 is testing the
missing zero distribution and contour estimates; Q846 attacks primitive
integrality independently.

## 153. Q830: an unconditional zero-free line, but only a conditional Rice theorem

After refetching the omitted Notion equation blocks, Q830 gives the exact
Golyshev--Zagier interpolation

`Acal(z)={}_4F_3(-z,-z,z+1,z+1;1,1,1;1)`.

On the reflection line `z=-1/2+it`, every summand is

`|(1/2+it)_k|^4/(k!)^4`,

so the series converges to a positive real number at least `1`.  Therefore

`Acal(-1/2+it) != 0`,

`|1/Acal(-1/2+it)| <= 1`.                           `(153.1)`

This is an unconditional and useful new lemma.

For the Rice kernel

`K_H(z)=Gamma(-z)Gamma(H+1)/Gamma(H+1-z)`,

the same line gives the Fourier representation

`K_H(-1/2+it)=FT[(1-e^(-u))^H e^(-u/2) 1_(u>0)](t)`.

The squared `L^2` norm of this function is exactly `1/(2H+1)`, and
the squared norm of its first derivative is `O(1/H)`.  The standard
weighted Cauchy--Schwarz/Plancherel inequality consequently gives

`int_R |K_H(-1/2+it)|dt=O(H^(-1/2))`.               `(153.2)`

Together with `(153.1)`, this bounds the terminal-line integral by
`O(H^(-1/2))` *if* the original finite Rice contour can first be moved to
that line.  The line bound by itself does not justify the displacement.
One still needs high-imaginary minimum-modulus estimates for `Acal`,
horizontal-boundary decay, and summability of every crossed zero
residue.  Q830 states these assumptions explicitly as its unproved
Theorem ZR.  Thus its later asymptotics are conditional, not a theorem
about `U_H`.

Conditionally on one simple strictly rightmost conjugate zero pair
`rho=alpha+i beta`, Q830's exact beta-mode algebra is sound:

`D_H(v(rho)+v(sigma))`

` =((rho-sigma)^2/H)v_H(rho)v_H(sigma)`.            `(153.3)`

For the conjugate pair this makes the leading Casoratian negative and
phase-free, of size `H^(2alpha-1)`.  The associated projective numerator
at `3H+1` is generically of size `H^(2-alpha)` and has a uniform
polynomial upper bound under the same contour and zero-gap hypotheses.
This would show that any exponential size of the primitive `a=1`
certificate comes from its Smith clearing multiplier, not its
Archimedean projective direction.

Two points remain under audit.  Q830 asserts weighted versions of
`(153.2)` for every fixed moment without displaying the derivative
calculation, and the infinite-residue/line cross terms used in the
Casoratian require precisely those bounds.  More importantly, neither
the proposed rightmost zero nor the absence of high-imaginary zeros to
its right has an interval/global certificate.  The severe numerical
mismatch in Section 152 is compatible with eventual dominance but gives
no support for using it at accessible heights.  A queued `dm4` follow-up
was routed as Q857 to `dm3`; it attacks the Barnes/high-imaginary step and the
ordinary-generating-function alternative directly.

## 154. Q817: bounded character order is harmless; the delivered sieve analysis is incomplete

Q817 correctly rederives the q=1 Lucas diagonal and the exact Mellin
coefficient identity

`A_r=-sum_(t in F_p^*) H_p(t)t^(-r) (mod p)`,

where `H_p(t)=sum_(m=0)^(p-1) A_m t^m` is the scalar Hasse--Witt
section of the canonical Apéry K3 family.  This distinguishes a moving
Mellin coefficient from a Hasse invariant evaluated at one fiber.

Its main exact arithmetic statement is the bounded-order slice lemma.
Writing

`d=(p-1)/gcd(p-1,n-p)`

for the order of the selected Kummer character, there is an integer
`1<=j<=d/2` with

`p=(dn+j)/(d+j)`.                                   `(154.1)`

Consequently characters of order at most `D` yield `O(D^2)` candidate
primes for each fixed `n`.  This is correct, but it is neither new nor
the strongest available count.  Q362 had already used

`g=gcd(p-1,n-1)`, `d=(p-1)/g`, `p=1+dg`, `g|(n-1)`

to obtain only `D tau(n-1)=D n^(o(1))` candidates with `d<=D`.
For every fixed `epsilon>0`, choosing `D=n^(1-epsilon)` makes their
weighted contribution `n^(1-epsilon+o(1))=o(n)`.  The genuinely hard
part can therefore be confined to

`d>n^(1-epsilon)`, `g<n^epsilon`,                    `(154.2)`

which is much sharper than the square-root cutoff obtained from
`(154.1)` alone.

This does not yet create an estimate for the high-order part.  A Kummer
sheaf can have bounded geometric conductor as its character order grows
inside one finite field, but the residual field, character, and selected
Mellin coefficient all move together here.  A vertical equidistribution
or large-sieve theorem over many characters of one field does not
control the single prescribed character in each different field.
Restriction-of-scalars rank growth is relevant to a characteristic-zero
compatible-system approach, but is not by itself an obstruction to a
finite-field trace estimate.

The delivered Q817 file is incomplete: it ends at Section 7, while its
executive verdict cites Sections 8--10 and a sufficient composite-modulus
local-limit theorem that do not occur in the answer.  Its unsupported
claim that the generic Mellin cohomology is rank two, pure of weight
three, symplectic, and has connected monodromy `SL_2` is also not
accepted without a primary citation and exact exceptional-character
hypotheses.  Q854 requests the missing sections, a citation audit, and
an explicit exponent calculation for the strongest applicable large
sieve.  It must incorporate Q362's stronger parametrization `(154.2)`;
until that returns, Q817 adds a useful geometric interpretation but no
new unconditional order-split gain.

## 155. Operational failure: `ask-gpt.py --help` is a submitted question

The polling wrapper has no ordinary `--help` branch.  Running
`ask-gpt.py --help` therefore submitted the literal text `--help` as
Q852 to `dm4`.  The task was explicitly completed through its own
bridge endpoint with marker `operator-cancelled-q852`; the later
all-connectors-failed notification is stale and must not trigger a retry.
This consumed one tab temporarily but did not alter any mathematical
state.  Future usage inspection must read the script source; real prompts
are passed as the sole positional argument.

The intended analytic follow-up became Q857 on `dm3`.  The newly free
`dm4` was then given Q859, an independent interpolation-capacity attack
on the integer values `A_j/(n-j)` at hit nodes.  Its precise purpose is
to compute whether Vandermonde/Hermite denominator gain can beat the
actual Apéry exponential growth constant, and, if not, to prove the
sharp saturation inequality.  The connector reported `dm4` degraded
after dispatch, so Q859 is watched for a late Notion drop and will not
be duplicated while the task remains pending.

## 156. Exact scan: the first multi-hit examples survive in the near-primitive regime

The canonical C++ scan was rerun through outer index `50000` and extended
to print

`g=gcd(p-1,n-1)`, `d=(p-1)/g`

for every maximum-multiplicity q=1 tuple.  The largest q=1 multiplicity
in this range is still three.  Several three-hit examples lie entirely
in one very small-kernel class:

- at `n=11576`, the three direct primes have `g=1` and orders
  `8892,9318,11436`;
- at `n=18444`, all three hits again have `g=1`;
- at `n=26164`, the three direct primes all have `g=3`;
- at `n=47066`, the three reflected primes all have `g=1`;
- at `n=47859`, all three hits have `g=2`.

Thus Q362's low-order estimate is a genuine reduction but does not
explain even the earliest multi-hit configurations: the hard
near-primitive stratum is already populated and can contain several
aligned zeros with the same small kernel.  Conversely, this is only
finite numerical evidence; it neither disproves an asymptotic horizontal
large sieve nor supplies a pointwise bound.  The exact output is
reproducible with

`g++ -O3 -std=c++17 q32_scan.cpp -o /tmp/q32_scan`

`/tmp/q32_scan 50000`.

## 157. Q850: the primitive `a=1` lattice is exactly one-dimensional

Q850 removes the remaining normalization ambiguity in Section 147.  Put

`L_H=lcm(A_0,...,A_H)`,

`X_H=L_H U_H`,

`Y_H=L_H H(U_H-U_(H-1))`,

`h_H=gcd(X_H,Y_H)`,

`z_H(x)=(Y_H-X_H x)/h_H`.

Then `z_H` is primitive linear, and the exact least multiplier making the
node quotients integral is

`m_H=lcm_(0<=s<=H) A_s/gcd(A_s,z_H(s))`.            `(157.1)`

The primitive integral Padé numerator is, up to sign,

`P_H^prim(x)=m_H z_H(x)`.                           `(157.2)`

This follows because the Pascal transformation between integer node
values and integer Newton coefficients is unimodular.  If a prime divided
every coefficient of the pair after the least multiplier was used, one
could divide the multiplier by that prime, contradicting minimality.
The local formula is

`v_l(m_H)=max_s max(0,v_l(A_s)-v_l(z_H(s)))`.        `(157.3)`

The formula was checked independently for every `2<=H<=12` against
`primitive_a_one_pair`; both primitive numerator coefficients agree up to
one global sign.

Q850's same-sign lemma is also valid.  If `U_H U_(H-1)>0`, the rational
root of `z_H` lies to the left of `H`; nodewise divisibility at `H-1` and
`H` then gives

`|P_H^prim(3H+1)|>=A_(H-1)`.                        `(157.4)`

This is a no-go result for obtaining a small integral certificate from the
fixed-`a=1` direction.  It is not a proof of the q=1 radical bound.
Q850's remaining proposed lower bound `log m_H>=cH` at sign changes is
another cross-index gcd problem.  Formula `(157.1)` is bankable; that
lower bound remains open and is not needed for the support audit below.

## 158. Q854: the Mellin large-sieve continuation confirms a quantifier barrier but loses the best order split

Q854 withdraws Q817's unsupported assertion that the generic Mellin
cohomology is automatically rank two, pure of weight three, symplectic,
and has connected monodromy `SL_2`.  A tame compact-support
Euler-characteristic calculation naturally starts at dimension six if
the endpoint invariants vanish; obtaining a rank-two primitive quotient
would require an explicit local-monodromy decomposition.  Moreover a
generic Kummer twist is paired with its inverse twist rather than being
self-dual by itself.  Thus none of the claimed symplectic package is
currently available.

Its theorem audit reaches the same negative conclusion as Section 154:
fixed-field character equidistribution, fixed-modulus trace-function prime
sums, and effective Chebotarev have the wrong quantifiers for one marked
character in each varying residual field and the same-characteristic
condition `p|A_(n-p)`.  Q854 isolates a sufficient new “marked Mellin
diagonal local-limit” estimate, but this is a name for the missing
power-saving theorem, not an existing result.

Q854 did not follow the requested strongest arithmetic reduction.  It
uses the weaker count `O(D^2)` and consequently stops at
`D<sqrt(n)/log^2(n)`.  The already proved parametrization

`p=1+dg`, `g|(n-1)`

gives only `D tau(n-1)=D n^(o(1))` candidates with `d<=D`, so every
`D=n^(1-epsilon)` is allowed.  The correct hard range remains `(154.2)`.
Accordingly Q854 contributes a useful hostile audit of the geometric
claims, but no improvement to the unconditional q=1 estimate.

## 159. The single-`a=1` candidate-window radical has unavoidable multi-zero pollution

A tempting sufficient statement after `(157.2)` was

`log W_H=o(H)`,

`W_H=rad_(2H<p<=3H+1)`

`     gcd(P_H^prim(3H+1),binom(3H+1,H+1))`.         `(159.1)`

Every direct target prime divides `W_H`, while the enormous prime factors
of `P_H^prim(3H+1)` outside the candidate window are irrelevant.  Exact
computation initially made `(159.1)` look plausible.  The following local
classification shows why that impression was misleading.

Fix a candidate prime `p>2H` and let

`Z_p(H)={s in [0,H]:p|A_s}`.

From `(157.3)`:

1. if `Z_p(H)` is empty, then `p` does not divide `m_H`;
2. if `Z_p(H)` contains two distinct nodes, then `p|m_H`;
3. if `Z_p(H)={r}` and `v_p(A_r)=1`, then `p` does not divide `m_H` and
   `z_H(r)=0 mod p`.

The second assertion is immediate: if `p` did not divide `m_H`, then
`z_H` would vanish at every node of `Z_p(H)`, impossible for a primitive
nonzero linear polynomial modulo `p`.  For the third, reduce the cleared
moments

`X_H=sum_s (-1)^s binom(H,s)L_H/A_s`,

`Y_H=sum_s (-1)^s binom(H,s)sL_H/A_s`

modulo `p`.  With a unique zero node, only its term survives, so
`Y_H=rX_H mod p` and `X_H` is a unit.  Higher `p`-adic multiplicity leaves
the exact residual condition in `(157.3)`.

Thus every prime with at least two prefix zeros is unconditional
good-prime pollution in the multiplier, independently of the projective
direction.  A unique simple zero is instead cancelled from the multiplier;
when it is the moving node `r=3H+1-p`, it reappears in the root factor
`z_H(3H+1)`, exactly as required for a target prime.  Root primes with no
prefix zero are a second, sporadic pollution source.

Two reproducible scans now separate these effects:

- `q32_a1_window_support.py` computes the exact local valuations of
  `X_H,Y_H,z_H,m_H` through `H=500`, without constructing the full
  degree-`H-1` denominator.
- `q32_a1_multiplier_window.cpp` scans only Apéry recurrences modulo
  candidate primes.  It gives an exact lower bound from primes having at
  least two prefix zeros and an upper bound from primes having at least
  one prefix zero.

At `H=50000`, the unavoidable multi-zero lower-bound weight divided by
`H` is `0.07427438645`, from 317 primes.  Over
`25000<=H<=50000`, its maximum ratio is `0.07430410809` at `H=49980`.
The one-zero-or-more upper ratio at `H=50000` is `0.3452990877`.
Meanwhile the exact moving target at that particular height is empty.

These finite computations do not disprove the asymptotic `(159.1)`.
They do overturn the sparse-support heuristic: the forced pollution
shows no decay through height `50000` and agrees with the previously
measured reflection-paired zero-orbit statistics.  Proving `(159.1)`
would require showing that primes with two lower-half zero orbits have
vanishing weighted density in every moving prime window, contrary to the
current data and stronger than the original alignment problem in a
different direction.  The single-`a=1` support bound should therefore be
treated as a likely false route.  A second certificate is needed to
remove this multiplier pollution; Q871 asks for exactly that separation.

## 160. A growing family of Padé numerators exactly absorbs every bounded zero set

The multi-zero obstruction has a clean general form.  Let
`(P_(H,a),Q_(H,a))` be the primitive integral interpolation pair with

`deg P_(H,a)<=a`, `deg Q_(H,a)<=H-a`,

`P_(H,a)(s)=A_s Q_(H,a)(s)` for `0<=s<=H`.

For a prime `p>H`, put

`Z_p(H)={s in [0,H]:p|A_s}`, `z_p(H)=|Z_p(H)|`.

Because every factorial in the Newton basis is a `p`-unit, reduction gives
ordinary polynomials of the same degree over `F_p`.  The following exact
zero-absorber lemma is then elementary:

`P_(H,a)=0 in F_p[x]` if and only if `z_p(H)>a`.     `(160.1)`

Indeed, `P` vanishes at every node of `Z_p(H)`, so more than `a` zeros
force `P=0`.  Conversely, if `P=0`, then the nonzero polynomial `Q` must
vanish at all `H+1-z_p(H)` complementary nodes.  Since
`deg Q<=H-a`, this is possible only when `z_p(H)>=a+1`.  The reduced pair
is nonzero because the integral coefficient vector is primitive.

If `z_p(H)<=a`, then the nonzero `P_(H,a)` has all `z_p(H)` prescribed
roots and at most `a-z_p(H)` extra interpolation-node roots.  In
particular,

`a=z_p(H) =>`

`{s in [0,H]:P_(H,a)(s)=0 mod p}=Z_p(H)`.           `(160.2)`

This yields an exact multi-certificate purification.  For

`G_(H,A)(n)=gcd_(0<=a<=A) P_(H,a)(n)`

and a candidate prime `p=n-j>2H`,

`p|G_(H,A)(n)`

if and only if

`j in Z_p(H)` or `z_p(H)>A`.                        `(160.3)`

If `z_p(H)<=A`, the member with `a=z_p(H)` has precisely the Apéry zero
nodes as roots, while every target node is a root of every member.  If
`z_p(H)>A`, `(160.1)` makes every numerator zero modulo `p`.

The previously proved uniform vertical bound
`|Z_p|<<p^(2/3)` therefore allows a sublinear family
`A=O(H^(2/3))` whose common candidate-prime support is exactly the direct
q=1 target.  This genuinely removes the positive-density multi-zero
pollution found in Section 159; it does not yet bound the resulting gcd.
The adjacent cross identities still carry the full binomial factor once
per target prime, so their universal height estimates saturate.  The new
opening is narrower: seek an `exp(o(H))` bound for the candidate-window
radical of the *whole growing-family determinantal divisor*, using the
codimension supplied by `(160.1)` rather than multiplying adjacent
identities.

The verifier `q32_pade_zero_absorber.py` checks `(160.1)--(160.3)` with
the exact primitive kernels for every `2<=H<=22` and every prime
`H<p<=3H+1`: 2,100 polynomial reductions and 1,071 growing-family support
tests pass.  The proof above is independent of the computation.

## 161. The zero-absorber-scale family gcd is tiny through height 40

The exact experiment `q32_pade_family_gcd.py` now computes

`G_(H,A)(3H+1)=gcd_(0<=a<=A) P_(H,a)(3H+1)`

at the first scale relevant to `(160.3)`,

`A=ceil(H^(2/3))`.                                  `(161.1)`

It uses the already verified primitive cofactor kernels, checks every
interpolation equation before taking a gcd, and independently checks the
candidate-window support predicted by `(160.3)`.  Selected exact outputs
are:

| `H` | `A` | bits of `G_(H,A)` | `log G_(H,A)/H` | factorization |
|---:|---:|---:|---:|:---|
| 10 | 5 | 5 | 0.321888 | `5^2` |
| 20 | 8 | 7 | 0.241416 | `5^3` |
| 25 | 9 | 13 | 0.359861 | `5^2*17*19` |
| 30 | 10 | 13 | 0.281308 | `5^3*37` |
| 35 | 11 | 10 | 0.183936 | `5^4` |
| 40 | 12 | 9 | 0.151302 | `5^2*17` |

This collapse is substantial.  At `H=40`, the constant-numerator member
alone has 3,666 bits, the gcd after `a=1` has 34 bits, after `a=6` it has
12 bits, and after `a=12` only 9 bits remain.  Across every
`2<=H<=40`, the final gcd has at most 15 bits in this computation.

There is an exact bordered-minor formulation behind the experiment.  Let

`M_(k,l)=binom(k,l) Delta^(k-l) A_l`

be the Newton multiplication matrix, put `b=H-a`, and let

`R_(H,a)=M_[a+1..H,0..b]`.

This is a `b` by `b+1` matrix.  If `delta_(H,a)` is the gcd of its maximal
minors, its signed cofactor vector divided by `delta_(H,a)` is precisely
the primitive denominator vector `q_(H,a)`.  Define the evaluation row

`w_l(n)=sum_(k=0)^H binom(n,k) M_(k,l)`, `0<=l<=b`.

The bottom Padé equations let one include all rows in this sum, and
cofactor expansion gives the exact identity

`P_(H,a)(n)=(-1)^b det([R_(H,a);w(n)])/delta_(H,a)`. `(161.2)`

Thus the new gcd is a gcd of *individually saturated bordered minors*.
The different saturation indices `delta_(H,a)` are the precise obstacle
to declaring it, without further work, a single ordinary determinantal
divisor.

The computation is only evidence.  Except for the small hit
`H=5, p=11`, the direct candidate window is empty through `H=40`; the
next direct hits occur at `H=74` and `H=80`.  Therefore the tiny values
above do not sample a dense target regime, and a pointwise asymptotic
cannot be inferred from them.  Hadamard bounds control each unnormalized
bordered determinant in the wrong direction and at quadratic-exponential
height.  Adjacent cross determinants still contain the full
`binom(3H+1,H+1)` carrier.  Generic many-integer gcd heuristics also give
no deterministic estimate.

The now sharply isolated sufficient lemma is

`log gcd_(0<=a<=ceil(H^(2/3)))`

`        |P_(H,a)(3H+1)| = o(H)`,                   `(161.3)`

or merely the same estimate for its `p>2H` radical.  Q885 asks for a
whole-family height/Smith argument for `(161.3)`.  Q886 asks first for a
single Fitting presentation that correctly incorporates the varying
indices in `(161.2)`, or a rigorous saturation counterexample.  Until one
of those steps is proved, `(161.3)` remains a promising numerical
strengthening, not an unconditional advance.

## 162. Growth and integer interpolation alone cannot bound the family gcd

There is a rigorous adversarial construction showing what a proof of
`(161.3)` must use.  Fix any `rho>1` and a rapidly separated sequence
`H_m=4^m`.  For every prime

`2H_m<p<=(5/2)H_m`

assign the distinct index

`j_(m,p)=3H_m+1-p`.

These indices lie in `[H_m/2,H_m]`, and the index blocks for distinct `m`
are disjoint.  Define a positive integer sequence `B_j` as follows.  At an
assigned index, take the least positive multiple of `p` not below
`rho^j`; at every other index take `ceil(rho^j)`.  Then

`B_j/rho^j -> 1`.                                  `(162.1)`

Indeed, the rounding error at an assigned index is at most
`p=O(H_m)`, whereas `j>=H_m/2`, so the relative error is exponentially
small.  The same construction works with the sharper model scale
`rho^j j^(-3/2)` by rounding that quantity to the next multiple of `p`.

For `n_m=3H_m+1`, every selected prime satisfies

`p=n_m-j_(m,p)` and `p|B_(j_(m,p))`.

Consequently its direct moving-prime weight obeys

`sum_(0<=j<=H_m, n_m-j prime, n_m-j|B_j) log(n_m-j)`

` >= theta((5/2)H_m)-theta(2H_m)`

` = (1/2+o(1))H_m`                                 `(162.2)`

by the prime number theorem.  Thus positivity, integrality, and even the
correct first-order exponential-plus-polynomial size of the sampled
sequence are compatible with a linear direct-prime weight along an
infinite subsequence.

The obstruction passes through every universal interpolation certificate:
if an integral Newton polynomial `P` satisfies

`P(s)=B_s Q(s)` for `0<=s<=H_m`,

then for a selected `p>2H_m`, degree `<p` and `n_m=j_(m,p) mod p` give

`P(n_m)=P(j_(m,p))=0 mod p`.

Hence every numerator in any Padé family retains all primes in `(162.2)`,
regardless of how its Smith presentation is organized.  This does not
disprove `(161.3)` for the Apéry sequence: the constructed `B_j` does not
satisfy the Apéry recurrence, reflection, or Dwork congruences.  It does
prove that a whole-family determinant argument based only on integer
interpolation and Archimedean growth cannot work.  A successful Fitting
estimate must insert a recurrence-specific restriction on the minors or a
cross-prime consequence of the Apéry differential equation.
