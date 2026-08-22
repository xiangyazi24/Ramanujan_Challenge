ANSWER Q2064 eae1ce2b

# Verdict: ACTUAL_SCOPED_SATURATION

There is a genuinely stronger actual-orbit statement than the generic CRT / endpoint-pin no-go, but it is a **saturation theorem**, not a first-hit contraction.

The key move is to translate all boundary polynomials to one common integer variable.  After that translation, the **entire actual Apéry prefix** enters through a rank-two coefficient lattice.  The first three genuine Apéry rows already generate this lattice with index dividing `24`.  Consequently every top prime `ell>=7` sees the full actual-prefix coefficient module as Smith-unimodular.  This rules out marked torsion coming from any linear prefix-aware Smith/Fitting construction, any Bezout combination of the boundary prefix, and any finite collection of linear finite-difference / WZ boundary operations.

For the two marked rows themselves there is one exact actual-orbit degeneracy determinant

```text
K_h = b_(h+1) - 5*b_h.
```

It is the determinant between the true initial Apéry state and the true `h`-th prefix state.  It gives a sharp dichotomy:

```text
if gcd(K_h,p*q)=1:
    the p- and q-support of the two hit equations is exactly diagonal;
    after endpoint/vertexwise saturation no mixed marked support remains.

if p | K_h:
    p*q | F_(s,0)(q).

if q | K_h:
    p*q | F_(s,h)(p).
```

Thus `K_h` is a real actual-prefix **singularity gate**.  However it is not a universal carrier: actual first hits occur abundantly in the generic branch `gcd(K_h,pq)=1`.  The standard exact examples `(p,q,m)=(139,181,200)` and the three pairs in row `321` are all first hits and all lie in that generic branch.

The full-prefix Smith calculation is even stronger: `K_h` is only one `2 x 2` minor of the actual prefix coefficient matrix.  The complete coefficient matrix has no `p`- or `q`-primary Smith torsion at all for top primes.  Hence adjoining more genuine prefix rows does not turn the `K_h` gate into a universal Smith carrier.

The natural mixed integer identity obtained from normalized lift quotients is exact but also saturated:

```text
u = F_(s,0)(q)/q,
v = F_(s,h)(p)/p,

K_h*A = p*v - b_h*q*u,
K_h*B = 5*p*v - b_(h+1)*q*u,
```

where

```text
A = N_s(q),
B = (q+1)^3*N_(s-1)(q+1).
```

When `K_h` is a unit at `p` and `q`, these are simply an invertible change of two coordinates.  They do not yield a new marked divisor after the two endpoint factors have been removed.

No sourced analytic theorem located in the current repository or in the targeted literature check gives a uniform saving for

```text
q | b_s,
p | b_(s+h),
q=p+h.
```

The repository itself labels its two-prime correlation bound as numerical only and notes that the naive two-prime sheaf formulation mixes residue characteristics.  The published fixed-prime results I checked do not bridge this gap.

So the bankable conclusion is:

```text
ACTUAL_SCOPED_SATURATION
```

for the whole **linear actual-prefix boundary / Smith / Fitting / finite-difference / WZ-boundary class**.  A genuinely nonlinear multi-time invariant, or a proved mixed-characteristic incidence estimate, remains open.

---

# 1. Repository and source audit

I read the current connected repository state before doing the algebra.

The `chatgpt-drop` branch head at the start of this audit was

```text
c1dbe771f9d10d0b7efda51ff4c155f51a7ad551
Replace Q2059 with adjacent-transfer Smith saturation audit
```

and I read the resulting

```text
drops/Q2059-d7657cf0.md
```

rather than reconstructing Q2059 from memory.  Its exact relevant conclusion is that the actual predecessor gives the new raw Lucas error

```text
L_p = b_(m-1) - 5*b_(m-p-1),
p | L_p,
```

but the actual predecessor is a marked-prime unit, the primitive two-coordinate state has no marked Smith torsion, and the extra raw prime is rank/height inflation rather than a larger primitive Smith exponent.

I also rechecked:

```text
SUBMIT/3.2/README.md
problems/3.2/CODEX_SPEC_wz_identity.md
problems/3.2/research/scripts/q32_endpoint_transfer_smith_audit.py
problems/3.2/research/scripts/q32_terminal_prefix_gcd_scan.sage
problems/3.2/research/scripts/q32_lower_terminal_prefix_audit.py
problems/3.2/research/scripts/q32_pade_projective_resultants.py
problems/3.2/THEOREMS_2026-08-01_campaign3_addendum.md
problems/3.2/campaign3_questions/answers/Q6665_zwin_section.tex
```

The first file gives the classical denominator-coordinate Lucas congruence used here.  The WZ spec gives the exact recurrence and confirms that the standard WZ certificate telescopes back to the same order-two Apéry recurrence.  The terminal-prefix scripts explicitly label their prefix-gcd output as diagnostic rather than a target-preserving theorem.  The Padé-resultant audit records that naive adjacent resultant gcds retain quadratic-size/polluted height.  None of those files already contains the first-hit contraction requested here.

The caller-local late proof file is still not connector-visible, so I use the exact normal form stated in Q2064 as authoritative where it is stronger than connector-visible state.

---

# 2. Exact first-hit normal form

Retain the question's notation:

```text
p < q,
q = p+h,
m = q+s = p+h+s,
0 < s < s+h < p,

q | b_s,
p | b_(s+h).
```

The Apéry recurrence is

```text
(n+1)^3*b_(n+1) = P(n)*b_n - n^3*b_(n-1),
P(n) = 34*n^3 + 51*n^2 + 27*n + 5,
(b_0,b_1) = (1,5).
```

The factorial-free continuants are

```text
N_0(a) = 0,
N_1(a) = 1,
N_(j+1)(a) = P(a+j)*N_j(a) - (a+j)^6*N_(j-1)(a).
```

The banked backward determinant theorem gives, for `1<=t<s`,

```text
ell | b_(q+t)
iff
ell | N_(s-t)(q+t),

ell in {p,q}.
```

Therefore first pair occurrence is exactly

```text
p*q does not divide N_(s-t)(q+t)
for every 1 <= t < s.
```

It will be useful to reindex this prefix by

```text
j = s-t,
C_j = N_j(q+s-j) = N_j(m-j).
```

Then the whole earlier first-hit packet is

```text
C_1, C_2, ..., C_(s-1),
```

with

```text
C_1 = N_1(m-1) = 1.
```

So any literal gcd/Bezout ideal of the **entire scalar backward-continuant prefix** is already the unit ideal.  The information in first occurrence is not a hidden common divisor of that packet; it is the open condition that no individual `C_j` carries both marked primes.

This is the first exact prefix-level reason ordinary prefix gcds cannot work.

---

# 3. Translate every boundary polynomial to one common packet

The question gives

```text
F_(s,r)(X)
  = N_s(X+r)*b_(r+1)
    - (X+r+1)^3*N_(s-1)(X+r+1)*b_r.
```

Do not compare `F_(s,0)(q)` and `F_(s,h)(p)` in two different variables.  Put

```text
Y = X+r
```

and define the translated actual-prefix boundary polynomial

```text
H_r(Y) = F_(s,r)(Y-r).
```

Then, identically over `Z[Y]`,

```text
H_r(Y)
 = b_(r+1)*N_s(Y)
   - b_r*(Y+1)^3*N_(s-1)(Y+1).
```

Define the common continuant packet

```text
A_s(Y) = N_s(Y),
B_s(Y) = (Y+1)^3*N_(s-1)(Y+1).
```

Then the entire actual boundary family factors as

```text
H_r(Y) = b_(r+1)*A_s(Y) - b_r*B_s(Y).
```

Equivalently, with the **actual Apéry coefficient row**

```text
v_r = (b_(r+1), -b_r),
```

we have

```text
H_r(Y) = v_r * (A_s(Y), B_s(Y))^T.
```

This is the first genuinely prefix-aware connected reduction: every boundary row in the actual prefix uses the same two polynomial coordinates.

At the first hit, both marked equations are evaluations at the **same integer `Y=q`**:

```text
H_0(q) = F_(s,0)(q),
H_h(q) = F_(s,h)(q-h) = F_(s,h)(p).
```

Hence

```text
q | H_0(q),
p | H_h(q).
```

The different residue characteristics remain, but the characteristic-zero packet is now genuinely common.

---

# 4. Actual-prefix 24-Bezout saturation theorem

This is the strongest new exact lemma of the audit.

The first four Apéry values are

```text
b_0 = 1,
b_1 = 5,
b_2 = 73,
b_3 = 1445.
```

Therefore the first three actual coefficient rows are

```text
v_0 = (5,    -1),
v_1 = (73,   -5),
v_2 = (1445, -73).
```

There are two exact integer Bezout identities:

```text
37*v_0 - 22*v_1 + v_2 = (24,0),
 3*v_0 - 20*v_1 + v_2 = (0,24).
```

Consequently, as polynomial identities in `Y`,

```text
24*A_s(Y)
  = 37*H_0(Y) - 22*H_1(Y) + H_2(Y),

24*B_s(Y)
  =  3*H_0(Y) - 20*H_1(Y) + H_2(Y).
```

This is not a fitted relation: it follows from the exact initial Apéry orbit and holds for every `s` and every `Y`.

## Theorem 4.1 — actual-prefix coefficient Smith saturation

Let `V_R` be the integer matrix whose rows are

```text
v_0, v_1, ..., v_R
```

for any `R>=2`.  Then for every prime `ell>=5`,

```text
V_R has full column rank over Z_(ell),
```

and both nonzero local Smith invariant exponents are zero.

### Proof

The first determinantal divisor is `1`, because the entry `-1` occurs in `v_0`.

For the second determinantal divisor it is enough to inspect two minors:

```text
det(v_0,v_1) = 48,
det(v_0,v_2) = 1080,
gcd(48,1080) = 24.
```

The gcd of **all** `2 x 2` minors divides the gcd of these two, hence divides `24`.  No prime `ell>=5` divides that determinantal divisor.

Equivalently, the explicit Bezout identities above show that the row lattice contains

```text
24*Z^2.
```

After localization at any `ell>=5`, `24` is a unit, so the actual prefix rows generate all of `Z_(ell)^2`.  QED.

### Consequence

Adding the *entire distinguished Apéry prefix* cannot create coefficient-side `p`- or `q`-primary Smith torsion for top primes.  This is an actual-orbit obstruction, not a free-vector model.

In particular, every actual boundary row can be recovered from the first three after inverting `24`.  Explicitly,

```text
24*H_r
 = (37*b_(r+1) - 3*b_r)*H_0
   + (-22*b_(r+1) + 20*b_r)*H_1
   + (b_(r+1)-b_r)*H_2.
```

Thus a prefix-aware linear presentation has characteristic-zero rank at most two in the boundary direction, and no marked prime can enter through its coefficient lattice.

---

# 5. The exact two-hit determinant: K_h

For the two marked coefficient rows `v_0` and `v_h`, the determinant is

```text
K_h
 = det(v_0,v_h)
 = b_(h+1) - 5*b_h.
```

Equivalently the two marked boundary values at `Y=q` are

```text
H_0(q) = 5*A - B,
H_h(q) = b_(h+1)*A - b_h*B,
```

where

```text
A = N_s(q),
B = (q+1)^3*N_(s-1)(q+1).
```

The coefficient matrix is

```text
C_h = [[5,       -1],
       [b_(h+1), -b_h]],

det(C_h) = K_h.
```

The exact row-elimination identities are

```text
H_h(q) - b_h*H_0(q) = K_h*A,

5*H_h(q) - b_(h+1)*H_0(q) = K_h*B.
```

This is the clean actual-orbit version of the proposed mixed boundary elimination.

---

# 6. The common packet is nonzero modulo both marked primes

The previous determinant is useful only after checking that `(A,B)` is not the zero vector modulo either marked prime.  Here this follows from the exact continuant transfer, with no synthetic-state assumption.

For `s>=2`, let

```text
T_j(a) = [[P(a+j), -(a+j)^6],
          [1,             0]],
```

and

```text
M_s(a) = T_(s-1)(a) * ... * T_1(a).
```

A direct continuant induction gives

```text
M_s(a)
 = [[N_s(a),       -(a+1)^6*N_(s-1)(a+1)],
    [N_(s-1)(a),  -(a+1)^6*N_(s-2)(a+1)]].
```

Also

```text
det M_s(a) = product_(j=1..s-1) (a+j)^6.
```

Take `a=q`.  For `ell=q`, every factor satisfies

```text
q+j == j mod q,
1 <= j <= s-1 < q.
```

For `ell=p`,

```text
q+j == h+j mod p,
1 <= h+j <= h+s-1 < p.
```

Thus the determinant is a unit modulo both `p` and `q`.  Also `q+1` is a unit modulo both primes.  If both

```text
A = N_s(q)
```

and

```text
B = (q+1)^3*N_(s-1)(q+1)
```

were zero modulo either marked prime, the first row of `M_s(q)` would vanish modulo that prime, contradicting the determinant-unit statement.

Therefore

```text
(A,B) != (0,0) mod p,
(A,B) != (0,0) mod q.
```

For `s=1`, this is immediate from `A=N_1(q)=1`.

---

# 7. Actual-prefix K_h gate theorem

We can now give the exact local Smith/support classification.

## Theorem 7.1 — generic branch is exactly diagonal

Assume

```text
gcd(K_h,p*q)=1.
```

Then

```text
p | H_h(q),
p does not divide H_0(q),

q | H_0(q),
q does not divide H_h(q).
```

### Proof

Modulo either marked prime, `C_h` is invertible because its determinant `K_h` is a unit.  The packet `(A,B)` is nonzero by Section 6.

At `p`, if `H_0(q)` also vanished, both rows of `C_h` would annihilate `(A,B)`, forcing `(A,B)=0`, contradiction.  Hence `H_0(q)` is a `p`-unit.

The argument at `q` is identical.  QED.

This means that in the generic branch the raw two-scalar presentation

```text
diag(H_0(q), H_h(q))
```

may well have a largest **global** Smith invariant containing `p*q`; for a model `diag(q*u,p*v)` that is exactly what ordinary Smith normalization does.  But that `pq` is merely the direct sum of the two vertexwise cyclic relations.  After saturating the forced endpoint factors `q` from the first scalar and `p` from the second, no marked factor remains.

So the generic branch fails condition 2 in Q2064 even if one packages it into one large Smith invariant.

## Theorem 7.2 — singular branch produces genuine mixed support

If

```text
p | K_h,
```

then

```text
p*q | H_0(q) = F_(s,0)(q).
```

If

```text
q | K_h,
```

then

```text
p*q | H_h(q) = F_(s,h)(p).
```

### Proof

First suppose `p|K_h`.  Since `h+1<p`, the two consecutive Apéry values `b_h,b_(h+1)` cannot both vanish modulo `p`.  If `p|b_h`, then `p|K_h` would force `p|b_(h+1)`, impossible.  Hence `b_h` is a `p`-unit.

The identity

```text
H_h(q) = b_h*H_0(q) + K_h*A
```

then gives `p|H_0(q)` from `p|H_h(q)`.  The original hit already gives `q|H_0(q)`, so `pq|H_0(q)`.

For `q|K_h`, the same identity and `q|H_0(q)` immediately give `q|H_h(q)`, while the original hit gives `p|H_h(q)`.  QED.

This is a genuine actual-prefix carrier on the singular subfamily: after dividing the ordinary `q`-endpoint content from `H_0`, an extra `p` remains, or vice versa.

But it is not universal.

---

# 8. K_h is nonzero and has only linear exponential height

The singularity gate is a real nonzero integer.

Define

```text
K_n = b_(n+1) - 5*b_n.
```

Then

```text
K_0 = 0,
K_1 = 48,
K_n > 0 for every n >= 1.
```

A short induction proves positivity.  If

```text
b_n >= 5*b_(n-1),
```

then the Apéry recurrence gives

```text
(n+1)^3*b_(n+1)
 >= (P(n) - n^3/5)*b_n.
```

The exact polynomial surplus is

```text
5*P(n) - n^3 - 25*(n+1)^3
 = 12*n*(12*n^2 + 15*n + 5) > 0
```

for `n>=1`.  Hence `b_(n+1)>5b_n`.  The base `b_1=5b_0` starts the induction.

The elementary binomial estimate

```text
b_n <= (n+1)*64^n
```

therefore yields

```text
0 < K_h < b_(h+1)
          <= (h+2)*64^(h+1),
```

so

```text
log K_h <= (h+1)*log 64 + log(h+2) = O(h).
```

This is much smaller than a generic boundary polynomial, but still not polylogarithmic for growing `h`.

---

# 9. Aggregate arithmetic of the K_h singular branch

For a fixed gap `h`, every singular first-hit pair has either

```text
p | K_h
```

or

```text
q | K_h.
```

All top primes are greater than `N`.  Therefore

```text
# {ell>N prime : ell | K_h}
 <= log(K_h)/log(N)
 = O(h/log N).
```

A prime divisor `p` determines at most one candidate pair `(p,p+h)`, and a prime divisor `q` determines at most one candidate pair `(q-h,q)`.  Hence

```text
# singular first-hit pairs at gap h
 <= 2*log(K_h)/log(N)
 = O(h/log N).
```

Summing over `h<=N` gives only

```text
O(N^2/log N).
```

This does **not** reach the first-hit support target from the corrected diagonal audit, which is on the scale

```text
o(N^2/log^2 N)
```

for unweighted pair support, equivalently

```text
o(N^2/log N)
```

at logarithmic prime weight.  The elementary `K_h` height argument is short by a logarithm even on the singular subfamily.

More importantly, the generic branch `gcd(K_h,pq)=1` remains and contains actual first hits, so even a perfect singular-branch estimate would not close the problem.

---

# 10. Entire actual prefix kills coefficient-side K_h torsion

The `K_h` gate is the determinant of the **two selected coefficient rows** `v_0,v_h`.  It is not a Smith factor of the whole actual prefix coefficient module.

Section 4 proves that the complete prefix matrix has second determinantal divisor dividing `24`.  Thus for any top prime `ell>=7`, even when

```text
ell | K_h,
```

some other true early Apéry row restores full coefficient rank modulo `ell`.

This distinction is important:

```text
K_h singularity
= degeneracy of the selected 2-row subpresentation,

not

a marked elementary divisor of the full actual-prefix coefficient lattice.
```

So an actual-prefix Smith/Fitting construction cannot obtain a universal marked divisor merely by stacking more true prefix rows.  Any new marked support must come from **evaluation-side arithmetic**, not from coefficient rank.

This is the promised compatible-family obstruction that preserves the distinguished Apéry orbit: the obstruction is the exact `24`-Bezout identity of the true first three Apéry rows.

No arbitrary boundary vector is introduced anywhere in this argument.

---

# 11. Normalized lift quotients: exact mixed identity, no contraction

Set

```text
u = H_0(q)/q = F_(s,0)(q)/q,
v = H_h(q)/p = F_(s,h)(p)/p.
```

These are integers at an actual hit.

Substitute

```text
H_0(q) = q*u,
H_h(q) = p*v
```

into the two row-elimination identities.  One obtains the exact mixed integer relations

```text
K_h*A = p*v - b_h*q*u,

K_h*B = 5*p*v - b_(h+1)*q*u.
```

The relation `q=p+h` changes these modulo `p` to

```text
K_h*A == -b_h*h*u           (mod p),
K_h*B == -b_(h+1)*h*u       (mod p),
```

and modulo `q` to the corresponding equations with `p==-h`.

When `K_h` is a marked unit, these equations are an invertible change between the packet coordinates and the normalized endpoint quotients.  There is no new divisibility.

The obvious universal integer

```text
J_triv = p*H_0(q) + q*H_h(q)
```

satisfies

```text
p*q | J_triv,
```

but after removing the separately vertexwise factors

```text
H_0(q)=q*u,
H_h(q)=p*v,
```

one has

```text
J_triv/(p*q) = u+v.
```

So this is exactly the forbidden endpoint-content packaging, not a primitive carrier.

The same applies to any determinant obtained by merely adjoining these two quotient coordinates: the generic coefficient change is `K_h`-unimodular at both marked primes.

---

# 12. Mixed finite differences and WZ boundary operations also saturate linearly

Let `L` be any integral linear operator on polynomials in `Y`: examples include a fixed finite difference, a fixed translated difference, or an endpoint-subtracted divided difference after exact divisibility has been established.

Because the actual prefix coefficients are constant with respect to `Y`,

```text
L(H_r)
 = b_(r+1)*L(A_s) - b_r*L(B_s).
```

Thus **the same actual coefficient matrix `V_R` survives every such operation**.

For a stack of finitely many linear operators `L_1,...,L_d`, the prefix/operator matrix factors through `d` copies of the same rank-two packet.  Locally at every top prime, the coefficient-side Smith factors are still units because `V_R` is already full over `Z_(ell)`.

Likewise a `2 x 2` determinant of two operated marked rows picks up the same actual determinant `K_h`:

```text
det([[L_1(H_0), L_2(H_0)],
     [L_1(H_h), L_2(H_h)]])

 = K_h
   * det([[L_1(A_s), L_2(A_s)],
          [L_1(B_s), L_2(B_s)]]).
```

So finite-difference Wronskians do not create a new marked factor in the generic branch.

The standard WZ certificate in `CODEX_SPEC_wz_identity.md` is linear and telescopes to the same order-two Apéry recurrence.  Its boundary dependence therefore remains in this two-dimensional prefix module.  No extra mixed-characteristic boundary equation is produced by applying that WZ telescope separately to the two hit rows.

**Scope:** this closes linear WZ-boundary / finite-difference constructions.  It does not claim to close a genuinely nonlinear multi-time WZ identity that is not generated by the separate order-two recurrence module.

---

# 13. Height audit for the actual boundary packet

The boundary packet itself is far too large to be a per-edge small-height carrier.

For positive `a`, the continuant positivity argument used elsewhere in the project gives

```text
0 < N_s(a)
  < product_(j=1..s-1) P(a+j).
```

In the current cell

```text
q+s = m <= 2N.
```

For `1<=x<=2N`, one may use the crude bound

```text
P(x) <= 936*N^3.
```

Hence

```text
log A = O(s*log N),
log B = O(s*log N).
```

Also

```text
log b_h = O(h)
```

from `b_h <= (h+1)64^h`.  Therefore

```text
log max(|H_0(q)|, |H_h(q)|)
 = O(s*log N + h).
```

There is no uniform polylogarithmic height here.  Multiplying or taking resultants of such values only makes the height worse unless a new cancellation theorem is proved.

This is consistent with the existing repository Padé-resultant audit, where adjacent resultant gcds remain large after primitive content stripping.

---

# 14. A further exact actual-orbit identity for K_n

Although `K_h` does not solve the first-hit problem, it has its own exact order-two holonomic recurrence.  This is potentially useful for a future singular-gate zero-count theorem.

Put

```text
K_n = b_(n+1) - 5*b_n,
R_n = 12*n*(12*n^2 + 15*n + 5).
```

The Apéry recurrence and the definition of `K_n` give first

```text
R_n*b_n
 = 5*(n+1)^3*K_n - n^3*K_(n-1).
```

Eliminating `b_n` between this equation and its shift by one yields, for `n>=1`,

```text
5*R_n*(n+2)^3*K_(n+1)
 = [ R_n*R_(n+1)
     + 25*R_(n+1)*(n+1)^3
     + R_n*(n+1)^3 ]*K_n
   - 5*R_(n+1)*n^3*K_(n-1).
```

with

```text
K_0 = 0,
K_1 = 48.
```

This is a new exact coordinate isolated by the present audit.  It does **not** by itself inherit the banked `ZWIN-RESTART` zero bound, because that theorem was proved for zeros of the original Apéry solution using the original restart polynomials.  Transferring the `2/3` zero-count to this derived recurrence requires a fresh restart/nonzero-polynomial proof.

A useful auxiliary target would be

```text
[K-ZWIN]
For every prime ell>=7 and H<ell,
#{1<=h<=H : ell | K_h} << H^(2/3).
```

If proved, the total `K_h`-singular top-pair support would be

```text
O( (N/log N)*N^(2/3) )
 = O(N^(5/3)/log N)
 = o(N^2/log^2 N).
```

So `[K-ZWIN]` would cleanly dispose of the singular gate.  It still would not touch the generic branch, which is the main residual.

---

# 15. The exact first uncovered object after saturation

The rank-two collapse has a clean projective formulation.

Define the actual continuant point

```text
U_s(Y) = [ A_s(Y) : B_s(Y) ] in P^1(Q).
```

The two hit equations at the common integer `Y=q` say

```text
U_s(q) == [1 : 5]                 (mod q),

U_s(q) == [b_h : b_(h+1)]         (mod p),

q = p+h.
```

Indeed

```text
H_0(q)=0 mod q
iff
5*A-B=0 mod q,
```

and

```text
H_h(q)=0 mod p
iff
b_(h+1)*A-b_h*B=0 mod p.
```

The determinant between the two target projective points is exactly `K_h`.

When `K_h` is a unit at both marked primes, the target points are transverse and the entire actual linear prefix module is locally unimodular.  Therefore the **first genuinely uncovered object** is not another linear Smith minor.  It is the mixed-characteristic incidence of one characteristic-zero projective packet with two different target points at two primes separated by `h`.

A sufficient next theorem can be stated without any synthetic states.

## [APERY-MIXED-PREFIX-INCIDENCE]

For the surviving q=1 first-hit geometry, prove

```text
sum_{p<q, q=p+h}
sum_s
  1_{ U_s(q) == [1:5] mod q }
  1_{ U_s(q) == [b_h:b_(h+1)] mod p }
  1_{ first-hit masks }

 = o(N^2/log^2 N),
```

or the corresponding weighted `o(N^2/log N)` statement.

The present theorem proves that any successful algebraic proof of this estimate must use a **nonlinear multi-time relation** or genuinely mixed-characteristic arithmetic.  Linear prefix/WZ/Smith algebra has been exhausted by the `24`-Bezout reduction.

---

# 16. Analytic alternative: exact character expansion, but no proved saving located

The root-set event has the exact additive-character expansion

```text
1_{p | b_(s+h)} * 1_{q | b_s}

 = 1/(p*q)
   * sum_(a mod p)
     sum_(c mod q)
       exp(2*pi*i*a*b_(s+h)/p)
       exp(2*pi*i*c*b_s/q),

q=p+h.
```

This identity is exact.  The difficulty is the nonzero-frequency part after summing over varying `p,q,s`.

A targeted literature check found the following relevant proved results, none of which supplies the required shifted two-prime saving:

1. M. Z. Garaev, F. Luca, I. E. Shparlinski, **Catalan and Apéry numbers in residue classes**, J. Combin. Theory Ser. A 113 (2006), 851--865, DOI `10.1016/j.jcta.2005.08.003`.  The Apéry application concerns residue-class/additive-basis properties for polynomial recurrences, not this cross-prime divisor correlation.

2. F. Luca, I. E. Shparlinski, **Arithmetic properties of Apéry numbers**, J. London Math. Soc. 78 (2008), 545--562, DOI `10.1112/jlms/jdn031`.  This gives density-one results on distinct/largest prime factors of individual Apéry numbers, not a uniform shifted pair estimate.

3. A. Malik, A. Straub, **Divisibility properties of sporadic Apéry-like numbers**, gives the Lucas-congruence framework used throughout the repository, but not mixed-prime shifted correlation.

4. X. Caruso, F. Fuernsinn, D. Vargas-Montoya, W. Zudilin, **Galois Groups of Apéry-like Series Modulo Primes**, Bull. Aust. Math. Soc. 114 (2026), 65--78, DOI `10.1017/S0004972725100932`.  This determines Galois groups of the reduction of the generating series **for a fixed prime p**.  It does not furnish a large sieve coupling two residue characteristics `p` and `p+h`.

The current repository's `SUBMIT/3.2/README.md` is explicit on the same point: its measured two-prime correlation is statistical/numerical, and the direct two-prime sheaf formulation mixes residue characteristics and does not type-check as a single fixed-field sheaf argument.

Therefore I do not bank an analytic contraction here.  The exact character expansion is available, but a proved uniform saving strong enough for Q2064 was not found.

---

# 17. Why this is stronger than the previous generic no-go

The previous endpoint-pin obstruction allowed one to worry that the failure came from treating two boundary states as unrelated formal variables.

That loophole is now closed for the linear prefix class.

Here:

- every coefficient is a genuine Apéry number from the distinguished orbit;
- all translated boundary polynomials share one genuine continuant packet;
- the first three true Apéry rows give an explicit `24`-Bezout inverse;
- the two hit rows have the exact actual determinant `K_h`;
- the packet is proved nonzero modulo both marked primes using the exact continuant transfer determinant;
- normalized lift quotients are incorporated exactly;
- finite-difference and ordinary WZ-boundary operations preserve the same coefficient module.

So this is an **actual-orbit compatible-family obstruction**, not a free-state model.

The theorem does not close nonlinear multi-time constructions, and it does not prove Problem 3.2.

---

# 18. Nonempty exact checker

The checker below uses only integer arithmetic.  It verifies:

- the Apéry recurrence;
- the exact first-hit examples `(139,181,200)` and all three marked pairs in row `321`;
- the boundary-polynomial hit equations;
- the translated common-packet identity;
- the `24`-Bezout prefix identities;
- the `K_h` elimination identities;
- the packet nonvanishing modulo both marked primes;
- the generic cross-unit theorem on the four exact examples;
- the normalized lift identities;
- the derived recurrence for `K_n`.

It is intentionally not executed here because Q2064 forbids Python/code-interpreter use.

```python
from math import gcd


def P(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery(M):
    b = [1, 5]
    for n in range(1, M):
        num = P(n)*b[n] - n**3*b[n-1]
        den = (n+1)**3
        assert num % den == 0
        b.append(num // den)
    return b


def N(j, a):
    if j == 0:
        return 0
    if j == 1:
        return 1
    n0, n1 = 0, 1
    for k in range(1, j):
        n0, n1 = n1, P(a+k)*n1 - (a+k)**6*n0
    return n1


def F(s, r, X, b):
    return (
        N(s, X+r)*b[r+1]
        - (X+r+1)**3*N(s-1, X+r+1)*b[r]
    )


def H(s, r, Y, b):
    return F(s, r, Y-r, b)


def det(u, v):
    return u[0]*v[1] - u[1]*v[0]


def kval(n, b):
    return b[n+1] - 5*b[n]


def Rcoef(n):
    return 12*n*(12*n*n + 15*n + 5)


def check_pair(p, q, m, b):
    assert p < q
    h = q-p
    s = m-q
    assert m == q+s == p+h+s
    assert 0 < s < s+h < p

    # Exact local hit and first occurrence.
    assert b[s] % q == 0
    assert b[s+h] % p == 0
    for t in range(1, s):
        assert not (b[t] % q == 0 and b[t+h] % p == 0)

    # Physical q=1 row is indeed common.
    assert b[m] % p == 0
    assert b[m] % q == 0

    # Question's boundary equations.
    f0 = F(s, 0, q, b)
    fh = F(s, h, p, b)
    assert f0 % q == 0
    assert fh % p == 0

    # Translate both to the same Y=q packet.
    Y = q
    A = N(s, Y)
    B = (Y+1)**3*N(s-1, Y+1)
    h0 = H(s, 0, Y, b)
    hh = H(s, h, Y, b)
    assert h0 == f0
    assert hh == fh
    assert h0 == 5*A - B
    assert hh == b[h+1]*A - b[h]*B

    # First three actual Apéry coefficient rows recover the packet.
    h1 = H(s, 1, Y, b)
    h2 = H(s, 2, Y, b)
    assert 24*A == 37*h0 - 22*h1 + h2
    assert 24*B == 3*h0 - 20*h1 + h2

    # Entire actual-prefix row is generated by the same rank-two packet.
    for r in range(0, s+h+1):
        hr = H(s, r, Y, b)
        rhs = (
            (37*b[r+1] - 3*b[r])*h0
            + (-22*b[r+1] + 20*b[r])*h1
            + (b[r+1] - b[r])*h2
        )
        assert 24*hr == rhs

    # Coefficient-lattice determinantal divisor is a top-prime unit.
    v0 = (b[1], -b[0])
    v1 = (b[2], -b[1])
    v2 = (b[3], -b[2])
    assert v0 == (5, -1)
    assert v1 == (73, -5)
    assert v2 == (1445, -73)
    assert det(v0, v1) == 48
    assert det(v0, v2) == 1080
    assert gcd(48, 1080) == 24
    assert 37*v0[0] - 22*v1[0] + v2[0] == 24
    assert 37*v0[1] - 22*v1[1] + v2[1] == 0
    assert 3*v0[0] - 20*v1[0] + v2[0] == 0
    assert 3*v0[1] - 20*v1[1] + v2[1] == 24

    # Exact actual determinant gate.
    K = kval(h, b)
    assert K > 0
    assert hh - b[h]*h0 == K*A
    assert 5*hh - b[h+1]*h0 == K*B

    # Common packet cannot vanish modulo either marked prime.
    for ell in (p, q):
        assert not (A % ell == 0 and B % ell == 0)

    # Generic branch: marked support is exactly diagonal.
    if K % p != 0:
        assert h0 % p != 0
    if K % q != 0:
        assert hh % q != 0

    # Singular branch implications, if a future example lands there.
    if K % p == 0:
        assert h0 % p == 0
        assert h0 % (p*q) == 0
    if K % q == 0:
        assert hh % q == 0
        assert hh % (p*q) == 0

    # Exact normalized-lift identities.
    assert h0 % q == 0
    assert hh % p == 0
    u = h0 // q
    v = hh // p
    assert K*A == p*v - b[h]*q*u
    assert K*B == 5*p*v - b[h+1]*q*u

    return {
        "p": p,
        "q": q,
        "m": m,
        "h": h,
        "s": s,
        "K_mod_p": K % p,
        "K_mod_q": K % q,
    }


def main():
    b = apery(322)

    # Exact first-hit examples.  These are not fitted into the proof;
    # they are regression tests of the generic branch.
    examples = [
        (139, 181, 200),
        (179, 193, 321),
        (193, 211, 321),
        (179, 211, 321),
    ]

    rows = [check_pair(p, q, m, b) for p, q, m in examples]

    # All four known first hits are genuinely generic K_h-unit examples.
    for row in rows:
        assert row["K_mod_p"] != 0
        assert row["K_mod_q"] != 0

    # Derived exact recurrence for K_n.
    K = [kval(n, b) for n in range(321)]
    assert K[0] == 0
    assert K[1] == 48
    for n in range(1, 319):
        Rn = Rcoef(n)
        Rnp = Rcoef(n+1)
        lhs = 5*Rn*(n+2)**3*K[n+1]
        rhs = (
            (Rn*Rnp + 25*Rnp*(n+1)**3 + Rn*(n+1)**3)*K[n]
            - 5*Rnp*n**3*K[n-1]
        )
        assert lhs == rhs

    print("Q2064_ACTUAL_PREFIX_CHECK=PASS")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
```

The four exact examples are nonempty hostile tests of the theorem's scope: they prove that the `K_h`-singular carrier cannot be promoted to a universal first-hit carrier by claiming that actual first hits force `p|K_h` or `q|K_h`.

---

# 19. Final ledger

The requested coordinates now have precise statuses.

## Full actual prefix / Bezout packet

```text
SATURATED.
```

The whole translated boundary prefix has rank two, and the first three actual Apéry rows give a `24`-Bezout inverse.  Top primes see no coefficient-side Smith torsion.

## q=p+h coupling of the two component equations

```text
EXACTLY RESOLVED.
```

It produces the actual determinant

```text
K_h = b_(h+1)-5*b_h
```

and no other linear coefficient obstruction.

## Normalized first-order quotients

```text
SATURATED IN THE GENERIC K_h-UNIT BRANCH.
```

They obey the exact mixed identities of Section 11, which are an invertible coordinate change when `K_h` is a marked unit.

## Mixed finite differences / linear WZ boundary packets

```text
SATURATED.
```

Every such operator preserves the same actual prefix coefficient lattice; `2 x 2` operated determinants only acquire the factor `K_h`.

## K_h-singular subfamily

```text
REAL CARRIER, INSUFFICIENT GLOBALLY.
```

A marked divisor of `K_h` transfers to the opposite endpoint scalar and creates genuine `pq` support, but elementary height gives only `O(N^2/log N)` support after summing gaps, and generic first hits remain.

## Analytic two-prime route

```text
UNRESOLVED.
```

The exact character expansion exists, but no sourced proved uniform cross-characteristic saving of the required strength was found.  Existing fixed-prime Apéry congruence/Galois/prime-factor results do not supply it.

## First uncovered theorem

The clean next target is a **nonlinear mixed-projective incidence theorem** for

```text
U_s(q) == [1:5] mod q,
U_s(q) == [b_h:b_(h+1)] mod p,
q=p+h,
```

with the first-hit masks retained.  A linear prefix-aware carrier cannot prove it: the actual prefix has already collapsed to a marked-unimodular rank-two module.

Therefore the correct Q2064 status is

```text
ACTUAL_SCOPED_SATURATION
```

and no claim of the full Problem 3.2 theorem is made.