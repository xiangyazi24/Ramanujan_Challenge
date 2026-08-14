# Q8377 — fixed-`g` top-strip Apéry/Jacobi packet

## Verdict

I do **not** obtain an unconditional proof of `K=o(L)`, nor a bound that closes the weighted 16-point packet gateway, from the recurrence, the level-6/Franel modular parametrization, the canonical first Gessel jet, or the currently banked horizontal Mellin/Jacobi theorems.

I do obtain a sharper and fully explicit terminal reduction for the actual packet.  Put

```text
M = m - 1 = g E,
X <= m < 2X.
```

After the usual folded top-strip branch inequalities, the two candidate maps are

```text
D: p_D(a) = m - g a          with 3a <= E,
R: p_R(a) = (m+1+g a)/2      with 3a < E and g(E+a) even.
```

They occupy disjoint prime intervals, are jointly injective, and satisfy

```text
p_D(a) | b_m  iff  p_D(a) | b_{g a},
p_R(a) | b_m  iff  p_R(a) | b_{g a}
```

for every candidate prime at least 7.  More importantly, the physical exponent is fixed across the packet:

```text
g a == +(m-1)  (mod p_D(a)-1),
g a == -(m-1)  (mod p_R(a)-1).
```

The exact Kummer orders are

```text
d_D(a) = (E-a)/gcd(E,a),
d_R(a) = (E+a)/gcd(E+a,2E).
```

Thus fixing `g` does **not** fix the character order or the cyclotomic coefficient field.  For example, when `E` is an odd prime, every admissible direct order exceeds `2E/3`, and every admissible reflected order is at least `E/2`.  The high-order sector can therefore contain the entire packet.

The decisive Apéry-specific simplification is also the obstruction:

```text
selected fixed-g packet
= congruence-restricted top-strip prime divisors of the one integer b_m.
```

The lower condition `p|b_{ga}` is an exact Lucas/reflection alias of `p|b_m`; it is not a second equation.  The recurrence gives one nonsingular zero state in each *different* characteristic.  The Franel square makes `b_{ga}` one cover-adapted coefficient, not an evaluation root.  The first Gessel jet is not forced to vanish: the actual fixed-`g` hit

```text
(m,g,E,a,p,s) = (16,5,3,1,11,5)
```

has `v_11(b_5)=1`, and the banked physical Racah/Gessel transport has first-primary valuation exactly one.  Hence no automatic `p^2` packet, derivative zero, or second selected prime ideal is available.

There is one unconditional packet-size estimate.  If `A` is any set of admissible parameters, `K` is the number of selected branch labels over `A`, and all candidates lie in the top strip, then

\[
 K\log\frac{m+1}{2}
 \le
 \min\left\{
   \log b_m,
   \sum_{a\in A}\log b_{ga}
 \right\}.
 \tag{V1}
\]

Using only the Apéry sum,

\[
 b_n\le(n+1)64^n,
\]

so

\[
 K
 \le
 \frac{
 \min\left\{
 6m\log2+\log(m+1),
 6g\log2\sum_{a\in A}a+\sum_{a\in A}\log(ga+1)
 \right\}
 }{\log((m+1)/2)}.
 \tag{V2}
\]

For an initial packet `A={1,...,A_0}` this is

\[
 K\ll \min\left\{\frac{X}{\log X},
                 \frac{gA_0^2}{\log X}+A_0\right\}.
 \tag{V3}
\]

Thus it proves `K=o(A_0)` only in the very short regime `gA_0=o(log X)`.  It gives no saving for the polynomial-size packets used by the 16-point gateway.

The exact magnitude comparison explains why no norm/height argument can repair this.  Every prime in the strict top strip divides the central binomial carrier,

\[
 \prod_{(m+1)/2<p\le m}p
 \mid \binom m{\lfloor m/2\rfloor}<2^m,
\]

whereas the single `k=m` Apéry summand gives

\[
 b_m\ge\binom{2m}{m}^2
     \ge \frac{16^m}{(2m+1)^2}.
\]

So the actual Apéry integer has ample exponential height to contain **every** candidate top-strip prime.  This does not say that it does; it proves that size, purity, conductor, a cyclotomic norm, and any other one-point height argument cannot force even one omission.

The smallest honest new horizontal statement is therefore the following weighted 16-point assertion.  Let `I_{m,g}` be the set of distinct admissible branch labels, let `epsilon_i=1` exactly when the corresponding candidate prime divides `b_{ga}`, and let `w_i>=0` be the gateway weights.  For

\[
 e_{16}(x_i)=\sum_{i_1<\cdots<i_{16}}x_{i_1}\cdots x_{i_{16}},
\]

one needs

\[
 \boxed{
 e_{16}(w_i\epsilon_i)
 =o\bigl(e_{16}(w_i)igr)
 }
 \tag{FGJ16}
\]

uniformly over the concentrated fixed-`g` packets consumed by the gateway.  This is strictly weaker than `K=o(L)`: a packet with at most fifteen selected labels has zero left side.  It is also genuinely horizontal—sixteen different defining characteristics occur—and it is not a Kummer norm or a local conductor statement.

No theorem currently banked in the repository proves `(FGJ16)`.  Proving it, or the stronger weighted one-point estimate `K=o(L)`, requires a canonical Apéry-specific correlation across the selected Teichmüller prime ideals as `p` varies.  That is the precise remaining gate.

---

## 1. Source state and ownership

The audit uses the connector-visible states

```text
main@734a5a84c1e4fd8703a811aadaa2b4c7f532b20e
chatgpt-drop@0d452109f1028cc0112cd8bf95da237efaa6e54f
```

and in particular:

```text
problems/3.2/ORACLE_COMM/chatgpt_q32_padic_horizontal_mellin.md
problems/3.2/ORACLE_COMM/chatgpt_q32_fixed_exponent_mellin_attack.md
problems/3.2/hasse_franel_descent.tex
problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_obstruction.md
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_obstruction.md
drops/Q8191-9171160c.md
drops/Q8193-0879e87d.md
```

No shared manuscript, doctrine, run log, or existing verifier is edited.  The owned Q8377 artifacts are

```text
problems/3.2/ORACLE_COMM/chatgpt_q32_fixedg_jacobi_report.md
problems/3.2/ORACLE_COMM/chatgpt_q32_fixedg_jacobi_verify.py
```

and the required delivery file.

---

# I. Exact packet algebra

## 2. Admissible branches

Put

\[
 M=m-1=gE,
 \qquad s_a=ga.
\]

The top-strip folded direct branch is

\[
 p_D(a)=m-s_a=g(E-a)+1.
 \tag{2.1}
\]

The residue of `m` modulo `p_D(a)` is `s_a`.  The folded condition is

\[
 2s_a\le p_D(a)-1,
\]

which is exactly

\[
 3a\le E.
 \tag{2.2}
\]

The reflected branch is

\[
 p_R(a)=\frac{m+1+s_a}{2}
       =\frac{g(E+a)+2}{2}.
 \tag{2.3}
\]

It exists as an integer exactly when `g(E+a)` is even.  Here

\[
 m=p_R(a)+(p_R(a)-1-s_a),
\]

and the strict folded reflected condition is

\[
 2s_a<p_R(a)-1,
\]

which is exactly

\[
 3a<E.
 \tag{2.4}
\]

The equality case is assigned to the direct branch, matching the branch convention in the fixed-exponent reduction.

### Proposition 2.1 — disjoint intervals and joint injectivity

For `a>=1`, every direct candidate satisfies

\[
 \frac{2m+1}{3}\le p_D(a)<m,
 \tag{2.5}
\]

while every reflected candidate satisfies

\[
 \frac{m+1}{2}<p_R(a)<\frac{2m+1}{3}.
 \tag{2.6}
\]

Both maps are injective in `a`, and their images are disjoint.

**Proof.**  From `3a<=E`,

\[
 p_D(a)=g(E-a)+1\ge\frac{2gE}{3}+1=\frac{2m+1}{3}.
\]

Since `a>=1`, `p_D(a)<m`.  From `3a<E`,

\[
 p_R(a)>\frac{gE+2}{2}=\frac{m+1}{2}
\]

and

\[
 p_R(a)<\frac{g(4E/3)+2}{2}=\frac{2m+1}{3}.
\]

The direct map has slope `-g`; the reflected map has slope `g/2` on its parity-compatible domain.  The disjoint intervals finish the proof. ∎

This removes a possible hidden multiplicity.  Sixteen distinct packet labels really do mean sixteen distinct rational primes and sixteen distinct defining characteristics.

## 3. Exact Apéry selection

### Theorem 3.1 — packet Lucas/reflection identity

Let a candidate prime `p>=7` occur on either branch.  Then

\[
 \boxed{p\mid b_m\iff p\mid b_{ga}.}
 \tag{3.1}
\]

**Direct proof.**  Since `m=p+ga` and both base-`p` digits are below `p`, Gessel--Lucas gives

\[
 b_m\equiv b_1b_{ga}=5b_{ga}\pmod p.
\]

The factor 5 is a unit for `p>=7`.

**Reflected proof.**  Here

\[
 m=p+(p-1-ga).
\]

Lucas and Apéry reflection give

\[
 b_m\equiv5b_{p-1-ga}\equiv5b_{ga}\pmod p.
\]

Again 5 is a unit. ∎

The prime 5 is a finite exceptional characteristic and belongs to none of the asymptotic packet issues.

### Corollary 3.2 — the packet is a restricted factorization of `b_m`

The selected direct primes are exactly

\[
 \left\{
 p:\ \frac{2m+1}{3}\le p<m,
 \ p\equiv1\pmod g,
 \ p\text{ prime},
 \ p\mid b_m
 \right\},
 \tag{3.2D}
\]

with the endpoint/parameter convention above.  The selected reflected primes are exactly

\[
 \left\{
 p:\ \frac{m+1}{2}<p<\frac{2m+1}{3},
 \ 2p\equiv m+1\pmod g,
 \ p\text{ prime},
 \ p\mid b_m
 \right\}.
 \tag{3.2R}
\]

Thus `p|b_{ga}` is not a second condition after the physical `m` is fixed.  It is the canonical lower-index explanation of the same prime factor of `b_m`.

---

# II. Fixed exponent and selected prime ideals

## 4. One exponent, varying orders

For the direct branch,

\[
 ga=M-(p_D(a)-1),
\]

and for the reflected branch,

\[
 ga=2(p_R(a)-1)-M.
\]

Hence

\[
 \boxed{
 ga\equiv M\pmod{p_D(a)-1},
 \qquad
 ga\equiv-M\pmod{p_R(a)-1}.
 }
 \tag{4.1}
\]

This is the exact fixed-exponent/half-orbit feature of the packet.

Let `omega_p` be the chosen Teichmüller generator and put

\[
 d_p=\operatorname{ord}(\omega_p^M)
    =\frac{p-1}{\gcd(p-1,M)}.
\]

### Proposition 4.1 — exact order formulas

On the two branches,

\[
 \boxed{
 d_D(a)=\frac{E-a}{\gcd(E,a)},
 \qquad
 d_R(a)=\frac{E+a}{\gcd(E+a,2E)}.
 }
 \tag{4.2}
\]

**Proof.**  Directly,

\[
 p_D(a)-1=g(E-a),
\]

so

\[
 \gcd(p_D(a)-1,M)=g\gcd(E-a,E)=g\gcd(E,a).
\]

For the reflected branch, `2(p_R(a)-1)=g(E+a)`.  Therefore

\[
 2\gcd(p_R(a)-1,gE)
 =g\gcd(E+a,2E),
\]

which gives the second formula. ∎

### Corollary 4.2 — fixing `g` does not bound the Kummer order

Suppose `E` is an odd prime.  For every `1<=a<E/3`,

\[
 d_D(a)=E-a>\frac{2E}{3}.
\]

Also `gcd(E+a,E)=1`, so `gcd(E+a,2E)` is 1 or 2 and

\[
 d_R(a)\ge\frac{E+a}{2}>\frac E2
\]

whenever the reflected candidate is integral.

Thus the entire packet may remain in the high-order sector.  The coefficient fields

\[
 \mathbb Q(\mu_{d_D(a)}),
 \qquad
 \mathbb Q(\mu_{d_R(a)})
\]

vary with `a` and have unbounded degree.  The congruence `p_D(a)=1 mod g` does not put the relevant Mellin traces in one fixed cyclotomic field: their minimal orders are the numbers in (4.2), not `g`.

The canonical Teichmüller prime ideal is therefore a moving prime in a moving field.  Complex conjugation accounts for the `+M/-M` half-orbit, but supplies no second independent prime ideal condition.  This is exactly the scope in which generic Jacobi/Kummer half-orbit models can saturate every candidate.

---

# III. What the special Apéry structures actually add

## 5. The Apéry sum: all minimum-slope terms are units

For every admissible label one has `2ga<=p-1`.  Hence for every `0<=k<=ga`,

\[
 ga+k<p.
\]

All factorials in

\[
 b_{ga}
 =\sum_{k=0}^{ga}
 \left(\binom{ga}{k}\binom{ga+k}{k}\right)^2
 \tag{5.1}
\]

are `p`-units.  Equivalently, the Gross--Koblitz/Morita carry polygon has exactly `ga+1` minimum-slope terms and no positive-slope tail in the strict folded range.  Therefore

\[
 p\mid b_{ga}
\]

is pure additive cancellation among `ga+1` nonzero units in `F_p`.

This is an Apéry-specific exact statement, but it goes in the wrong direction for a slope sieve: the number of competing units grows with `a`.  There is no unique-minimum exclusion, no bounded number of Jacobi summands, and no valuation gap.

## 6. The recurrence: one zero state in one characteristic

Write

\[
 P(n)=34n^3+51n^2+27n+5
\]

and

\[
 \binom{b_{n+1}}{b_n}
 =T_n\binom{b_n}{b_{n-1}},
 \qquad
 T_n=
 \begin{pmatrix}
 P(n)/(n+1)^3&-n^3/(n+1)^3\\
 1&0
 \end{pmatrix}.
 \tag{6.1}
\]

Its determinant is

\[
 \det T_n=\frac{n^3}{(n+1)^3}.
 \tag{6.2}
\]

At a packet prime `p>2ga`, every coefficient needed to propagate through the selected lower index is a `p`-unit.  If `p|b_{ga}`, the recurrence state is one projective zero line, while `b_{ga-1}` and `b_{ga+1}` are nonzero modulo `p`; otherwise backward propagation would force `b_0=0`.

Joint injectivity from Proposition 2.1 means no other packet label is forced in this same characteristic.  Thus a fixed-gap continuant, a Wronskian, or a resultant between two decimated rows receives only one forced zero.  For two direct labels,

\[
 p_D(a)-p_D(b)=g(b-a),
\]

and for two reflected labels,

\[
 2(p_R(a)-p_R(b))=g(a-b),
\]

but these are relations between *different characteristics*.  They do not turn `p_D(a)|b_{ga}` into `p_D(a)|b_{gb}`.  This is precisely the mixed-characteristic obstruction isolated in Q8345.

The `g`-step decimation `c_a=b_{ga}` is indeed holonomic: multiplying the matrices (6.1) over one block gives a rational two-dimensional transfer with polynomially clearable denominators.  What is missing is not a recurrence in `a`; it is a theorem controlling divisibility of its `a`th value by a different affine prime `p(a)` at every step.

## 7. Franel/Hasse modular parametrization: coordinate, not point

The repository proves

\[
 K_p(x)^2
 =\sum_{r=0}^{p-1}b_r\Psi_{p,r}(x),
 \qquad
 \Psi_{p,r}(x)=x^r(1-8x)^r(1+x)^{p-1-r},
 \tag{7.1}
\]

and the `Psi_{p,r}` are a triangular basis.  Thus

\[
 [K_p^2]_{\Psi_{p,r}}=b_r.
 \tag{7.2}
\]

The selected condition is one **coordinate zero** of the Hasse square.  It is not a zero of `K_p`, of `A_p`, or of the level-6 modular parameter at a point labeled by `r`.

The exact actual-data counterexample is

```text
p = 11, s = 5,
11 | b_5,
tau(5) = 8 mod 11,
A_11(8) = 5 mod 11.
```

Root-of-unity filtering on the multiples of `g` does not change the type of the statement.  When `mu_g` is available it projects (7.1) to the coordinates with `g|r`; it still does not convert the vanishing of one coordinate into an evaluation root.  Moreover both the Hasse polynomial and the cover-adapted basis vary with `p`.

So the special modular parametrization supplies a precise coordinate model, but no selected-prime-ideal restriction across the packet.

## 8. Canonical first Gessel jet: no automatic second primary condition

The fixed-`g` packet itself contains the banked first-primary counterexample

\[
 (m,g,E,a,p,s)=(16,5,3,1,11,5).
 \tag{8.1}
\]

Indeed

\[
 p=m-ga=11,
 \qquad
 3a=E,
 \qquad
 b_5=819005=11\cdot74455,
 \qquad
 v_{11}(b_5)=1.
\]

The physical triangular relation is

\[
 \tau(16)-\tau(5)=242=22\cdot11.
\]

Q8336 verifies for the canonical integral physical Racah/Gessel transport that the transported value is `77 mod 121`, hence also has valuation exactly one.  Q8193 proves more generally that the first root/jet lattice aliases the complete target product by CRT; a nonzero first-jet defect is divisible by, and at least as large as, that old product.

Therefore the zeroth Gessel/Lucas condition does not force the first divided jet to vanish.  Any use of the canonical jet must prove a genuinely new horizontal estimate for its first-primary conormal classes; merely adjoining the jet to each label neither raises the prime depth nor compresses the product.

For a fixed set of sixteen packet primes, any fixed number `h` of local `p`-adic digits has CRT modulus at most

\[
 \left(\prod_{j=1}^{16}p_j\right)^h=X^{O_h(1)},
\]

whose logarithm is `O_h(log X)`.  This is negligible beside the characteristic-zero height `log b_m=Theta(X)`.  Thus a fixed-depth 16-point Gessel packet cannot contradict the size of the canonical integer; it needs a nonzero low-height identity coupling the sixteen characteristics.

---

# IV. The strongest unconditional count and its stopping point

## 9. Product divisibility

Let `I` be any set of admissible prime labels, and let `S` be its selected subset.  All candidate primes are distinct.  Put

\[
 P_S=\prod_{i\in S}p_i,
 \qquad K=|S|.
\]

By Theorem 3.1,

\[
 \boxed{P_S\mid b_m.}
 \tag{9.1}
\]

If `A` is the set of parameter values occurring in `I`, then all selected primes attached to one `a` are distinct divisors of `b_{ga}`.  Therefore

\[
 \boxed{P_S\mid\prod_{a\in A}b_{ga}.}
 \tag{9.2}
\]

Since every candidate prime is greater than `(m+1)/2`, (9.1)--(9.2) give (V1).

The elementary estimate

\[
 b_n\le(n+1)64^n
 \tag{9.3}
\]

follows from

\[
 \binom nk\le2^n,
 \qquad
 \binom{n+k}{k}\le4^n.
\]

Substituting (9.3) proves (V2) and (V3).

This is a real selected-prime restriction.  It is useful only when the sum of the lower indices is sublogarithmic relative to the packet length.  It does not reach a concentrated packet whose parameter length is a fixed power of `X`.

## 10. Exact height slack

Let

\[
 P_{\rm top}(m)=\prod_{(m+1)/2<p\le m}p.
\]

Every such prime occurs once in the numerator of

\[
 \binom m{\lfloor m/2\rfloor}
\]

and not in either denominator factorial.  Hence

\[
 P_{\rm top}(m)
 \mid\binom m{\lfloor m/2\rfloor}<2^m.
 \tag{10.1}
\]

On the other hand, the `k=m` summand in the Apéry sum is

\[
 \binom{2m}{m}^2.
\]

Since the largest coefficient in `(1+1)^{2m}` is at least the average,

\[
 \binom{2m}{m}\ge\frac{4^m}{2m+1},
\]

and therefore

\[
 b_m\ge\frac{16^m}{(2m+1)^2}.
 \tag{10.2}
\]

Equations (10.1)--(10.2) are an exact death certificate for pure size arguments.  The common Apéry value has exponentially more room than the whole strict top-strip primorial.  A cyclotomic norm enlarges the available height further because its degree is `phi(d_p)`.

---

# V. Smallest surviving horizontal theorem

## 11. Weighted packet notation

Let

\[
 \mathcal I_{m,g}
 =\{(D,a)\text{ admissible prime}\}
  \sqcup
  \{(R,a)\text{ admissible prime}\}.
\]

For `i=(sigma,a)` write

\[
 p_i=p_\sigma(a),
 \qquad s_i=ga,
 \qquad
 \epsilon_i=\mathbf1_{p_i\mid b_{s_i}}.
\]

By Theorem 3.1, `epsilon_i` is also `1_{p_i|b_m}`.  Give the packet nonnegative gateway weights `w_i` and put

\[
 L_w=\sum_iw_i,
 \qquad
 K_w=\sum_iw_i\epsilon_i.
\]

The strong desired statement is

\[
 K_w=o(L_w).
 \tag{11.1}
\]

Nothing proved above yields (11.1) at gateway scale.

## 12. Fixed-`g` Jacobi 16-point dispersion

Define the weighted elementary symmetric sum

\[
 e_{16}(x_i)
 =\sum_{i_1<\cdots<i_{16}}
   x_{i_1}\cdots x_{i_{16}}.
\]

### `FGJ16` — smallest gateway-specific new statement

Uniformly for every concentrated fixed-`g` packet admitted by the top-strip decomposition,

\[
 \boxed{
 e_{16}(w_i\epsilon_i)
 =o\bigl(e_{16}(w_i)\bigr).
 }
 \tag{FGJ16}
\]

This is the exact selected-prime-ideal restriction needed at the first surviving arity:

- it uses the canonical Apéry events, not a generic Kummer model;
- it keeps the direct/reflected branch, parity, primality, and folding masks;
- its sixteen labels are sixteen distinct defining characteristics by Proposition 2.1;
- it is strictly weaker than `K=o(L)`, because up to fifteen selected labels are invisible;
- it directly matches the weighted 16-point gateway rather than passing through an unsaturated norm.

In Mellin language, `(FGJ16)` says that the sixteen selected Teichmüller prime ideals cannot simultaneously divide the sixteen canonical fixed-exponent Mellin traces with the generic frequency predicted by the ambient packet.  The exact rational reduction says the same thing without moving coefficient fields:

\[
 e_{16}\left(
 w_i\mathbf1_{p_i\mid b_{ga_i}}
 \right)
 =o(e_{16}(w_i)).
\]

No one-characteristic recurrence identity, no all-character equidistribution theorem, and no fixed-depth Gessel jet controls this expression.  A proof must introduce a genuinely horizontal canonical relation among the sixteen Apéry coefficient cancellations.

A stronger but simpler sufficient theorem is the one-point estimate (11.1).  A two-point decorrelation theorem would also imply `(FGJ16)` under the usual weight regularity, but it is not presently available and is stronger than the gateway asks for.

---

## 13. Final status

### Proved

1. Exact branch ranges, parity condition, and joint injectivity.
2. Exact Lucas/reflection equivalence `p|b_m iff p|b_{ga}`.
3. Exact fixed-exponent congruences and Kummer-order formulas.
4. Growing-order obstruction even for fixed `g`.
5. Apéry-specific unit-cancellation normal form with `ga+1` minimum terms.
6. Nonsingular recurrence-state obstruction to same-characteristic resultants.
7. Franel coordinate/evaluation type mismatch and the actual `(p,s)=(11,5)` counterexample.
8. Fixed-`g` simple-primary counterexample to an automatic Gessel-jet lift.
9. Product bounds (V1)--(V3), including the only unconditional short-packet saving.
10. Exact exponential height slack showing that the common `b_m` can accommodate the full top-strip primorial.

### Not proved

Neither `K=o(L)` nor `(FGJ16)`.

### Smallest new arithmetic input

`(FGJ16)`: weighted 16-point horizontal anti-concentration for the actual canonical fixed-`g` Apéry packet.  It is not another Kummer norm.  It is the first statement that couples the selected defining-prime ideals rather than studying each local Mellin/Jacobi object separately.

The sibling verifier `chatgpt_q32_fixedg_jacobi_verify.py` checks every finite algebraic identity used above, the exact branch partition, the Lucas/reflection equivalence, the order formulas, product divisibility, the binomial height slack, and the `(16,5,3,1,11,5)` / `(11,5)` regressions.