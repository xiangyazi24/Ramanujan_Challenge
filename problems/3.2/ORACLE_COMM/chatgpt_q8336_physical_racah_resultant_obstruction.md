# Q8336 — physical-`n` Racah/resultant obstruction for the triangular core

## Verdict

I do **not** obtain an unconditional proof of `T_n=o(n)`, and I do not have a known theorem whose verified hypotheses imply a strictly smaller positive target. I do obtain option **(3)** in the question: a theorem-level obstruction to the most natural “retain the common physical `n` and take a resultant/height certificate” route.

The obstruction is exact and has three layers.

1. **Physical-point transport is an exact alias of the old congruence.**  There is a canonical integral Racah/Newton polynomial `Rhat_s(Y)` such that

   ```text
   Rhat_s(tau(s)) = (s!)^2 b_s.
   ```

   For every `n,s`, polynomial division gives

   ```text
   Rhat_s(tau(n))
     = (s!)^2 b_s
       + (tau(n)-tau(s)) Qhat_s(tau(n)),       Qhat_s in Z[Y].
   ```

   Consequently

   ```text
   gcd(Rhat_s(tau(n)), tau(n)-tau(s))
   = gcd((s!)^2 b_s, tau(n)-tau(s)).           (A)
   ```

   On the range relevant to `T_n`, where `p>2s`, the factor `s!` is a `p`-unit. Thus replacing `p|b_s` by `p|Rhat_s(tau(n))` while keeping `p|tau(n)-tau(s)` gives **exactly the same prime incidence**, not a smaller one.

2. **A target prime belongs to only one folded row.**  If `2s<p`, `2t<p`, and

   ```text
   tau(s) == tau(n) == tau(t)  (mod p),
   ```

   then `s=t`. Therefore a pairwise or fixed-depth resultant between different transported rows does not receive two forced zeros from one `T_n` prime. This is the precise “one coincidence per characteristic” obstruction in the present triangular coordinates.

3. **The adjacent row is provably a unit, not a hidden second equation.**  The monic Racah rows satisfy

   ```text
   P_{j+1}(Y)
     = (Y + (j^3+(j+1)^3)/(2(2j+1))) P_j(Y)
       - j^6/(4(4j^2-1)) P_{j-1}(Y),
   ```

   and hence

   ```text
   Res(P_s,P_{s-1})
     = (-1)^(s(s-1)/2)
       product_{1<=j<s} (j^6/(4(4j^2-1)))^j.  (B)
   ```

   Every prime divisor of the numerator or denominator on the right is at most `2s-1`. Hence for every prime `p>2s`, `P_s` and `P_{s-1}` have no common root modulo `p`. At a `T_n` hit,

   ```text
   p | R_s(tau(n)),
   p ∤ R_{s-1}(tau(n)).                         (C)
   ```

   Thus the obvious two-row Wronskian/resultant companion loses the target rather than doubling it.

This is a genuine obstruction to a concrete proposed route. It is not another radical/gcd rewrite: it proves that the canonical common-physical-point lift has **no additional local algebraic content** and that the first natural second row is actually transverse.

There is a parallel obstruction to the level-6/Hasse idea. The exact Franel pullback in the repository makes `b_s` a **coordinate/Mellin coefficient** of the Hasse square, not an evaluation value. The congruence `tau(n)==tau(s) mod p` transports evaluation points; it does not transport coefficient coordinates. The finite identity

```text
p=11, s=5:
11 | b_5,
tau(5) = 8 (mod 11),
A_11(8) = 5 (mod 11)
```

is an exact counterexample to the nearby stronger claim that a coefficient zero becomes an evaluation root at the triangular point.

Finally, even the first-primary strengthening is false automatically. The triple

```text
(n,s,p)=(16,5,11)
```

is a genuine `T_16` incidence, but

```text
Rhat_5(tau(16)) = 77 (mod 11^2),
```

so its transported value has `v_11=1`, not `>=2`. The joint equations do not manufacture a second `p`-adic digit.

The current top-strip four-prime rectangle gateway is therefore **not replaced by a smaller positive theorem here**. This result sits strictly *below* that gateway: it proves that a one-prime/one-row or adjacent-two-row common-`n` resultant cannot be the missing shortcut. A genuine four-prime rectangle couples different characteristics and escapes exactly the one-characteristic uniqueness statement above. So the rectangle gateway remains the first surviving type of positive cross-prime input; the present result is a scoped no-go, not a competing sufficient theorem.

---

## 1. Source state and scope

The connector-visible repository state used for this audit is

```text
main@734a5a84c1e4fd8703a811aadaa2b4c7f532b20e
```

and the current `chatgpt-drop` branch. I read in particular:

```text
problems/3.2/hasse_franel_descent.tex
drops/Q8058-912cda27.md
drops/Q8191-9171160c.md
drops/Q8193-0879e87d.md
drops/Q8239-a15a4d3f.md
```

The displayed definition of `T_n` and the proved reduction

```text
0 <= log R_n - T_n << n^(2/3)
```

are taken as the authoritative current interface from Q8336. I do not re-prove that reduction here and I do not use it as a claimed new result.

No shared TeX, `DOCTRINE.md`, or `RUN_LOG_P32.md` is edited. The only new owned files are this report and

```text
problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_verify.py
```

on `chatgpt-drop`.

---

# I. PROVED

## 2. The triangular condition has one folded root

Write

\[
\tau(u)=u(u+1).
\]

Then exactly

\[
\tau(n)-\tau(s)
=(n-s)(n+s+1).
\tag{2.1}
\]

For an odd prime `p`,

\[
\tau(n)\equiv\tau(s)\pmod p
\iff
s\equiv n\pmod p
\quad\text{or}\quad
s\equiv-1-n\pmod p.
\tag{2.2}
\]

If additionally `2s<p`, there is a unique representative in the folded interval. If `r=n mod p`, it is

\[
\boxed{
\rho_p(n)=\min(r,p-1-r).
}
\tag{2.3}
\]

Indeed the two roots of `X(X+1)=tau(n)` are `r` and `p-1-r`; they sum to `p-1` and exactly one lies in `0<=s<p/2`, except at the fixed midpoint where they coincide.

A form that will be used below is the following uniqueness lemma.

### Lemma 2.1 — one prime, one folded row

Let `p` be odd and let `s,t` be nonnegative integers satisfying `2s<p` and `2t<p`. If

\[
p\mid\tau(n)-\tau(s),
\qquad
p\mid\tau(n)-\tau(t),
\]

then `s=t`.

**Proof.** Subtraction gives

\[
p\mid\tau(s)-\tau(t)=(s-t)(s+t+1).
\]

We have `|s-t|<p`. Also `s+t+1<=p`; equality can occur only when
`s=t=(p-1)/2`. Hence either factor can be divisible by `p` only when
`s=t`. ∎

Thus a `T_n` prime cannot be made into a same-characteristic two-row collision merely by retaining the common physical `n`.

### Quotient form

Equation (2.1) also gives the two exact branches

\[
n-s=qp
\quad\text{or}\quad
n+s+1=qp,
\tag{2.4}
\]

with `q>=1`. Since `p^3>n^2`, `p>n^(2/3)`. Hence the quotient is on the cube-root scale; this is useful bookkeeping but does not add a second Apéry equation.

---

## 3. The integral Racah/Newton row

Put

\[
\lambda_j=j(j+1)=\tau(j),
\qquad
U_{j,k}=\binom jk\binom{j+k}{k}.
\tag{3.1}
\]

Define the Newton basis

\[
\phi_k(Y)
=\frac{\prod_{a=0}^{k-1}(Y-\lambda_a)}{(k!)^2},
\qquad \phi_0=1,
\tag{3.2}
\]

and the all-zero Racah/Wilson row

\[
R_s(Y)=\sum_{k=0}^s U_{s,k}\phi_k(Y).
\tag{3.3}
\]

At a triangular node `lambda_j`,

\[
\prod_{a=0}^{k-1}(\lambda_j-\lambda_a)
=\prod_{a=0}^{k-1}(j-a)(j+a+1),
\]

so

\[
\phi_k(\lambda_j)=
\begin{cases}
U_{j,k},&k\le j,\\
0,&k>j.
\end{cases}
\tag{3.4}
\]

Therefore

\[
\boxed{
R_s(\lambda_j)
=\sum_{k=0}^{\min(s,j)}U_{s,k}U_{j,k}
=R_j(\lambda_s).
}
\tag{3.5}
\]

At the diagonal node,

\[
\boxed{
R_s(\lambda_s)=\sum_{k=0}^sU_{s,k}^2=b_s.
}
\tag{3.6}
\]

The denominators in (3.2) are harmless for a prime `p>s`, but for the resultant obstruction it is cleaner to clear them integrally. Define

\[
\widehat R_s(Y):=(s!)^2R_s(Y).
\tag{3.7}
\]

Explicitly,

\[
\boxed{
\widehat R_s(Y)
=\sum_{k=0}^s
U_{s,k}\left(\frac{s!}{k!}\right)^2
\prod_{a=0}^{k-1}(Y-\lambda_a)
\in\mathbf Z[Y].
}
\tag{3.8}
\]

and

\[
\widehat R_s(\lambda_s)=(s!)^2b_s.
\tag{3.9}
\]

This is the canonical integral object for a physical-`n` resultant attack.

---

## 4. Main obstruction theorem: physical-point transport collapses exactly

### Theorem 4.1 — integral transport/remainder identity

For every `s>=0` there is a polynomial `Qhat_s(Y) in Z[Y]` such that for every integer `n`,

\[
\boxed{
\widehat R_s(\lambda_n)
=(s!)^2b_s
 +(\lambda_n-\lambda_s)\widehat Q_s(\lambda_n).
}
\tag{4.1}
\]

Consequently

\[
\boxed{
\gcd\!\left(\widehat R_s(\lambda_n),\lambda_n-\lambda_s\right)
=
\gcd\!\left((s!)^2b_s,\lambda_n-\lambda_s\right).
}
\tag{4.2}
\]

For every prime `p>s`,

\[
\boxed{
\begin{aligned}
p\mid\lambda_n-\lambda_s
\quad\Longrightarrow\quad
&p\mid\widehat R_s(\lambda_n)
\iff p\mid b_s,\\
&p^e\mid\gcd(\widehat R_s(\lambda_n),\lambda_n-\lambda_s)
\iff
p^e\mid\gcd(b_s,\lambda_n-\lambda_s)
\end{aligned}}
\tag{4.3}
\]

for every exponent `e>=1` for which `p` is a unit on `s!`.

**Proof.** Since `Rhat_s` is an integer polynomial, the polynomial remainder theorem over `Z[Y]` gives

\[
\widehat R_s(Y)-\widehat R_s(\lambda_s)
=(Y-\lambda_s)\widehat Q_s(Y)
\]

with `Qhat_s in Z[Y]`. Insert (3.9) and then `Y=lambda_n` to get (4.1). The elementary identity

\[
\gcd(A+BC,B)=\gcd(A,B)
\]

gives (4.2). If `p>s`, then `p` does not divide `s!`, so multiplying `b_s` by `(s!)^2` does not alter its `p`-adic valuation. ∎

### Corollary 4.2 — exact resultant collapse

Using the convention `Res_Y(Y-a,F(Y))=F(a)`,

\[
\boxed{
\operatorname{Res}_Y(Y-\lambda_s,\widehat R_s(Y))=(s!)^2b_s,
}
\tag{4.4}
\]

\[
\boxed{
\operatorname{Res}_Y(Y-\lambda_n,\widehat R_s(Y))=\widehat R_s(\lambda_n),
}
\tag{4.5}
\]

and their difference is exactly a multiple of the triangular constraint:

\[
\boxed{
\operatorname{Res}_Y(Y-\lambda_n,\widehat R_s)
-
\operatorname{Res}_Y(Y-\lambda_s,\widehat R_s)
=(\lambda_n-\lambda_s)\widehat Q_s(\lambda_n).
}
\tag{4.6}
\]

Thus a resultant that “retains the common physical `n`” has not produced a new numerator. Modulo the very equation used to transport the node, it is the old Apéry coefficient.

### Corollary 4.3 — application to `T_n`

Every prime counted in `T_n` satisfies `p>2s>s`. Hence (4.3) applies without any denominator caveat. The replacement

\[
(p\mid b_s,\ p\mid\tau(n)-\tau(s))
\quad\rightsquigarrow\quad
(p\mid\widehat R_s(\tau(n)),\ p\mid\tau(n)-\tau(s))
\]

is an exact equivalence prime by prime, including multiplicities in their common gcd. It cannot lower `T_n`.

This is stronger than saying “the same primes seem to occur”: it is an equality of integer gcds before any asymptotic estimate.

---

## 5. The adjacent row is transverse at every `T_n` prime

The rows `R_j` satisfy the denominator-free recurrence

\[
(j+1)^3R_{j+1}(Y)
=\left(j^3+(j+1)^3+2(2j+1)Y\right)R_j(Y)
-j^3R_{j-1}(Y).
\tag{5.1}
\]

Let `P_j=R_j/kappa_j` be monic. Comparing leading coefficients in (5.1) gives the monic three-term recurrence

\[
\boxed{
P_{j+1}(Y)
=
\left(Y+\frac{j^3+(j+1)^3}{2(2j+1)}\right)P_j(Y)
-
\beta_jP_{j-1}(Y),
}
\tag{5.2}
\]

where

\[
\boxed{
\beta_j=\frac{j^6}{4(4j^2-1)}.
}
\tag{5.3}
\]

Reduce (5.2) modulo `P_j`. If `alpha` ranges over the `j` roots of monic `P_j`, then

\[
P_{j+1}(\alpha)=-\beta_jP_{j-1}(\alpha).
\]

Taking the product over those roots gives

\[
\operatorname{Res}(P_{j+1},P_j)
=(-\beta_j)^j\operatorname{Res}(P_j,P_{j-1}).
\]

Starting from `Res(P_1,P_0)=1` yields

\[
\boxed{
\operatorname{Res}(P_s,P_{s-1})
=(-1)^{s(s-1)/2}
\prod_{j=1}^{s-1}
\left(\frac{j^6}{4(4j^2-1)}\right)^j.
}
\tag{5.4}
\]

Every prime occurring in the numerator or denominator of (5.4) is at most `2s-1`. Therefore, for an odd prime `p>2s`, both monic rows are `p`-integral and their resultant is a `p`-adic unit. They have no common root over `F_p`.

### Theorem 5.1 — no adjacent companion at a `T_n` hit

If

\[
p>2s,
\qquad
p\mid b_s,
\qquad
p\mid\tau(n)-\tau(s),
\]

then

\[
\boxed{
p\mid R_s(\tau(n)),\qquad p\nmid R_{s-1}(\tau(n)).}
\tag{5.5}
\]

**Proof.** The first statement is Theorem 4.1. If the second failed, then modulo `p` the common value `Y=tau(n)=tau(s)` would be a common root of `P_s` and `P_{s-1}`, contradicting the `p`-unit resultant (5.4). ∎

This is the decisive extra obstruction beyond the tautological remainder identity. A two-row Wronskian, adjacent subresultant, or Christoffel-Darboux numerator does not give a second target-divisible scalar here: the adjacent row is forced to be nonzero.

---

## 6. Why fixed-depth cross-row resultants do not collect `T_n` primes

There are two different ways one might try to go beyond the single row.

### 6.1 Same prime, two folded labels

This is impossible by Lemma 2.1. If a `T_n` prime is attached to `s`, the triangular equation cannot attach the same prime to another `t` with `2t<p`.

Thus a pairwise resultant `Res(R_s,R_t)` for `s!=t` receives no **forced** divisibility from the `T_n` data. Any such divisibility would be an accidental additional Apéry zero, not a consequence of the triangular core.

### 6.2 Adjacent or fixed local stencil around the selected row

The first adjacent row is already a unit by Theorem 5.1. Therefore no fixed local determinant whose target mechanism requires both `R_s` and `R_{s-1}` to vanish can capture a `T_n` prime.

This is exactly why retaining the common physical `n` does not break the old one-primary wall: the physical equation gives one zero in one characteristic, and the canonical contiguous equation is transverse.

---

## 7. A small unconditional pruning: the reflection midpoint is negligible

Although it is not needed for Theorem 5.1, the central folded case can be removed at logarithmic cost.

If `p=2s+1`, then

\[
4\tau(s)+1=(2s+1)^2\equiv0\pmod p.
\]

If also `p|tau(n)-tau(s)`, then

\[
(2n+1)^2=4\tau(n)+1\equiv0\pmod p,
\]

so

\[
\boxed{p\mid2n+1.}
\tag{7.1}
\]

For `p^3>n^2`, every such prime satisfies `p>n^(2/3)`. For all sufficiently large `n`, two distinct such primes would have product greater than `2n+1`; hence there is at most one central prime. Its contribution is at most `log(2n+1)=O(log n)`.

Therefore one may, at the cost `O(log n)`, work only with strict nonmidpoint rows `p>2s+1`. This removes the fixed-point degeneration but does not change the hard core.

---

# II. REFUTED ROUTES

## 8. REFUTED: transport a level-6/Hasse evaluation zero along `tau`

The repository proves the exact Franel pullback

\[
A_p(a)=\sum_{m=0}^{p-1}b_ma^m,
\qquad
K_p(x)=\sum_{m=0}^{p-1}f_mx^m,
\qquad
f_m=\sum_{k=0}^m\binom mk^3,
\]

\[
\boxed{
K_p(x)^2
=(1+x)^{p-1}A_p\!\left(\frac{x(1-8x)}{1+x}\right)
=\sum_{m=0}^{p-1}b_m\Psi_{p,m}(x),
}
\tag{8.1}
\]

with

\[
\Psi_{p,m}(x)=x^m(1-8x)^m(1+x)^{p-1-m}.
\tag{8.2}
\]

The polynomials `Psi_{p,m}` form a triangular basis, so

\[
[K_p^2]_{\Psi_{p,m}}=b_m.
\tag{8.3}
\]

Thus `p|b_s` is a **coordinate zero** of the Hasse square. It is not the assertion that `A_p` or `K_p` vanishes at an evaluation point canonically labeled by `s`.

The triangular congruence

\[
\tau(n)\equiv\tau(s)\pmod p
\]

is an equality of evaluation labels. It cannot transport a coordinate zero without an additional theorem converting that coordinate functional to evaluation at a distinguished point.

### Exact finite falsification

For `p=11`, the Apéry row is

\[
(b_0,\ldots,b_{10})
\equiv
(1,5,7,4,1,0,1,4,7,5,1)\pmod{11}.
\tag{8.4}
\]

Hence `11|b_5`. But

\[
\tau(5)=30\equiv8\pmod{11},
\]

and direct evaluation gives

\[
\boxed{A_{11}(8)\equiv5\not\equiv0\pmod{11}.}
\tag{8.5}
\]

So the nearby stronger claim

```text
b_s == 0 (mod p)  =>  A_p(tau(s)) == 0 (mod p)
```

is false on actual Apéry data.

This does not prove that level-6 geometry is useless. It proves that the simple “transport the modular zero from `s` to `n` because `tau(s)=tau(n)`” argument is ill-typed: the zero lives in coefficient/Mellin space, while `tau` transports evaluation points.

---

## 9. REFUTED: the joint equations automatically give a second `p`-adic digit

A possible escape from Theorem 4.1 would be an automatic strengthening

```text
p | b_s and p | tau(n)-tau(s)
    => p^2 | Rhat_s(tau(n)).
```

This is false even inside the exact `T_n` range.

Take

\[
(n,s,p)=(16,5,11).
\]

Then

\[
H=\lfloor16^{1/3}\rfloor=2,
\qquad
2<5\le\frac{15}{2},
\]

\[
11^3=1331>16^2,
\qquad
2s=10<11\le16,
\]

\[
\tau(16)-\tau(5)=272-30=242=22\cdot11,
\]

and

\[
b_5=819005,
\qquad
v_{11}(b_5)=1.
\]

Thus this is a genuine `T_16` incidence. The exact Racah transport gives

\[
\boxed{
\widehat R_5(\tau(16))\equiv77\pmod{121},
\qquad
v_{11}(\widehat R_5(\tau(16)))=1.
}
\tag{9.1}
\]

So there is no free first-primary lift hidden in the triangular transport.

---

# III. CONDITIONAL / SURVIVING DIRECTION

## 10. What remains after the obstruction

The theorem above removes only a class of candidate proofs. It does **not** bound `T_n` by itself.

A successful theorem must introduce information that is absent from the one-characteristic pair

\[
(p\mid b_s,\ p\mid\tau(n)-\tau(s)).
\]

There are three ways that could happen in principle:

1. a genuinely cross-prime estimate coupling several distinct `T_n` primes at the same physical `n`;
2. a first-primary invariant that is not forced by the mod-`p` Racah row and whose aggregate height is smaller than the target prime mass;
3. a bilinear/large-sieve theorem for the actual defining-characteristic coefficient-zero packets, uniform in the moving characteristic and in the quotient branch.

I do not prove any of these here.

### Comparison with the current top-strip four-prime rectangle gateway

The present theorem is **not** a smaller sufficient theorem than the four-prime rectangle gateway. It is a negative theorem showing why one might have hoped to avoid that gateway but cannot do so with the common-`n` Racah resultant.

The distinction is structural:

- the failed route uses one characteristic and one folded label, optionally a fixed local stencil of adjacent labels;
- Lemma 2.1 says that one characteristic supplies only one folded `T_n` row;
- Theorem 5.1 says the first adjacent row is a unit;
- a four-prime rectangle, by definition, brings in several distinct characteristics and therefore is not subject to same-prime folded-label uniqueness.

So the rectangle gateway remains genuinely beyond this obstruction. If it can be proved, it supplies new cross-prime information. Replacing it by `Rhat_s(tau(n))`, an adjacent Wronskian, or a level-6 evaluation resultant does not.

This is the correct size comparison: **the obstruction is earlier than the rectangle, not stronger than it as a closing theorem.**

---

# IV. Status ledger

## PROVED

1. `tau(n)-tau(s)=(n-s)(n+s+1)` and uniqueness of the folded row `s<p/2`.
2. The exact Racah/Newton interpolation formulas (3.3)–(3.6).
3. The integral scaling `Rhat_s in Z[Y]`.
4. The exact transport/remainder identity (4.1).
5. The exact gcd equality (4.2), hence prime-by-prime equivalence (4.3) on every `T_n` prime.
6. The linear-factor resultant identities (4.4)–(4.6).
7. The monic Racah recurrence and adjacent resultant product (5.4).
8. Adjacent-row nonvanishing at every `T_n` prime.
9. At most logarithmic contribution from the reflection midpoint.
10. The level-6 Franel pullback is coefficient/coordinate data, not evaluation-root data.

## CONDITIONAL / OPEN

1. `T_n=o(n)` remains open in this report.
2. No known large-sieve, modular-parametrization, or resultant theorem is claimed to supply the required pointwise varying-`p` estimate.
3. The four-prime rectangle gateway remains a surviving positive cross-prime target.
4. A genuinely new first-primary or mixed-characteristic carrier could evade this no-go; none is constructed here.

## REFUTED

1. **Refuted:** transporting `b_s=0 mod p` to `Rhat_s(tau(n))=0 mod p` creates a new independent numerator. It is exactly equivalent by (4.2).
2. **Refuted:** the same `T_n` prime can be forced into two folded transported rows. Lemma 2.1 forbids this.
3. **Refuted:** the adjacent Racah row supplies a second target-divisible equation. The resultant product makes it a `p`-unit.
4. **Refuted:** a level-6/Hasse coefficient zero is an evaluation zero at the triangular point. Equation (8.5) is an actual Apéry counterexample.
5. **Refuted:** the joint triangular equations automatically lift the transported value to `p^2`. Equation (9.1) is a genuine `T_16` counterexample.

---

## 11. Verifier

The owned verifier is

```text
problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_verify.py
```

It uses only the Python standard library and exact integer/rational arithmetic. It checks:

- the integral `Rhat_s` construction and diagonal Apéry identity;
- the exact polynomial remainder and gcd identities for a finite regression grid;
- self-duality at triangular nodes;
- uniqueness of the folded triangular root for all tested primes/indices;
- the monic three-term recurrence and the Sylvester determinant formula for adjacent resultants;
- every actual `T_n` event through `n<=100` against transported divisibility and adjacent-row nonvanishing;
- the exact `p=11,s=5` Hasse coefficient/evaluation counterexample;
- the genuine `(n,s,p)=(16,5,11)` `p^2` counterexample;
- finite instances of the Franel pullback identity.

Its final success line is

```text
Q8336_PHYSICAL_RACAH_OBSTRUCTION=PASS
```

The finite checks are falsification/regression only; none is used as a proof of an asymptotic statement.
