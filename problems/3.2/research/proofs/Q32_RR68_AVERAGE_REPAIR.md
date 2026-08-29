# Q32 RR68 average repair

## Verdict

The proposed repair is proved at the staged RR68 level, under the literal
RR68c hypotheses and the recorded physical range `1 <= B <<_v N`.
The four requested estimates are

```text
w_k <= N^(-3+o(1)) z_k m_k,
m_k <= N^(1+o(1))/B,
sum_k z_k <= N^(8/5+o(1)),
W(theta) = sum_k w_k <= N^(-2/5+o(1))/B.
```

The decisive overlap is from a fixed source prime `p` to the physical
integer `k`, not from a zero label to `k`.  For fixed `theta` and `p`,
the congruence fixes `k mod b`; the physical interval has length `O_v(B)`
and `b ~_v B`, so it contains only `O_v(1)` lifts.  The factor `Z(p)`
already counts every actual zero label carried by `p`.  It must not be
inserted a second time when the exact-`k` packets are overlapped.

Together with the staged positive first-moment ledger, the repair gives

```text
sum_theta W(theta)^2 <= N^(2/5+o(1)) M/B,
M = min(B, D/N^(2/5)).
```

The literal host sum then gives

```text
E_mean <= F M K N^(o(1)),
E_cent <= N^(3/5+o(1)) K (B M F)^(1/2).
```

Thus the staged local ridge `sigma < 1/40` is restored: its downstream
optimizer receives exactly the two estimates it was designed to consume.
The staged opposite high-corner minimax remains

```text
min_ell max(5/2 - ell/4, 7/3 + ell)
  = 37/15 at ell = 2/15.
```

This is still above the staged fully unconditional target
`12/5 = 36/15` by `1/15`.  The repair therefore does not close Problem
3.2.

The separate third-moment sufficient premise is not proved.  The actual
pointwise estimate and AF1 give only

```text
sum_{r ~ N} Z(r)^3 <= N^(43/15+o(1)).
```

Nothing below proves `N^(29/15-6 sigma+o(1))` or
`N^(2-6 sigma+o(1))`.

---

## 1. Source and dependency boundary

This is a research proof, not an authoritative insertion.

The connector-visible pins used for the audit were

```text
repository: xiangyazi24/Ramanujan_Challenge
canonical main: c5d932b66ce5e4f1657b587215d290ae7a13018b
chatgpt-drop before this write: b012211d9ff90db90e03726213ec40dc840bb6e0
canonical oracleA_result.tex blob on main:
  f0bfccc441bda22c17658c5586fbb2e6c6431238
```

The task's source audit states that the authoritative AF/RC/RR chain in
the current canonical working source reaches only roughly RR48; RR49--RR68
is staged research.  The configured `chatgpt-drop` branch is itself a
staging/drop branch and does not supply a connected authoritative
`problems/3.2/oracleA_result.tex` surface containing RR68c and its
consumers.  Therefore this change does not edit either
`oracleA_result.tex` or `proof.tex`.  Section 10 supplies a future
insertion block to use only after the staged definitions have been
promoted.

### 1.1 Canonical analytic inputs

Only the following current estimates are used:

```text
(A0)  Z(p) <= p^(2/3+o(1)).
(A1)  sum_{p ~ N} Z(p) <= N^(8/5+o(1)).
(A2)  sum_{p ~ N} |u_p| Z(p) <= N^(3/5+o(1)).
(A3)  sum_{p ~ N} |u_p| Z(p)^2 <= N^(6/5+o(1)).
(A4)  |u_p| <= N^(-1+o(1)).
```

The exact-`k` overlap uses (A1), not the false pointwise assertion
`Z(p) <= N^(3/5)`.

### 1.2 Literal staged RR68 inputs

Fix one unrefined radial datum

```text
theta = (q, eta, b, c0, eps),
C = b eta + eps q c0 = b n + rho,
b ~_v B.
```

For every literal RR68c occurrence, the two source primes `p,p'` share
one genuine integer `k` with

```text
|k| <= C_v B,
gcd(k,b) = gcd(rho,b) = 1,
p k  == -rho (mod b),
p' k == -rho (mod b).
```

The source-prime shell has fixed compact support `p,p' ~_v N`, and the
staged physical range has `1 <= B <<_v N`.  All source, sign, primitive,
actual-zero, chart, radial, and lift masks remain present.  In a positive
upper bound they may only be discarded as `0`--`1` deletions.

For an exact integer `k`, let `P_k = P_k(theta)` be the actual restricted
source-prime set.  Put

```text
m_k = |P_k|,
z_k = sum_{p in P_k} Z(p),
w_k = |u_q| C_k,
W(theta) = sum_k w_k.
```

The literal staged factorization is

```text
(RR68-F)
C_k <= N^(o(1))
      (sum_{p in P_k} |u_p| Z(p))
      (sum_{p' in P_k} |u_p'|).
```

The staged positive first-moment ledger is

```text
(RR68-L1)
sum_theta W(theta) <= N^(4/5+o(1)) M,
M = min(B, D/N^(2/5)).
```

The host consumer is

```text
(RR-host-I)
I_r <= F Z(r)/r sum_theta W(theta)
       + sqrt(F Z(r)(r-Z(r))) B
         (sum_theta W(theta)^2)^(1/2),

(RR-host-E)
E_nonper <= N^(-1+o(1)) K
            sum_{r ~ N} |u_r| Z(r) I_r.
```

---

## 2. Marked-occurrence multiplicities

Three different multiplicities occur.

### 2.1 Actual zero labels on the first source-prime factor

For fixed `p`, all surviving actual zero labels are majorized by `Z(p)`.
This is already the factor in

```text
sum_{p in P_k} |u_p| Z(p).
```

It is also already the factor in `z_k`.  No later overlap may multiply
by another `Z(p)`.

### 2.2 Marks on the second source-prime factor

The second factor in (RR68-F) carries no additional zero-count weight.
For fixed `theta` and exact `k`, the finite sign, chart, primitive-branch,
and physical-lift multiplicity is already absorbed into `N^(o(1))`.
Therefore the second factor is bounded by

```text
N^(o(1)) sum_{p' in P_k} |u_p'|,
```

not by a second copy of `z_k`.

### 2.3 One prime occurring in several exact-k packets

Define

```text
nu_theta(p) = #{k : p in P_k}.
```

Section 5 proves `nu_theta(p) <<_v 1`.  Thus

```text
sum_k z_k = sum_{p ~ N} Z(p) nu_theta(p).
```

This is prime-to-`k` overlap.  It is independent of the number of zero
labels after the single factor `Z(p)` has been included.

---

## 3. Exact-k carrier bound

### Proposition 3.1

For every exact physical `k`,

```text
w_k <= N^(-3+o(1)) z_k m_k.
```

### Proof

By (A4),

```text
sum_{p in P_k} |u_p| Z(p) <= N^(-1+o(1)) z_k,
sum_{p' in P_k} |u_p'|    <= N^(-1+o(1)) m_k.
```

Substitution in (RR68-F) gives

```text
C_k <= N^(-2+o(1)) z_k m_k.
```

The fixed source prime `q` in `theta` contributes its weight once, and
`|u_q| <= N^(-1+o(1))`.  Hence

```text
w_k = |u_q| C_k <= N^(-3+o(1)) z_k m_k.
```

No cancellation is used and every omitted literal mask is a deletion
from a nonnegative factor.  QED.

---

## 4. Size of one exact-k source-prime packet

### Proposition 4.1

Uniformly in physical `k`,

```text
m_k <= N^(1+o(1))/B.
```

Without using the active range `B << N`, the robust statement is

```text
m_k <<_v 1 + N/B.
```

### Proof

Since `gcd(k,b)=1`, the congruence gives one residue class

```text
p == -rho k^(-1) (mod b).
```

The compact source shell has length `O_v(N)`.  One residue class modulo
`b` has at most `1+O_v(N/b)` integer representatives, and the prime
representatives form a subset.  Since `b ~_v B`, this is `O_v(1+N/B)`.
In the staged physical range `B <<_v N`, the `1` is absorbed by `N/B`.
QED.

### Boundary warning

If a future RR version allows `B>N`, it must retain `1+N/B`; a nonempty
packet cannot satisfy a bound strictly smaller than one.  The current
physical range is load-bearing for the shorter displayed form.

---

## 5. Fixed-prime overlap across exact k

### Proposition 5.1

For fixed `theta` and fixed shell prime `p`,

```text
#{k : p in P_k} <<_v 1.
```

Consequently

```text
sum_k z_k <= N^(8/5+o(1)).
```

### Proof

First `gcd(p,b)=1`.  If `d` divides both `p` and `b`, then `d` divides
`p k`; RR68c implies that `d` divides `rho`; and `gcd(rho,b)=1` forces
`d=1`.

The congruence therefore fixes `k` modulo `b`:

```text
k == -rho p^(-1) (mod b).
```

Every physical representative lies in `[-C_v B,C_v B]`, an interval of
length `O_v(B)`.  Since `b ~_v B`, this interval contains only `O_v(1)`
integers in the fixed residue class.  All literal masks can only delete
such representatives.

Now expand the overlap without introducing a new mark:

```text
sum_k z_k
 = sum_k sum_{p in P_k} Z(p)
 = sum_{p ~ N} Z(p) #{k : p in P_k}
 <<_v sum_{p ~ N} Z(p)
 <= N^(8/5+o(1))
```

by (A1).  The zero-label multiplicity is the displayed `Z(p)` and is not
repeated.  QED.

### Corollary 5.2: the common k is essential

Because `gcd(k,b)=1`, the two source primes satisfy

```text
p == p' == -rho k^(-1) (mod b).
```

If the two occurrences had unrelated integers `k,k'`, they would not
belong to one common packet `P_k`, and the factorization could not be
bounded by `z_k m_k` in this form.

---

## 6. Pointwise repaired radial mass

### Theorem 6.1

For every fixed literal `theta`,

```text
W(theta) <= N^(-2/5+o(1))/B.
```

### Proof

All terms are nonnegative.  Sections 3--5 give

```text
W(theta)
 <= N^(-3+o(1)) sum_k z_k m_k
 <= N^(-3+o(1)) (max_k m_k) sum_k z_k
 <= N^(-3+o(1)) (N/B) N^(8/5)
 =  N^(-2/5+o(1))/B.
```

The exact exponent arithmetic is

```text
-3 + 1 + 8/5 = -2/5.
```

The proof uses the first zero moment once and nowhere uses a pointwise
`3/5` estimate.  QED.

---

## 7. Repaired radial square

### Corollary 7.1

Assuming (RR68-L1),

```text
sum_theta W(theta)^2 <= N^(2/5+o(1)) M/B.
```

### Proof

By positivity,

```text
sum_theta W(theta)^2
 <= (max_theta W(theta)) sum_theta W(theta)
 <= (N^(-2/5+o(1))/B) (N^(4/5+o(1)) M)
 =  N^(2/5+o(1)) M/B.
```

The exponent check is `-2/5+4/5=2/5`.  QED.

---

## 8. Literal host summation

Write

```text
S1 = sum_theta W(theta),
S2 = sum_theta W(theta)^2.
```

Then

```text
S1 <= N^(4/5+o(1)) M,
S2 <= N^(2/5+o(1)) M/B.
```

The sum over `r` must be performed before any pointwise zero majorant.

### 8.1 Mean term

From (RR-host-I),

```text
I_r^mean <= F Z(r) S1/r.
```

Therefore

```text
E_mean
 <= N^(-1+o(1)) K F S1
    sum_{r ~ N} |u_r| Z(r)^2/r.
```

Since `r ~ N`, (A3) gives

```text
sum |u_r| Z(r)^2/r <= N^(-1+o(1)) N^(6/5).
```

Using `S1 <= N^(4/5+o(1)) M`,

```text
E_mean <= F M K N^(o(1)).
```

The complete exponent check is

```text
-1 -1 +4/5 +6/5 = 0.
```

The weighted second zero moment is the correct input.

### 8.2 Centered term

The centered part is

```text
I_r^cent
 <= sqrt(F Z(r)(r-Z(r))) B S2^(1/2).
```

Thus

```text
E_cent
 <= N^(-1+o(1)) K sqrt(F) B S2^(1/2)
    sum_{r ~ N} |u_r| Z(r)^(3/2) sqrt(r-Z(r)).
```

The tautological bound is `sqrt(r-Z(r)) <= N^(1/2+o(1))`.  Weighted
Cauchy--Schwarz gives

```text
sum |u_r| Z(r)^(3/2)
 <= (sum |u_r| Z(r))^(1/2)
    (sum |u_r| Z(r)^2)^(1/2)
 <= N^((3/5+6/5)/2+o(1))
 =  N^(9/10+o(1)).
```

Also

```text
S2^(1/2) <= N^(1/5+o(1)) (M/B)^(1/2).
```

Hence

```text
E_cent
 <= N^(-1+1/2+9/10+1/5+o(1))
    K sqrt(F) B (M/B)^(1/2)
 =  N^(3/5+o(1)) K (B M F)^(1/2).
```

The exponent check is

```text
-1 + 1/2 + 9/10 + 1/5 = 3/5.
```

No pointwise `Z(r) <= N^(3/5)` is used.  The only pointwise step is
`r-Z(r) <= r`.

---

## 9. Downstream exponent status

### 9.1 Local ridge

The average repair itself is sigma-free.  Its complete downstream output
is exactly the intended pair

```text
E_mean <= F M K,
E_cent <= N^(3/5+o(1)) K (B M F)^(1/2).
```

Therefore the staged local optimizer receives no worsened exponent and
its strict range `sigma<1/40` is restored.  This is a staged dependency
statement, not an authoritative theorem insertion, because the sigma
consumer belongs to RR49--RR68 research.

### 9.2 Opposite high corner

The staged opposite-corner exponent lines are

```text
H1(ell) = 5/2 - ell/4,
H2(ell) = 7/3 + ell.
```

Their difference is `1/6-5 ell/4`, so they cross at `ell=2/15`.
Both then equal `37/15`.  Hence

```text
min_ell max(H1(ell),H2(ell)) = 37/15.
```

The staged unconditional target is `12/5=36/15`; the excess is `1/15`.
The repair restores the local ridge but leaves this corner open.

### 9.3 Third moment remains absent

Using the actual pointwise estimate and AF1,

```text
sum Z(r)^3
 <= (max Z(r)) sum Z(r)^2
 <= N^(2/3+o(1)) N^(11/5+o(1))
 =  N^(43/15+o(1)).
```

The identity is `2/3+11/5=43/15`.  This does not prove either staged
sufficient premise

```text
sum Z(r)^3 <= N^(29/15-6 sigma+o(1)),
```

or

```text
sum Z(r)^3 <= N^(2-6 sigma+o(1)).
```

### 9.4 Smallest remaining residual

At the staged endpoint, the smallest explicitly quantified residual is
the opposite high-corner excess `N^(1/15+o(1))` over the fully
unconditional target.  A genuinely stronger third zero moment is one
sufficient route, but it is not implied by AF1 and the pointwise `2/3`
bound.

```text
RR68 common-k average repair:       PROVED
staged local sigma<1/40 ridge:      RESTORED
staged high-corner 37/15:           UNCHANGED
fully unconditional target 12/5:   NOT REACHED
third-moment sufficient premise:    UNPROVED
full Problem 3.2:                   OPEN
```

---

## 10. Future authoritative insertion block

Do not insert this block until RR49--RR68, the exact-`k` packet `P_k`,
`C_k,w_k,W(theta)`, and the host consumer have been promoted into the
authoritative source.

```latex
\begin{lemma}[Average repair for the RR68 common radial integer]
\label{lem:oracleA-rr68-average-repair}
Fix one literal unrefined radial datum
\[
 \theta=(q,\eta,b,c_0,\varepsilon),\qquad
 b\asymp_v B,\qquad 1\le B\ll_vN,
\]
and write
\[
 b\eta+\varepsilon q c_0=bn+\rho.
\]
Assume that every literal RR68c occurrence belonging to an exact
integer $k$ satisfies
\[
 |k|\ll_vB,\qquad (k,b)=(\rho,b)=1,
\]
and, for its two source primes,
\[
 pk\equiv p'k\equiv-\rho\pmod b.
\]
Let $\mathcal P_k$ be the actual restricted source-prime set, and put
\[
 m_k=|\mathcal P_k|,\qquad
 z_k=\sum_{p\in\mathcal P_k}Z(p),\qquad
 w_k=|u_q|C_k.
\]
Suppose the literal positive carrier factorization is
\[
 C_k\ll N^{o(1)}
 \left(\sum_{p\in\mathcal P_k}|u_p|Z(p)\right)
 \left(\sum_{p'\in\mathcal P_k}|u_{p'}|\right).
\]
Then
\[
 w_k\ll N^{-3+o(1)}z_km_k,\qquad
 m_k\ll_v {N^{1+o(1)}\over B},
\]
\[
 \sum_kz_k\ll N^{8/5+o(1)},
 \qquad
 W(\theta):=\sum_kw_k
 \ll {N^{-2/5+o(1)}\over B}.
\]
All literal source, sign, primitive, actual-zero, chart and lift masks
may be retained.
\end{lemma}

\begin{proof}
The source weights give
\[
 C_k\ll N^{-2+o(1)}z_km_k,
 \qquad w_k\ll N^{-3+o(1)}z_km_k.
\]
For fixed $k$, the congruence and $(k,b)=1$ put every source prime in
one residue class modulo $b$, whence
$m_k\ll1+N/b\ll N/B$.  Conversely, for fixed $p$, the congruence and
$(\rho,b)=1$ first imply $(p,b)=1$ and then fix
$k\equiv-\rho\bar p\pmod b$.  The physical interval $|k|\ll_vB$ has
length $O_v(B)$ and $b\asymp_vB$, so it contains $O_v(1)$ such lifts.
Therefore
\[
 \sum_kz_k
 \ll_v\sum_{p\asymp N}Z(p)
 \ll N^{8/5+o(1)}.
\]
Combining these estimates proves the result.  The factor $Z(p)$ already
counts all zero labels for $p$ and is not repeated in the $k$-overlap.
\end{proof}

\begin{corollary}[Repaired RR68 radial square]
\label{cor:oracleA-rr68-average-square}
If
\[
 \sum_\theta W(\theta)\ll N^{4/5+o(1)}M,
 \qquad M=\min\!\left(B,{D\over N^{2/5}}\right),
\]
then
\[
 \sum_\theta W(\theta)^2
 \ll {N^{2/5+o(1)}M\over B}.
\]
Consequently the literal RR host inequalities give
\[
 E_{\rm mean}\ll FMK,
 \qquad
 E_{\rm cent}\ll N^{3/5+o(1)}K(BMF)^{1/2}.
\]
\end{corollary}
```

---

## 11. Computation and verification scope

The companion script

```text
problems/3.2/research/scripts/q32_rr68_average_repair_verify.py
```

checks exact rational exponents and exhaustively tests the elementary
residue/lift statements for small moduli, physical intervals, and prime
shells.  That computation is not the asymptotic proof.  The asymptotic
proof is Sections 3--8 above.

The intended commands are

```bash
git diff --check
python3 problems/3.2/research/scripts/q32_rr68_average_repair_verify.py
cd problems/3.2
latexmk -pdf -interaction=nonstopmode -halt-on-error proof.tex
```

Neither TeX file is edited or newly included, so the TeX source graph is
unchanged.  A successful `latexmk` run checks the pre-existing canonical
document, not the staged RR68 insertion block.
