# Q32 RR68 average repair -- corrected source audit

## Verdict

The **common-`k` average repair survives**.  For one literal unrefined
radial datum `theta`, the exact-`k` packets still satisfy

```text
w_k <= N^(-3+o(1)) z_k m_k,
m_k <= N^(1+o(1))/B,
sum_k z_k <= N^(8/5+o(1)),
W(theta) = sum_k w_k <= N^(-2/5+o(1))/B.
```

Consequently, with the staged positive first-carrier ledger,

```text
sum_theta W(theta)^2 <= N^(2/5+o(1)) M/B,
M = min(B, D/N^(2/5)).
```

The repair does not use a false pointwise estimate `Z(p) <= N^(3/5)`.
Its load-bearing overlap is the bounded number of physical integer lifts
of one residue class for `k`; the zero-label multiplicity is already the
single factor `Z(p)` inside `z_k`.

Two downstream bookkeeping statements in the previous version were
wrong and are corrected here:

```text
RR68av high-corner competitors:
  H_-(ell) = 13/5 - ell,
  H_+(ell) = 12/5 + ell/2.

min_ell max(H_-(ell), H_+(ell))
  = 37/15 at ell = 2/15.

fully-unconditional staged target:
  T(sigma) = 11/5 - 2 sigma.

high-corner residual:
  37/15 - T(sigma) = 4/15 + 2 sigma.
```

Thus the local `sigma<1/40` RR68av ridge is restored, but the opposite
high corner remains well above the fully-unconditional target.  Full
Problem 3.2 remains open.

---

## 0. Source and dependency boundary

This file is a research proof.  It is not an authoritative TeX
insertion.

The source pins audited for the repository artifact are

```text
repository: xiangyazi24/Ramanujan_Challenge
canonical main: c5d932b66ce5e4f1657b587215d290ae7a13018b
canonical oracleA_result.tex blob:
  f0bfccc441bda22c17658c5586fbb2e6c6431238
artifact base: a9f8e0cdaae1ef7c94d77b50d549c6b846e347db
independent audit parent: 2ef90cbc892edb911977d94deca3c6e83f2677c6
```

The connector-readable canonical `oracleA_result.tex` does not contain
the staged RR49--RR68 chain.  In particular, it does not license
promoting an RR68av statement into authoritative TeX, and it does not
license the previous shorthand

```text
sum_{all shell primes r} |u_r| Z(r)^2 <= N^(6/5+o(1))
```

as a global theorem.

The host-local staging path

```text
/Users/huangx/.openclaw/workspace/staged_patches/p32_q5349_rr.tex
```

is outside the GitHub connector surface, and the task forbids a local
shell/sandbox route.  Accordingly, this correction does not pretend to
have read unavailable bytes: it uses the literal RR68av formulas and
moment scope supplied in Q5567, checks their algebra against the
connector-readable canonical source and the two committed artifacts,
and keeps every staged conclusion inside this research file.

Neither `problems/3.2/oracleA_result.tex` nor
`problems/3.2/proof.tex` is edited.

### 0.1 Inputs used by the common-k repair

The common-`k` part uses only

```text
(A0) |u_p| <= N^(-1+o(1)) on the active source shell;
(A1) sum_{p ~ N} Z(p) <= N^(8/5+o(1));
(A2) the literal RR68c common-k congruences and physical range;
(A3) the literal positive factorization of C_k;
(A4) the staged positive L1 ledger for W(theta).
```

It does **not** use a weighted second zero moment.

### 0.2 Exact scope of the weighted second-zero input

Let `R_cen` denote the exact masked central host-`r` layer consumed by
RR68av.  The safe staged input is the factorial moment on that layer:

```text
(CF2)
sum_{r in R_cen} |u_r| Z(r)(Z(r)-1)
  <= N^(6/5+o(1)).
```

The diagonal is supplied by the weighted first moment, restricted by
deletion to the same layer:

```text
(CF1)
sum_{r in R_cen} |u_r| Z(r)
  <= N^(3/5+o(1)).
```

Since, identically for every nonnegative integer `Z`,

```text
Z^2 = Z(Z-1) + Z,
```

one obtains the exact RR68av consumer bound

```text
(CZ2)
sum_{r in R_cen} |u_r| Z(r)^2
 <= N^(6/5+o(1)) + N^(3/5+o(1))
 <= N^(6/5+o(1)).
```

This is a **central-layer consequence**, not a new global theorem about
all shell primes.  Every use below is explicitly over `R_cen`.

---

## 1. Literal RR68c data and multiplicities

Fix one literal unrefined radial datum

```text
theta = (q, eta, b, c0, eps),
C = b eta + eps q c0 = b n + rho,
b ~_v B,
1 <= B <<_v N.
```

For every literal RR68c occurrence, the two source primes `p,p'` share
one genuine integer `k` satisfying

```text
|k| <<_v B,
gcd(k,b) = gcd(rho,b) = 1,
p k  == -rho (mod b),
p' k == -rho (mod b).
```

For this exact integer `k`, let `P_k=P_k(theta)` be the actual restricted
**set** of source primes.  Define

```text
m_k = |P_k|,
z_k = sum_{p in P_k} Z(p),
w_k = |u_q| C_k,
W(theta) = sum_k w_k.
```

The literal positive factorization is

```text
(RR68-F)
C_k <= N^(o(1))
      (sum_{p in P_k} |u_p| Z(p))
      (sum_{p' in P_k} |u_p'|).
```

The marked-occurrence convention is therefore fixed:

1. The first source-prime factor pays for all surviving actual zero
   labels by the single factor `Z(p)`.
2. The second source-prime factor carries no second zero label.
3. Source, sign, primitive, actual-zero, chart, radial, and physical-lift
   masks are `0`--`1` deletions, or bounded literal chart multiplicities
   already absorbed in `N^(o(1))`.
4. `P_k` is a support set.  Replacing it by an occurrence multiset would
   require a separate multiplicity theorem and is not done here.

---

## 2. Exact-k carrier bound

By `|u_p| <= N^(-1+o(1))` on the source shell,

```text
sum_{p in P_k} |u_p| Z(p)
  <= N^(-1+o(1)) z_k,

sum_{p' in P_k} |u_p'|
  <= N^(-1+o(1)) m_k.
```

Substitution into (RR68-F) gives

```text
C_k <= N^(-2+o(1)) z_k m_k.
```

The fixed source prime `q` in `theta` contributes its weight once:

```text
w_k = |u_q| C_k
    <= N^(-3+o(1)) z_k m_k.
```

No cancellation is used.

---

## 3. Fixed k: size of the source-prime packet

Because `gcd(k,b)=1`, the congruence fixes one class for the source
prime:

```text
p == -rho k^(-1) (mod b).
```

The source shell has length `O_v(N)`.  One class modulo `b` has at most
`O_v(1+N/b)` integer representatives, and the prime representatives are
a subset.  Since `b ~_v B`,

```text
m_k <<_v 1 + N/B.
```

In the active physical range `B <<_v N`, this becomes

```text
m_k <= N^(1+o(1))/B.
```

The robust `1+N/B` form must be retained if a future version allows
`B>N`.

---

## 4. Fixed p: bounded physical k lifts

First the congruence forces `gcd(p,b)=1`.  Indeed, a common divisor of
`p` and `b` divides `pk`, hence divides `rho`; `gcd(rho,b)=1` then makes
that divisor one.

Thus fixed `theta` and fixed source prime `p` determine one class

```text
k == -rho p^(-1) (mod b).
```

The physical interval `|k| <<_v B` has length `O_v(B)`, while
`b ~_v B`.  It contains only `O_v(1)` representatives of that class.
Consequently

```text
nu_theta(p) := #{k : p in P_k} <<_v 1.
```

Expanding the marked overlap,

```text
sum_k z_k
 = sum_k sum_{p in P_k} Z(p)
 = sum_{p ~ N} Z(p) nu_theta(p)
 <<_v sum_{p ~ N} Z(p)
 <= N^(8/5+o(1)).
```

There is no extra zero-label multiplier: `Z(p)` is already the weight in
`z_k`.  This is the load-bearing average repair.

The fact that the two source primes share the same integer `k` is also
load-bearing.  It gives

```text
p == p' == -rho k^(-1) (mod b).
```

Unrelated integers `k,k'` would not define the single packet appearing
in (RR68-F).

---

## 5. Fixed-theta mass and the radial square

By positivity and Sections 2--4,

```text
W(theta)
 <= N^(-3+o(1)) sum_k z_k m_k
 <= N^(-3+o(1)) (max_k m_k) sum_k z_k
 <= N^(-3+o(1)) (N/B) N^(8/5)
 =  N^(-2/5+o(1))/B.
```

The exact exponent identity is

```text
-3 + 1 + 8/5 = -2/5.
```

Now use the staged positive first-carrier ledger

```text
(RR68-L1)
sum_theta W(theta) <= N^(4/5+o(1)) M,
M = min(B, D/N^(2/5)).
```

Then

```text
sum_theta W(theta)^2
 <= (max_theta W(theta)) sum_theta W(theta)
 <= (N^(-2/5+o(1))/B) (N^(4/5+o(1)) M)
 =  N^(2/5+o(1)) M/B.
```

The exponent identity is `-2/5+4/5=2/5`.

---

## 6. Literal central-host summation

Put

```text
S1 = sum_theta W(theta),
S2 = sum_theta W(theta)^2.
```

The staged host inequalities, with their actual central mask retained,
are

```text
I_r <= F Z(r)/r S1
       + sqrt(F Z(r)(r-Z(r))) B S2^(1/2),

E_nonper <= N^(-1+o(1)) K
            sum_{r in R_cen} |u_r| Z(r) I_r.
```

### 6.1 Mean term

Summing `r` first,

```text
E_mean
 <= N^(-1+o(1)) K F S1
    sum_{r in R_cen} |u_r| Z(r)^2/r.
```

Since `r ~ N` on the central layer, (CZ2) gives

```text
sum_{r in R_cen} |u_r| Z(r)^2/r
 <= N^(-1+o(1)) N^(6/5).
```

Using `S1 <= N^(4/5+o(1)) M`,

```text
E_mean <= F M K N^(o(1)).
```

The exponent ledger is

```text
-1 -1 + 4/5 + 6/5 = 0.
```

### 6.2 Centered term

The centered contribution is

```text
E_cent
 <= N^(-1+o(1)) K sqrt(F) B S2^(1/2)
    sum_{r in R_cen}
      |u_r| Z(r)^(3/2) sqrt(r-Z(r)).
```

Use `sqrt(r-Z(r)) <= N^(1/2+o(1))`.  Weighted
Cauchy--Schwarz on the same central layer gives

```text
sum_{r in R_cen} |u_r| Z(r)^(3/2)
 <= (sum_{r in R_cen} |u_r| Z(r))^(1/2)
    (sum_{r in R_cen} |u_r| Z(r)^2)^(1/2)
 <= N^((3/5+6/5)/2+o(1))
 =  N^(9/10+o(1)).
```

Also

```text
S2^(1/2) <= N^(1/5+o(1)) (M/B)^(1/2).
```

Therefore

```text
E_cent
 <= N^(-1+1/2+9/10+1/5+o(1))
    K sqrt(F) B (M/B)^(1/2)
 =  N^(3/5+o(1)) K (B M F)^(1/2).
```

The exponent identity is

```text
-1 + 1/2 + 9/10 + 1/5 = 3/5.
```

The only second-zero estimate used in either host term is (CZ2), whose
scope is exactly `R_cen` and whose proof is factorial plus diagonal.

---

## 7. Corrected RR68av high-corner ledger

The literal competitors are

```text
H_-(ell) = 13/5 - ell,
H_+(ell) = 12/5 + ell/2.
```

They cross when

```text
13/5 - ell = 12/5 + ell/2,
1/5 = 3 ell/2,
ell = 2/15.
```

At this point,

```text
H_-(2/15) = 13/5 - 2/15 = 37/15,
H_+(2/15) = 12/5 + 1/15 = 37/15.
```

The first function decreases and the second increases.  Hence, on the
staged admissible interval containing `2/15`,

```text
min_ell max(H_-(ell), H_+(ell)) = 37/15.
```

The relevant fully-unconditional target in this ledger is

```text
T(sigma) = 11/5 - 2 sigma = 33/15 - 2 sigma.
```

Thus the exact residual is

```text
37/15 - T(sigma)
 = 37/15 - (33/15 - 2 sigma)
 = 4/15 + 2 sigma.
```

It is **not** `1/15`.  For the local range `0<=sigma<1/40`, the residual
is at least `4/15` and approaches `19/60` at the upper endpoint.

The average repair is sigma-free: it restores the host inputs that the
local optimizer consumes.  It does not remove this opposite-corner gap.

---

## 8. Third moment and final status

The separate cubic sufficient premise remains unproved.  With the actual
pointwise estimate and AF1,

```text
sum_{r ~ N} Z(r)^3
 <= (max_{r ~ N} Z(r)) sum_{r ~ N} Z(r)^2
 <= N^(2/3+o(1)) N^(11/5+o(1))
 =  N^(43/15+o(1)).
```

The exponent identity is `2/3+11/5=43/15`.  This does not prove either
stronger staged premise

```text
sum Z(r)^3 <= N^(29/15-6 sigma+o(1))
```

or

```text
sum Z(r)^3 <= N^(2-6 sigma+o(1)).
```

Final classification:

```text
(i) common-k average repair:                    PROVED
(ii) fixed-theta W(theta) bound:                PROVED
(iii) central host E_mean/E_cent ledger:        RESTORED
      using central factorial + diagonal only
(iv) local sigma<1/40 ridge:                    RESTORED
(v) literal high-corner minimax:                37/15
(vi) target:                                    11/5-2 sigma
(vii) exact high-corner residual:               4/15+2 sigma
(viii) stronger cubic moment premise:           UNPROVED
(ix) full Problem 3.2:                          OPEN
```

The smallest explicit residual in the corrected staged high-corner
ledger is therefore `4/15+2 sigma`, not `1/15`.

---

## 9. Future insertion block -- do not insert yet

The following block is only a template for a future authoritative source
in which RR49--RR68, `R_cen`, `P_k`, `C_k`, `w_k`, and `W(theta)` have
already been defined.

```latex
\begin{lemma}[RR68 common-integer average repair]
\label{lem:oracleA-rr68-average-repair}
Fix one literal RR68 datum
\[
 \theta=(q,\eta,b,c_0,\varepsilon),\qquad
 b\asymp_v B,\qquad 1\le B\ll_vN,
\]
and write
\[
 b\eta+\varepsilon q c_0=bn+\rho.
\]
Suppose every exact occurrence has one common integer $k$ satisfying
\[
 |k|\ll_vB,\quad (k,b)=(\rho,b)=1,
 \quad pk\equiv p'k\equiv-\rho\pmod b,
\]
and suppose the literal factorization is
\[
 C_k\ll N^{o(1)}
 \left(\sum_{p\in P_k}|u_p|Z(p)\right)
 \left(\sum_{p'\in P_k}|u_{p'}|\right).
\]
With $m_k=|P_k|$, $z_k=\sum_{p\in P_k}Z(p)$, and
$w_k=|u_q|C_k$, one has
\[
 w_k\ll N^{-3+o(1)}z_km_k,\qquad
 m_k\ll_vN^{1+o(1)}/B,\qquad
 \sum_kz_k\ll N^{8/5+o(1)},
\]
and hence
\[
 W(\theta)=\sum_kw_k\ll_v N^{-2/5+o(1)}/B.
\]
\end{lemma}

\begin{remark}[Scope of the RR68 host second moment]
\label{rem:oracleA-rr68-central-z2}
The host layer uses only
\[
 \sum_{r\in R_{\rm cen}}|u_r|Z(r)(Z(r)-1)
 \ll N^{6/5+o(1)}
\]
together with
\[
 \sum_{r\in R_{\rm cen}}|u_r|Z(r)
 \ll N^{3/5+o(1)}.
\]
Thus $Z(r)^2=Z(r)(Z(r)-1)+Z(r)$ supplies the required
central-layer square moment.  No global all-shell-prime assertion is
being made.
\end{remark}
```

No part of this block is inserted by Q5567.

---

## 10. Role of the verifier

`research/scripts/q32_rr68_average_repair_verify.py` checks exact rational
exponent arithmetic, the factorial-plus-diagonal identity, the corrected
RR68av high-corner and residual arithmetic, and exhaustive small
congruence/lift incidence models.

The computation is not the asymptotic proof.  The asymptotic content is
the interval-counting and moment argument written above.
