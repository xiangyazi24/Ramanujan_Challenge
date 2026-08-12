# Q7696 — finite-characteristic multiplicity of transverse Apéry common rows

## Status

I do **not** have a proof of a uniform `O(1)` bound, and I found no counterexample to the sharp conjecture

\[
\boxed{|Z_{\rm common}(p)|\le 1.}
\]

The new exact Sage scan through every prime `p <= 20000` finds exactly the same two common rows as before,

\[
\boxed{(p,r)=(17,13),\quad (2237,492),}
\]

and no prime has multiplicity two.  This is stronger finite evidence, not a proof.

The main structural conclusion is sharper than “we still need a second equation.”  There is only **one** marked/Green shooting equation after an ordinary Hasse return.  The continuant marked return and the Green-period condition are two coordinate expressions for that same scalar equation.  The actual missing theorem is a zero-normalized level-six Eichler endpoint-avoidance statement for the canonical source and canonical boundary line.

The verifier is

```text
problems/3.2/research/scripts/q7696_common_row_multiplicity.sage
```

on audit branch

```text
q7696-common-row-6cb4d9fd
```

and was run exactly by GitHub Actions run `31570501901` with `Q7696_PMAX=20000`.

---

## 1. Normalization and the exact Green/Casoratian identity

Write

\[
P(n)=34n^3+51n^2+27n+5=(2n+1)(17n^2+17n+5).
\]

For `n<p` the homogeneous Apéry recurrence is

\[
n^3b_n-P(n-1)b_{n-1}+(n-1)^3b_{n-2}=0,
\]

with `b_0=1,b_1=5`.

For the actual level-six source

\[
F(t)=\sum_{n\ge0}b_nt^n,
\qquad
\Delta(t)=1-34t+t^2,
\]

\[
G(t)=\frac1{F(t)^2\sqrt{\Delta(t)}}=\sum_{n\ge0}g_nt^n,
\]

take the canonically normalized inhomogeneous companion

\[
\kappa_0=0,\qquad \kappa_1=-36,
\]

\[
n^3\kappa_n-P(n-1)\kappa_{n-1}+(n-1)^3\kappa_{n-2}=-5g_n.
\]

Define

\[
W_n=n^3(b_n\kappa_{n-1}-b_{n-1}\kappa_n).
\]

Multiply the inhomogeneous recurrence by `b_{n-1}`, the homogeneous recurrence by `\kappa_{n-1}`, and subtract.  The middle terms cancel and give

\[
W_n-W_{n-1}=5g_nb_{n-1}.
\]

At `n=1`, `W_1=36`, while `g_1=7`, so

\[
-W_n=-1-5\sum_{m=1}^n g_mb_{m-1}.
\]

Therefore the exact bridge is

\[
\boxed{
\Xi_n
=n^3(b_{n-1}\kappa_n-b_n\kappa_{n-1}).
}
\tag{1.1}
\]

At a Hasse zero `b_n=0`, consecutive homogeneous zeros are impossible for `0<n<p`, since all recurrence leading/trailing factors are units.  Hence `b_{n-1}` is a unit and

\[
\boxed{b_n=0\quad\Longrightarrow\quad
\Xi_n=0\iff \kappa_n=0.}
\tag{1.2}
\]

This is checked for every relevant row in the new verifier.

---

## 2. Correct Duhamel shooting law: only one endpoint equation

Put

\[
v_n=\binom{\kappa_n}{\kappa_{n-1}}.
\]

Then

\[
v_{n+1}=T_n v_n+u_n,
\]

where

\[
T_n=
\begin{pmatrix}
P(n)/(n+1)^3 & -n^3/(n+1)^3\\
1&0
\end{pmatrix},
\qquad
u_n=\binom{-5g_{n+1}/(n+1)^3}{0}.
\]

For `r<s<p`, with

\[
\Phi(s,r)=T_{s-1}\cdots T_r,
\]

Duhamel gives

\[
\boxed{
v_s=\Phi(s,r)v_r+
\sum_{j=r}^{s-1}\Phi(s,j+1)u_j.}
\tag{2.1}
\]

Now suppose `b_r=b_s=0`.  The homogeneous return says the appropriate off-diagonal transfer coefficient vanishes; in the cleared continuant language of the marked-return analysis this is

\[
N_h(r)=0,\qquad h=s-r,
\]

while backward recurrence gives

\[
\boxed{N_{h-1}(r)\ne0\pmod p.}
\tag{2.2}
\]

Thus after the ordinary return, the terminal marked condition is a single affine condition with unit coefficient on the final forcing value `g_s`.

The Green form is obtained immediately from (1.1): for Hasse-zero endpoints,

\[
\boxed{
s^3b_{s-1}\kappa_s-r^3b_{r-1}\kappa_r
=-5\sum_{m=r+1}^{s}g_mb_{m-1}.}
\tag{2.3}
\]

If the first endpoint is common, then `\kappa_r=0`; since `s^3b_{s-1}` is a unit,

\[
\boxed{
\kappa_s=0
\iff
\sum_{m=r+1}^{s}g_mb_{m-1}=0.}
\tag{2.4}
\]

Equivalently,

\[
\boxed{
g_s=-b_{s-1}^{-1}
\sum_{m=r+1}^{s-1}g_mb_{m-1}.}
\tag{2.5}
\]

The cleared marked-continuant equation is

\[
M_h(r)=
\sum_{j=1}^{h-1}
 g_{r+j+1}N_j(r)
 \prod_{q=j+1}^{h-1}(r+q)^3=0,
\]

whose coefficient of `g_s` is exactly `N_{h-1}(r)`, a unit.  Equations (2.5) and `M_h(r)=0` are the **same shooting equation**, after multiplying by the unit relating `N_{h-1}(r)` to the homogeneous endpoint coefficient.  They must not be counted as two independent codimension-one constraints.

The verifier checks (2.3) and the equivalence (2.4) for **every pair of Hasse-zero rows** for every prime in the scan.

---

## 3. Why a free-source recurrence cannot prove multiplicity one

For a formal source `s_n`, keep the canonical initial line and solve

\[
n^3x_n-P(n-1)x_{n-1}+(n-1)^3x_{n-2}=-5s_n.
\]

Changing only `s_r` by `\delta` changes the same-row value by

\[
\Delta x_r=-\frac5{r^3}\delta.
\]

For `p>=7` this coefficient is a unit.  Perturbations at later rows do not change earlier values.  Hence values at any prescribed increasing list of Hasse-zero rows are triangularly controllable by source perturbations supported at those rows.

So no proof based only on

```text
- the second-order transfer recurrence,
- the homogeneous return N_h=0,
- unit denominators,
- the formal Green identity,
- projective dynamics,
- or a generic inhomogeneous recurrence
```

can prove `|Z_common(p)|<=1`.

This is not hypothetical: the Q7684 audit produced the exact `p=19`, rows `8,10` source perturbation forcing both marked values to zero while leaving the homogeneous Apéry orbit unchanged.

---

## 4. The actual source and boundary normalization are both essential

The source is not arbitrary.  On the level-six parametrization,

\[
t(\tau)=
\left(\frac{\eta(\tau)\eta(6\tau)}
{\eta(2\tau)\eta(3\tau)}\right)^{12},
\]

\[
E(\tau)=
\frac{\eta(2\tau)^7\eta(3\tau)^7}
{\eta(\tau)^5\eta(6\tau)^5}
=F(t(\tau)).
\]

For `K(t)=\sum\kappa_nt^n`, the canonical extension is the normalized triple Eichler antiderivative of

\[
\mathcal M_4=
\frac{-3E_4(\tau)+4E_4(2\tau)-9E_4(3\tau)+108E_4(6\tau)}{20},
\]

in the form

\[
\frac{K(t(q))}{E(q)}=\mathcal E_{\mathcal M_4}(q).
\]

There is also an exact first-block characteristic-`p` description.  Let

\[
A_p(t)=\sum_{j=0}^{p-1}b_jt^j.
\]

Then

\[
F(t)=A_p(t)F(t^p)
\]

and

\[
\boxed{
G(t)=R_p(t)G(t^p),\qquad
R_p(t)=\frac{\Delta(t)^{(p-1)/2}}{A_p(t)^2}.}
\tag{4.1}
\]

Thus for `0<=m<p`,

\[
\boxed{g_m=[t^m]R_p(t).}
\tag{4.2}
\]

The new scan computes precisely this canonical `G=1/(F^2\sqrt\Delta)` coefficient sequence and the canonical boundary line `\kappa_0=0,\kappa_1=-36`; it never substitutes an arbitrary forcing sequence.

The boundary line is genuinely additional arithmetic data.  If `a_n` is the second homogeneous solution with `a_0=0,a_1=6`, then

\[
a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}.
\]

Replacing `\kappa` by `\kappa^{(x)}=\kappa+xa` keeps **the same source** but changes

\[
\Xi_n^{(x)}=\Xi_n+6x.
\]

At `p=41`, the actual source has

\[
b_{10}=b_{30}=0,\qquad \Xi_{10}=\Xi_{30}=7.
\]

Taking `x=-7/6` would move both to zero.  Therefore neither the recurrence alone nor even the exact source without its canonical Eichler normalization can prove zero-fiber simplicity.

---

## 5. Hasse reciprocity: exact, useful, but not a uniqueness theorem

For the Apéry integers one has the classical finite sum

\[
b_n=\sum_{k=0}^n
\binom nk^2\binom{n+k}{k}^2.
\]

For `0<=n<p`, reduction modulo `p` gives

\[
\binom{p-1-n}{k}
\equiv(-1)^k\binom{n+k}{k},
\]

and

\[
\binom{p-1-n+k}{k}
\equiv(-1)^k\binom nk,
\]

with the terms `k>n` vanishing modulo `p`.  Hence

\[
\boxed{b_{p-1-n}\equiv b_n\pmod p.}
\tag{5.1}
\]

The same reflection is compatible with the recurrence because

\[
P(-n-1)=-P(n).
\]

The verifier asserts (5.1) for every row and every prime `7<=p<=20000`.

Consequences:

- Hasse-zero rows occur in reflected pairs `r <-> p-1-r`, except possibly the midpoint.
- `b_{p-1}=b_0=1`, so a common row never occurs at `p-1`.
- Hasse reciprocity does **not** imply any corresponding injectivity of `\Xi`.

The guard prime is decisive:

\[
p=41,\quad 30=40-10,
\]

and the reflected Hasse pair has the same actual Eichler value `7`.  Thus any argument asserting `\Xi_r\ne\Xi_{p-1-r}` on Hasse zeros is false.

More generally, Hasse reciprocity is a symmetry in the coefficient/Mellin index.  It does not turn those indices into points of a bounded-degree divisor on the geometric `t`-line.

---

## 6. What the second-order / Green-kernel viewpoint does and does not give

The Apéry Picard-Fuchs equation is the familiar third-order K3 equation, with a rank-two/symmetric-square origin in the modular description.  The reciprocal source

\[
G=\frac1{F^2\sqrt\Delta}
\]

has exactly the reduction-of-order/Green-kernel shape expected from that rank-two structure.  This explains why, after an ordinary homogeneous return, only one affine marked coordinate remains: the Duhamel defect is one-dimensional.

But this does **not** give a finite-field Sturm theorem in the row index.

- Real sign/interlacing arguments have no order structure modulo `p`.
- The row `r` is a coefficient index, not a geometric parameter value.
- The free-source controllability countermodel preserves the homogeneous transfer operator.
- The `p=41` actual-source collision shows that even the canonical Green coordinate can repeat on distinct Hasse-zero rows.

So the second-order viewpoint explains the shooting rank correctly, but no present second-order theorem rules out the special value zero twice.

---

## 7. Marked returns and resultants

For a fixed gap `h=s-r`, the ordinary return is the explicit continuant equation

\[
N_h(r)=0.
\]

The marked return is linear in the final source coefficient:

\[
M_h(r)=C_{r,h}+N_{h-1}(r)g_s,
\]

and `N_{h-1}(r)` is a unit at an ordinary return.  Hence the marked condition cuts the actual forcing orbit transversely in the terminal source coordinate.

This is the correct use of a resultant: for fixed/short `h`, eliminate the ordinary return and the **single** marked equation.  It is not legitimate to use the Green-period equation as a third independent condition.

For the actual source, (4.2) rewrites the endpoint condition as

\[
[t^s]R_p(t)=
-b_{s-1}^{-1}\sum_{m=r+1}^{s-1}
 b_{m-1}[t^m]R_p(t).
\tag{7.1}
\]

The current obstruction is that the reduced rational source `R_p` has Cartier/denominator complexity growing linearly with `p` (the earlier audit gives denominator degree between `p-1` and `2p-2`).  Thus current unit-root theory does not turn (7.1) into a bounded-degree polynomial equation in the row index.  Short-gap continuant resultants remain useful for aggregate estimates; they do not give a uniform pointwise multiplicity bound for arbitrary gaps.

---

## 8. What two actual common rows would force

Assume, toward a contradiction, that

\[
0<r<s<p,\qquad b_r=b_s=\Xi_r=\Xi_s=0.
\]

Then all of the following hold, but they are **one collision condition viewed in different coordinates**:

1. `N_h(r)=0` for `h=s-r`.
2. `N_{h-1}(r)` is a unit.
3. `M_h(r)=0`.
4. The canonical Green period vanishes:
   \[
   \sum_{m=r+1}^{s}g_mb_{m-1}=0.
   \]
5. The terminal canonical source coefficient equals the unique shooting value (2.5).
6. Equivalently, the normalized level-six Eichler companion hits `\kappa_r=\kappa_s=0` on two Hasse-zero rows.

There is no contradiction in this list from formal recurrence algebra.  The `p=41` example proves that item 4 can occur for the actual source between distinct Hasse zeros; it fails to be a common-zero counterexample only because the shared endpoint value is `7`, not `0`.  The homogeneous-boundary shift described in §4 can even move that repeated value to zero without changing the source.

Therefore an actual proof must use the **specific zero normalization of the canonical Eichler extension**, not merely source reciprocity or a vanishing Green period.

---

## 9. Exact scan through p <= 20000

The new verifier performs, for every prime `7<=p<=20000`:

```text
1. exact F_p Apéry recurrence for b_0,...,b_{p-1};
2. exact Hasse reflection check b_r=b_{p-1-r};
3. exact G=1/(F^2 sqrt(Delta)) coefficients through degree p-1;
4. exact Xi prefix sums;
5. exact canonical kappa recurrence with (0,-36) boundary line;
6. Xi/Casoratian identity at every row;
7. Xi_r=0 <-> kappa_r=0 at every Hasse zero;
8. exact Green/Duhamel endpoint identity for every Hasse-zero pair;
9. multiplicity and repeated-Xi-value census;
10. mandatory p=41 noninjectivity guard.
```

GitHub Actions run `31570501901` completed successfully.  Its theorem-relevant output is

```text
Q7696_SCAN_PMAX 20000
PRIMES_WITH_HASSE_ZEROS 903
MAX_HASSE_ZERO_COUNT 8
MAX_HASSE_ZERO_COUNT_PRIMES_FIRST20 [3727, 6841, 13463, 16451]
COMMON_RECORDS [(17, 13), (2237, 492)]
COMMON_BY_PRIME [(17, (13,)), (2237, (492,))]
COMMON_MULTIPLICITY_HISTOGRAM [(0, 2257), (1, 2)]
MAX_COMMON_MULTIPLICITY 1
REPEATED_COMMON_PRIMES []
XI_COLLISION_PRIME_COUNT 1
FIRST_XI_COLLISIONS [(41, [(10, 30, 7)])]
GUARD_P41 (10, 7) (30, 7)
Q7696_EXACT_SCAN=PASS
```

The histogram above covers the 2259 primes from `7` through `20000`.  The omitted small primes are immediate:

- `p=2`: `b_1=1 mod 2`, so no common row.
- `p=3`: `b_1=2`, `b_2=73=1 mod 3`, so no Hasse row.
- `p=5`: `Xi_r=-1 mod 5` for every `r`, since the defining correction is multiplied by `5`; hence no common row.

Thus among **all 2262 primes `p<=20000`**, exactly two have a common row and none has two.

A useful stress test is that the Hasse-zero multiplicity itself reaches eight, at `p=3727,6841,13463,16451`.  So the observed common-row simplicity is not an artifact of every Hasse zero set having size at most two.

The scan is finite data and proves nothing for primes above `20000`.

To rerun:

```bash
Q7696_PMAX=20000 sage problems/3.2/research/scripts/q7696_common_row_multiplicity.sage
```

---

## 10. Sharp conjecture and the exact missing theorem

The data give no reason to weaken the pointwise conjecture.  The sharp formulation is still

\[
\boxed{|Z_{\rm common}(p)|\le1\quad(p\ge7).}
\tag{10.1}
\]

The smallest source-specific theorem sufficient for (10.1) is:

### EIS-2RET / zero-normalized endpoint avoidance

For every prime `p>=7` and `0<r<s<p`,

\[
\boxed{
 b_r=b_s=0,\quad \Xi_r=0
 \quad\Longrightarrow\quad
 \sum_{m=r+1}^{s}b_{m-1}g_m\ne0\pmod p.}
\tag{10.2}
\]

By (2.3), this is exactly `\Xi_s\ne0`.  It says **nothing** about Green periods between arbitrary Hasse-zero rows, so it survives the `p=41` guard.

In first-block Hasse/Cartier form it is

\[
\boxed{
[t^s]R_p(t)\ne
-b_{s-1}^{-1}\sum_{m=r+1}^{s-1}
 b_{m-1}[t^m]R_p(t),}
\tag{10.3}
\]

under `b_r=b_s=\Xi_r=0`.

This is where the canonical modular/Eichler normalization has to enter.  Neither Hasse reflection, the second-order transfer, nor marked-return transversality proves (10.2).

---

## 11. What is actually proved about multiplicity today

Without a new arithmetic theorem such as (10.2), the current formal machinery does **not** prove a uniform constant or subpower pointwise multiplicity bound.

There is a trivial linear bound because homogeneous zeros cannot be consecutive:

\[
|Z_{\rm common}(p)|
\le \#\{0<r<p:b_r=0\}
\le \frac{p-1}{2}.
\]

Hasse reflection organizes those zeros but does not improve this to `p^{o(1)}` for the coefficient-zero set.

So the status should not be overstated:

```text
proved structurally: exact one-dimensional shooting law;
proved structurally: Hasse reflection and terminal coefficient is a unit;
proved computationally to p<=20000: max common multiplicity = 1;
not proved: absolute O(1) multiplicity;
not proved: p^{o(1)} pointwise multiplicity;
not proved: <=1 for all p.
```

If a future counterexample to `<=1` appears, a single isolated double return would still be compatible with an absolute `O(1)` theorem.  The next natural hierarchy would be

\[
|Z_{\rm common}(p)|\le C
\]

for an absolute `C`, then only if that also fails,

\[
|Z_{\rm common}(p)|\ll_\varepsilon p^\varepsilon.
\]

At present there is no evidence forcing this retreat: the exact scan found no double return at all.

For the dyadic repeated-prime application, a weaker aggregate far-gap theorem can suffice even without pointwise `O(1)`, but that is a different statement from the multiplicity question here.

---

## 12. Recommended next attack

The most focused next target is not another generic recurrence resultant.  It is a finite-characteristic identity for the **zero-normalized Eichler extension** that distinguishes the value `0` from the repeated nonzero value in the `p=41` example.

A useful successful theorem would have one of these forms:

```text
(A) a Fricke/Hasse reflection formula for the canonical Eichler coordinate
    whose correction term is explicit and nonzero whenever one endpoint is 0;

(B) a Cartier relation for the pair (A_p, K_p) showing that two Hasse-coefficient
    zeros with K-coordinate zero force a forbidden global Eichler period;

(C) an endpoint resultant in which the canonical boundary normalization
    contributes a fixed nonzero factor that cannot vanish at two ordinary
    returns.
```

Any candidate theorem must be tested against

```text
p=41, r=10, s=30, Xi_r=Xi_s=7,
```

and against the source-preserving homogeneous shift `Xi -> Xi+6x`.  A theorem invariant under that shift cannot prove common-zero uniqueness.

The sharp mathematical bottleneck is therefore:

\[
\boxed{
\text{Why can the canonical level-six Eichler extension hit the Hasse-zero
coefficient locus at value }0\text{ at most once?}}
\]

The present audit reduces the problem to that question, verifies all algebraic bridges exactly, and raises the exact finite search bound to `p=20000` without finding a counterexample.
