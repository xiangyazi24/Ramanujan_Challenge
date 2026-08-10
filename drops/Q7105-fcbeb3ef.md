ANSWER Q7105 fcbeb3ef

# R1 tactical audit: local algebra of the resultant pencil

## 1. Setup and the exact invariant

Let

\[
R(T)=\operatorname{Res}_x(F,G+TJ)\in \mathbb Z[T],
\]

and write

\[
C=\operatorname{cont}_T(R(T)).
\]

The relevant local object is not the length of `A/(G,J)` in general. The exact statement is a Fitting-ideal statement.

Put

\[
A_p=\mathbb Z_p[x]/(F),\qquad M=A_p/(G,J).
\]

Assume first that `F` is monic. Multiplication by `G+TJ` is an endomorphism

\[
m_{G+TJ}:A_p\otimes_{\mathbb Z_p}\mathbb Q_p[T]\to A_p\otimes_{\mathbb Z_p}\mathbb Q_p[T].
\]

Its determinant is the resultant (up to the harmless unit coming from the monic convention). The constant content of the determinant is the zeroth Fitting ideal of the cokernel after removing the free `T`-part:

\[
(C) = \operatorname{Fitt}_0\bigl(A_p[T]/(G+TJ)\bigr)\cap \mathbb Z_p.
\]

Equivalently,

\[
 v_p(C)=\operatorname{length}_{\mathbb Z_p}(\operatorname{Tor}_{\mathbb Z_p}(A_p/(G,J),\mathbb Z_p))
\]

only under an additional flatness/transversality hypothesis. Without that hypothesis, length of `A_p/(G,J)` is not the invariant measured by `C`.

The always valid statement is:

\[
 v_p(C)=\min_{t\in\mathbb Z_p}v_p(\operatorname{Res}_x(F,G+tJ)).
\]

Indeed the content of a polynomial over a DVR is the minimum valuation of its values on the DVR (after the usual Gauss lemma normalization).

---

## 2. The correct inequality

Let `\bar A=\mathbb F_p[x]/(\bar F)`. If `\bar G` and `\bar J` have `t` distinct common roots in an algebraic closure of `\mathbb F_p`, then

\[
v_p(C)\ge t.
\]

This is the statement used in the paper. It is a consequence of the local DVR factors: each distinct common geometric point contributes at least one factor of `p` to every specialization `G+TJ`.

The stronger general statement is:

For the factorization after passing to a finite unramified extension,

\[
A_p\otimes \mathcal O_K\simeq\prod_i A_i,
\]

one has

\[
v_p(C)=\sum_i \min_{t\in\mathcal O_K}v_i(G_i+tJ_i),
\]

where `v_i` is the valuation on the corresponding DVR factor. Thus the contribution is the minimum simultaneous valuation of the two elements `(G_i,J_i)` in each local branch.

It is a branchwise minimum valuation, not a global quotient length.

---

## 3. Why `length(A/(G,J))` is false

The tempting identity

\[
v_p(C)=\operatorname{length}_{\mathbb Z_p}(A_p/(G,J))
\]

fails in several independent ways.

### Counterexample 1: repeated roots

Take

\[
F=x^2,
\quad G=x,
\quad J=x.
\]

Then

\[
A=\mathbb Z_p[x]/(x^2),\qquad A/(G,J)=A/(x)\cong\mathbb Z_p,
\]

so the length is infinite (not a finite torsion length).

But

\[
\operatorname{Res}(x^2,x+Tx)=0
\]

identically. The content is not a finite length invariant.

### Counterexample 2: nilpotent multiplicity is invisible

Take

\[
F=x^2-p,
\quad G=x,
\quad J=1.
\]

Then

\[
\operatorname{Res}(x^2-p,x+T)=T^2-p.
\]

Hence

\[
C=1,
\qquad v_p(C)=0.
\]

But

\[
A/(G,J)=0
\]

because `J=1`. The quotient length gives no information.

### Counterexample 3: nonmonic scaling

Let

\[
F=p x-1.
\]

The algebra is zero after reduction modulo `p`, but

\[
\operatorname{Res}(px-1,G+TJ)
\]

contains the leading coefficient contribution of `F`. Thus the resultant content depends on the chosen presentation unless the leading coefficient is saturated away.

Therefore degree padding and nonmonic models must be normalized before interpreting the resultant.

---

## 4. What is true after the paper's structural saturation

After removing structural primes and saturating the leading coefficient, the finite flat monic model is the correct one. In that setting:

* simple common roots contribute exactly one;
* repeated roots contribute the ramification-weighted minimum valuation;
* unrelated roots contribute zero;
* the invariant is the zeroth Fitting valuation of the pencil cokernel.

Thus the paper's root-count lower bound is a coarse shadow of the exact local formula.

---

## 5. Effect of the continuant Bezout transformation

The paper compares

\[
J_1=N_{b+c}(x+a)
\]

with the direct pencil

\[
J_2=N_c(x+a+b).
\]

The continuant identity gives a Bezout relation of the form

\[
U(x)J_1+V(x)J_2=\Delta(x),
\]

where `\Delta` is the explicit boundary factor coming from the continuant determinant.

Consequently the local ideals satisfy

\[
(G,J_1)\subseteq (G,J_2,\Delta),
\]

and the reverse inclusion holds after inverting `\Delta`.

Therefore the Fitting valuations differ only at primes dividing the norm of `\Delta`:

\[
 v_p(C_1)-v_p(C_2)=0
\]

for every non-structural prime `p` with

\[
p\nmid N(\Delta).
\]

At structural primes the difference is exactly the valuation of the corresponding determinant factor. This is the explicit correction factor; there is no hidden large mass.

---

## 6. Consequence for compression to `H^{3+o(1)}`

The local-algebra reformulation is useful because it identifies the correct invariant:

\[
\text{aligned mass}=\sum_p v_p(\operatorname{Fitt}_0(M_p)).
\]

However it does not by itself yield a compression theorem. The obstruction is that the Fitting valuation is a sum over local branches of a minimum valuation. It is not a single global gcd or a polynomial gcd.

In particular:

* scalar gcds still charge unrelated branches;
* adjoining all cut-edge values still has cubic cost;
* the local algebra does not remove the need to identify which primes have simultaneous vanishing of `G` and `J`.

The possible route to `H^{3+o(1)}` is therefore not replacing the invariant by `length(A/(G,J))`, which is false, but exploiting the Fitting description together with the sparse support of exceptional primes and the continuant invariance away from explicitly controlled structural factors.

The exact conclusion is:

> `v_p(C)` is a zeroth-Fitting/minimum-valuation invariant of the pencil. It dominates the number of distinct common geometric roots modulo `p`, but it is neither equal to nor bounded above by `length_{Z_p}(A/(G,J))` without strong flatness and separability hypotheses. The continuant transformation preserves it outside the explicitly computable boundary determinant primes, so any further compression must come from distribution of exceptional local branches, not from a quotient-length replacement.
