ANSWER Q1270 3fe4aeaa

# Hostile audit: primitive normalization in the bordered-minor formula

## 1. Verdict on the previous cofactor statement

The previous statement is not correct as a normalization theorem.  The maximal-minor construction gives a **projective kernel vector** over the fraction field, but it does not by itself give the primitive combined Padé pair.

There are two separate content operations:

1. clearing denominators in the rational kernel computation;
2. making the resulting numerator/denominator pair primitive after applying the Apéry transfer matrix.

They cannot be collapsed into one gcd of rational maximal minors.

Let

\[
N=U^{-1}D^{-1}U,
\qquad C_a=N[H-a+1:H,0:a].
\]

The vector of signed maximal minors of `C_a` is a generator of

\[
\ker_{\mathbb Q}(C_a),
\]

but the phrase "divide by the gcd of the rational minors" is meaningless until an integral lattice has been chosen. Different denominator clearings produce different integral representatives.

Moreover the challenge object is not merely a primitive kernel vector. It is the primitive pair

\[
(p_a,q_a),\qquad q_a=Np_a,
\]

with

\[
\gcd(p_{a,0},\ldots,p_{a,a},q_{a,0},\ldots,q_{a,H-a})=1.
\]

A primitive integer kernel vector for `C_a` can acquire a common factor after multiplication by `N`; conversely a denominator clearing chosen on `q_a` can force a compensating rescaling of `p_a`.

Therefore the bordered determinant formula must include all content factors.

---

## 2. Integral formulation

Use instead

\[
M=U^{-1}DU.
\]

Let

\[
T_a=M[a+1:H,0:H-a].
\]

This is an integral matrix of size

\[
(H-a)\times(H-a+1).
\]

The denominator vector is an integer kernel vector:

\[
T_a q_a=0.
\]

For each column index `j`, define the signed maximal minors

\[
\Delta_{a,j}=(-1)^j\det(T_a^{\widehat j}),
\]

where the hat means deletion of column `j`.

The determinantal divisor is

\[
\delta_a=\gcd_j(\Delta_{a,j}).
\]

Then

\[
\widetilde q_{a,j}=\Delta_{a,j}/\delta_a
\]

is the primitive integral kernel vector, up to a global sign.

However this is still only the primitive denominator vector.  The numerator produced by the transfer matrix is

\[
\widetilde p_a=M\widetilde q_a.
\]

Its first `H-a` rows vanish, so the support is in rows `0,\ldots,a`.

Define the combined content

\[
\gamma_a=\gcd\bigl(\widetilde p_{a,0},\ldots,\widetilde p_{a,a},
\widetilde q_{a,0},\ldots,\widetilde q_{a,H-a}\bigr).
\]

The actual primitive pair is therefore

\[
\boxed{
 p_a=\widetilde p_a/\gamma_a,
 \qquad
 q_a=\widetilde q_a/\gamma_a .
}
\]

The scalar `gamma_a` is not generally equal to one and is not visible from the kernel minors alone.

---

## 3. Exact bordered-minor quotient

Let `B_{a,k}` be the bordered determinant obtained by adjoining the evaluation column corresponding to the basis vector for

\[
\Phi_k(x(x+1))=\binom{x+k}{2k}.
\]

The unnormalized numerator coefficients are

\[
\widetilde p_{a,k}=(M\Delta_a)_k,
\]

where

\[
\Delta_a=(\Delta_{a,0},\ldots,\Delta_{a,H-a}).
\]

Equivalently, by Cramer expansion,

\[
\widetilde p_{a,k}=\pm B_{a,k}.
\]

The evaluation is therefore not simply a bordered determinant divided by a rational gcd.  The exact expression is

\[
\boxed{
E_a(x)=
\frac{1}{\gamma_a\delta_a}
\sum_{k=0}^{a}
(\pm B_{a,k})\binom{x+k}{2k}.
}
\]

The denominator vector has the analogous formula

\[
\boxed{
q_{a,j}=\frac{\Delta_{a,j}}{\delta_a\gamma_a}.
}
\]

All signs are fixed only up to the simultaneous sign change

\[
(p_a,q_a)\mapsto(-p_a,-q_a).
\]

The invariant object is the primitive projective class.

---

## 4. What survives denominator clearing?

For the rational matrix `C_a`, a row-wise clearing or global clearing changes the integer lattice. Hence the gcd of rational maximal minors is not invariant.

For the integral matrix `M`, the following are invariant:

* the rank of `T_a`;
* the one-dimensional rational kernel line;
* the determinantal divisor
  \[
  \delta_a=\Delta_{H-a}(T_a),
  \]
  the gcd of all maximal minors of `T_a`;
* the Smith invariant factors of `T_a`.

The following are not invariant:

* the raw maximal minors after arbitrary denominator multiplication;
* the primitive kernel vector before fixing an integral lattice;
* the numerator scale before the combined content reduction.

---

## 5. Small cases `a=0,1,2`

### Case `a=0`

`T_0` has shape `H x (H+1)`.  The kernel vector is the signed maximal-minor vector

\[
\Delta_0.
\]

After

\[
q_0=\Delta_0/(\delta_0\gamma_0),
\]

we have

\[
p_0=Mq_0,
\]

with support only in row zero:

\[
p_0=(p_{0,0},0,\ldots,0).
\]

Hence

\[
E_0(x)=p_{0,0}.
\]

The missing factor in the old formula is exactly `gamma_0`.

### Case `a=1`

Now

\[
T_1=M[2:H,0:H-1].
\]

The signed maximal minors give the projective denominator direction.  Multiplying by `M` gives

\[
\widetilde p_1=(\widetilde p_{1,0},\widetilde p_{1,1},0,\ldots,0).
\]

The primitive pair is

\[
(p_1,q_1)=
(\widetilde p_1,\widetilde q_1)/\gamma_1.
\]

The evaluation is

\[
E_1(x)=
\frac{\widetilde p_{1,0}+\widetilde p_{1,1}\binom{x+1}{2}}
{\gamma_1}.
\]

### Case `a=2`

Similarly,

\[
T_2=M[3:H,0:H-2],
\]

and

\[
\widetilde p_2=(\widetilde p_{2,0},\widetilde p_{2,1},\widetilde p_{2,2},0,\ldots,0).
\]

The exact normalized evaluation is

\[
E_2(x)=
\frac{
\widetilde p_{2,0}
+\widetilde p_{2,1}\binom{x+1}{2}
+\widetilde p_{2,2}\binom{x+2}{4}}
{\gamma_2}.
\]

These examples show the structural defect: the maximal minors determine the numerator only before the final pair-content reduction.

---

## 6. Adjacent `a` values and Jacobi/Smith relations

After the corrected normalization, adjacent values `a` and `a+1` are not related by a simple raw Jacobi identity on primitive coefficients.  Jacobi identities apply to minors before content reduction.

The correct statement is:

* the unnormalized vectors `Delta_a` and `Delta_{a+1}` satisfy Plucker/Jacobi relations inside the common determinantal lattice;
* the Smith data of the corresponding bordered blocks control the possible common factors;
* the primitive Padé pairs require dividing by the independent contents `delta_a gamma_a` and `delta_{a+1} gamma_{a+1}`.

Therefore an adjacent recurrence may exist at the level of determinantal representatives, but any claimed identity between primitive pairs must include the content quotients.

In particular, there is no usable adjacent-`a` Jacobi relation of the form

\[
p_{a+1}=J(p_a)
\]

without tracking these Smith/determinantal divisors.

The correct frontier statement is thus: the bordered-minor construction provides a canonical projective representative; the primitive Apéry Padé pair is obtained only after the two-stage normalization

\[
\text{maximal minors}
\longrightarrow
\frac{1}{\delta_a}
\longrightarrow
\frac{1}{\gamma_a}
\longrightarrow
(p_a,q_a).
\]

The omitted second content factor is the normalization defect in the previous claim.
