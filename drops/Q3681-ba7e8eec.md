ANSWER Q3681 ba7e8eec

# P3.2: smallest honest nonlocal Smith/Fitting candidate

## Verdict

I do **not** see a proved universal Smith column with the requested selection sensitivity yet. The smallest mathematically natural object worth pursuing is a **global Lucas-Bockstein orbit, compressed by a Green/Casoratian (Wronskian) pairing**. It is genuinely nonlocal and can be presented by an integer band matrix before the hits are selected. But the Green identity by itself does **not** make the other endpoint's selectedness visible modulo the first endpoint.

There is a sharp algebraic reason. If `p != q`, then modulo `p` the integer `q` is a unit. Therefore merely adjoining the selected quotient

\[
h_q=b_{r_q}/q
\]

through the equation `q h_q=b_{r_q}` cannot create a new `p`-rank defect: over `F_p` that variable can be eliminated. A successful construction must prove a genuinely **cross-prime identity among universal integer coefficients** whose reduction modulo `p` changes when `q|b_{r_q}`. That is the non-tautological arithmetic input still missing.

For the genuine test edge

\[
(m,p,q)=(2932,439,443)
\]

we have

\[
r_p=2932-6\cdot439=298,\qquad
r_q=2932-6\cdot443=274,
\]

so the residue gap is

\[
r_p-r_q=24=6(q-p).
\]

I would use this edge to test the two exact constructions below. The first is a minimal transport/Wronskian sanity check and demonstrates the obstruction. The second is the smallest genuinely nonlocal candidate I think is still alive.

---

## 1. Exact barrier: a selected quotient does not survive as selectedness modulo the other prime

Suppose a proposed universal integer presentation introduces a variable `h` with a row

\[
q h-y=0.
\]

After tensoring with `F_p`, where `p != q`, `q` is invertible, hence

\[
h=q^{-1}y.
\]

Thus elementary row/column elimination over `F_p` removes this quotient variable. The mod-`p` presentation has the same rank information as the presentation in which `h` is substituted by `q^{-1}y`.

Consequences:

* It is fine for a **kernel vector after selection** to contain `h_q=b_{r_q}/q`.
* It is **not** fine to argue that the mere integrality of `h_q` creates the `p`-kernel. Modulo `p`, the same coordinate `q^{-1}b_{r_q}` exists whether or not `q|b_{r_q}` over `Z`.
* Putting `h_q`, `b_m/q`, or a remainder test for `b_{r_q} mod q` into the **matrix coefficients** is selected-set-dependent (or simply computes the selection predicate in advance) and is circular for the intended carrier.

So the desired mechanism cannot be `q`-division alone. It has to be an identity of the form

\[
\text{(both endpoint hits)}\quad\Longrightarrow\quad
p\mid \mathcal G_{p,q}
\]

for a **universally prescribed integer** `\mathcal G_{p,q}` (or a universally prescribed column/minor), with the implication proved from Apéry arithmetic and not from inserting `b_{r_q}/q` by hand.

---

## 2. Formula A: the exact 24-step transport bridge

This is the smallest exact nonlocal bridge between the two residue states. It is useful as a diagnostic, but it does **not** solve the selection-sensitivity problem.

Write

\[
P(n)=34n^3+51n^2+27n+5,
\]

and define the integral transfer matrix

\[
R_n=
\begin{pmatrix}
P(n)&-n^3\\
(n+1)^3&0
\end{pmatrix}.
\]

For

\[
s_n=\binom{b_n}{b_{n-1}},
\]

the Apéry recurrence is exactly

\[
R_n s_n=(n+1)^3s_{n+1}.
\tag{2.1}
\]

For `u<v`, put

\[
T_{u,v}=R_{v-1}\cdots R_u
 =\begin{pmatrix}A&B\\ C&E\end{pmatrix},
\qquad
D_{u,v}=\prod_{n=u}^{v-1}(n+1)^3
       =\left(\frac{v!}{u!}\right)^3.
\]

Then

\[
T_{u,v}s_u=D_{u,v}s_v.
\tag{2.2}
\]

At the test edge take `u=274`, `v=298`. This is only a product of 24 explicit `2 x 2` integer matrices, so it is an easy exact computation.

If both hits occur, write

\[
b_{274}=443h_q,\qquad b_{298}=439h_p.
\]

The first row of (2.2) gives the exact integer identity

\[
443A h_q+B b_{273}=439D_{274,298}h_p.
\tag{2.3}
\]

Hence modulo 439,

\[
\begin{pmatrix}443A&B\end{pmatrix}
\binom{h_q}{b_{273}}=0.
\tag{2.4}
\]

So the obvious right-kernel vector is

\[
\boxed{\kappa_{443\to439}=\binom{h_q}{b_{273}}.}
\]

There is a symmetric exact formula without inverting `T`. Since

\[
\det R_n=n^3(n+1)^3,
\]

we have

\[
\det(T_{u,v})s_u=D_{u,v}\operatorname{adj}(T_{u,v})s_v.
\]

Its first coordinate at `(u,v)=(274,298)` is

\[
\det(T)\,443h_q
 =D\bigl(439E h_p-Bb_{297}\bigr).
\tag{2.5}
\]

All factors in `D` are units modulo 443, so modulo 443

\[
\begin{pmatrix}439E&-B\end{pmatrix}
\binom{h_p}{b_{297}}=0,
\tag{2.6}
\]

with right-kernel vector

\[
\boxed{\kappa_{439\to443}=\binom{h_p}{b_{297}}.}
\]

### Why this is not the desired Smith column

Modulo 439, equation (2.3) is simply

\[
A b_{274}+B b_{273}=0
\]

because `b_298=0 mod 439`. Since 443 is a unit modulo 439, the coordinate `h_q` in (2.4) can always be replaced by

\[
443^{-1}b_{274}\pmod{439},
\]

even if `443` does **not** divide `b_274` in `Z`.

Thus the mod-439 kernel relation is caused by the **439-hit at the right endpoint**, not by the 443-hit at the left endpoint. The symmetric statement has the same defect. This is exactly the two-state/one-Smith-slot obstruction in another guise.

Also note that

\[
\det T_{u,v}=\prod_{n=u}^{v-1}n^3(n+1)^3
\]

is pure structural content. Factors coming from this determinant, from `D`, or from the explicit primes appearing as scalar multipliers are **raw factors**, not selected incidence. They must not be charged as if they recorded a second endpoint hit.

So Formula A is a useful exact control computation, but not a solution.

---

## 3. Formula B: a full Lucas-Bockstein orbit and its Green/Wronskian class

This is the smallest construction I would now test seriously.

Use the standard quotient-six Lucas divisibility already underlying the divided-Lucas work: for a prime `ell` and `0<=t<ell`,

\[
b_{6\ell+t}\equiv b_6b_t\pmod\ell.
\]

Therefore the **full first Bockstein orbit**

\[
\boxed{
Y_\ell(t)=\frac{b_{6\ell+t}-b_6b_t}{\ell}
}
\tag{3.1}
\]

is a universally defined integer for every admissible `t`. This is not a vertex-local jet: retain the whole orbit on a macroscopic interval before taking any endpoint functional.

At the test edge the common-`m` boundary values are

\[
443Y_{443}(274)=b_{2932}-b_6b_{274},
\tag{3.2}
\]

\[
439Y_{439}(298)=b_{2932}-b_6b_{298}.
\tag{3.3}
\]

Subtracting gives the exact cross identity

\[
439Y_{439}(298)-443Y_{443}(274)
 =b_6(b_{274}-b_{298}).
\tag{3.4}
\]

If both hits occur, then

\[
439\bigl(Y_{439}(298)+b_6h_p\bigr)
=443\bigl(Y_{443}(274)+b_6h_q\bigr)
=b_{2932}.
\tag{3.5}
\]

Equation (3.5) is useful, but it is still **not** selected-sensitive modulo the other prime: over `F_439`, `443` is a unit and `h_q` can again be eliminated. This explains conceptually why merely taking another local divided-Lucas row cannot fix the problem.

The reason to keep the **whole** orbit `Y_ell(t)` is that it satisfies an exact inhomogeneous recurrence, and the forcing has a canonical Green/Casoratian class.

### 3.1 Exact inhomogeneous recurrence

Let

\[
(Lx)_t=(t+1)^3x_{t+1}-P(t)x_t+t^3x_{t-1}.
\]

Put `X_ell(t)=b_{6ell+t}`. Define

\[
A^+_\ell(t)
 =\frac{(t+1+6\ell)^3-(t+1)^3}{\ell}
 =18(t+1)^2+108\ell(t+1)+216\ell^2,
\]

\[
A^-_\ell(t)
 =\frac{(t+6\ell)^3-t^3}{\ell}
 =18t^2+108\ell t+216\ell^2,
\]

and

\[
B_\ell(t)
 =\frac{P(t+6\ell)-P(t)}{\ell}
 =612t^2+612t+162
   +\ell(3672t+1836)+7344\ell^2.
\]

Direct subtraction of the shifted and unshifted Apéry recurrences gives

\[
\boxed{
(LY_\ell)_t
 =F_\ell(t)
 =-A^+_\ell(t)X_\ell(t+1)
   +B_\ell(t)X_\ell(t)
   -A^-_\ell(t)X_\ell(t-1).
}
\tag{3.6}
\]

This identity is exact over `Z`.

For a low-height presentation one should **not** store the huge `X_ell(t)` as unrelated coefficients. Substitute

\[
X_\ell(t)=b_6b_t+\ell Y_\ell(t)
\]

back into the shifted recurrence. Equivalently,

\[
\boxed{
\begin{aligned}
&(t+1+6\ell)^3Y_\ell(t+1)
-P(t+6\ell)Y_\ell(t)
+(t+6\ell)^3Y_\ell(t-1)\\
&\qquad=-b_6\bigl(
 A^+_\ell(t)b_{t+1}
-B_\ell(t)b_t
+A^-_\ell(t)b_{t-1}
\bigr).
\end{aligned}}
\tag{3.7}
\]

Now every recurrence coefficient is polynomial-size in `t,ell`; the only Apéry orbit is the **shared** base orbit `b_t`. This is the right form for an integral band/Fitting presentation.

### 3.2 Exact Green/Wronskian functional

Define

\[
W_t^{(\ell)}
 =b_tY_\ell(t+1)-b_{t+1}Y_\ell(t).
\tag{3.8}
\]

Using `Lb=0` and (3.6), one gets the exact discrete Green identity

\[
\boxed{
(t+1)^3W_t^{(\ell)}-t^3W_{t-1}^{(\ell)}
=b_tF_\ell(t).
}
\tag{3.9}
\]

Hence for every `u<v`,

\[
\boxed{
v^3W_{v-1}^{(\ell)}-u^3W_{u-1}^{(\ell)}
 =\sum_{t=u}^{v-1}b_tF_\ell(t).
}
\tag{3.10}
\]

For the test edge define the exact integer

\[
\boxed{
G_{443}
:=298^3W_{297}^{(443)}-274^3W_{273}^{(443)}
 =\sum_{t=274}^{297}b_tF_{443}(t).
}
\tag{3.11}
\]

This is a 24-step **nonlocal** quantity, not a local divided jet. It is immediately suitable for exact integer computation at `(2932,439,443)`.

There is also a canonical cross-Bockstein Wronskian. Put

\[
\Omega_t^{(p,q)}
=Y_p(t)Y_q(t+1)-Y_p(t+1)Y_q(t).
\]

Because `LY_p=F_p` and `LY_q=F_q`, the Lagrange identity gives

\[
\boxed{
(t+1)^3\Omega_t^{(p,q)}-t^3\Omega_{t-1}^{(p,q)}
=Y_p(t)F_q(t)-Y_q(t)F_p(t).
}
\tag{3.12}
\]

Thus at the test edge

\[
\boxed{
\begin{aligned}
H_{439,443}
&:=298^3\Omega_{297}^{(439,443)}
   -274^3\Omega_{273}^{(439,443)}\\
&=\sum_{t=274}^{297}
 \bigl(Y_{439}(t)F_{443}(t)-Y_{443}(t)F_{439}(t)\bigr).
\end{aligned}}
\tag{3.13}
\]

This is probably the **smallest genuinely two-endpoint nonlocal scalar** I would compute next. It couples the two full quotient-six Bockstein orbits before taking a boundary determinant.

But I emphasize: (3.9)--(3.13) are identities for **all** admissible `p,q`. I do not see an identity in the current material proving

\[
439\mid G_{443},\qquad
439\mid H_{439,443},
\]

or the symmetric 443-divisibility **because of the other hit**. Such a divisibility must be checked, then proved symbolically if it holds. It cannot be inferred from the Green identity alone.

A failure of primitive `G_443 mod 439` or `H_439,443 mod 439` at the genuine edge would kill that particular boundary functional immediately, without killing the whole Bockstein-orbit strategy.

---

## 4. The corresponding universal Fitting presentation

The natural matrix is the presentation of the **coupled base orbit plus Bockstein source orbit**, not a 2-column endpoint-state matrix.

For a fixed macroscopic interval `I=[u,v]`, start with variables for the base states `b_t` and, for each **prescribed candidate label** `ell` in the block, variables for the Bockstein orbit `Y_ell(t)`. Put in rows for

1. the ordinary Apéry recurrence for `b_t`; and
2. the coupled shifted recurrence (3.7) for `Y_ell(t)`.

All candidate `ell` are included before knowing whether `ell|b_{r_ell}`. No selected quotient occurs in the matrix coefficients.

After eliminating interior variables, the Schur/Fitting boundary columns are exactly transfer/Green objects of the type (3.10)--(3.13). This is the mathematically natural place where a selected-neighbor column **could** live.

If `q` is selected, an integral augmented solution vector may contain

\[
h_q=b_{r_q}/q.
\]

For the test edge the naive augmented right-kernel vector would have the shape

\[
\widetilde\kappa_{443\to439}
=\bigl((Y_{443}(t))_{t\in I},\ h_q,\ 1\bigr)^T
\tag{4.1}
\]

(together with whatever shared base-orbit coordinates are kept in the presentation). The recurrence rows annihilate the `Y_443` part exactly, and the boundary row `443h_q-b_274=0` is integral when the 443-hit occurs.

However, **(4.1) is not yet a selected Smith vector**: modulo 439 the boundary equation solves `h_q=443^{-1}b_274` for every candidate 443, selected or not. If the only extra row is this quotient row, the construction is tautological.

What is needed is an additional **universal Green/Fitting boundary row** `Lambda_{p,q}` such that

\[
\Lambda_{p,q}\widetilde\kappa_{q\to p}\equiv0\pmod p
\]

is a consequence of the two selected hits and is **not** already a consequence of the `p`-hit alone after eliminating `h_q`. Formulae (3.11) and (3.13) are concrete first candidates for the coefficient appearing in such a boundary row, but there is presently no proved selected-divisibility identity for them.

This gives a clean design contract: **the matrix can be universal; the kernel vector may depend on the selected hit; but after eliminating all coordinates multiplied by units modulo `p`, a genuinely new `p`-relation must remain.**

---

## 5. What exact arithmetic theorem is still missing?

There are really two independent missing theorems.

### A. Cross-prime selectedness transduction

One needs a universal primitive boundary/Fitting class `C_I(p,q)` built from the full orbit/Green data such that, for a common quotient-six collision,

\[
p\mid b_{r_p},\quad q\mid b_{r_q}
\quad\Longrightarrow\quad
C_I(p,q)\equiv0\pmod p,
\tag{5.1}
\]

where (5.1) is **not** true merely from `p|b_{r_p}` and is not obtained by writing `b_{r_q}=qh_q` and then dividing by the unit `q mod p`.

Equivalently, after localizing the presentation at `p`, the `q`-selected condition must make a new universal coefficient/minor vanish. That is a genuine cross-characteristic statement.

The most plausible source is a **second, global Bockstein/Frobenius identity** for the entire quotient-six orbit, not another local jet. Schematically, one would need an extra divisibility of a Green/Casoratian numerator after the first universal Lucas division. I do not claim such an identity exists; proving it is exactly the research problem.

A useful falsification protocol is:

* compute (3.11) and (3.13) at `(2932,439,443)`;
* remove every obvious structural content factor from the transfer/Green presentation;
* test whether a genuinely **extra** factor 439 (respectively 443) remains;
* then repeat at several genuine and non-genuine pairs before spending effort on a symbolic proof.

A raw factor 439 or 443 that comes from `D`, `det T`, an explicit `p`/`q` multiplier, or a row content inserted by clearing denominators is not evidence.

### B. Neighbor independence/full spark

Even if each selected `q` gives one valid mod-`p` kernel vector, repeated neighbors only yield Smith multiplicity if those vectors are independent. The Q3667 audit already isolates this point: distinct gaps/residues do not imply vector independence.

For candidate neighbors `q_1,...,q_rho` assigned to the same endpoint `p`, one still needs a nonzero coordinate minor

\[
\det\bigl[\kappa_{q_1}\ \cdots\ \kappa_{q_\rho}\bigr]_{S}
\not\equiv0\pmod p
\tag{5.2}
\]

for some canonical coordinate set `S`, uniformly in the block. A Vandermonde-looking parametrization is not enough unless an exact determinant identity proves (5.2).

So even a successful cross-Bockstein identity solves only the **activation** problem; it does not automatically solve the **multiplicity** problem.

---

## 6. Height audit

There is one encouraging point and one serious danger.

### Encouraging: the coupled recurrence can be low-height

The 24-step transport entries in Formula A have logarithmic size

\[
O((v-u)\log T).
\]

For a fixed-`J` macroscopic interval of length `O_J(T)`, a single continuant/band determinant with polynomial-size recurrence coefficients has log-height on the order of

\[
O_J(T\log T),
\]

which is safely

\[
o_J(T^2/\log T).
\]

Likewise, the coupled formulation (3.7) uses only polynomial-size coefficients in `t` and `ell`; it avoids putting the exponentially large shifted Apéry values into the matrix as unrelated coefficients.

This is the main reason the Bockstein-orbit/Fitting route is not obviously dead on height.

### Danger: duplicating a whole orbit for every candidate prime is fatal or near-fatal

If one takes a separate length-`Theta(T)` orbit block for each of `Theta(T/log T)` candidate neighbors, the presentation has essentially quadratic total size. Likewise, if one keeps one raw source column per candidate whose entries have `log|b_n|=Theta(T)`, a maximal determinant using many such columns naturally reaches the `T^2/log T` scale.

Therefore a viable construction must **share the recurrence operator and eliminate the interior**, leaving only a small boundary/Fitting column per candidate, with a determinant-collapse identity that prevents the heights of those columns from adding naively.

This is where the Green/Wronskian form matters: it gives a canonical Schur complement of a band system. But one still has to prove the primitive maximal-minor gcd of the resulting fixed-`J` block is subcritical. Entrywise bounds alone will not do it.

Also, all predictable transfer content

\[
\prod n^3(n+1)^3,
\]

factorial cubes, denominator-clearing factors, and explicit scalar `p,q` factors must be stripped before interpreting the primitive determinantal divisor. Otherwise the apparent Smith mass is structural rather than selected.

---

## 7. Recommended exact computation at `(2932,439,443)`

I would do the following before designing a larger matrix.

1. **Transport control.** Compute
   \[
   T_{274,298}=R_{297}\cdots R_{274}
   \]
   and verify (2.3)--(2.6) exactly. This confirms the expected tautological kernel and records the structural gcd/content that must later be removed.

2. **Full Bockstein orbits.** Compute
   \[
   Y_{443}(t)=\frac{b_{2658+t}-b_6b_t}{443},\qquad
   Y_{439}(t)=\frac{b_{2634+t}-b_6b_t}{439}
   \]
   for `t=273,...,298` (all these indices are below both primes where needed for the Lucas division).

3. **Verify the forcing identity.** Check (3.6), or preferably the lower-height coupled form (3.7), at every `t=274,...,297`.

4. **Compute the primitive Green classes** `G_443`, `G_439` and the cross class `H_439,443` from (3.11)--(3.13). Factor only after dividing out deterministic row/transfer content. The key question is whether the **other endpoint prime** occurs as an extra factor.

5. If an extra factor survives at this genuine edge, test the same primitive class on non-hit pairs with comparable `(p,q,m)`. A factor that appears for every `q` once the `p`-hit is imposed is still structural and gives no selection sensitivity.

6. Only if step 5 survives should one build a many-neighbor block and test full-spark minors. Distinct `d=q-p< p/6` or distinct residues are only distinct labels, not an independence proof.

---

## Bottom line

The smallest construction I would keep alive is:

> **full quotient-six Lucas-Bockstein orbit -> exact inhomogeneous recurrence -> Green/Casoratian Schur complement -> primitive Fitting column.**

It is nonlocal, universal before selection, and has a plausible route to subcritical height because the coupled recurrence (3.7) has only polynomial-size coefficients.

But it is **not yet a solution**. The exact missing identity is a cross-characteristic selectedness transducer: after all `q`-unit quotient variables are eliminated modulo `p`, the `q`-hit must still force a new universal Green/Fitting coefficient to vanish modulo `p`. Formulae (3.11) and (3.13) are natural exact test objects for that phenomenon; no such divisibility follows from the recurrence identities written above.

And even if that activation theorem is found, a separate full-spark/minor identity is still required to turn multiple selected neighbors into repeated `p`-corank.
