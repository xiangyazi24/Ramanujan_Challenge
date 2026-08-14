# Q8373 — full physical top-strip selector: exact saturation, branch audit, and first nonlinear survivor

## Verdict

The canonical selector

\[
C=\sum_{i=1}^{16}\frac{R}{p_i}b_{h_i},
\qquad R=\prod_{i=1}^{16}p_i,
\]

is an exact all-target carrier, but it is maximally saturated.

For every one of the sixteen actual target primes,

\[
\boxed{p_i\mid C.}
\]

Indeed all foreign terms contain `p_i` in their coefficient and the own term contains `p_i` in `b_{h_i}`. Since the primes are distinct,

\[
\boxed{R\mid C,\qquad
\frac CR=Q:=\sum_i\frac{b_{h_i}}{p_i}\in\mathbf Z_{>0}.}
\tag{V1}
\]

More precisely,

\[
\boxed{v_{p_i}(C)=1+v_{p_i}(Q),}
\tag{V2}
\]

and

\[
p_i\mid Q
\iff
\frac{b_{h_i}}{p_i}
 +\sum_{j\ne i}b_{h_j}p_j^{-1}
\equiv0\pmod{p_i}.
\tag{V3}
\]

The cross-unit conditions make every foreign summand in (V3) a unit, but do not decide whether their sum cancels. Thus the packet forces exactly one copy of every target prime in the selector; no second copy is forced.

The quotient `Q=C/R` is certainly nonzero—indeed every summand is a positive integer and `Q>=16`—but its height is controlled by the largest **absolute folded row**, not by the four Boolean gaps. If

\[
H=\max_i h_i,
\qquad X<p_i\le2X,
\]

then

\[
\boxed{
\frac{b_H}{2X}\le Q\le\frac{16b_H}{X},
}
\tag{V4}
\]

and the elementary Apéry bounds give

\[
\boxed{
H\log16-O(\log X+\log(H+1))
\le\log Q
\le H\log64+O(\log(H+1)).
}
\tag{V5}
\]

So `C` is not a small nonzero multiple of `R`; it is a multiple at least `16R`, normally exponentially larger.

The same verdict holds for the canonical mixed Boolean finite difference. On either branch put

\[
q_n=\begin{cases}
 b_n/(m-n),&p=m-n\quad\text{(direct)},\\[2mm]
 2b_n/(m+1+n),&2p=m+1+n\quad\text{(reflected)}.
\end{cases}
\]

At a packet node, `q_{h_i}=b_{h_i}/p_i` is an integer. For a Boolean cube

\[
h_\epsilon=h_0+\epsilon_1\delta_1+\cdots+\epsilon_4\delta_4,
\]

define

\[
F_\square=
\sum_{\epsilon\in\{0,1\}^4}
(-1)^{4-|\epsilon|}q_{h_\epsilon}.
\]

The exact growth `b_{n+1}>=5b_n` implies

\[
\boxed{
F_\square\ge\frac34q_H>0
\quad\text{on the direct branch},
}
\tag{V6D}
\]

and

\[
\boxed{
F_\square\ge\frac13q_H>0
\quad\text{on the reflected branch}.
}
\tag{V6R}
\]

Hence the canonical Boolean boundary is nonzero, but again has logarithmic height `Omega(H)`. Dividing by a product of the four small gaps, or even by the full sixteen-node Vandermonde, removes only logarithmic gap height. It does not remove the exponential absolute-row height.

There is an exact saturation theorem behind this failure. Let `x_i` range over integral rows satisfying

```text
p_i | x_i,
p_j does not divide x_i  for j != i,
```

and let `L_A=sum_i A_i x_i` be an integer linear selector whose divisibility by every `p_j` is to follow universally from those zero/unit incidences. Then necessarily and sufficiently

\[
\boxed{
\frac{R}{p_i}\mid A_i\qquad(1\le i\le16).
}
\tag{V7}
\]

Thus the universal coefficient module is exactly

\[
\bigoplus_i (R/p_i)\mathbf Z\,e_i,
\]

of index

\[
\prod_i(R/p_i)=R^{15}.
\]

Every universally target-divisible linear selector therefore has the form

\[
R\sum_i c_i\frac{b_{h_i}}{p_i}.
\]

Finite differences and cleared Lagrange weights merely alter the integers `c_i`; they cannot remove the one universal `R`. If one divides out `R`, the surviving boundary can be nonzero but has no forced target divisor.

This strictly classifies the selector relative to Q8360. The Q8360 post-`R` Boolean weight

\[
V_i=\prod_{j\ne i}(p_j-p_i)
\]

satisfies

\[
R/p_i\equiv V_i\pmod{p_i}.
\]

The present selector multiplies this small local coefficient by `b_{h_i}=0` in the defining field, so it erases the Boolean derivative geometry. Q8360 instead divides the common physical value `b_m` by `R` first and uses the first `p_i`-adic digit:

\[
B V_i\equiv5z_i+\Gamma_i\pmod{p_i},
\qquad
B=b_m/R,
\quad z_i=b_{h_i}/p_i.
\tag{BRJ}
\]

That is the first interface not killed by the linear selector theorem. It still does not close TP16: `B` and the sixteen local first digits remain.

The fixed-`g` packet also changes the Q8360 height ledger. Write `m=1+gt` and `h_i=ga_i`. Then

\[
\begin{array}{ll}
\text{direct:}&p_i=1+g(t-a_i),\\
\text{reflected:}&2p_i=2+g(t+a_i).
\end{array}
\]

For odd shell primes, `gcd(g,p_i)=1` on both branches, while

\[
\boxed{
V_i^{\rm dir}=(-g)^{15}\prod_{j\ne i}(a_j-a_i),
}
\tag{V8D}
\]

and

\[
\boxed{
2^{15}V_i^{\rm ref}=g^{15}\prod_{j\ne i}(a_j-a_i).
}
\tag{V8R}
\]

If the four Boolean gaps in the `a`-coordinate have total diameter `D_a`, then

\[
|V_i^{\rm dir}|\le(gD_a)^{15},
\qquad
|V_i^{\rm ref}|\le(gD_a/2)^{15},
\]

and the full Vandermonde height contains `240 log g`. Since the proved concentration fixes `g` but does not make it `X^{o(1)}`, the BRJ coefficient height is not controlled by the four normalized Boolean gaps alone.

Therefore I do **not** obtain a target-divisor/height contradiction sufficient for the weighted full-packet TP16 gateway.

The smallest objects not covered by the linear saturation theorem are quadratic. If

\[
L_i=5z_i+\Gamma_i
\]

(or `L_i=5z_i+10D_{h_i}` on the proved direct chart), then the first common-`B` eliminating face minor is

\[
\boxed{
\mathcal M_{ij}=V_jL_i-V_iL_j.
}
\tag{V9}
\]

It has small **coefficient** height, but no divisibility by `p_i` or `p_j` follows: the two BRJ equations live in different residue fields. Conversely, if

\[
E_i=BV_i-L_i,
\]

then `p_i|E_i`, so the smallest nonlinear target-divisible carrier is

\[
\boxed{E_iE_j,\qquad p_ip_j\mid E_iE_j.}
\tag{V10}
\]

It retains `B` and the local digits and has no height compression. Thus (V9) is the smallest `B`-eliminating survivor, while (V10) is the smallest target-divisible nonlinear survivor. A new cross-characteristic theorem must connect these two properties; neither the canonical selector nor one-`R` Lagrange interpolation does so.

---

## 0. Source state and owned files

The connector-visible default branch used here is

```text
main@734a5a84c1e4fd8703a811aadaa2b4c7f532b20e
```

and the inspected `chatgpt-drop` predecessor is

```text
0d452109f1028cc0112cd8bf95da237efaa6e54f
```

(the Q8345 report commit). The same-project inputs used are:

```text
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_obstruction.md
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_verify.py
problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_obstruction.md
problems/3.2/ORACLE_COMM/chatgpt_q32_scale_sensitive_crossrow.md
```

and the connector-visible Q8360 report `Q8360 61c4588c`, especially its post-`R` Boolean jet. The user-supplied full-packet theorem—one fixed `g|m-1`, one direct/reflected branch, four matchings, and sixteen actual zero/unit nodes—is taken as the authoritative newest combinatorial input.

No shared TeX file is edited. The owned files are

```text
problems/3.2/ORACLE_COMM/chatgpt_q32_fullpacket_selector.md
problems/3.2/ORACLE_COMM/chatgpt_q32_fullpacket_selector_verify.py
```

plus the required delivery copy in `drops/`.

I did not find a materialized numerical sixteen-node packet in the visible repository state. The verifier therefore accepts such a packet as JSON and also locks two actual lower-arity branch regressions already present in the project.

---

# 1. Exact branch algebra

Let the sixteen folded rows be

\[
h_i=ga_i,
\qquad g\mid m-1,
\]

and assume the target primes are distinct odd shell primes.

## 1.1 Direct branch

Here

\[
p_i=m-h_i.
\tag{1.1D}
\]

Thus

\[
m=p_i+h_i.
\]

Gessel/Lucas gives

\[
b_m=b_{p_i+h_i}\equiv5b_{h_i}\pmod{p_i}.
\tag{1.2D}
\]

Therefore every target `p_i|b_{h_i}` also divides `b_m`.

Writing `m=1+gt`,

\[
p_i=1+g(t-a_i),
\]

so

\[
p_i\equiv1\pmod g,
\qquad \gcd(p_i,g)=1.
\tag{1.3D}
\]

## 1.2 Reflected branch

Here

\[
2p_i=m+1+h_i.
\tag{1.1R}
\]

Equivalently,

\[
m=p_i+(p_i-1-h_i).
\]

The folded reflection congruence and Gessel/Lucas give

\[
b_m
\equiv5b_{p_i-1-h_i}
\equiv5b_{h_i}\pmod{p_i}.
\tag{1.2R}
\]

Thus the reflected target also satisfies `p_i|b_m`.

Since

\[
2p_i=2+g(t+a_i),
\]

any common divisor of `p_i` and `g` divides `2`; the shell prime is odd, hence

\[
\gcd(p_i,g)=1.
\tag{1.3R}
\]

In both branches the sixteen primes divide the same actual integer `b_m`, so

\[
R\mid b_m,
\qquad B=b_m/R\in\mathbf Z.
\tag{1.4}
\]

---

# 2. Canonical selector theorem

Put

\[
x_i=b_{h_i},
\qquad q_i=x_i/p_i\in\mathbf Z_{>0}.
\]

Then

\[
C=\sum_i(R/p_i)x_i
 =R\sum_iq_i.
\tag{2.1}
\]

This proves (V1). It also proves the exact valuation formula (V2), because `R` is squarefree.

For the quotient criterion, reduce

\[
Q=\sum_jq_j
\]

modulo `p_i`. The own quotient `q_i` is defined integrally, while for `j!=i`

\[
q_j=b_{h_j}p_j^{-1}\pmod{p_i}.
\]

This gives (V3). Isolation says the foreign terms are nonzero; it does not prohibit a sixteen-term cancellation.

Two immediate consequences are worth keeping separate.

1. `C` is a nonzero target-divisible integer, but it is never smaller than `R`:

   \[
   C\ge16R.
   \]

2. `Q=C/R` is a nonzero primitive post-`R` object, but no `p_i|Q` is forced.

That is the complete linear divisor ledger.

---

# 3. Exact height of `C/R`

Let `H=max_i h_i`, and choose `k` with `h_k=H`. Since all Apéry numbers are positive and increasing,

\[
Q\ge b_H/p_k\ge b_H/(2X).
\]

Also

\[
Q\le16b_H/X.
\]

This proves (V4).

For the lower Apéry bound, retain the central summand in

\[
b_H=\sum_{r=0}^H\binom Hr^2\binom{H+r}{r}^2.
\]

The term `r=H` is `binom(2H,H)^2`, and

\[
\binom{2H}{H}\ge\frac{4^H}{2H+1}.
\]

Hence

\[
b_H\ge\frac{16^H}{(2H+1)^2}.
\tag{3.1}
\]

For the upper bound, every summand is at most `64^H`, so

\[
b_H\le(H+1)64^H.
\tag{3.2}
\]

Combining (3.1), (3.2), and (V4) proves (V5).

This is fatal to a gap-height contradiction. The selector has inherited the absolute Apéry row height before any packet topology is used.

---

# 4. Canonical finite difference: nonzero but still tall

## 4.1 Exact growth lemma

The defining sum shows `b_n>b_{n-1}` for `n>=1`. The recurrence is

\[
(n+1)^3b_{n+1}
=(34n^3+51n^2+27n+5)b_n-n^3b_{n-1}.
\]

Using `b_{n-1}<=b_n`,

\[
(n+1)^3b_{n+1}
\ge(33n^3+51n^2+27n+5)b_n.
\]

But

\[
33n^3+51n^2+27n+5-5(n+1)^3
=28n^3+36n^2+12n\ge0.
\]

Therefore

\[
\boxed{b_{n+1}\ge5b_n.}
\tag{4.1}
\]

## 4.2 Direct branch

Define the positive rational sequence

\[
q_n^{\rm dir}=\frac{b_n}{m-n}.
\]

For every `n` before the largest target row,

\[
\frac{q_{n+1}^{\rm dir}}{q_n^{\rm dir}}
=\frac{b_{n+1}}{b_n}\frac{m-n}{m-n-1}
\ge5.
\tag{4.2D}
\]

Thus

\[
\sum_{n<H}q_n^{\rm dir}\le\frac14q_H^{\rm dir}.
\tag{4.3D}
\]

The top Boolean vertex has sign `+1`, so every other signed term is bounded in absolute value by the left side of (4.3D). This proves (V6D).

## 4.3 Reflected branch

Define

\[
q_n^{\rm ref}=\frac{2b_n}{m+1+n}.
\]

Then

\[
\frac{q_{n+1}^{\rm ref}}{q_n^{\rm ref}}
\ge5\frac{m+1+n}{m+2+n}
\ge\frac52.
\tag{4.2R}
\]

Therefore

\[
\sum_{n<H}q_n^{\rm ref}\le\frac23q_H^{\rm ref},
\tag{4.3R}
\]

which proves (V6R).

## 4.4 Primitive divided boundaries

Let `D_gap` be any nonzero integer built from the four edge gaps (or from their fixed finite products) and assume

\[
\gcd(D_{\rm gap},R)=1.
\tag{4.4}
\]

This is automatic when its prime factors lie below every shell prime. The reduced numerator of

\[
\frac{R F_\square}{D_{\rm gap}}
\]

is exactly

\[
\boxed{
R\frac{F_\square}{\gcd(F_\square,D_{\rm gap})}.
}
\tag{4.5}
\]

It is nonzero and still divisible by all sixteen target primes. Moreover

\[
\left|\operatorname{num}\frac{RF_\square}{D_{\rm gap}}\right|
\ge R\frac{F_\square}{|D_{\rm gap}|}.
\tag{4.6}
\]

For the full node Vandermonde

\[
\Delta_h=\prod_{i<j}|h_i-h_j|,
\]

there are `120` factors. If the physical diameter is `D_h`, then

\[
\Delta_h\le D_h^{120}.
\tag{4.7}
\]

Thus even maximal Lagrange denominator clearing subtracts at most `120 log D_h` from a boundary whose logarithm is `Omega(H)`.

If one first divides by `R`, the reduced numerator of `F_square/D_gap` is still nonzero by (V6D)/(V6R), but it has no forced target-prime divisor. This is the exact pre-`R`/post-`R` dichotomy.

---

# 5. Fixed-`g` height in the BRJ weights

For `i!=j`, the branch formulas give

\[
p_j-p_i=-g(a_j-a_i)
\quad\text{(direct)},
\tag{5.1D}
\]

and

\[
2(p_j-p_i)=g(a_j-a_i)
\quad\text{(reflected)}.
\tag{5.1R}
\]

Multiplying over the other fifteen vertices gives (V8D) and (V8R).

If

\[
a_\epsilon=a_0+\epsilon_1d_1+\cdots+\epsilon_4d_4,
\qquad D_a=d_1+\cdots+d_4,
\]

then

\[
|a_\eta-a_\epsilon|\le D_a.
\]

Hence

\[
|V_\epsilon^{\rm dir}|\le(gD_a)^{15},
\qquad
|V_\epsilon^{\rm ref}|\le(gD_a/2)^{15}.
\tag{5.2}
\]

The product over the sixteen vertices is the squared Vandermonde:

\[
\prod_i|V_i^{\rm dir}|
=g^{240}\prod_{i<j}|a_i-a_j|^2,
\tag{5.3D}
\]

\[
2^{240}\prod_i|V_i^{\rm ref}|
=g^{240}\prod_{i<j}|a_i-a_j|^2.
\tag{5.3R}
\]

Thus the correct coefficient-height statement is

\[
\sum_i\log|V_i|
\le240\log g+240\log D_a+O(1),
\]

not `O(log D_a)` unless a separate theorem makes `g` small. The hypothesis `g|m-1` makes the branch arithmetic clean and ensures `g` is a target-prime unit; it does not bound `g` away from the physical scale.

---

# 6. Universal linear-selector saturation theorem

### Theorem 6.1

Let `p_1,...,p_N` be distinct odd primes, `R=prod_i p_i`, and let `A_i in Z`. The following are equivalent.

1. For every tuple of integers `x_i` satisfying

   ```text
   p_i | x_i,
   p_j does not divide x_i for j != i,
   ```

   every prime `p_j` divides `sum_i A_i x_i`.

2. For every `i`,

   \[
   R/p_i\mid A_i.
   \]

### Proof

The reverse implication is immediate. For the forward direction fix `j`. Modulo `p_j`, the own row vanishes and the foreign residues may vary independently through `F_{p_j}^*`. Vary only the `i`th foreign residue between two distinct units while holding the others fixed. The difference of the two universally zero sums is

\[
A_i(u-v)=0\pmod{p_j}.
\]

Choose `u-v!=0`; hence `p_j|A_i`. This holds for every `j!=i`, so `R/p_i|A_i`. QED.

For `N=16`, the coefficient lattice has index

\[
\prod_iR/p_i=R^{15}.
\]

The theorem is stronger than checking the displayed selector. It classifies **every** linear construction whose target divisibility is derived only from the complete zero/unit packet.

A rational Lagrange form does not escape. After clearing a gap-only denominator, either its integer coefficients fail condition 2, in which case target divisibility is not forced, or they satisfy condition 2, in which case the form contains the universal factor `R` and is exactly of the saturated shape.

---

# 7. Comparison with the one-`R` Lagrange obstruction

The earlier moving-target/Lagrange constructions place all defining primes in one characteristic-zero content or one product `R`. They are exact support encodings, but their candidate content has the same height scale as the target radical.

The present selector is the finite packet version of that obstruction:

- its coefficients are the CRT idempotent numerators `R/p_i`;
- its target divisibility is termwise;
- after division by `R`, the local quotients remain unrelated across characteristics;
- clearing a divided-difference denominator restores the same coefficient lattice.

The four matchings and Boolean topology do not enter the proof of `R|C`. Therefore `C` cannot exploit the extra packet geometry that created the sixteen nodes. It is a one-`R` support encoding, not a transverse boundary.

---

# 8. Comparison with Q8360 BRJ

Define

\[
V_i=\prod_{j\ne i}(p_j-p_i),
\qquad
B=b_m/R,
\qquad
z_i=b_{h_i}/p_i.
\]

Also take the canonical integral lift

\[
\widehat\Gamma_i=rac{b_m-5b_{h_i}}{p_i}\in\mathbf Z.
\tag{8.1}
\]

It exists on both branches by (1.2D)/(1.2R). Then

\[
5z_i+\widehat\Gamma_i=b_m/p_i=B(R/p_i).
\tag{8.2}
\]

Consequently the integral BRJ residual is exactly

\[
\boxed{
B V_i-5z_i-\widehat\Gamma_i
=B\left(V_i-R/p_i\right).
}
\tag{8.3}
\]

Since `R/p_i congruent V_i mod p_i`, (8.3) is divisible by `p_i`. This shows precisely what part of BRJ1 is universal coefficient algebra.

On the proved direct `p^2` chart, the additional arithmetic statement is

\[
\widehat\Gamma_i\equiv10D_{h_i}\pmod{p_i},
\]

so

\[
BV_i\equiv5z_i+10D_{h_i}\pmod{p_i}.
\tag{8.4D}
\]

That substitution is not contained in the selector theorem. It is the genuine direct Gessel first-jet input.

On the reflected branch, (8.3) remains exact, but one must retain `Gamma_i` or prove the appropriate reflected canonical lift. Reusing (8.4D) without such a theorem is invalid. The obvious integral lift (8.1) has the height of `b_m/p_i`, so it is not a small reflected boundary.

Thus Q8360 sits strictly after the selector saturation: it preserves `V_i` only by retaining a first local digit. The remaining obstruction is exactly the one Q8360 stated—cross-characteristic elimination of `B,z_i` (and the relevant canonical defects) before CRT.

---

# 9. Smallest surviving nonlinear objects

Write schematically

\[
BV_i\equiv L_i\pmod{p_i},
\]

where

\[
L_i=5z_i+\Gamma_i
\]

and, on the direct chart,

\[
L_i=5z_i+10D_{h_i}.
\]

## 9.1 Smallest `B`-eliminating object

For two vertices define

\[
\mathcal M_{ij}=V_jL_i-V_iL_j.
\]

This is quadratic in the packet rows/digits and cancels the common `B` formally. Its coefficients are built from the `V` weights, so after the fixed-`g` correction their height is controlled by `g` and the Boolean gaps.

However, modulo `p_i`, the equation at vertex `j` gives no information. Indeed

\[
\mathcal M_{ij}
\equiv V_i(BV_j-L_j)\pmod{p_i},
\]

and the bracket is constrained only modulo `p_j`. Thus neither `p_i` nor `p_j` is forced to divide `M_ij`.

## 9.2 Smallest target-divisible nonlinear object

Set

\[
E_i=BV_i-L_i.
\]

Then `p_i|E_i`, and therefore

\[
p_ip_j\mid E_iE_j.
\]

This is the first nonlinear carrier forcing two target characteristics without inserting `R/p_i` directly into each linear coefficient. But it retains `B` and the local digits, and its height is not gap-controlled.

The exact frontier is therefore:

- `M_ij`: eliminates `B`, no target divisor;
- `E_iE_j`: has target divisors, no height elimination.

A TP16 breakthrough requires a new arithmetic theorem joining those two properties across defining characteristics. Larger determinants formed without such a theorem merely iterate this dichotomy.

---

# 10. Consequence for the weighted TP16 gateway

A pointwise target-divisor/height contradiction would need a nonzero integer divisible by at least one (and preferably many) shell primes while having logarithmic height strictly below the corresponding prime product.

The canonical selector fails in the strongest possible way:

\[
C=RQ,
\qquad Q\ge16.
\]

The canonical Boolean divided boundary also fails:

- before removing `R`, its primitive numerator remains `R`-divisible but has additional height `Omega(H)`;
- after removing `R`, it is nonzero but has no forced target divisor.

The fixed-`g` normalization does not repair Q8360’s coefficient ledger because the weights retain `g^{15}` per vertex. The one-`R` Lagrange route is exactly the linear saturation theorem. The first quadratic face minor loses divisibility, and the first quadratic divisible product loses height.

Therefore the selector and its finite-difference variants do not prove the weighted full-packet TP16 gateway.

A sufficient next theorem would have to assert, for a positive proportion of actual packet faces or for the full four-cube, that a `B`-eliminating minor such as `M_ij` (or a higher face analogue) is nonzero and divisible by a specified target prime, with aggregate height controlled after the exact `g` factor is included. No currently banked selector, Lagrange, or BRJ identity supplies that cross-characteristic divisibility.

---

# 11. Regressions and verifier scope

The owned verifier locks two actual lower-arity packet regressions.

## Direct

```text
m=39,
g=2,
h=8,
p=31,
```

with

```text
p=m-h,
g | m-1,
p | b_h,
p | b_m.
```

## Reflected

The Q8345 actual rows give

```text
m=321,
g=4,
h=(36,64,100),
p=(179,193,211),
```

and

```text
2p=m+1+h,
g | m-1,
p | b_h,
p | b_m.
```

These are not represented as a sixteen-node witness. They test the branch formulas, selector factorization, exact valuation ledger, fixed-`g` Vandermonde scaling, and BRJ normalization on actual Apéry data.

For a materialized H16 witness, the verifier accepts JSON containing `m`, `g`, `branch`, and either the sixteen `h` values or `a0` plus four Boolean gaps. It recomputes the actual Apéry integers, verifies zero/unit isolation when requested, checks `C=RQ`, the primitive divided-boundary formulas, the Boolean nonvanishing lower bound, and the BRJ residuals.

The verifier is standard-library and exact. It was not executed in this response because the question explicitly forbids Python/code-interpreter execution; it is supplied for repository-side running.

---

# Final ruling

**PROVED:** all sixteen target primes divide the canonical selector; `C/R` is a positive integer; the exact valuation and cancellation criterion are (V2)–(V3); `C/R` has absolute-row exponential height; the canonical direct and reflected Boolean boundaries are nonzero but equally tall; any gap-only primitive denominator coprime to `R` leaves the pre-`R` numerator divisible by `R`; and every universally target-divisible linear selector lies in the coefficient lattice of Theorem 6.1.

**REFUTED as a TP16 route:** a small nonzero `C/R`; a primitive Boolean/Lagrange boundary whose height depends only on the four normalized gaps; reuse of the direct `10D_h` formula on the reflected branch; and the claim that the fixed-`g` BRJ weights are independent of `m`—they contain `g^{15}`.

**SMALLEST SURVIVORS:** the quadratic `B`-eliminating minor `M_ij`, which has no forced target divisor, and the quadratic target-divisible product `E_iE_j`, which has no height compression.

Hence the full-packet selector saturates exactly at one `R`. It does not close the weighted TP16 gateway; the remaining missing input is a genuinely cross-characteristic divisibility theorem for a post-`R`, `B`-eliminating nonlinear face or four-cube boundary.