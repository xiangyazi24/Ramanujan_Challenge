# Q7735 — finite-part transport for the Apéry rational function

## Verdict

Let

\[
F_n(X):=R_n(X)^2,
\qquad
R_n(X)=\frac{\prod_{u=1}^n(X+u)}{\prod_{u=0}^n(X-u)}.
\]

The rational-function route closes all of the local work requested in Q7735:

1. the displayed harmonic formula for `e_n` is recovered from the *global* finite-part trace, with an explicit warning that it is **not** a pole-by-pole identification;
2. the exact old-pole finite-part transformations under `n -> n+1` and `n -> n+2` follow from one universal multiplier lemma;
3. all new-pole terms are explicit differential boundary operators;
4. the central-binomial square comes canonically from the new double-pole coefficient;
5. the complete recurrence reduces to one explicit, harmonic-free trace identity in
   \(\mathbf Q[X]/\prod_{a=0}^n(X-a)\).

I do **not** have a proof of that last trace identity for symbolic `n`.  The
quotient-ring and Hermite formulations expose it cleanly, but a direct bounded-degree
collapse does not occur: in the exact Sage audit the reduced trace integrand has
degree `n` for every `1 <= n <= 8`, while the canonical Hermite correction has
degree `2n+1`.  Thus declaring the recurrence proved from the present rational
manipulation would merely hide the hard summation inside quotient reduction.

The exact missing identity is boxed in §8 below.  It contains no harmonic
numbers and is substantially smaller than a harmonic WZ certificate.

Verifier:

`problems/3.2/research/scripts/q7735_finite_part_trace_verify.sage`

The successful Sage 10.6 run ends with `Q7735_SYMBOLIC_VERIFY PASS`.

---

## 1. First normalization: the harmonic summand is not a local finite part

Write the exact partial fraction expansion as

\[
F_n(X)=\sum_{b=0}^n T_{n,b}
\left(\frac1{(X-b)^2}+\frac{2q_{n,b}}{X-b}\right).
\]

Let

\[
C_{n,a}:=\operatorname{FP}_{X=a}F_n(X)
\]

be the regular finite part at the pole `a`.  Reading the finite part from the
*other* partial fractions gives

\[
C_{n,a}=\sum_{\substack{b=0\\b\ne a}}^n
T_{n,b}\left(\frac1{(a-b)^2}+\frac{2q_{n,b}}{a-b}\right).
\tag{1}
\]

Therefore

\[
\begin{aligned}
-\sum_{a=0}^n C_{n,a}
&=-\sum_{b=0}^n T_{n,b}
\left(
\sum_{a\ne b}\frac1{(a-b)^2}
+2q_{n,b}\sum_{a\ne b}\frac1{a-b}
\right)\\
&=\sum_{b=0}^nT_{n,b}\left[
-H_b^{(2)}-H_{n-b}^{(2)}
+2q_{n,b}(H_b-H_{n-b})
\right].
\end{aligned}
\]

This proves the stated equivalence

\[
\boxed{e_n=-\sum_{a=0}^n C_{n,a}.}
\tag{2}
\]

But (2) is only a **trace identity**.  The summand displayed in the definition
of `e_n` is not `-C_{n,a}` at the same `a`.

The smallest exact regression is `n=1`:

\[
(-C_{1,0},-C_{1,1})=(-8,-5),
\]

whereas the two displayed harmonic summands are

\[
(-5,-8).
\]

They have the same sum `-13`, but they are swapped.  Any proof that transports
the displayed summand pole-by-pole is therefore false.  Everything below
transports the genuine finite parts `C_{n,a}`.

---

## 2. Universal finite-part multiplication lemma

Suppose near `X=a`, with `z=X-a`,

\[
F(X)=\frac{T}{z^2}+\frac{2qT}{z}+C+O(z),
\]

and `A` is analytic at `a`:

\[
A(X)=A(a)+A'(a)z+\frac12A''(a)z^2+O(z^3).
\]

Multiplying and reading the constant term gives

\[
\boxed{
\operatorname{FP}_a(AF)
=A(a)C+2A'(a)qT+\frac12A''(a)T.}
\tag{3}
\]

Equivalently, if

\[
\varepsilon_{n,a}:=-C_{n,a},
\]

then

\[
\boxed{
\varepsilon(AF)
=A(a)\varepsilon-2A'(a)qT-\frac12A''(a)T.}
\tag{4}
\]

No harmonic manipulation occurs in (3) or (4).

---

## 3. Exact `n -> n+1` transformation

Put

\[
c=n+1,
\qquad
A_c(X):=\left(\frac{X+c}{X-c}\right)^2.
\]

Then

\[
F_{n+1}=A_cF_n.
\]

At every old pole `0 <= a <= n`,

\[
\boxed{
\varepsilon_{n+1,a}
=A_c(a)\varepsilon_{n,a}
-2A_c'(a)q_{n,a}T_{n,a}
-\frac12A_c''(a)T_{n,a}.}
\tag{5}
\]

The derivatives are completely rational.  Define

\[
L_c(X):=\frac{4c}{c^2-X^2}.
\]

Then

\[
\frac{A_c'}{A_c}=L_c,
\qquad
\frac{A_c''}{A_c}
=\frac{8c(2c+X)}{(c^2-X^2)^2}.
\]

Thus (5) may be written as

\[
\boxed{
\varepsilon_{n+1,a}
=A_c(a)\left[
\varepsilon_{n,a}
-\frac{8c}{c^2-a^2}q_{n,a}T_{n,a}
-\frac{4c(2c+a)}{(c^2-a^2)^2}T_{n,a}
\right].}
\tag{6}
\]

Also

\[
A_c(a)=\frac{T_{n+1,a}}{T_{n,a}},
\]

so (6) has exactly the expected principal-part normalization.

---

## 4. Exact `n -> n+2` transformation

Put

\[
d=n+2=c+1,
\qquad
B_{c,d}(X):=A_c(X)A_d(X).
\]

Then

\[
F_{n+2}=B_{c,d}F_n.
\]

At every old pole,

\[
\boxed{
\varepsilon_{n+2,a}
=B_{c,d}(a)\varepsilon_{n,a}
-2B_{c,d}'(a)q_{n,a}T_{n,a}
-\frac12B_{c,d}''(a)T_{n,a}.}
\tag{7}
\]

If

\[
L(X)=L_c(X)+L_d(X),
\]

then

\[
\frac{B'}B=L,
\]

and

\[
\frac{B''}B=L^2
+\frac{8cX}{(c^2-X^2)^2}
+\frac{8dX}{(d^2-X^2)^2}.
\]

So the two-step transport is still one rational multiplier calculation; no
iteration of harmonic identities is required.

---

## 5. New poles and the boundary differential operator

The old-pole lemma does not apply to the new pole `X=c`.  Here the multiplier
itself creates the pole.  If `H` is analytic at `s`, write `z=X-s`.  Since

\[
A_s(X)=\frac{(2s+z)^2}{z^2}
=\frac{4s^2}{z^2}+\frac{4s}{z}+1,
\]

we obtain

\[
\operatorname{FP}_{X=s}(A_sH)
=H(s)+4sH'(s)+2s^2H''(s).
\]

Define the boundary functional

\[
\boxed{
\beta_s(H):=H(s)+4sH'(s)+2s^2H''(s).}
\tag{8}
\]

Then the three new-pole finite parts appearing over two steps are exactly

\[
\boxed{
C_{n+1,c}=\beta_c(F_n),}
\tag{9}
\]

\[
\boxed{
C_{n+2,c}=\beta_c(A_dF_n),
\qquad
C_{n+2,d}=\beta_d(A_cF_n).}
\tag{10}
\]

These formulas are the complete endpoint transformation.  There is no hidden
summation.

---

## 6. Where the central-binomial square comes from

The double-pole coefficient created at `X=c` is

\[
4c^2F_n(c).
\]

But

\[
R_n(c)
=\frac{(2n+1)!/(n+1)!}{(n+1)!}
=\frac1{n+1}\binom{2n+1}{n}
=\frac1c\binom{2n+1}{n}.
\]

Hence

\[
\boxed{
c^2F_n(c)=\binom{2n+1}{n}^2,}
\tag{11}
\]

and therefore

\[
\boxed{
T_{n+1,n+1}=4\binom{2n+1}{n}^2.}
\tag{12}
\]

So the central-binomial-square scale is not guessed from the recurrence.  It is
forced by the first newly created double pole.  At the next new pole one gets
`4 binom(2n+3,n+1)^2`, and

\[
\binom{2n+3}{n+1}
=\frac{2(2c+1)}{c+1}\binom{2n+1}{n},
\]

so every two-step boundary contribution is again a rational multiple of the
same square in (11).

---

## 7. Quotient-ring representation of the finite-part trace

Set

\[
D_n(X)=\prod_{a=0}^n(X-a),
\qquad
U_n(X)=\prod_{u=1}^n(X+u),
\]

so `R_n=U_n/D_n`, and let

\[
\mathcal A_n:=\mathbf Q[X]/(D_n).
\]

Because the roots `0,...,n` are simple, `D_n'` is a unit in `A_n`; because
`U_n(a) != 0` for every root `a`, `U_n` is also a unit.  All quotients below
are therefore legitimate classes in `A_n`.

At a root `a`, put

\[
r_a(X):=(X-a)R_n(X)=\frac{U_n(X)}{D_n(X)/(X-a)}.
\]

Then

\[
T_{n,a}=r_a(a)^2
=\left(\frac{U_n}{D_n'}\right)^2\!(a).
\]

The logarithmic derivative is

\[
q_{n,a}
=\left(
\frac{U_n'}{U_n}-\frac12\frac{D_n''}{D_n'}
\right)(a).
\]

Define in `A_n`

\[
\mathsf T_n
:=\left(\frac{U_n}{D_n'}\right)^2,
\tag{13}
\]

\[
\mathsf Q_n
:=\frac{U_n'}{U_n}-\frac12\frac{D_n''}{D_n'}.
\tag{14}
\]

Since

\[
C_{n,a}=\frac12(r_a^2)''(a),
\]

a direct logarithmic differentiation gives the class

\[
\boxed{
\begin{aligned}
\mathsf C_n:=\mathsf T_n\bigg[&
\frac{U_n''}{U_n}
+\left(\frac{U_n'}{U_n}\right)^2
-2\frac{U_n'}{U_n}\frac{D_n''}{D_n'}\\
&+\frac34\left(\frac{D_n''}{D_n'}\right)^2
-\frac13\frac{D_n'''}{D_n'}
\bigg].
\end{aligned}}
\tag{15}
\]

Its value at `a` is exactly `C_{n,a}`.  Since

\[
\mathcal A_n\cong\prod_{a=0}^n\mathbf Q,
\]

the algebra trace of multiplication by a class is the sum of its evaluations.
Therefore

\[
\boxed{
e_n=-\operatorname{Tr}_{\mathcal A_n/\mathbf Q}(\mathsf C_n).}
\tag{16}
\]

This is a completely harmonic-free representation of the finite-part trace.

---

## 8. Exact recurrence reduction and the one missing trace identity

Let the Apéry polynomial be

\[
P(t)=34t^3+51t^2+27t+5,
\]

and define

\[
\boxed{
K_n(X):=c^3-P(c)A_c(X)+d^3A_c(X)A_d(X).}
\tag{17}
\]

At every old pole, linearity of (3) gives

\[
\operatorname{FP}_a(K_nF_n)
=K_n(a)C_{n,a}
+2K_n'(a)q_{n,a}T_{n,a}
+\frac12K_n''(a)T_{n,a}.
\]

Hence define the old-pole trace

\[
\boxed{
\mathscr T_n:=
\operatorname{Tr}_{\mathcal A_n/\mathbf Q}
\left(
K_n\mathsf C_n
+2K_n'\mathsf Q_n\mathsf T_n
+\frac12K_n''\mathsf T_n
\right).}
\tag{18}
\]

All denominators of `K_n,K_n',K_n''` are powers of `X-c` and `X-d`, hence are
units modulo `D_n`; (18) is an honest quotient-ring trace.

Splitting the recurrence into old poles and the new poles (9)--(10) gives the
**exact identity**

\[
\boxed{
\begin{aligned}
&c^3e_n-P(c)e_{n+1}+d^3e_{n+2}\\
&\qquad=-\mathscr T_n
+P(c)\beta_c(F_n)
-d^3\bigl[\beta_c(A_dF_n)+\beta_d(A_cF_n)\bigr].
\end{aligned}}
\tag{19}
\]

No recurrence fitting and no WZ step has entered.

Now

\[
287n^2+813n+578=287c^2+239c+52,
\]

and by (11)

\[
\binom{2n+1}{n}^2=c^2F_n(c).
\]

Therefore the desired recurrence is **equivalent** to the single trace identity

\[
\boxed{
\begin{aligned}
\mathscr T_n
={}&P(c)\beta_c(F_n)
-d^3\bigl[\beta_c(A_dF_n)+\beta_d(A_cF_n)\bigr]\\
&+\frac{c^2(287c^2+239c+52)}{d}\,F_n(c).
\end{aligned}}
\tag{TRACE}
\]

This is the exact remaining theorem.  It contains:

- no harmonic numbers;
- no `a`-indexed hypergeometric summand;
- one trace in the semisimple algebra `Q[X]/D_n`;
- three explicit boundary differential evaluations;
- the central square represented simply as `c^2 F_n(c)`.

A proof of `(TRACE)` would finish Q7735 immediately via (19).  Conversely, the
stated recurrence and the already proved transport imply `(TRACE)`, so nothing
is lost in the reduction.

---

## 9. Hermite interpolation reformulation of the same missing trace

There is a canonical finite-part extractor depending only on the pole
polynomial `D_n`.

Let `Pi_n` be the unique polynomial of degree `< 2(n+1)` satisfying, for every
root `a=0,...,n`,

\[
\Pi_n(a)=\frac{D_n''(a)}{2D_n'(a)},
\]

\[
\Pi_n'(a)=
\frac{D_n'''(a)}{3D_n'(a)}
-\frac{D_n''(a)^2}{4D_n'(a)^2}.
\]

Set

\[
W_n(X)=\frac{D_n'(X)}{D_n(X)}-\Pi_n(X).
\]

The two Hermite jets are exactly what is needed to get

\[
\boxed{
W_n(X)=\frac1{X-a}+O((X-a)^2)
\quad(X\to a).}
\tag{20}
\]

Therefore, for every rational `H` having at most a double pole at the roots of
`D_n`,

\[
\boxed{
\operatorname{FP}_{X=a}H
=\operatorname{Res}_{X=a}(H W_n).}
\tag{21}
\]

The Hermite condition is equivalently the polynomial Riccati congruence

\[
\boxed{
D_n''-2\Pi_nD_n'
+D_n(\Pi_n^2-\Pi_n')\equiv0\pmod{D_n^2}.}
\tag{22}
\]

Taking `H=K_nF_n` and using the global residue theorem yields

\[
\boxed{
\mathscr T_n
=-\operatorname{Res}_{X=c}(K_nF_nW_n)
-\operatorname{Res}_{X=d}(K_nF_nW_n)
-\operatorname{Res}_{X=\infty}(K_nF_nW_n).}
\tag{23}
\]

Thus `(TRACE)` can equivalently be attacked by evaluating two endpoint residues
and one infinity residue of a single Hermite kernel.

This does not yet close the proof.  The exact computation shows why: the
Hermite polynomial is not low degree.  In the audit it has degrees

\[
3,5,7,\ldots,17
\]

for `n=1,...,8`, i.e. `2n+1` in every tested nontrivial case.  Consequently the
infinity residue in (23) retains global information rather than collapsing to
a fixed number of leading coefficients.

Likewise the quotient representative of the integrand in (18) has degree
exactly

\[
1,2,3,\ldots,8
\]

for `n=1,...,8`.  These are exact finite computations, not a theorem for all
`n`, but they decisively rule out claiming that the observed instances already
exhibit an `O(1)`-degree trace collapse.

---

## 10. Additional exact checks

The multiplier itself has the compact value at infinity

\[
\boxed{
K_n(\infty)=-4(2c+1)^3.}
\tag{24}
\]

The verifier proves (24) symbolically in independent variables `c,X`.

The successful remote Sage 10.6 audit checked:

- the exact partial-fraction finite-part trace against the displayed harmonic
  definition;
- the `n=1` local-normalization regression `(-8,-5) != (-5,-8)`;
- the quotient class (15) at every pole;
- the Hermite jets (20), residue extraction (21), and congruence (22);
- every old-pole transformation (5) and (7);
- every new-pole boundary formula (8)--(10);
- `T_{n+1,n+1}=4 binom(2n+1,n)^2`;
- the complete recurrence exactly through `n=10` in the main pass;
- quotient/Hermite degree diagnostics through `n=8`.

The output of the successful main pass included

```text
SYMBOLIC_MULTIPLIER PASS
EXACT_N_RANGE 0 10
E_PREFIX [0, -13, -905/4, -167965/36, -15533105/144, -9676971997/3600]
Q7735_SYMBOLIC_VERIFY PASS
```

The degree-audit pass also ended in `Q7735_SYMBOLIC_VERIFY PASS`.

This verification is not presented as a proof of `(TRACE)` for arbitrary `n`;
it is a guard against sign, endpoint, finite-part, and normalization errors in
the reduction.

---

## 11. Sharp stopping point

The requested rational-function strategy succeeds up to one precise theorem:
`(TRACE)`.

What has been removed compared with the original harmonic recurrence is
substantial:

- the local harmonic finite-part formula is gone;
- the two `n`-shifts are handled by derivatives of one rational multiplier;
- all endpoint terms are explicit evaluations of `F_n` and its first two
  derivatives;
- the central-binomial-square factor is forced by the new pole;
- the remaining global object is one algebra trace of an explicit rational
  class.

What has **not** been removed is the genuinely global trace evaluation.  The
quotient representative continues to fill the whole degree-`n` algebra in the
exact tests, and the Hermite extractor has linear degree.  A proof of `(TRACE)`
therefore still needs a nontrivial global identity — perhaps a structured
falling-factorial trace evaluation or a new relation for the Hermite kernel —
but it no longer needs a giant harmonic WZ certificate.
