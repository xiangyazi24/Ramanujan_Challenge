# Intermediate-range cube arithmetic: canonical determinant audit and scoped saturation theorem

## Verdict

I pinned this audit to repository `main` at `734a5a84c1e4fd8703a811aadaa2b4c7f532b20e` and did not modify shared TeX.

I do **not** obtain a target-selective nonzero integer of height controlled only by the cube span.  There is a very natural canonical candidate—the eight-node Vandermonde/barycentric determinant—and its exact factorization explains the obstruction:

\[
\det\bigl(1,p_i,p_i^2,\ldots,p_i^6,BV_i\bigr)_{i=1}^8
=-8B\prod_{i<j}(p_j-p_i),
\]

where
\[
V_i=\prod_{j\ne i}(p_j-p_i),\qquad B=\frac{b_m}{\prod_i p_i}.
\]

The entire span-controlled factor is the Vandermonde.  The remaining factor is exactly the uncontrolled global quotient `B`.  Thus ordinary interpolation does not manufacture a small target scalar; it isolates the free quotient.

More generally, after localizing at the 56 directed cross units and the 28 short-transfer units, **finite-depth local lifts and finite short transfers saturate componentwise under CRT**.  At the `p_i` component there is only one defining zero, at `h_i`; the other seven target rows are units.  A first divided digit `z_i=b_{h_i}/p_i mod p_i` is a free Hensel coordinate, and a triangular `p_i^2` lift solves for that digit rather than imposing a relation among different primes.  Since

\[
\mathbf Z/(\prod_i p_i^k)\mathbf Z\simeq\prod_i\mathbf Z/p_i^k\mathbf Z,
\]

a single integral boundary state or a single integer `B` carries arbitrary independent local residues unless one imports an additional global height/arithmetic restriction.  Consequently no boundary determinant/resultant obtained solely from the stated local package contracts, after unit saturation and elimination, to a new nonzero equation in `(m,a,u,d,e,f)`.

This is a **scoped** no-go.  It does not rule out a genuinely global Apéry identity, modular/Hecke relation, moving-index geometric realization, or a cross-prime theorem that couples the first divided digits before CRT.

The weighted aggregate is more plausible than a pointwise cell bound.  If the number of active cells is

\[
C_X\ll \frac{X\log^2X}{L}
\]

and the forced aggregate cube count is

\[
F_X\asymp \frac{L^{15}}{X^{14}\log^{14}X},
\]

then a pointwise argument would need the very strong scale

\[
C_{\rm cube}(\mathfrak c)=o\!\left(
\frac{L^{16}}{X^{15}\log^{16}X}
\right)
\]

uniformly in the cell.  By contrast the second-moment target

\[
\sum_{\mathfrak c} C_{\rm cube}(\mathfrak c)^2
=o\!\left(
\frac{L^{31}}{X^{29}\log^{30}X}
\right)
\]

already implies `sum C_cube=o(F_X)` by Cauchy.  Even weaker in correlation order, if all cube sides are at most `H_*`, a cube has a canonical square face and at most `O(H_*)` choices for its third side, so the first genuinely new correlation theorem can be a **four-point square aggregate**

\[
\sum_{\mathfrak c}Q_{\mathfrak c}(H_*)
=o\!\left(
\frac{L^{15}}{X^{14}H_*\log^{14}X}
\right).
\]

That square theorem is the smallest remaining target I recommend: pairs do not see the commuting face and are explicitly neutralized here by the directed cross-unit masks.

---

## 1. Exact fixed-cell geometry, both orientations

Let `eps in {0,1}` denote orientation, with `eps=0` direct and `eps=1` reflected.  Put

\[
c=a+\varepsilon,\qquad \tau=2\varepsilon-1,\qquad M=m+\varepsilon.
\]

Then both formulas in the question are

\[
\boxed{p(h)=\frac{M+\tau h}{c}.}
\tag{1.1}
\]

Indeed:

* direct: `c=a`, `tau=-1`, `M=m`, so `p=(m-h)/a`;
* reflected: `c=a+1`, `tau=+1`, `M=m+1`, so `p=(m+h+1)/(a+1)`.

Write the Boolean vertices as

\[
h_S=u+\sum_{j\in S}\delta_j,
\qquad (\delta_1,\delta_2,\delta_3)=(d,e,f),
\qquad S\subseteq\{1,2,3\}.
\]

If all eight vertices lie in one quotient/orientation cell, every numerator in (1.1) is divisible by `c`.  Taking differences shows

\[
\boxed{c\mid d,\quad c\mid e,\quad c\mid f.}
\tag{1.2}
\]

Set

\[
\Delta_j=\delta_j/c,
\qquad P=(M+\tau u)/c.
\]

Then the eight defining primes themselves form an affine Boolean cube:

\[
\boxed{p_S=P+\tau\sum_{j\in S}\Delta_j.}
\tag{1.3}
\]

In particular

\[
|p_T-p_S|=\frac{|h_T-h_S|}{c}
\le \frac{d+e+f}{c}.
\tag{1.4}
\]

This identity is purely quotient geometry.  It is present before any Apéry zero is imposed, so a determinant that vanishes only because the eight `p_S` form an affine cube is not target-selective.

The reflected chart introduces no new algebraic issue: it changes `tau` from `-1` to `+1` and `c` from `a` to `a+1`.  The `O(1)` reflection exceptions from the reduction must simply be discarded before applying the unit-local argument below.

---

## 2. What the local `p^2` lift does and does not provide

The repository contains the exact direct quotient-one lift

\[
b_{p+r}\equiv5b_r+10pD_r\pmod{p^2},
\tag{2.1}
\]

with the harmonic-binomial correction `D_r`.  The audit below does **not** silently promote (2.1) to a formula for every ordinary quotient `a`.  For a general fixed `(a,eps)` I only use a local lift if it has actually been proved in that cell, and I write its first divided equation abstractly as

\[
\boxed{B V_i\equiv \beta_i z_i+\Gamma_i\pmod{p_i},}
\tag{2.2}
\]

where

\[
R=\prod_{i=1}^8p_i,
\quad B=b_m/R,
\quad z_i=b_{h_i}/p_i,
\quad V_i=\prod_{j\ne i}(p_j-p_i),
\tag{2.3}
\]

and `beta_i` is a unit on the legal nonsingular branch.  In the direct quotient-one case, `beta_i=5` and `Gamma_i=10D_{h_i}`.

The key congruence behind (2.2) is independent of the Apéry lift.  Since

\[
\frac{R}{p_i}=\prod_{j\ne i}p_j
\equiv\prod_{j\ne i}(p_j-p_i)=V_i\pmod{p_i},
\tag{2.4}
\]

we have

\[
\frac{b_m}{p_i}=B\frac{R}{p_i}\equiv BV_i\pmod{p_i}.
\]

Because the defining primes are distinct, `V_i` is a `p_i`-unit.  Equation (2.2) therefore has the opposite effect from the desired elimination: if `beta_i` is a unit, it **solves uniquely for the free divided digit**

\[
z_i\equiv \beta_i^{-1}(BV_i-\Gamma_i)\pmod{p_i}.
\tag{2.5}
\]

It does not constrain `B`.

This also explains why the smallness of `V_i` is not enough.  One cannot discard `z_i`, `Gamma_i`, or `B`; they are part of the literal lifted equation.

---

## 3. The canonical Vandermonde boundary and its exact failure

Let

\[
F(T)=\prod_{i=1}^8(T-p_i).
\]

Then

\[
F'(p_i)=\prod_{j\ne i}(p_i-p_j)=-V_i,
\tag{3.1}
\]

because there are seven sign reversals.  For arbitrary values `y_i`, the coefficient of `T^7` in the unique degree-at-most-seven interpolant through `(p_i,y_i)` is

\[
[T^7]Q(T)=\sum_{i=1}^8\frac{y_i}{F'(p_i)}.
\tag{3.2}
\]

Taking `y_i=BV_i` gives the exact identity

\[
\boxed{[T^7]Q=-8B.}
\tag{3.3}
\]

Equivalently, with the row order used to define

\[
\Delta(p)=\prod_{i<j}(p_j-p_i),
\]

we obtain

\[
\boxed{
\det
\begin{pmatrix}
1&p_1&\cdots&p_1^6&BV_1\\
\vdots&\vdots&&\vdots&\vdots\\
1&p_8&\cdots&p_8^6&BV_8
\end{pmatrix}
=-8B\Delta(p).
}
\tag{3.4}
\]

Let `H=d+e+f`.  From (1.4),

\[
|V_i|\le(H/c)^7,
\qquad
|\Delta(p)|\le(H/c)^{28}.
\tag{3.5}
\]

So (3.4) achieves the desired span height **only for the geometric factor**.  The target information leaves behind `B`, whose Archimedean size is not controlled by `H`.

The same calculation kills ordinary face interpolation in a more conceptual way.  Eight arbitrary labels on a Boolean 3-cube have an eight-dimensional multilinear interpolation space.  Unless one proves a new rank/degree drop, interpolation has exactly enough freedom to absorb all eight target values.  The local `p^2` equations do not give such a rank drop.

Finally, one cannot insert the eight congruences (2.2) into a single determinant as if they held in one field.  Row `i` is valid modulo `p_i`, not modulo every `p_j`.  Writing an integral lift

\[
y_i=BV_i+p_it_i
\tag{3.6}
\]

introduces eight independent integers `t_i`; the resulting determinant acquires the corresponding free cofactor terms.  Those are exactly the CRT degrees of freedom that the informal cross-prime determinant suppresses.

---

## 4. Scoped saturation theorem

The following is the precise no-go I can justify from the current package.

### Theorem (finite local/transfer/CRT saturation)

Fix one nonexceptional cell `(m,a,eps)` and one nondegenerate eight-vertex cube.  Assume:

1. the eight defining primes `p_i` are distinct;
2. at the `p_i` component, `b_{h_i}=0 mod p_i` and all seven cross values `b_{h_j}`, `j!=i`, are units (the 56 directed cross-unit masks);
3. every short transfer used inside the cube is invertible at the relevant component (the 28 transfer-unit masks);
4. every extra finite-depth lift used in the argument is **Hensel-triangular**: after the lower digits have been fixed, each new scalar equation contains a fresh divided local digit with unit coefficient.  The proved first divided equation (2.2) is of this form;
5. different primes are coupled only through shared integral boundary variables, `B`, and ordinary algebraic interpolation/determinants.  No extra cross-prime modular/Hecke identity or Archimedean height restriction is imported.

Let `G_geom` be the ideal generated by the quotient geometry (1.1)-(1.3).  Form the polynomial/congruence system consisting of the finite transfer equations, the eight target zeros, all finite triangular lift equations, and the shared boundary/quotient variables.  Localize at the product of all cross-unit, transfer-unit, and triangular leading coefficients.  Then elimination of the boundary variables, divided digits, finite correction states, lift multipliers, and `B` contracts to **no new target equation**:

\[
\boxed{I_{\rm local}^{\rm sat}\cap
\mathbf Z[m,a,u,d,e,f]=G_{\rm geom}}
\tag{4.1}
\]

within this formal local model.

In particular, this package alone cannot produce a nonzero target-selective scalar depending only on `(a,eps,d,e,f)` whose height is controlled by `d+e+f`.

### Proof

Work first at one component `p_i`.  A short transfer is an invertible linear map on the two-dimensional Apéry state.  The condition `b_{h_i}=0 mod p_i` is one nonzero linear functional on that state.  By the cross-unit hypotheses, the seven other evaluation functionals are nonzero on the actual kernel line and hence remain units after localization.

Lift the boundary state from modulo `p_i` to modulo `p_i^2`.  Replacing it by

\[
s_i\longmapsto s_i+p_i t_i
\]

does not change any mod-`p_i` target or unit mask.  Because the target evaluation functional is nonzero, its value on `t_i` can be prescribed arbitrarily; therefore the divided digit

\[
z_i=b_{h_i}/p_i\pmod{p_i}
\]

is a free local coordinate.  Equation (2.2), or any triangular analogue, then determines that fresh digit for every prescribed residue of `B`.  Inductively the same statement holds for any fixed number of triangular lift levels.

Now take all eight primes.  For every fixed `k`, CRT gives

\[
\mathbf Z/R^k\mathbf Z
\simeq\prod_{i=1}^8\mathbf Z/p_i^k\mathbf Z.
\tag{4.2}
\]

Thus a *single* integral boundary vector modulo `R^k` is equivalent to eight arbitrary local boundary vectors, and a *single* integer `B mod R^k` is equivalent to eight arbitrary local residues.  The common-source notation does not remove a degree of freedom at the congruence level.

Therefore every local component can be solved independently for every geometric point in the localized open set.  Existential elimination can only recover identities already imposed on the geometry.  Ordinary interpolation does not change this conclusion: full-degree interpolation is an isomorphism between the eight nodal values and eight coefficients.  This proves (4.1).

### Scope

The theorem is intentionally not universal.  It does **not** cover a new identity that couples two defining characteristics before CRT, a modular/Hecke relation, a global integral Fitting object with an independently proved small height, or a new geometric realization of the moving index.  It also does not say that fixing the true Apéry boundary state at index zero leaves it arbitrary.  Rather: **the finite local package has no access to the global restriction that selects that boundary state.**  Importing the true boundary data is a genuinely global arithmetic input, and its natural integer carriers need not have short-span height.

---

## 5. Why short transfer and pair resultants do not escape the theorem

The repository's exact gap kernel has transfer degree `3(G-1)` for a gap `G`, and fixed-pair resultants have logarithmic height growing like the product of the two gap lengths.  That machinery is useful when the **same prime** divides two endpoint values.

Here the cube masks are deliberately the opposite.  At component `p_i`,

\[
p_i\mid b_{h_i},
\qquad
p_i\nmid b_{h_j}\quad(j\ne i).
\tag{5.1}
\]

Hence no directed pair inside the cube supplies two zeros for one prime.  A gap resultant `Res(N_g(...),...)` therefore receives no forced `p_i` divisor from the cube.  The 28 transfer-unit masks confirm that the short words stay on the invertible branch; they are open conditions, not additional vanishing equations.

A Fitting/determinant construction can still be new **only if it uses many rows before localizing and proves an integral small-height minor that is not rendered surjective by CRT**.  That is exactly the kind of genuinely global construction not covered by the theorem.  No such minor follows from the current local transfer identities alone.

---

## 6. Direct versus reflected chart

Nothing in the saturation proof depends on the sign of the affine prime cube.  The two charts are:

\[
\begin{array}{c|c|c|c}
&c&\tau&p_S\\ \hline
\text{direct}&a&-1&P-\sum_{j\in S}\delta_j/a\\
\text{reflected}&a+1&+1&P+\sum_{j\in S}\delta_j/(a+1).
\end{array}
\]

The only chart-specific arithmetic is the precise local lift replacing (2.2).  I use the explicit `5,10D` formula only where it is proved (direct quotient one).  For a general fixed ordinary quotient, the no-go needs only the weaker structural fact that any available first-divided lift is triangular with a unit coefficient; if that fact is not proved, then that lift cannot be used at all, and the local package is weaker still.

The finite set of reflection exceptions from the cube reduction should be excised before applying the theorem; absorbing them into the final `O(1)` loss is legitimate and does not alter the aggregate exponents.

---

## 7. Pointwise versus weighted aggregate

Let `C_cube(c)` denote the number of retained actual cubes in cell `c`.  The reduction supplies

\[
\#\{\text{active cells}\}\le C_X
\ll \frac{X\log^2X}{L}
\tag{7.1}
\]

and forces

\[
\sum_c C_{\rm cube}(c)\ge c_0F_X,
\qquad
F_X=\frac{L^{15}}{X^{14}\log^{14}X}.
\tag{7.2}
\]

To contradict (7.2) by a uniform pointwise estimate one would need

\[
C_{\rm cube}(c)=o(F_X/C_X)
=o\!\left(\frac{L^{16}}{X^{15}\log^{16}X}\right).
\tag{7.3}
\]

That is much stronger than necessary.  Cauchy gives

\[
\sum_c C_{\rm cube}(c)
\le C_X^{1/2}
\left(\sum_c C_{\rm cube}(c)^2\right)^{1/2}.
\]

Hence it is enough to prove

\[
\boxed{
\sum_c C_{\rm cube}(c)^2
=o\!\left(\frac{F_X^2}{C_X}\right)
=o\!\left(
\frac{L^{31}}{X^{29}\log^{30}X}
\right).
}
\tag{7.4}
\]

This route can exploit cancellation/rarity across cells and is not forced to defeat a single pathological cell.

There is also a lower-order projection.  Suppose the exact cube reduction gives `d,e,f<=H_*` (or only `d+e+f<=H_*`; the same crude bound works).  Order the side lengths canonically and project each cube to its first two directions.  For a fixed square there are at most `O(H_*)` integral choices for the third side (in fact at most `H_*/c` inside one quotient cell).  Therefore

\[
\sum_c C_{\rm cube}(c)
\ll H_*\sum_c Q_c(H_*),
\tag{7.5}
\]

where `Q_c(H_*)` counts actual retained four-vertex square faces with the literal quotient/folded masks.  A sufficient theorem is

\[
\boxed{
\sum_c Q_c(H_*)
=o\!\left(
\frac{L^{15}}{X^{14}H_*\log^{14}X}
\right).
}
\tag{IC4}
\]

`IC4` is the weakest new correlation I would attack next.  A two-point estimate could also imply the cube bound after paying `O(H_*^2)` completions, but the current directed cross-unit structure means pairs see no same-prime double zero and hence no new local arithmetic.  A square is the first face on which additive commutation can support a genuinely cross-prime dispersion theorem.

---

## 8. What would actually escape the saturation theorem

Any successful arithmetic closure must add at least one of the following kinds of input.

1. **Cross-prime divided-lift rigidity.**  A relation involving the first divided digits/corrections at two or more distinct defining primes that is not equivalent, under CRT, to separate local equations.  It must remain nontrivial after the unit variables are eliminated.
2. **A global small-height Fitting minor.**  Construct the many-row banded system first, then prove a primitive minor is nonzero and has logarithmic height `o(log R)` (or whatever the exact forced-prime product demands) after all legitimate small factors are saturated.  Merely observing small Vandermonde coefficients is insufficient if `B` or divided digits remain.
3. **A weighted defining-characteristic dispersion theorem.**  Prove `IC4`, (7.4), or a stronger mixed-cell estimate directly for the actual Apéry rows.  This is a cross-prime theorem, not a pointwise Deligne statement at an auxiliary characteristic.
4. **A global modular/Hecke or moving-index geometric identity.**  Such an identity can evade the direct-transfer degree obstruction, but it must control divisibility in the defining characteristic, not merely give an auxiliary-ell trace bound.

The existing fixed-rank differential compression of the gap kernel and the affine codimension-two ambient model make such routes plausible targets, but neither supplies the defining-characteristic mixed-cell estimate by itself.

---

## 9. Exact status classification

### Proved in this audit

* the unified direct/reflected affine-prime-cube formula (1.1)-(1.4);
* the barycentric congruence `R/p_i = V_i mod p_i`;
* the exact canonical interpolation determinant (3.4), with span factor at most `(H/c)^28` and surviving quotient `B`;
* the componentwise Hensel/CRT saturation theorem under the explicitly stated finite-triangular hypotheses;
* failure of same-prime pair resultants on a cube satisfying all directed cross-unit masks;
* exact aggregate exponent bookkeeping (7.3)-(7.4);
* the sufficient four-point square target `IC4`.

### Not proved

* a general-`a` `p^2` lift beyond formulas already banked for that quotient/orientation;
* a universal no-go against modular, Hecke, toric/Kummer, or genuinely global boundary constructions;
* `IC4` or the second-moment estimate (7.4);
* an asymptotic upper bound for actual cube counts strong enough to contradict the forced aggregate.

Accordingly this is a **negative local-arithmetic result plus a sharpened remaining target**, not a proof of the P3.2 gateway.

---

## 10. Verifier

The owned standard-library verifier is

`problems/3.2/ORACLE_COMM/chatgpt_q32_intermediate_cube_arithmetic_verify.py`.

It checks, using exact integer/rational arithmetic only:

* the unified direct/reflected quotient formula and the divisibility of side lengths by `a+eps`;
* an explicit eight-prime affine cube and all barycentric weights;
* `R/p_i = V_i mod p_i` for every node;
* the interpolation leading coefficient `-8B` and the determinant identity `det=-8B*Vandermonde` by Bareiss elimination;
* full eight-degree Boolean interpolation (no hidden face constraint at ordinary interpolation degree);
* CRT freedom for an arbitrary tuple of local `B_i` and arbitrary local boundary states modulo `p_i^2`;
* solvability of a model triangular first-divided lift for every prescribed local `B_i`;
* the pointwise and second-moment exponent bookkeeping.

The synthetic residues in that verifier are **only algebraic witnesses for CRT/interpolation surjectivity**.  They are not presented as Apéry rows or as arithmetic counterexamples.

Per the delivery constraint for this task, I did not invoke Python/Sage or a sandbox to execute the verifier in this turn.