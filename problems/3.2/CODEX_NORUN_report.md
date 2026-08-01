# [NO-RUN] terminal report: the proposed certificate targets the wrong orbit

## Verdict

**REFUTED-AS-STATED.**  The claimed companion action

\[
\xi_{n+1}=M_n\xi_n,\qquad \xi_n=(b_n:c_n),
\]

is false.  Consequently the advertised collapse

\[
r,r+1\in Z_d,\ c_r\ne0
\quad\Longrightarrow\quad
\xi_r=(\alpha_r-\alpha_{r+d}:\beta_r-\beta_{r+d})
\]

is also false on the live Apéry orbit.  The machine certificate gives an
exact rational counterexample to the action and an exact finite-field
counterexample to the collapse.  The requested uniform no-three-run theorem,
its mixed-branch constants, and its Strike2 corollary therefore cannot be
reported as proved.

This is a structural failure, not a failed factorization or a numerical
anomaly.  Following the specification's stop rule, the verifier terminates
with `FINAL GATE: FAIL` rather than enumerating branches of the wrong
dynamical system.

## 1. The two projective objects that were conflated

Put

\[
u_n=(b_n,c_n).
\]

Because both coordinates solve the scalar Apéry recurrence,

\[
u_{n+1}=\alpha_nu_n-\beta_nu_{n-1}.
\]

Thus the companion matrix has the valid matrix identity

\[
\begin{pmatrix}u_n\\u_{n+1}\end{pmatrix}
=M_n\begin{pmatrix}u_{n-1}\\u_n\end{pmatrix},
\qquad
M_n=\begin{pmatrix}0&1\\-\beta_n&\alpha_n\end{pmatrix}.
\]

It propagates the two-time state `(u_(n-1),u_n)` of **one scalar
solution**.  It does not propagate the row `u_n=(b_n,c_n)` formed from
**two different solutions**.  The collision set in the specification is

\[
Z_d=\{r:u_{r+d}\parallel u_r\},
\]

so it belongs to the latter, row-projective orbit.

The mismatch is visible over the rationals at the first step.  With

\[
(b_1,c_1)=(5,6),\quad (b_2,c_2)=(73,351/4),\quad
(\alpha_1,\beta_1)=(117/8,1/8),
\]

one has

\[
M_1(b_1,c_1)^T=(6,697/8)^T,
\]

and

\[
\det\!\begin{pmatrix}6&73\\697/8&351/4\end{pmatrix}
=-46669/8\ne0.
\]

Hence even projective equality fails.

## 2. Exact live counterexample to the collapse

Take

\[
p=997,\qquad d=182,\qquad r=248.
\]

All indices lie in the nonwrapping window `1,...,p-2`.  Direct recurrence
evaluation with `c_0=0,c_1=6` gives

| index | `(b_n,c_n)` | affine key `b_n/c_n` |
|---:|:---|---:|
| 248 | `(717,994)` | 758 |
| 430 | `(588,106)` | 758 |
| 249 | `(190,566)` | 409 |
| 431 | `(763,153)` | 409 |

Therefore `248,249 in Z_182`, and `c_248=994 != 0`.  But

\[
\alpha_{248}-\alpha_{430}=384,\qquad
\beta_{248}-\beta_{430}=960\pmod {997},
\]

so the predicted point has affine key

\[
384/960=798\pmod {997},
\]

whereas `xi_248` has key `758`.

There is a second symmetric witness at `r=565`:

\[
\xi_{565}=409,qquad
(\alpha_{565}-\alpha_{747})/(\beta_{565}-\beta_{747})=165
\pmod {997}.
\]

The verifier checks both collision equalities at each start before checking
the failed prediction.  Thus this does not depend on the gap-polynomial
implementation or on a normalization of `N_d`.  Replacing `c_1=6` by
`c_1=1` only rescales the second solution and leaves every collision set
unchanged; it does not repair the companion action.

## 3. What the symbolic scratch computation actually proves

The advertised polynomial is algebraically correct for a projective point
`x` that really does obey `x_(n+1)=M_n x_n`.  Set

\[
A=\alpha_r-\alpha_{r+d},\qquad
B=\beta_r-\beta_{r+d},
\]

and similarly `A_1,B_1` with `r` replaced by `r+1`.  In the chart `B != 0`,

\[
\det(M_rx,M_{r+d}x)=0
\]

forces `x=(A:B)`.  The next determinant is

\[
B B_1-(-\beta_rA+\alpha_rB)A_1.
\]

After cancellation its denominator is

\[
(r+1)^3(r+2)^3(r+d+1)^3(r+d+2)^3,
\]

and its numerator is

\[
F(r,d)=-24d^2G(r,d),
\]

where

\[
\begin{aligned}
G={}&108r^8+864r^7+2763r^6+4482r^5+3849r^4+1644r^3+276r^2\\
&+d(432r^7+3024r^6+8289r^5+11205r^4+7698r^3+2466r^2+276r)\\
&+d^2(648r^6+3888r^5+8826r^4+9384r^3+4719r^2+1038r+92)\\
&+d^3(432r^5+2160r^4+3837r^3+2871r^2+870r+108)\\
&+d^4(108r^4+432r^3+537r^2+210r+32).
\end{aligned}
\]

Thus `deg_r F=8`, `deg_d F=6`, and the total degree after removing the
denominator is 10.  The verifier records and checks all 32 nonzero
coefficients.  This is a valid symbolic certificate for the companion-state
problem, but it supplies no condition on the collision sets `Z_d` defined
from `(b_n:c_n)`.

For the correct row orbit, if a three-run gives scalars

\[
u_{r+d+j}=\lambda_j u_{r+j}\qquad(j=0,1,2),
\]

comparison of the two endpoint recurrences yields

\[
\lambda_2\alpha_{r+1}=\lambda_1\alpha_{r+d+1},\qquad
\lambda_2\beta_{r+1}=\lambda_0\beta_{r+d+1}.
\]

The missing ratio `lambda_1/lambda_0` is orbit-dependent.  The false
companion action is exactly what removed this datum in the proposed proof.

## 4. Independent claims that do verify

### [GAP-ONE-CASORATIAN] — PROVED

Let

\[
W_r=\det(u_r,u_{r+1}).
\]

The recurrence gives

\[
W_r=\beta_rW_{r-1},\qquad W_0=6,
\]

and hence the exact telescoping identity

\[
W_r=6\prod_{j=1}^r\left(\frac j{j+1}\right)^3
=\frac6{(r+1)^3}.
\]

For `p>3` and `0<=r<=p-2` this is nonzero.  In particular two consecutive
`c`-values cannot both vanish.  The verifier checks this identity at every
legal index modulo 997 as well as symbolically through the recurrence.

### [COMPANION-DEGENERACY] — PROVED, but not applicable to `Z_d`

The equality `beta_r=beta_(r+d)` is equivalent to

\[
[r(r+d+1)]^3=[(r+d)(r+1)]^3.
\]

Writing

\[
r(r+d+1)=z(r+d)(r+1),\qquad z^3=1,
\]

the branch `z=1` reduces exactly to `d=0`, so it is empty for
`1<=d<p`.  Each nontrivial cube root gives a quadratic equation in `r`.
There are at most two such roots `z`, hence at most four `r` values per
fixed `d`.  Independently, the cleared numerator of
`alpha_r-alpha_(r+d)` has degree four in `r` and leading coefficient
`-51d`; for `p>17` it is nonzero.  Therefore the simultaneous companion
degeneracy locus has the explicit bound

\[
\#\{r:\alpha_r=\alpha_{r+d},\ \beta_r=\beta_{r+d}\}\le4.
\]

This bound is correct but does not fix the orbit mismatch.

## 5. Branch enumeration and numerical survey

The generic and mixed `c_(r+j)=0` branches requested in the specification
all use the implication `xi_(n+1)=M_n xi_n`.  Enumerating them after the
action gate fails would enumerate a different state-space problem.  In
particular, the mandated numerical check of every live two-run against
`xi_r=v(r,d)` already fails at the two witnesses above.  The verifier
therefore marks both the mixed-branch enumeration and the forty-prime survey
as `SKIPPED` and exits with status 1.

The earlier observation that the four tested primes had no three-runs is
compatible with this report.  It remains empirical evidence; the displayed
degree-eight polynomial does not prove it for the actual orbit.

## 6. Strike2 consequence

**NOT ESTABLISHED.**  If a genuine Apéry theorem bounded the total number of
three-run starts over `d<=D` by `cD`, it would kill precisely the interval
clause of the primitive-padding construction in `CODEX_STRIKE2_report.md`,
Section 6.2, equations (6.8)--(6.9): in the block for `d_j`, pairing the
first `d_j` positions with the last `d_j` positions makes `Z_(d_j)` contain
an entire interval of length `d_j` (and a reflected copy).  Those blocks
contribute at least `d_j-2` to

\[
\sum_d(\operatorname{maxrun}(Z_d)-2)_+,
\]

and summing over `t_q` selected gaps is of order `t_qD_q`, not `D_q`.

That conditional implication is sound, but the needed Apéry premise was not
proved here.  Even such a corrected bounded-run theorem would not exclude
scattered padding in which the gaps vary from position to position and no
single `Z_d` contains a long interval.

## Inventory

| Item | Status | Meaning |
|:---|:---|:---|
| Exact `F=-24d^2G` certificate | **PROVED-SYMBOLIC** | Valid only for a companion-state orbit |
| Gap-one Casoratian `6/(r+1)^3` | **PROVED-all-r** | Valid for the two-solution Apéry rows |
| Companion degeneracy bound `<=4` | **PROVED for p>17** | Correct but irrelevant to the failed action |
| Collapse `xi_r=v(r,d)` | **REFUTED** | Counterexample `p=997,d=182,r=248` |
| Uniform three-run bound | **NOT PROVED** | Its proposed implication uses the refuted collapse |
| Strike2 interval-padding exclusion | **NOT PROVED** | Would follow only from a corrected run theorem |

Reproduction:

```text
python3 CODEX_NORUN_verify.py
...
FINAL GATE: FAIL -- CODEX_SPEC_laststand_norun.md conflates two projective orbits
```
