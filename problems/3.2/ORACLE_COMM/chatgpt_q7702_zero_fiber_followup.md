# Q7702 — zero-fiber endpoint avoidance after Q7696

## Verdict

I do **not** have a proof of

\[
\#\{0<r<p:b_r=\kappa_r=0\}\le 1
\]

or even a new unconditional `O(1)` multiplicity bound, and I found no actual-source counterexample.  The sharp conjecture remains compatible with the exact scan through `p<=20000` from Q7696.

What can be proved now is a substantially sharper obstruction theorem about every natural "second global constraint" I could make precise.

1. The normalized second homogeneous Apéry solution is **also Hasse-palindromic**.  Thus reflection acts trivially on the entire two-dimensional homogeneous solution space, not just on `b`.
2. Nevertheless a two-return monodromy/adjoint determinant does **not** give two endpoint equations.  The space of homogeneous Dirichlet test solutions for a return pair is one-dimensional, so every such determinant is a scalar multiple of the already-known Green/Duhamel period.
3. The second variation coordinate `Phi` is genuinely independent, but two common endpoints do not constrain it.  At the actual-source example `p=41`, the `b`-weighted Green period from `10` to `30` is zero while the `u`-weighted period is `6 != 0`.
4. The canonical truncated Eichler companion has a genuine global Cartier defect:
   \[
   \mathcal C_p(LK_{<p})=-\Xi_{p-1}t.
   \]
   This is a real extension-class residue, not a restatement of an internal Green period.
5. Even the **full two-coordinate terminal variation vector** `(Xi_{p-1},Phi_{p-1})` cannot be the missing invariant by itself.  For any two Hasse zeros `r<s`, formal source changes at `r,s` can force both rows common, and changes at the last two rows `p-2,p-1` can restore both terminal periods exactly.  The response determinant is the proved unit `25/8`.
6. The actual source does satisfy the rational reciprocity
   \[
   R_p(1/t)=t^{p-1}R_p(t),\qquad
   R_p=\Delta^{(p-1)/2}/A_p^2,
   \]
   but this is reciprocity between the expansions at `0` and `infinity`, not a reflection law for the finite coefficient vector.  The obvious induced coefficient and affine-`Xi` partner laws are false; `p=181` already gives two Hasse reflection pairs with different `Xi`-pair sums.

Therefore the remaining theorem is now sharply localized: it must use **coefficient-level rigidity of the actual level-six Eichler source**, not merely the Picard--Fuchs operator, homogeneous Hasse reciprocity, a two-return transfer determinant, or finitely many terminal Frobenius/period coordinates.

The exact missing implication is still the zero-normalized source-specific endpoint avoidance

\[
\boxed{
 b_r=b_s=0,\quad \Xi_r=0,\quad 0<r<s<p
 \Longrightarrow
 \sum_{m=r+1}^s b_{m-1}[t^m]
 \frac{\Delta(t)^{(p-1)/2}}{A_p(t)^2}\ne0.
}
\tag{EIS-2RET}
\]

No argument below uses the disproved reflection-depth law at `(2237,492)`.

---

## 1. Unit ledger and the endpoint `p=5`

For `p=5`,
\[
\Xi_n=-1-5\sum_{m=1}^n g_mb_{m-1}=-1\ne0\pmod5,
\]
so the common set is empty.  Hence all divisions by `5` below are made only for `p>=7`.

For `p>=7` and `0<n<p`, `n`, `n^3`, and every product of such indices are units.  If `b_n=0`, then `b_{n-1}` and `b_{n+1}` are units: two consecutive zero states would propagate backward to `b_0=0`, impossible.  Hasse reciprocity gives
\[
b_{p-1}=b_0=1,\qquad b_{p-2}=b_1=5,
\]
so every Hasse zero lies at most at `p-3`.

Let `u` be the normalized second homogeneous solution
\[
u_0=0,\qquad u_1=1.
\]
Its Casoratian with `b` is
\[
\boxed{
 n^3(b_{n-1}u_n-b_nu_{n-1})=1
}
\tag{1.1}
\]
for `1<=n<p`.  This follows at `n=1`, and the two homogeneous recurrences propagate the same value from `n-1` to `n`.

At a Hasse zero `b_r=0`, (1.1) becomes
\[
r^3b_{r-1}u_r=1,
\]
so both `b_{r-1}` and `u_r` are units.

---

## 2. New theorem: full homogeneous Hasse reciprocity

### Theorem 2.1

For every prime `p>=5` and `0<=r<=p-1`,
\[
\boxed{u_{p-1-r}=u_r.}
\tag{2.1}
\]
Consequently reflection `J(y)_r=y_{p-1-r}` is the identity on the full two-dimensional homogeneous Apéry solution space over `F_p`.

### Proof

Use
\[
P(-x-1)=-P(x).
\tag{2.2}
\]
If `y` is homogeneous and `z_r=y_{p-1-r}`, put `m=p-1-r`.  The forward recurrence for `y` at `m` is
\[
(m+1)^3y_{m+1}-P(m)y_m+m^3y_{m-1}=0.
\]
Modulo `p`, `m+1=-r`, `m=-(r+1)`, and (2.2) turns the negative of this equation into
\[
(r+1)^3z_{r+1}-P(r)z_r+r^3z_{r-1}=0.
\]
Thus `J` preserves the homogeneous solution space.

The already-proved Hasse reciprocity gives `Jb=b`.  Since `(b,u)` is a basis by (1.1), write
\[
Ju=A b+B u.
\]
At rows `0,1`,
\[
A=u_{p-1},\qquad B=u_{p-2}-5u_{p-1}.
\]
Apply (1.1) at `n=p-1`.  Using `b_{p-2}=5`, `b_{p-1}=1`, and `(p-1)^3=-1`,
\[
-\bigl(5u_{p-1}-u_{p-2}\bigr)=1,
\]
so `B=1`.  Hence `Ju=A b+u`.  Applying `J` again and using `J^2=1`,
\[
u=J^2u=2A b+u.
\]
Because `p` is odd and `b` is nonzero, `A=0`.  Therefore `Ju=u`. ∎

This is stronger than the reflection of `b`, but it does **not** imply a reflected-depth law for common rows.

---

## 3. The exact rank-two variation frame

Define
\[
\Phi_n=5\sum_{m=1}^n g_m u_{m-1},\qquad \Phi_0=0.
\tag{3.1}
\]
Together with
\[
\Xi_n=-1-5\sum_{m=1}^n g_mb_{m-1},
\]
the canonical companion satisfies
\[
\boxed{
\kappa_n=\Xi_nu_n+\Phi_nb_n.
}
\tag{3.2}
\]
Equivalently,
\[
\binom{\kappa_n}{\kappa_{n-1}}
=
\begin{pmatrix}u_n&b_n\\u_{n-1}&b_{n-1}\end{pmatrix}
\binom{\Xi_n}{\Phi_n},
\tag{3.3}
\]
whose determinant is `1/n^3`, a unit for `n<p`.

At a Hasse zero `b_r=0`,
\[
\kappa_r=\Xi_r u_r.
\tag{3.4}
\]
Since `u_r` is a unit, commonness is exactly `Xi_r=0`, but **Phi_r is completely invisible to the endpoint equation**.  This is the fundamental reason the rank-two homogeneous frame does not automatically give two equations at a common row.

The increments are
\[
\Xi_n-\Xi_{n-1}=-5g_nb_{n-1},
\qquad
\Phi_n-\Phi_{n-1}=5g_nu_{n-1}.
\tag{3.5}
\]
Thus there really is a second global period, but common endpoints do not force it to vanish.

---

## 4. Exact no-go: every homogeneous two-return determinant is Green/Duhamel

For any homogeneous `y`, define the Lagrange concomitant
\[
W_n(y,\kappa)
=n^3(y_{n-1}\kappa_n-y_n\kappa_{n-1}).
\tag{4.1}
\]
Subtracting the homogeneous equation for `y`, multiplied by `kappa_{n-1}`, from the inhomogeneous equation for `kappa`, multiplied by `y_{n-1}`, gives
\[
\boxed{
W_n(y,\kappa)-W_{n-1}(y,\kappa)
=-5g_ny_{n-1}.
}
\tag{4.2}
\]

### Theorem 4.1 — Dirichlet test space is rank one

Assume `0<r<s<p` and `b_r=b_s=0`.  Then
\[
\{y:\ y\text{ homogeneous},\ y_r=y_s=0\}=\mathbf F_p b.
\tag{4.3}
\]

Indeed, `y_r=0` leaves only the value `y_{r+1}` free; since `b_{r+1}` is a unit, `y=(y_{r+1}/b_{r+1})b` everywhere.  The condition at `s` is then automatic because `b_s=0`.

Now suppose in addition `kappa_r=kappa_s=0`.  For every nonzero Dirichlet test `y=lambda b`, both boundary concomitants vanish, so telescoping (4.2) gives
\[
0=-5\lambda\sum_{m=r+1}^s g_mb_{m-1}.
\tag{4.4}
\]
Since `5` and `lambda` are units,
\[
\boxed{
\sum_{m=r+1}^s g_mb_{m-1}=0.
}
\tag{4.5}
\]
This is exactly `Xi_s-Xi_r=0`, hence exactly the Green/Duhamel shooting condition already isolated in Q7696.

Therefore any two-return monodromy determinant or adjoint pairing whose endpoint covector is a homogeneous solution vanishing at both endpoints is **one-dimensional** and cannot provide an independent second equation.  In continuant coordinates, the old fact `N_{h-1}(r)!=0` is precisely the statement that this unique terminal covector is nondegenerate.

Taking `y=u` in (4.2) does produce the independent period in (3.5), but `u_r,u_s` are units rather than zero.  At common endpoints its boundary values are `-Phi_r` and `-Phi_s`, not zero, so there is no second vanishing equation.

---

## 5. A genuine global Cartier residue for the canonical extension

Let
\[
\theta=t\frac d{dt},\qquad
L=\theta^3-tP(\theta)+t^2(\theta+1)^3,
\tag{5.1}
\]
and truncate
\[
K_{<p}(t)=\sum_{n=0}^{p-1}\kappa_nt^n,
\qquad
G_{<p}(t)=\sum_{n=0}^{p-1}g_nt^n.
\]
The row `n=1` has the canonical boundary defect
\[
-36=-5g_1-1
\]
because `g_1=7`.  Rows `2,...,p-1` satisfy the inhomogeneous recurrence.  At degree `p`, the truncation boundary is
\[
-P(p-1)\kappa_{p-1}+(p-1)^3\kappa_{p-2}
=5\kappa_{p-1}-\kappa_{p-2}.
\]
Using Hasse reciprocity `b_{p-2}=5,b_{p-1}=1`,
\[
\Xi_{p-1}
=(p-1)^3(5\kappa_{p-1}-\kappa_{p-2})
=-(5\kappa_{p-1}-\kappa_{p-2}).
\]
The degree `p+1` boundary term is `p^3 kappa_{p-1}=0`.  Hence:

### Theorem 5.1 — extension-class defect

\[
\boxed{
LK_{<p}
=-5(G_{<p}-1)-t-\Xi_{p-1}t^p
\quad\text{in }\mathbf F_p[t].
}
\tag{5.2}
\]
If `C_p(sum a_nt^n)=sum a_{pn}t^n` is Cartier extraction, then
\[
\boxed{
C_p(LK_{<p})=-\Xi_{p-1}t.
}
\tag{5.3}
\]

This is a genuine global extension residue.  It is not the internal Green period between `r` and `s`.

However, the existing strict-triangle Cartier defect for the homogeneous Green kernel reads only the Hasse coefficient `b_{a-1}`.  Thus the two available simple Cartier defects see, respectively, a local Hasse condition and a single global extension residue; neither couples two internal common rows.

---

## 6. Strong no-go: even both terminal variation periods can be preserved

The previous scalar defect might suggest adjoining the second terminal period `Phi_{p-1}`.  That still does not close the problem at the level of abstract source data.

### Theorem 6.1 — two-period terminal compensation

Let `p>=7` and let `r<s` be two Hasse zeros.  Start from any source sequence `s_n` in
\[
n^3x_n-P(n-1)x_{n-1}+(n-1)^3x_{n-2}=-5s_n
\]
with the canonical initial line `x_0=0,x_1=-36`.  There is a source perturbation supported at
\[
\{r,s,p-2,p-1\}
\]
which makes `x_r=x_s=0` while leaving both terminal variation coordinates `(Xi_{p-1},Phi_{p-1})` unchanged.

### Proof

Because `b_{p-2}=5` and `b_{p-1}=1`, every Hasse zero is at most `p-3`; hence the two terminal corrections occur strictly after `s`.

A source change `delta` at row `j` changes `x_j`, with earlier rows fixed, by
\[
\Delta x_j=-\frac5{j^3}\delta.
\]
The coefficient is a unit.  Therefore choose the perturbation at `r` to kill `x_r`, propagate, then choose the perturbation at `s` to kill the updated `x_s`.

A source perturbation `delta_j` changes the terminal variation vector by
\[
\binom{\Delta\Xi_{p-1}}{\Delta\Phi_{p-1}}
=5\delta_j\binom{-b_{j-1}}{u_{j-1}}.
\tag{6.1}
\]
For the last two rows the response matrix is
\[
M=5
\begin{pmatrix}
-b_{p-3}&-b_{p-2}\\
u_{p-3}&u_{p-2}
\end{pmatrix}.
\tag{6.2}
\]
By (1.1) at `n=p-2`,
\[
(p-2)^3(b_{p-3}u_{p-2}-b_{p-2}u_{p-3})=1.
\]
Since `(p-2)^3=-8`, the bracket is `-1/8`, and therefore
\[
\boxed{\det M=25/8,}
\tag{6.3}
\]
a unit for every `p>=7`.  Thus the two terminal perturbations are uniquely solvable so as to restore both terminal coordinates.  Being later than `s`, they do not disturb the two forced rows. ∎

This theorem is an obstruction, **not** an actual-source counterexample.  It says that a proof using only the Picard--Fuchs recurrence, the canonical initial line, and the terminal rank-two period vector cannot distinguish the canonical Eichler source from a source with two common rows.  Some coefficient-level identity special to
\[
G=1/(F^2\sqrt\Delta)
\]
must be used.

The executed exact example is `p=19`, Hasse zeros `8,10`.  Starting from the canonical source, the forcing changes are
\[
\delta_8=14,\qquad\delta_{10}=10,
\]
and the terminal compensators are
\[
\delta_{17}=14,\qquad\delta_{18}=8.
\]
Both forced rows are common and the original terminal periods
\[
(\Xi_{18},\Phi_{18})=(12,17)
\]
are restored exactly.

---

## 7. Actual-source reciprocity: exact theorem, exact failure of the tempting consequence

For the first characteristic-`p` block let
\[
A_p(t)=\sum_{j=0}^{p-1}b_jt^j,
\qquad
R_p(t)=\frac{\Delta(t)^{(p-1)/2}}{A_p(t)^2}.
\]
The actual source satisfies `g_m=[t^m]R_p` for `m<p`.  Because both `A_p` and `Delta^{(p-1)/2}` are reciprocal,
\[
\boxed{R_p(1/t)=t^{p-1}R_p(t).}
\tag{7.1}
\]
This is a real global symmetry of the canonical source.

But it relates the Taylor expansion at `0` to the Laurent expansion at `infinity`; it does **not** imply a reflection identity among the first `p` Taylor coefficients.  In particular the tempting law
\[
g_{p-m}b_m=g_mb_{m-1}
\tag{7.2}
\]
is false.  Already at `p=7,m=1`, the two sides are `6` and `0`.

The induced affine partner law
\[
\Xi_r+\Xi_{p-1-r}=\text{constant}
\]
is also false.  More strongly, at the single prime `p=181`, all four rows
\[
19,47,133,161
\]
are Hasse zeros, but the two reflection-pair sums are
\[
\Xi_{19}+\Xi_{161}=121,
\qquad
\Xi_{47}+\Xi_{133}=10.
\]
Thus even restricting to the Hasse-zero locus does not rescue a constant reciprocal-partner relation.

This explains why a `t`-plane residue/norm theorem is not automatically a coefficient-index theorem: `Xi_r` is a prefix of the **Hadamard product** `g_m b_{m-1}`, not a root or residue of `R_p(t)` as a rational function of `t`.

---

## 8. The independent second period and the `p=41` guard

At `p=41`, the actual canonical source has
\[
b_{10}=b_{30}=0,
\qquad
\Xi_{10}=\Xi_{30}=7.
\]
Therefore
\[
\sum_{m=11}^{30}g_mb_{m-1}=0.
\]
The consolidated verifier computes, exactly,
\[
\boxed{
\sum_{m=11}^{30}g_mu_{m-1}=6\ne0,
\qquad
\Phi_{30}-\Phi_{10}=30.
}
\tag{8.1}
\]
So the second rank-two period is genuinely independent of the Green period.  But because common endpoints constrain only `Xi`, not `Phi`, this independence does not prove uniqueness.

---

## 9. Targeted finite experiments

The consolidated verifier is

```text
problems/3.2/research/scripts/q7702_zero_fiber_followup.sage
```

on branch

```text
q7702-zero-fiber-635f6233
```

GitHub Actions run `31572862600`, job `94038403530`, completed successfully with final line

```text
Q7702_ZERO_FIBER_FOLLOWUP=PASS
```

The mechanism scan through every prime `7<=p<=5000` found:

```text
RECIPROCAL_Q_FAILURE_COUNT 666
AFFINE_XI_PARTNER_FAILURE_COUNT 666
P181_PAIR_SUMS (121, 10)
P41_B_PERIOD 0 P41_U_PERIOD 6 P41_PHI_DIFF 30
RESIDUE_ZERO_FIRST20 []
RESIDUE_ONE_FIRST20 [(7, (), 1)]
RESIDUE_AT_COMMON_PRIMES [(17, (13,), 7, 7), (2237, (492,), 1023, 1438)]
P19_FORCE_AND_RESTORE ([(8, 14), (10, 10)], (14, 8), 12, 17)
```

This is not a repeat of the old uniqueness count: it tests the candidate reciprocal law, the full homogeneous reflection theorem, the exact extension defect, the second period, and terminal-period compensation.

The same verifier then tested the frame and extension-defect identities at the first twelve primes **strictly above 20000**:

```text
(20011, 0, (), 721, 2105)
(20021, 2, (), 19648, 286)
(20023, 2, (), 16463, 7796)
(20029, 0, (), 16010, 19023)
(20047, 0, (), 4967, 7240)
(20051, 0, (), 2883, 10821)
(20063, 0, (), 20012, 12238)
(20071, 0, (), 4230, 19911)
(20089, 0, (), 6901, 6516)
(20101, 0, (), 13320, 17094)
(20107, 2, (), 9503, 215)
(20113, 2, (), 3271, 11778)
```

Each tuple is `(p, number_of_Hasse_zeros, common_rows, Xi_{p-1}, Phi_{p-1})`.  The absence of a new common row in these twelve primes is incidental finite data and is **not** used as evidence for a proof.

---

## 10. Strongest surviving theorem and exact missing implication

No `O(1)` multiplicity theorem follows from the present argument.  The elementary unconditional bound remains only linear (`<= (p-3)/2` for `p>=7`, from no consecutive Hasse zeros and the two nonzero terminal coefficients).

The new unconditional content is structural:

- full homogeneous Hasse reciprocity, Theorem 2.1;
- rank-one Dirichlet/monodromy collapse, Theorem 4.1;
- the canonical extension Cartier residue, Theorem 5.1;
- terminal two-period compensation, Theorem 6.1;
- exact counterexamples to simple reciprocal partner laws.

Put
\[
H_p(n)=\sum_{m=1}^n b_{m-1}[t^m]R_p(t).
\]
For `p>=7`,
\[
\Xi_n=0\iff H_p(n)=-1/5.
\]
Thus the sharp missing theorem can be stated as the special-level Hadamard/Eichler assertion
\[
\boxed{
 b_r=b_s=0,\quad H_p(r)=-1/5,\quad r<s
 \Longrightarrow H_p(s)\ne-1/5.
}
\tag{10.1}
\]
Equivalently, it is `(EIS-2RET)` above.

The experiments and no-go theorems say what a successful proof must add: a defining-characteristic identity that constrains the **Hadamard-prefix/Mellin coefficient index** of the canonical Eichler source.  Root geometry of `A_p(t)`, ordinary rational reciprocity of `R_p(t)`, the homogeneous transfer monodromy, and finitely many global period residues do not supply that implication.

A plausible next target is therefore not another determinant.  It is a modular-symbol/Cartier statement for the canonical Eichler extension after taking the coefficient-index (Hadamard/Mellin) transform, with enough rigidity to make the distinguished level `-1/5` simple on the Hasse-zero set.  That exact implication remains unproved here.
