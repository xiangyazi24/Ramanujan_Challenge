ANSWER Q1402 debadaed

# Singleton one-wrap Smith layer: exact calculation

## Verdict

There is a clean closed form for the part of the first divided Schur layer that controls the **primitive numerator after saturation**.

Put

\[
r^\vee=p-1-r,\qquad K=p+r,\qquad
T_s:=U_{s,1}=\binom{s+1}{2}=\frac{s(s+1)}2,
\]

and write all first divided quantities in \(\mathbf F_p\) as

\[
\beta:=\frac{b_r}{p}\pmod p,\qquad D:=D_r\pmod p,
\qquad a_r:=\frac{2r+1}{2}\pmod p.
\]

Because \(r<(p-1)/2\), \(a_r\ne0\).

The three exceptional Apéry rows have divided digits

\[
\mu_0=\beta,\qquad
\mu_-=\beta-2D,\qquad
\mu_+=5(\beta+2D),
\tag{0.1}
\]

at \(r,r^\vee,K\), respectively.  The two denominator-only modular kernel directions give the \(2\times2\) minor

\[
M_{\rm den}=
\begin{pmatrix}
2D&-(\beta-2D)\\
-4\beta-10D&5(\beta+2D)
\end{pmatrix},
\tag{0.2}
\]

whose determinant is

\[
\boxed{
\det M_{\rm den}
=-4\Delta_r,
\qquad
\Delta_r:=\beta^2-2\beta D_r-10D_r^2\pmod p.
}
\tag{0.3}
\]

This is the decisive scalar.

**Exact conclusion.** The following are equivalent:

1. \(\Delta_r\ne0\pmod p\);
2. the first divided Schur map has rank two and its one-dimensional kernel has nonzero numerator projection;
3. the \(p\)-primary maximal determinantal divisor has valuation exactly two and, after dividing the cofactor kernel by that \(p^2\), the primitive numerator is nonzero modulo \(p\).

Under these equivalent conditions,

\[
P_{\rm prim}(s)\equiv u\,(T_s-T_r)\pmod p,
\qquad u\in\mathbf F_p^\times.
\tag{0.4}
\]

In the notation of the question this is precisely

\[
\boxed{
P_{\rm prim}(Y)
\equiv u\bigl(\Phi_1(Y)-\Phi_1(\tau(r))\bigr)\pmod p.
}
\tag{0.5}
\]

Thus its unique folded root is \(r\).

There is, however, one important correction to the literal wording of the question.  **Bare rank two alone is not characterized by \(\Delta_r\).**  If \(\Delta_r=0\), the full \(2\times3\) block can still have rank two, but then its unique kernel is denominator-only and the saturated numerator is zero modulo \(p\).  Bare rank requires one additional lower-half scalar \(\chi_r\), written explicitly in §7 below.  Consequently there is no rank-two criterion depending only on \((r,\beta,D_r)\) unless one proves an additional identity eliminating \(\chi_r\).  For the phenomenon actually observed in the computation—**rank two plus a surviving primitive numerator**—\(\Delta_r\ne0\) is necessary and sufficient.

Finally, the currently banked Apéry recurrence and no-consecutive-zero facts do **not** prove \(\Delta_r\ne0\).  The exact extra theorem needed is the first-jet transversality statement

\[
\boxed{
Z_p^{\rm fold}=\{r\},\ r<(p-1)/2
\quad\Longrightarrow\quad
\beta^2-2\beta D_r-10D_r^2\not\equiv0\pmod p.
}
\tag{0.6}
\]

I found no actual Apéry counterexample in the repository material.  The computation stated in the question implies, by (0.3)--(0.5), that (0.6) holds for every singleton case already scanned.  But (0.6) is a new arithmetic input; it is not a consequence of the presently banked recurrence/no-adjacency statements.

For reference, the complete result can be summarized as:

```text
three exceptional rows:       r, p-1-r, p+r
first divided Apéry digits:   beta, beta-2D, 5(beta+2D)
triangular first displacement: -a_r, +a_r
mod-p row defect:              2
mod-p right-kernel dimension:  3
p-primary divisor if Delta!=0: p^2 exactly
primitive numerator slope:    unit * (-4 Delta)
Delta:                         beta^2 - 2 beta D - 10 D^2
```

---

## 1. Why there are exactly three exceptional rows

The full zero set below \(p\) is the reflected pair

\[
r,\qquad r^\vee=p-1-r,
\]

because the folded zero set is the singleton \(\{r\}\) and \(r<(p-1)/2\).

For the wrapped rows write \(s=p+u\), \(0\le u\le r\).  Apéry--Lucas gives

\[
b_{p+u}\equiv b_1b_u=5b_u\pmod p.
\]

Since the only folded zero at or below \(r\) is \(r\), the only wrapped zero row is \(u=r\), namely

\[
K=p+r.
\]

Therefore among the rows \(0,\ldots,K\), exactly

\[
\boxed{r,\quad r^\vee=p-1-r,\quad K=p+r}
\tag{1.1}
\]

have \(b_s\equiv0\pmod p\).

At these three rows the denominator block vanishes modulo \(p\).  Their numerator rows are identical because

\[
T_r\equiv T_{r^\vee}\equiv T_K\pmod p.
\tag{1.2}
\]

Thus every exceptional row is

\[
[1,T_r\mid0,\ldots,0]
\pmod p.
\]

---

## 2. Exact mod-\(p\) defect: two

Let \(S\) be the set of the \(K-2\) nonexceptional rows.  On \(S\), every \(b_s\) is a \(p\)-unit.  Consider the denominator columns indexed by the same labels \(k=s\in S\).  Since

\[
U_{s,k}=0\quad(k>s),
\qquad U_{s,s}=1,
\]

the submatrix

\[
(U_{s,k})_{s,k\in S}
\]

is triangular with diagonal one.  Hence the denominator block on the nonexceptional rows has rank exactly

\[
K-2.
\tag{2.1}
\]

Those \(K-2\) nonexceptional rows are linearly independent after projection to the denominator coordinates.  Adding one exceptional row increases the full row rank by one, because an exceptional row has zero denominator part and no nontrivial combination of the independent nonexceptional rows can have zero denominator part.

Therefore

\[
\operatorname{rank}_{\mathbf F_p}(A)=K-1.
\tag{2.2}
\]

Since \(A\) has \(K+1\) rows and \(K+2\) columns,

\[
\boxed{
\dim\operatorname{coker}_{\rm left}(A\bmod p)=2,
\qquad
\dim\ker(A\bmod p)=3.
}
\tag{2.3}
\]

Two particularly convenient left-defect functionals are simply

\[
\ell_-:=e_{r^\vee}-e_r,
\qquad
\ell_+:=e_K-e_r.
\tag{2.4}
\]

They annihilate \(A\bmod p\), are independent, and hence form a basis of the left defect.

On the right, an element of \(\ker(A\bmod p)\) has numerator

\[
N_s=c_0+c_1T_s.
\]

The exceptional row forces

\[
c_0+c_1T_r=0.
\]

Thus the numerator projection of the modular kernel is one-dimensional.  Write

\[
\lambda:=c_1,
\qquad
N_s=\lambda(T_s-T_r).
\tag{2.5}
\]

The remaining two modular kernel dimensions are denominator-only.

This proves the observed mod-\(p\) row defect two without any Smith assumption.

---

## 3. The three first Apéry digits

The wrapped digit is the given first-jet formula with \(N=1\):

\[
b_{p+r}
\equiv b_1(b_r+2pD_r)
=5b_r+10pD_r
\pmod{p^2}.
\tag{3.1}
\]

Hence

\[
\frac{b_K}{p}\equiv5(\beta+2D_r)\pmod p.
\tag{3.2}
\]

The reflected digit is also explicit.  Recall

\[
b_n=\sum_{k=0}^n
\binom nk^2\binom{n+k}{k}^2.
\]

For \(n=p-1-r\), all terms with \(k>r\) are divisible by \(p^2\): the factor \(\binom{p-1-r+k}{k}\) is divisible by \(p\), and it occurs squared.  For \(0\le k\le r\), since \(2r<p-1\), all harmonic denominators below are \(p\)-units and

\[
\binom{p-1-r}{k}
\equiv
(-1)^k\binom{r+k}{k}
\left(1-p(H_{r+k}-H_r)\right)
\pmod{p^2},
\tag{3.3}
\]

\[
\binom{p-1-r+k}{k}
\equiv
(-1)^k\binom rk
\left(1-p(H_r-H_{r-k})\right)
\pmod{p^2}.
\tag{3.4}
\]

Squaring and multiplying gives

\[
\binom{p-1-r}{k}^2
\binom{p-1-r+k}{k}^2
\equiv
\binom rk^2\binom{r+k}{k}^2
\left(1-2p(H_{r+k}-H_{r-k})\right)
\pmod{p^2}.
\]

Summing over \(k\) gives the exact reflected first jet

\[
\boxed{
 b_{p-1-r}\equiv b_r-2pD_r\pmod{p^2}.
}
\tag{3.5}
\]

Consequently the three divided digits are indeed

\[
\boxed{
\mu_0:=\frac{b_r}{p}=\beta,
\quad
\mu_-:=\frac{b_{r^\vee}}p=\beta-2D,
\quad
\mu_+:=\frac{b_K}p=5(\beta+2D)
}
\quad(\bmod p).
\tag{3.6}
\]

No assumption that \(\beta\) is nonzero is permitted.  In fact the repository's `q5731_cartier_first_digit_audit.py` records the useful warning example \(p=17,r=3\), where \(17^2\mid b_3\), so \(\beta=0\).

---

## 4. The triangular first displacement

The three exceptional triangular coordinates agree modulo \(p\), but their divided differences do not.  Direct calculation gives

\[
T_{r^\vee}-T_r
=\frac{p(p-2r-1)}2,
\]

\[
T_K-T_r
=\frac{p(p+2r+1)}2.
\]

Hence

\[
\boxed{
\frac{T_{r^\vee}-T_r}{p}\equiv-a_r,
\qquad
\frac{T_K-T_r}{p}\equiv+a_r,
\qquad
a_r=\frac{2r+1}{2}\ne0
}
\pmod p.
\tag{4.1}
\]

This is the only first derivative of the numerator carrier that enters the two defect rows.

---

## 5. The two explicit first-divided Schur equations

Let \(x\in\ker(A\bmod p)\).  Write its numerator as in (2.5) and let

\[
q_0:=Q(r),\qquad q_-:=Q(r^\vee),\qquad q_+:=Q(K)
\pmod p
\]

be the three exceptional values of its denominator polynomial in the triangular basis.

Apply the left defect rows (2.4) to an integral lift of \(x\), divide by \(p\), and reduce modulo \(p\).  Using (3.6) and (4.1) gives

\[
\boxed{
-a_r\lambda-\mu_-q_-+\mu_0q_0=0,
}
\tag{5.1}
\]

\[
\boxed{
+a_r\lambda-\mu_+q_++\mu_0q_0=0.
}
\tag{5.2}
\]

Equivalently,

\[
\boxed{
2\mu_0q_0-\mu_-q_--\mu_+q_+=0,
}
\tag{5.3}
\]

\[
\boxed{
2a_r\lambda+\mu_-q_--\mu_+q_+=0.
}
\tag{5.4}
\]

Equations (5.1)--(5.2) are the requested symbolic first divided Schur layer in an invariant form.  They are not a schematic factorization: every coefficient is explicit in \(r,\beta,D_r\).

---

## 6. The denominator-only plane and the decisive \(2\times2\) minor

It remains to identify exactly which triples \((q_0,q_-,q_+)\) occur for the two denominator-only modular kernel directions.

Let

\[
V_K=(U_{s,k})_{0\le s\le K,\ 0\le k\le K-1}.
\]

The columns of \(V_K\) satisfy one exact integral row relation.  One convenient normalization is

\[
L_s=(-1)^{K-s}
\left[
\binom{2K}{K-s}-\binom{2K}{K-s-1}
\right],
\tag{6.1}
\]

for which

\[
\sum_{s=0}^K L_sU_{s,k}=0,
\qquad0\le k\le K-1.
\tag{6.2}
\]

One way to see this is to write

\[
U_{s,k}
=\frac{2^k}{(2k)!}
\prod_{j=0}^{k-1}(T_s-T_j).
\]

Thus the columns span the degree-\(<K\) polynomials in the distinct rational nodes \(T_0,\ldots,T_K\), and (6.1) is the integer-scaled barycentric relation.  The difference-of-binomial form is useful modulo \(p\).

For \(K=p+r\), Lucas reduction gives

\[
L_r\equiv-2,
\qquad
L_{r^\vee}\equiv1,
\qquad
L_K\equiv1
\pmod p.
\tag{6.3}
\]

A denominator-only kernel vector has \(Q(s)=0\) at every nonexceptional row, because \(b_s\) is then a unit.  Reducing (6.2) therefore gives

\[
\boxed{-2q_0+q_-+q_+=0.}
\tag{6.4}
\]

The evaluation map from the two-dimensional denominator kernel to these three values is injective: if all three exceptional values also vanish, then \(Q(s)=0\) at every row, and the first \(K\) rows of \(V_K\) are triangular with unit diagonal, forcing \(Q=0\).  Hence (6.4) is exactly the two-dimensional image plane.

Choose the basis

\[
h_0=(1,1,1),
\qquad
h_1=(0,1,-1).
\tag{6.5}
\]

Substitution in (5.1)--(5.2), with \(\lambda=0\), gives (0.2):

\[
M_{\rm den}=
\begin{pmatrix}
\mu_0-\mu_-&-\mu_-\\
\mu_0-\mu_+&\mu_+
\end{pmatrix}
=
\begin{pmatrix}
2D&-(\beta-2D)\\
-4\beta-10D&5(\beta+2D)
\end{pmatrix}.
\tag{6.6}
\]

Its determinant is

\[
\begin{aligned}
\det M_{\rm den}
&=\mu_0(\mu_-+\mu_+)-2\mu_-\mu_+\\
&=\beta(6\beta+8D)
 -10(\beta-2D)(\beta+2D)\\
&=-4\beta^2+8\beta D+40D^2\\
&=-4(\beta^2-2\beta D-10D^2).
\end{aligned}
\tag{6.7}
\]

This proves (0.3).

---

## 7. The complete \(2\times3\) block and the bare-rank caveat

For completeness, here is the full \(2\times3\) block, including the one quantity that disappears from the desired numerator criterion.

The Lucas reduction of (6.1) is slightly more precise.  For \(0\le u\le r\), put

\[
w_u=(-1)^{r-u}
\left[
\binom{2r}{r-u}-\binom{2r}{r-u-1}
\right].
\tag{7.1}
\]

Then

\[
L_u\equiv-2w_u,
\qquad
L_{p-1-u}\equiv w_u,
\qquad
L_{p+u}\equiv w_u
\pmod p,
\tag{7.2}
\]

and all remaining \(L_s\) vanish modulo \(p\).  Note \(w_r=1\), recovering (6.3).

Take a modular kernel vector with \(\lambda=1\).  At a nonexceptional triple \(u,p-1-u,p+u\) with \(u<r\), the kernel equation \(N_s=b_sQ(s)\), reflection, and Apéry--Lucas give

\[
Q(u)=\frac{T_u-T_r}{b_u},
\quad
Q(p-1-u)=\frac{T_u-T_r}{b_u},
\quad
Q(p+u)=\frac{T_u-T_r}{5b_u}
\pmod p.
\tag{7.3}
\]

The singleton hypothesis makes every \(b_u\), \(u<r\), a unit.  Inserting (7.3) into the row relation yields the invariant exceptional combination

\[
\boxed{
\chi_r:=-2q_0+q_-+q_+
=\frac45\sum_{u=0}^{r-1}
 w_u\frac{T_u-T_r}{b_u}
\pmod p.
}
\tag{7.4}
\]

Adding a denominator-only kernel vector does not change \(\chi_r\), by (6.4).  Conversely, using the two denominator directions one can choose the \(\lambda=1\) representative with

\[
(q_0,q_-,q_+)=(0,0,\chi_r).
\tag{7.5}
\]

In the right-kernel basis \((h_0,h_1,v)\) given by (6.5) and (7.5), and the left basis \((\ell_-,\ell_+)\), the first divided Schur block is therefore

\[
\boxed{
B_r=
\begin{pmatrix}
2D&-(\beta-2D)&-a_r\\
-4\beta-10D&5(\beta+2D)&a_r-5(\beta+2D)\chi_r
\end{pmatrix}.
}
\tag{7.6}
\]

This is an actual \(2\times3\) matrix, not only a defect-space schematic.

Its three \(2\times2\) minors are

\[
\delta_{12}=-4\Delta_r,
\tag{7.7}
\]

\[
\delta_{13}
=-2(\beta+2D)\bigl(2a_r+5D\chi_r\bigr),
\tag{7.8}
\]

\[
\delta_{23}
=4a_r(\beta+3D)
 +5(\beta^2-4D^2)\chi_r.
\tag{7.9}
\]

Define

\[
\Theta_r:=4a_r(\beta+3D)
 +5(\beta^2-4D^2)\chi_r.
\tag{7.10}
\]

If \(\Delta_r=0\) and \((\beta,D)\ne(0,0)\), the first two columns of (7.6) have rank one, and the second column is nonzero.  Hence \(\delta_{23}\) alone tests whether the third column leaves that line.  If \((\beta,D)=(0,0)\), both denominator columns vanish and (7.6) has rank one because \(a_r\ne0\).  Thus the exact bare-rank criterion is

\[
\boxed{
\operatorname{rank}B_r=2
\quad\Longleftrightarrow\quad
\Delta_r\ne0\ \text{or}\ \Theta_r\ne0.
}
\tag{7.11}
\]

This is the promised caveat: bare rank depends on \(\chi_r\) when \(\Delta_r=0\).  Therefore a statement that **bare rank two is equivalent to a scalar involving only \(\beta,D_r\)** is false unless one adds a new identity that eliminates \(\chi_r\).

If one literally insists on a single field element rather than the pair in (7.11), Fermat indicators can package the logical OR as

\[
\mathcal R_r
=1-(1-\Delta_r^{p-1})(1-\Theta_r^{p-1}),
\tag{7.12}
\]

for which \(\mathcal R_r=1\) exactly when \(B_r\) has rank two.  This is algebraically correct but conceptually less informative than (7.11).

The good news is that \(\chi_r\) is irrelevant to the primitive-numerator theorem.

---

## 8. Why \(\Delta_r\) is exactly the primitive numerator coordinate

For any \(2\times3\) rank-two matrix, a kernel generator is given by its signed \(2\times2\) minors.  In the column order of (7.6), the coordinate multiplying the third basis vector \(v\)—the only right-kernel basis vector with nonzero numerator—is

\[
\delta_{12}=-4\Delta_r.
\tag{8.1}
\]

The first two basis vectors are denominator-only.  Therefore, whenever \(B_r\) has rank two,

\[
\text{numerator projection of }\ker B_r\ne0
\quad\Longleftrightarrow\quad
\Delta_r\ne0.
\tag{8.2}
\]

Moreover \(\Delta_r\ne0\) already makes the first two columns independent, so it automatically gives rank two.  Conversely, if \(\Delta_r=0\) but \(\Theta_r\ne0\), then \(B_r\) still has rank two, but the pre-existing one-dimensional kernel of its denominator restriction is the entire kernel; hence the unique first-layer survivor is denominator-only.

Thus

\[
\boxed{
\bigl(\operatorname{rank}B_r=2\bigr)
\ \&\ 
\bigl(\pi_{\rm num}\ker B_r\ne0\bigr)
\quad\Longleftrightarrow\quad
\Delta_r\ne0.
}
\tag{8.3}
\]

This is exactly the empirical package in the question.

---

## 9. Determinantal divisor and saturation

Now work over \(\mathbf Z_p\).  Since \(A\bmod p\) has rank \(K-1\), unit row and column operations put it in the form

\[
\begin{pmatrix}
I_{K-1}&0\\
0&0_{2\times3}
\end{pmatrix}
\pmod p.
\]

Choose the last two rows to be lifts of \(\ell_-,\ell_+\) and the last three columns to be lifts of the right-kernel basis used above.  Clearing the off-diagonal blocks with the unit pivots gives, modulo \(p^2\),

\[
A\sim
\begin{pmatrix}
I_{K-1}&0\\
0&pB_r
\end{pmatrix}
\pmod{p^2}.
\tag{9.1}
\]

This is the usual divided Schur reduction, but here \(B_r\) has already been computed explicitly in (7.6).

If \(\Delta_r\ne0\), then \(B_r\) has rank two.  Consequently the \(p\)-primary Smith valuations of the nonzero invariant factors of \(A\) are

\[
\underbrace{0,\ldots,0}_{K-1\text{ times}},1,1.
\tag{9.2}
\]

Hence the top determinantal divisor satisfies

\[
\boxed{
v_p\bigl(\Delta_{K+1}(A)\bigr)=2.
}
\tag{9.3}
\]

For a full-row-rank \((K+1)\times(K+2)\) integer matrix, the signed maximal minors form a kernel vector, and their gcd is the top determinantal divisor.  Thus every raw cofactor coordinate is divisible by \(p^2\), at least one is not divisible by \(p^3\), and dividing by the primitive gcd removes exactly \(p^2\) at this prime.

After that division, reduction modulo \(p\) gives the signed-minor kernel of \(B_r\).  By (8.1), its numerator slope is a \(p\)-unit multiple of

\[
-4\Delta_r.
\]

Therefore the primitive numerator is

\[
P_{\rm prim}(s)
\equiv u(T_s-T_r)\pmod p,
\qquad u\ne0.
\tag{9.4}
\]

This is the requested determinantal-divisor tracking.  It is stronger than saying merely that the modular numerator projection is one-dimensional: it identifies the precise first nonzero cofactor digit after the \(p^2\) Smith content is removed.

If \(\Delta_r=0\) and \(B_r\) nevertheless has rank two because \(\Theta_r\ne0\), (9.3) still holds, but the primitive cofactor vector reduces to the denominator-only kernel direction.  In that case every primitive numerator coefficient remains divisible by \(p\).  This is exactly why bare Schur rank is not enough.

---

## 10. Translation to \(\Phi_1\) and the folded root

The numerator columns of the original matrix are \(U_{s,0}=1\) and \(U_{s,1}=T_s\).  Thus (9.4) says

\[
(c_0,c_1)\equiv u(-T_r,1)\pmod p.
\]

Whatever normalization of the project variable \(Y\) is used, evaluation at row \(s\) sends the degree-one carrier to

\[
\Phi_1(\tau(s))=U_{s,1}=T_s.
\]

Hence

\[
P_{\rm prim}(Y)
\equiv
u\left(\Phi_1(Y)-\Phi_1(\tau(r))\right)
\pmod p,
\quad \nu\in\mathbf F_p^\times.
\tag{10.1}
\]

The reflected row \(r^\vee\) and the wrapped row \(p+r\) have the same triangular coordinate modulo \(p\); folding identifies them with \(r\).  Therefore the unique folded root is exactly \(r\).

---

## 11. Do the Apéry recurrence and no-consecutive-zero theorem prove \(\Delta_r\ne0\)?

No, not from the currently banked facts.

At a mod-\(p\) zero, the ordinary Apéry recurrence only gives the neighbor relation

\[
(r+1)^3b_{r+1}+r^3b_{r-1}\equiv0\pmod p,
\tag{11.1}
\]

because the middle term \(P(r)b_r\) vanishes.  No-consecutive-zero says that the two neighboring values are units; it does not involve either divided digit \(\beta=b_r/p\) or the first jet \(D_r\).

Lifting the recurrence modulo \(p^2\) introduces the first \(p\)-digits of the two neighboring values.  Schematically it is one linear relation of the form

\[
(r+1)^3\dot b_{r+1}+r^3\dot b_{r-1}
\equiv P(r)\beta
\pmod p,
\tag{11.2}
\]

so it does not isolate \(\beta\).

Similarly, differentiating the lower-half first-jet recurrence gives the linear inhomogeneous relation

\[
\begin{aligned}
2\bigl((r+1)^3D_{r+1}-P(r)D_r+r^3D_{r-1}\bigr)
={}&P'(r)b_r\\
&-3(r+1)^2b_{r+1}-3r^2b_{r-1}.
\end{aligned}
\tag{11.3}
\]

At \(b_r\equiv0\pmod p\), this relates \(D_r\) to the *neighboring jets* \(D_{r\pm1}\) and neighboring Apéry values.  Again it gives no quadratic exclusion for the pair \((\beta,D_r)\).

The Smith transversality condition is instead

\[
\Delta_r=\beta^2-2\beta D_r-10D_r^2\ne0.
\tag{11.4}
\]

That is a genuinely quadratic relation between two first divided coordinates.  It is not a reformulation of (11.1), no-adjacency, or either of the linear lifted recurrences above.

A useful way to see its arithmetic shape is

\[
\Delta_r
=(\beta-(1+\sqrt{11})D_r)
 (\beta-(1-\sqrt{11})D_r)
\tag{11.5}
\]

in a quadratic extension.  Thus, if \(D_r\ne0\), failure would force

\[
\frac\beta{D_r}=1\pm\sqrt{11}.
\tag{11.6}
\]

In particular, for primes with \(11\) a quadratic nonresidue, a failure can only occur through the deeper simultaneous condition \(\beta=D_r=0\).  The current recurrence/no-consecutive package does not rule that out either.

The example \((p,r)=(17,3)\) is also instructive: \(\beta=0\) because \(17^2\mid b_3\).  Therefore any attempted proof based on simplicity of the zero \(b_r/p\ne0\) is already false.  The observed successful Smith layer at this prime must come from the \(D_r\)-term, exactly as (11.4) predicts.

So the missing theorem is not another modular-zero theorem.  It is precisely the first-jet transversality statement (0.6).

---

## 12. Exact theorem that can be banked conditionally

Here is the clean theorem statement suggested by the calculation.

### Theorem (singleton one-wrap primitive Smith, conditional only on jet transversality)

Let \(p\ge7\) be prime.  Assume the folded Apéry zero set is \(\{r\}\) with \(0\le r<(p-1)/2\), let \(K=p+r\), and let \(A=A_{K,1}\).  Define

\[
\beta=b_r/p\pmod p,
\qquad
\Delta_r=\beta^2-2\beta D_r-10D_r^2\pmod p.
\]

If \(\Delta_r\ne0\), then:

1. \(\operatorname{rank}(A\bmod p)=K-1\), so the mod-\(p\) row defect is two;
2. the first divided \(2\times3\) Schur layer has rank two;
3. the top determinantal divisor has exactly \(p^2\) as its \(p\)-part;
4. the primitive saturated kernel has nonzero numerator projection modulo \(p\);
5. that numerator is a nonzero scalar multiple of
   \[
   \Phi_1(Y)-\Phi_1(\tau(r));
   \]
6. hence its unique folded root is \(r\).

Conversely, if the first divided layer has rank two and the saturated primitive numerator is nonzero modulo \(p\), then \(\Delta_r\ne0\).

Thus the empirical singleton theorem is **equivalent** to the single arithmetic nonvanishing condition \(\Delta_r\ne0\).

### What is still open

To remove the conditional hypothesis one must prove

\[
\beta^2-2\beta D_r-10D_r^2\ne0\pmod p
\]

for every singleton folded Apéry zero.  I do not see a derivation from the recurrence and no-consecutive-zero theorem alone, and the repository material I inspected contains no such first-jet nonvanishing theorem.

If one only wants bare rank two, the exact extra lower-half invariant is \(\chi_r\) from (7.4), and the criterion is (7.11).  This is the precise additional hypothesis/data omitted by a formulation that asks bare rank to depend only on \(\beta,D_r\).

---

## 13. Bottom line for Q1402

The symbolic calculation does reach the determinantal target, but it exposes one new arithmetic lemma rather than closing it from the existing zero-set facts.

\[
\boxed{
\Delta_r=\left(\frac{b_r}{p}\right)^2
-2\left(\frac{b_r}{p}\right)D_r
-10D_r^2
\pmod p
}
\]

is the exact singleton Smith transversality invariant.

* \(\Delta_r\ne0\) is necessary and sufficient for the observed **rank-two first layer + nonzero primitive numerator** package.
* Under \(\Delta_r\ne0\), the \(p\)-part of the maximal determinantal divisor is exactly \(p^2\), and saturation leaves
  \[
  P_{\rm prim}(Y)\equiv
  u\bigl(\Phi_1(Y)-\Phi_1(\tau(r))\bigr)\pmod p.
  \]
* The existing Apéry recurrence and no-consecutive-zero facts do not prove \(\Delta_r\ne0\); the missing theorem is exactly the first-jet transversality (0.6).
* There is no known actual prime/residue counterexample in the scanned data.  If one is eventually found, it must satisfy
  \[
  \beta^2-2\beta D_r-10D_r^2=0,
  \]
  and then the desired primitive-numerator conclusion fails at the first divided layer.
* If \(\Delta_r=0\) but the bare \(2\times3\) layer still has rank two, that rank is supplied by \(\Theta_r\) in (7.10), and the unique saturated first-layer survivor is denominator-only.  This is the exact distinction that a rank-only formulation misses.
