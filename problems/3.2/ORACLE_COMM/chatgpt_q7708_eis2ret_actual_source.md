# Q7708 — EIS-2RET on the actual level-six source

## Verdict

I do **not** prove

\[
\#\{0<r<p:b_r=\Xi_r=0\}\le1
\]

or a new uniform `O(1)` multiplicity bound, and I found no actual-source counterexample.  The exact FLINT scan now present on `main` is stronger than the bound quoted in the prompt: `chatgpt_q32_transverse_common_fast_scan.md` exhausts every prime `p<=100000` and every `0<r<p` and still finds only

\[
(17,13),\qquad(2237,492).
\]

That is finite evidence only.

The proof-first outcome is a sharper localization of EIS-2RET on the **actual level-six source**.

1. `g_m` has an exact level-six modular residue formula, an explicit finite eta-product coefficient formula, and an exact Franel/elliptic first-block Cartier formula.
2. The homogeneous Faber family is diagonal for Cartier:
   \[
   \boxed{\mathcal C_p\!\left(\Omega_m\frac{dq}{q}\right)
   =b_m\Psi\frac{dq}{q}},\qquad \Omega_m=\Psi t^{-m},\ 0\le m<p.
   \]
   Thus **every Hasse row is killed before `Xi_m` is visible**.  Ordinary `U_p`/Cartier cannot be the second constraint.
3. There is an exact source-specific two-return Bol/Eichler cocycle using the fixed canonical Eichler boundary.  Its pairing reduces *exactly* to the known Green/Duhamel shooting law, so it is not a second equation.
4. The inhomogeneous reflection cocycle which landed concurrently on current `main`,
   \[
   C=J\kappa-\kappa,
   \qquad \mathcal L_rC=5(g_{p-r}+g_{r+1}),
   \]
   also fails to give a second return equation.  For a generic common pair its endpoint values are unconstrained reflected `kappa` values; for a reflected pair, where its endpoints do vanish, its shooting law is exactly `2` times the original Green law.
5. For Hasse zeros `r<s`, the actual Green period is one coefficient against a **contiguous reflected block of the Apéry Hasse polynomial**:
   \[
   \boxed{
   S_{r,s}:=\sum_{m=r+1}^{s}b_{m-1}g_m
   =[T^{p-1}]R_p(T)C_{r,s}(T),
   }
   \]
   \[
   R_p(T)=\frac{\Delta(T)^{(p-1)/2}}{A_p(T)^2},
   \qquad
   C_{r,s}(T)=T^{p-s-1}\sum_{j=0}^{s-r-2}b_{s-1-j}T^j.
   \]
   The two endpoint coefficients of the inner block are the units `b_{s-1}` and `b_{r+1}`.

Hence the smallest exact missing statement is

\[
\boxed{
 b_r=b_s=0,\quad \Xi_r=0,\quad 0<r<s<p
 \Longrightarrow
 [T^{p-1}]R_p(T)C_{r,s}(T)\ne0.
}
\tag{RB-EIS}
\]

The hypothesis `Xi_r=0` is indispensable: at the actual prime `p=41`,

\[
b_{10}=b_{30}=0,\qquad \Xi_{10}=\Xi_{30}=7,
\]

and the same block coefficient is exactly zero.

The strongest surviving Frobenius route is **extension-valued Cartier**.  If `K(t)=sum kappa_n t^n` and `mathcal E=K/E`, then

\[
\boxed{
\mathcal C_p\!\left(\mathcal E\Omega_m\frac{dq}{q}\right)
=\mathcal C_p\!\left(K(t)t^{-m}\frac{dt}{t}\right)
=\sum_{j\ge0}\kappa_{m+pj}t^j\frac{dt}{t}.
}
\tag{EC_m}
\]

A common row kills only the constant term of `(EC_m)`, not the whole section.  Q7702's finite-terminal-period compensation does not rule out rigidity of this full coefficient-index family.  A Dwork/Frobenius law for the **inhomogeneous extension** `(E,K)`, rather than another homogeneous determinant, is now the cleanest unresolved route.

No reflected-depth law is used anywhere below.

---

## 1. Endpoint and unit ledger

Fix `p>=7`.  Every `n` with `0<n<p` and hence every `n^3` is a unit in `F_p`.

If `b_n=0`, then the adjacent homogeneous values are nonzero whenever they occur.  Two consecutive zeros in a nonzero second-order solution would propagate backwards to `b_0=0`, contradicting `b_0=1`.

Hasse reciprocity is

\[
\boxed{b_{p-1-j}=b_j.}
\tag{1.1}
\]

Therefore

\[
b_{p-1}=1,\qquad b_{p-2}=5,
\]

so every Hasse zero lies in `1,...,p-3`.  If `r<s` are Hasse zeros, then

\[
r^3b_{r-1},\quad s^3b_{s-1},\quad b_{r+1},\quad b_{s-1}
\]

are all units.  These are exactly the endpoint factors used below.

---

## 2. The actual level-six source as a modular residue

Use

\[
t(\tau)=\left(\frac{\eta(\tau)\eta(6\tau)}
{\eta(2\tau)\eta(3\tau)}\right)^{12},
\]

\[
E(\tau)=\frac{\eta(2\tau)^7\eta(3\tau)^7}
{\eta(\tau)^5\eta(6\tau)^5}=F(t(\tau)),
\]

and

\[
D_q=q\frac d{dq},\qquad H=D_q\log t,\qquad \Psi=EH.
\]

The exact level-six identity is

\[
\boxed{H^2=E^2\Delta(t)},\qquad \Delta(t)=1-34t+t^2.
\tag{2.1}
\]

At `q=0`, `H/E=1+O(q)`, so the selected square root gives

\[
H=E\sqrt\Delta,
\qquad
\boxed{\Psi=E^2\sqrt\Delta=\frac1{g(t)}}.
\tag{2.2}
\]

Hence

\[
g(t)\frac{dt}{t}
=g(t)H\frac{dq}{q}
=\frac1E\frac{dq}{q}.
\]

Residue extraction yields the first exact answer to the coefficient-index question:

### Theorem 2.1 — modular residue formula

For every `m>=0`,

\[
\boxed{
 g_m=[t^m]g(t)
 =\operatorname{CT}_q\frac{t^{-m}}E.
}
\tag{2.3}
\]

This is a modular residue indexed by the coefficient label `m`; it is not evaluation of a modular form at the field element `m`.

Expanding the eta quotient gives

\[
\frac{t^{-m}}E
=q^{-m}\prod_{n\ge1}
\frac{(1-q^{2n})^{12m-7}(1-q^{3n})^{12m-7}}
{(1-q^n)^{12m-5}(1-q^{6n})^{12m-5}}.
\]

Therefore

\[
\boxed{
 g_m=[q^m]\prod_{1\le n\le m}
\frac{(1-q^{2n})^{12m-7}(1-q^{3n})^{12m-7}}
{(1-q^n)^{12m-5}(1-q^{6n})^{12m-5}}.
}
\tag{2.4}
\]

Only `n<=m` can contribute to `[q^m]`, so (2.4) is a finite exact eta-product coefficient formula before reduction modulo `p`.

---

## 3. Franel/elliptic pullback and a finite-field Cartier formula

Let

\[
f_n=\sum_{k=0}^n\binom nk^3,
\qquad h(x)=\sum_{n\ge0}f_nx^n,
\]

and

\[
\phi(x)=\frac{x(1-8x)}{1+x}.
\]

The Apéry-Franel pullback used in the repository is

\[
F(\phi(x))=(1+x)h(x)^2.
\tag{3.1}
\]

Also

\[
\Delta(\phi(x))
=\left(\frac{1-16x-8x^2}{1+x}\right)^2,
\quad
\phi'(x)=\frac{1-16x-8x^2}{(1+x)^2}.
\tag{3.2}
\]

Thus

\[
\boxed{
 g(\phi(x))\phi'(x)
 =\frac1{(1+x)^3h(x)^4}.
}
\tag{3.3}
\]

Taking the residue of `g(t)t^{-m-1}dt` after the substitution `t=phi(x)` gives

\[
\boxed{
 g_m=[x^m](1+x)^{m-2}(1-8x)^{-m-1}h(x)^{-4}.
}
\tag{3.4}
\]

This is an explicit Franel-hypergeometric expression for the actual source.

For the first characteristic-`p` block set

\[
H_p(x)=\sum_{j=0}^{p-1}f_jx^j.
\]

Franel Lucas/Dwork gives

\[
h(x)=H_p(x)h(x^p)\pmod p.
\tag{3.5}
\]

For `m<p`, `h(x^p)^{-4}=1+O(x^p)`, so

\[
\boxed{
 g_m=[x^m](1+x)^{m-2}(1-8x)^{-m-1}H_p(x)^{-4}.
}
\tag{3.6}
\]

There is a cleaner first-block Frobenius form.  Put

\[
W(x)=\frac1{(1+x)^3h(x)^4},
\qquad
\rho_p(x)=\frac{(1+x)^{3p-3}}{H_p(x)^4}.
\]

Equation (3.5) implies

\[
W(x)=\rho_p(x)W(x^p).
\tag{3.7}
\]

Using `phi(x^p)=phi(x)^p`,

\[
W(x)\phi(x)^{-m-1}
=\rho_p(x)\phi(x)^{p-m-1}
\left(\frac W\phi\right)(x^p).
\]

Since `(W/phi)(y)=y^{-1}+O(1)` and the first factor is regular for `0<=m<p`, the only possible contribution to the `x^{-1}` residue is its `x^{p-1}` coefficient times `x^{-p}`.  Hence:

### Theorem 3.1 — first-block Franel-Cartier source formula

\[
\boxed{
 g_m=[x^{p-1}]\rho_p(x)\phi(x)^{p-m-1},
\qquad 0\le m<p.
}
\tag{3.8}
\]

The repository's exact rational pullback

\[
A_p(\phi(x))=(1+x)^{1-p}H_p(x)^2,
\qquad A_p(t)=\sum_{j=0}^{p-1}b_jt^j,
\]

also gives

\[
\rho_p(x)=\frac{(1+x)^{p-1}}{A_p(\phi(x))^2}.
\tag{3.9}
\]

Equations (3.6)-(3.8) keep the coefficient index separate from the elliptic fiber variable `x`.

---

## 4. Homogeneous Cartier diagonalization

Define the weight-four Faber family

\[
\Omega_m=\Psi t^{-m}.
\]

As a logarithmic differential,

\[
\Omega_m\frac{dq}{q}
=E(t)t^{-m}\frac{dt}{t}.
\tag{4.1}
\]

The Apéry Dwork factorization is

\[
E(t)=A_p(t)E(t^p),
\qquad A_p(t)=\sum_{j=0}^{p-1}b_jt^j.
\tag{4.2}
\]

Cartier on logarithmic differentials is coordinate invariant and semilinear:

\[
\mathcal C_p(f(t^p)\omega)=f(t)\mathcal C_p(\omega).
\]

Therefore

\[
\mathcal C_p\!\left(E(t)t^{-m}\frac{dt}{t}\right)
=E(t)\mathcal C_p\!\left(A_p(t)t^{-m}\frac{dt}{t}\right).
\]

For `0<=m<p`, the exponents `j-m`, `0<=j<p`, lie strictly between `-p` and `p`.  The only multiple of `p` is zero, at `j=m`.  Hence

\[
\mathcal C_p\!\left(A_p(t)t^{-m}\frac{dt}{t}\right)
=b_m\frac{dt}{t}.
\]

Thus:

### Theorem 4.1 — Faber-Cartier eigenbasis

\[
\boxed{
\mathcal C_p\!\left(\Omega_m\frac{dq}{q}\right)
=b_m\Psi\frac{dq}{q},
\qquad 0\le m<p.
}
\tag{4.3}
\]

Equivalently `U_p Omega_m=b_m Psi` in q-expansion notation.

If `b_m=0`, then the **entire** homogeneous Cartier image is zero.  Thus every Hasse row, common or not, is indistinguishable after this projection.  Iterating homogeneous Cartier cannot recover `Xi_m` once the first image is zero.

A dedicated q-expansion audit checked (4.3) for every row at `p=7,11,17,19,41` through four extracted `U_p` coefficients; GitHub Actions run `31575567966`, job `94046785277`, had no failures.  This computation checks the theorem but is not its proof.

---

## 5. Fixed Eichler boundary and the level-six Bol identity

Let

\[
\theta=t\frac d{dt},
\]

\[
L=\theta^3-t(34\theta^3+51\theta^2+27\theta+5)
+t^2(\theta+1)^3.
\tag{5.1}
\]

Equivalently,

\[
L=\Delta\theta^3+(3t^2-51t)\theta^2
 +(3t^2-27t)\theta+(t^2-5t).
\]

### Lemma 5.1 — gauge/Bol identity

For every Laurent series `f(t)`,

\[
\boxed{
D_q^3\!\left(\frac fE\right)=\Psi Lf.
}
\tag{5.2}
\]

#### Proof

Write `e=theta E/E`.  Directly expanding `D_q=H theta` gives

\[
E D_q^3E^{-1}
=H^3(\theta-e)^3
+3H^2(\theta H)(\theta-e)^2
+H\theta(H\theta H)(\theta-e).
\tag{5.3}
\]

Use `H^2=E^2 Delta`.  The `theta^3` coefficient becomes `E^2H Delta`; differentiating that identity gives the `theta^2` coefficient, while the remaining two coefficient identities reduce to `LE=0` and its logarithmic derivative.  Hence

\[
E D_q^3E^{-1}=E^2HL,
\]

which is equivalent to (5.2). ∎

For

\[
\lambda_m=\frac{t^{-m}}E,
\]

use `P(-m)=-P(m-1)` to obtain

\[
\boxed{
D_q^3\lambda_m
=-m^3\Omega_m+P(m-1)\Omega_{m-1}
-(m-1)^3\Omega_{m-2}.
}
\tag{5.4}
\]

Now define

\[
\mathcal M_4=
\frac{-3E_4(\tau)+4E_4(2\tau)-9E_4(3\tau)+108E_4(6\tau)}{20}.
\]

The exact level-six identity is

\[
\boxed{\mathcal M_4=(5-t)\Psi.}
\tag{5.5}
\]

One modular proof is short.  The cusp orders of `t` at the four denominator cusps `c=1,2,3,6` are

\[
(1,-1,-1,1),
\]

and those of `E` are

\[
(0,1,1,0).
\]

The logarithmic derivative `H` is holomorphic of weight two, so `(5-t)Psi` is holomorphic of weight four on `Gamma_0(6)`.  The index is `12`, hence the weight-four Sturm bound is `4`.  Both sides have expansion

\[
5-36q-276q^2-1116q^3-2196q^4+O(q^5),
\]

proving (5.5).  The Q7708 verifier checks farther through `O(q^14)`.

The fixed Eichler normalization is

\[
\mathcal E=\frac KE,
\qquad
D_q^3\mathcal E=\mathcal M_4-5.
\tag{5.6}
\]

Applying (5.2) to `f=K` and then (5.5) gives

\[
LK=\frac{\mathcal M_4-5}{\Psi}=5-t-5g.
\tag{5.7}
\]

Thus the coefficient of `t` is `-1-5g_1=-36`, while every coefficient from row `2` onward is `-5g_n`.  The modular identity reproduces exactly the canonical inhomogeneous recurrence and its one-unit left boundary defect.

Finally, for `m>=1`,

\[
\boxed{
\operatorname{CT}_q(\mathcal M_4\lambda_m)=-\delta_{m,1}.
}
\tag{5.8}
\]

Indeed, `M4 lambda_m=(5-t)Ht^{-m}`, and

\[
\operatorname{CT}_q(Ht^k)
=\operatorname{Res}_t t^k\frac{dt}{t}=\delta_{k,0}.
\]

So the simplest fixed-boundary modular-symbol functional is identically zero on all Faber labels `m>=2`.

---

## 6. A canonical two-return Bol/Eichler cocycle collapses to Green

Let `0<r<s<p` be Hasse zeros and define

\[
Y_{r,s}
=\frac1E\sum_{m=r+1}^{s}b_{m-1}t^{-m}
=\sum_{m=r+1}^{s}b_{m-1}\lambda_m.
\tag{6.1}
\]

Insert (5.4).  Every interior `Omega_k` coefficient is

\[
-k^3b_{k-1}+P(k)b_k-(k+1)^3b_{k+1}=0.
\]

At the left endpoint, `b_r=0` gives

\[
-(r+1)^3b_{r+1}=r^3b_{r-1}.
\]

The right endpoint is direct.  Hence:

### Theorem 6.1 — two-return Bol bridge

\[
\boxed{
D_q^3Y_{r,s}
=r^3b_{r-1}\Omega_r-s^3b_{s-1}\Omega_s.
}
\tag{6.2}
\]

Both endpoint coefficients are units by Section 1.

This is a genuine source-specific two-return cocycle and is not one of Q7702's arbitrary homogeneous adjoints.  But its **canonical Eichler pairing is not independent**.

First,

\[
\operatorname{CT}_q(\mathcal E\Omega_n)=\kappa_n.
\tag{6.3}
\]

Thus pairing (6.2) with `mathcal E` gives the endpoint expression

\[
r^3b_{r-1}\kappa_r-s^3b_{s-1}\kappa_s.
\]

On the other hand, three integrations by parts in the constant-term pairing give

\[
\operatorname{CT}(\mathcal E D_q^3Y)
=-\operatorname{CT}((\mathcal M_4-5)Y).
\]

All labels in `Y` are at least `2`, so (5.8) kills the `M4` term.  Equation (2.3) then yields

\[
\boxed{
 r^3b_{r-1}\kappa_r-s^3b_{s-1}\kappa_s
 =5\sum_{m=r+1}^sb_{m-1}g_m.
}
\tag{6.4}
\]

At Hasse endpoints the left side is `Xi_r-Xi_s`.  Therefore (6.4) is exactly

\[
\Xi_s-\Xi_r=-5\sum_{m=r+1}^sb_{m-1}g_m,
\]

the existing Green/Duhamel shooting equation.  This canonical Bol/Eichler construction must be rejected as a second independent constraint.

---

## 7. Current-main inhomogeneous reflection cocycle

While this audit was running, current `main` added

```text
problems/3.2/research/scripts/q7709_inhomogeneous_hasse_reflection.sage
```

and an exact audit of the actual source.  Its algebra can be stated cleanly here.

Put `N=p-1` and

\[
(Jy)_r=y_{N-r}.
\]

For the row operator

\[
\mathcal L_r y
=(r+1)^3y_{r+1}-P(r)y_r+r^3y_{r-1},
\]

`P(-r-1)=-P(r)` gives, on interior rows,

\[
\boxed{\mathcal L J=-J\mathcal L.}
\tag{7.1}
\]

Q7702's stronger homogeneous reflection theorem says the full two-dimensional homogeneous kernel is pointwise `J`-fixed.  In particular there is no anti-invariant homogeneous direction hidden here.

For the **canonical** `kappa`, define

\[
C=J\kappa-\kappa.
\tag{7.2}
\]

Since

\[
\mathcal L_r\kappa=-5g_{r+1},
\]

(7.1) gives the exact inhomogeneous cocycle

\[
\boxed{
\mathcal L_rC=5\bigl(g_{p-r}+g_{r+1}\bigr).
}
\tag{7.3}
\]

This uses both the actual source and the fixed canonical boundary through `kappa`.

### Why (7.3) does not give a general second return equation

If `r` is common, then `kappa_r=0`, but

\[
C_r=\kappa_{N-r},
\]

which is not forced to vanish.  The exact current-main run gives, for example,

\[
p=17,\ r=13:\quad C_{13}=\kappa_3=2,
\]

and

\[
p=2237,\ r=492:\quad C_{492}=\kappa_{1744}=220.
\]

Thus two hypothetical common rows `r<s` do not make the endpoints of `C` zero unless their reflected values happen to be common too.

### The only automatic zero-endpoint case collapses to Green

Suppose specifically

\[
s=N-r=p-1-r
\]

and both `r,s` are common.  Then

\[
C_r=\kappa_s-\kappa_r=0,
\qquad
C_s=\kappa_r-\kappa_s=0.
\]

The Green shooting condition for (7.3) is therefore

\[
\sum_{u=r}^{s-1}b_u\bigl(g_{p-u}+g_{u+1}\bigr)=0.
\tag{7.4}
\]

The second half equals the original period

\[
S_{r,s}=\sum_{m=r+1}^{s}b_{m-1}g_m.
\]

For the first half set `m=p-u`.  Since `s=p-1-r`, Hasse reciprocity gives

\[
\sum_{u=r}^{s-1}b_ug_{p-u}
=\sum_{m=r+2}^{s+1}b_{m-1}g_m.
\]

The difference from `S_{r,s}` is

\[
b_sg_{s+1}-b_rg_{r+1}=0.
\]

Hence (7.4) is exactly

\[
\boxed{2S_{r,s}=0.}
\tag{7.5}
\]

Since `p>=7`, `2` is a unit, so this is again the original Green condition and supplies no second equation.

The current-main Q7709 run (`31575925992`, job `94047902597`) also shows that the exact reciprocal-tail source in (7.3) is generically dense in the tested primes: support `15/15` at `p=17`, `39/39` at `p=41`, and `2235/2235` at `p=2237`.  This density is finite evidence about complexity, not a theorem of noncompressibility.

Thus the natural inhomogeneous reflection cocycle is useful structure but does not close EIS-2RET.

---

## 8. Reflected Hasse-block form of the actual return period

For the first characteristic-`p` source block put

\[
R_p(T)=\frac{\Delta(T)^{(p-1)/2}}{A_p(T)^2},
\qquad g_m=[T^m]R_p(T),\quad m<p.
\tag{8.1}
\]

For Hasse zeros `r<s`, define

\[
S_{r,s}=\sum_{m=r+1}^sb_{m-1}g_m
\tag{8.2}
\]

and

\[
C_{r,s}(T)=\sum_{m=r+1}^{s}b_{m-1}T^{p-1-m}.
\tag{8.3}
\]

Coefficient convolution gives

\[
\boxed{
S_{r,s}=[T^{p-1}]R_p(T)C_{r,s}(T).
}
\tag{8.4}
\]

Use Hasse reciprocity with `k=p-m`:

\[
C_{r,s}(T)
=\sum_{k=p-s}^{p-r-1}b_kT^{k-1}.
\tag{8.5}
\]

Because `b_r=b_s=0`, the reflected coefficients immediately outside the block are also zero:

\[
b_{p-s-1}=b_s=0,
\qquad
b_{p-r-1}=b_r=0.
\]

Factoring the lower exponent gives

\[
\boxed{
C_{r,s}(T)
=T^{p-s-1}B_{r,s}(T),
\qquad
B_{r,s}(T)=\sum_{j=0}^{s-r-2}b_{s-1-j}T^j.
}
\tag{8.6}
\]

The omitted next coefficient would be `b_r=0`.  The two active endpoints are

\[
B_{r,s}(0)=b_{s-1}\ne0,
\qquad
[T^{s-r-2}]B_{r,s}=b_{r+1}\ne0.
\tag{8.7}
\]

Thus

\[
\boxed{\deg B_{r,s}=s-r-2}
\]

with unit constant and leading coefficients.

Combining (3.8) with (8.2) gives the equivalent Franel/elliptic coefficient

\[
\boxed{
S_{r,s}
=[x^{p-1}]\rho_p(x)
\sum_{m=r+1}^{s}b_{m-1}\phi(x)^{p-m-1}.
}
\tag{8.8}
\]

And (8.4) is the residue

\[
S_{r,s}
=\operatorname{Res}_{T=0}
\frac{R_p(T)C_{r,s}(T)}{T^p}\,dT.
\tag{8.9}
\]

The differential has zero residue at infinity, hence

\[
\boxed{
S_{r,s}
=-\sum_{\alpha:A_p(\alpha)=0}
\operatorname{Res}_{T=\alpha}
\frac{\Delta(T)^{(p-1)/2}C_{r,s}(T)}
{T^pA_p(T)^2}\,dT.
}
\tag{8.10}
\]

Via the Franel pullback, (8.10) is a supersingular fourth-order residue/jet trace against `H_p(x)^4`.  This is an exact global source pairing, but the residue theorem alone gives no nonvanishing; `p=41` is an actual cancellation example for a Hasse pair.

---

## 9. EIS-2RET is exactly reflected-block nonvanishing

If `r` is common and `s>r` is another Hasse zero, then `Xi_r=0` and (6.4) gives

\[
\Xi_s=-5S_{r,s}.
\]

Because `5` is a unit,

\[
\boxed{s\text{ is common}\iff S_{r,s}=0.}
\tag{9.1}
\]

Therefore uniqueness is equivalent to the following source-specific statement.

### `(RB-EIS)` — distinguished-Eichler reflected-block nonvanishing

For every prime `p>=7` and every `0<r<s<p`,

\[
\boxed{
 b_r=b_s=0,\quad \Xi_r=0
 \Longrightarrow
 [T^{p-1}]R_p(T)C_{r,s}(T)\ne0.
}
\tag{9.2}
\]

Using (8.6), this is equivalently

\[
\boxed{
[T^s]R_p(T)B_{r,s}(T)\ne0.
}
\tag{9.3}
\]

Using (8.8), it is equivalently

\[
\boxed{
[x^{p-1}]\rho_p(x)
\sum_{m=r+1}^{s}b_{m-1}\phi(x)^{p-m-1}\ne0.
}
\tag{9.4}
\]

I do not know how to prove (9.2)-(9.4) uniformly.  This is the precise missing implication.

It is stronger localization than simply restating Duhamel: the source is the actual Cartier factor `R_p` (or `rho_p`), and the return interval is a contiguous reflected Hasse block with two proved unit endpoints.

---

## 10. Why `Xi_r=0` is essential

At `p=41`, the actual canonical source has

\[
b_{10}=b_{30}=0,
\qquad
\Xi_{10}=\Xi_{30}=7.
\]

Therefore

\[
\boxed{S_{10,30}=0.}
\]

The Q7708 verifier checks simultaneously

\[
[T^{40}]R_{41}(T)C_{10,30}(T)=0,
\qquad
\deg B_{10,30}=18,
\]

with unit endpoints.

Thus all statements obtained by deleting the distinguished Eichler condition are false, including

```text
Every Hasse pair has nonzero Green/Cartier block pairing.
Every reflected Hasse block pairs nontrivially with R_p.
Hasse/p-curvature geometry alone prevents two-return cancellation.
```

The fixed boundary `Xi_r=0` is genuinely additional arithmetic information.

---

## 11. Targeted exact computations

The committed Q7708 verifier is

```text
problems/3.2/research/scripts/q7708_eis2ret_verify.sage
```

GitHub Actions run `31575840359`, job `94047634657`, completed successfully with

```text
RETURN_PAIR 41 10 30 XI 7 7 S 0 BLOCK_DEG 18
RETURN_PAIR 2237 492 1744 XI 0 1314 S 632 BLOCK_DEG 1250
RETURN_PAIR 2237 492 2213 XI 0 1763 S 1437 BLOCK_DEG 1719
P2237_LATER_SHOOTING 632 1437
PULLBACK_CARTIER 7 (0, 1, 2, 5, 6)
PULLBACK_CARTIER 17 (0, 1, 2, 3, 13, 15, 16)
PULLBACK_CARTIER 41 (0, 1, 2, 10, 30, 39, 40)
PULLBACK_CARTIER 181 (0, 1, 2, 19, 47, 133, 161, 179, 180)
LEVEL6_QSERIES 14 M4_OVER_PSI=5-t CT_G_OK
Q7708_EIS2RET_VERIFY=PASS
```

For `p=2237`, the exact Hasse-zero set is

\[
\{23,492,1744,2213\},
\]

and

\[
(\Xi_{23},\Xi_{492},\Xi_{1744},\Xi_{2213})
=(367,0,1314,1763).
\]

Thus the two later actual Hasse returns from the known common row satisfy

\[
\boxed{
S_{492,1744}=632,
\qquad
S_{492,2213}=1437
\pmod{2237},
}
\]

both nonzero.

The separate homogeneous Cartier audit, run `31575567966` / job `94046785277`, checked every coefficient label at `p=7,11,17,19,41` through four `U_p` coefficients and found no failure of (4.3).

The concurrent current-main Q7709 audit, run `31575925992` / job `94047902597`, independently verifies the actual-source reflection cocycle (7.3) at `p=17,41,181,2237` and records the nonzero reflected endpoint values quoted in Section 7.

All of these are finite checks of exact identities or finite examples; none is promoted to a proof of uniqueness.

---

## 12. Scoped no-go ledger

The following natural mechanisms are now rigorously insufficient in their straightforward form.

### 12.1 Homogeneous two-return adjoints

Q7702: the two-endpoint Dirichlet test space is rank one, so every homogeneous adjoint determinant is the Green period.

### 12.2 Canonical linear Bol/Eichler bridge

Section 6 is genuinely source-specific and uses the fixed Eichler boundary, but its canonical pairing is exactly Green/Duhamel.

### 12.3 One-step or iterated homogeneous Cartier

Theorem 4.1 kills the whole Faber form at every Hasse zero.  `Xi_m` has already been discarded.

### 12.4 Fixed `M4` modular-symbol functional

Equation (5.8) is `-delta_{m,1}`.  It contains no high-label information.

### 12.5 Inhomogeneous reflection cocycle

For a generic common pair its endpoint values are unconstrained reflected `kappa` values.  For a reflected common pair its zero-endpoint shooting condition is exactly `2S_{r,s}=0`, hence Green again.

### 12.6 Bare rational reciprocity or residue theorem

Q7702 already refutes naive coefficient reflection.  The refined residue trace (8.10) is exact but can cancel; `p=41` proves that cancellation occurs on an actual Hasse pair.

These no-go statements do **not** rule out an inhomogeneous Frobenius extension matrix, a nonlinear modular-symbol/intersection pairing, or a coefficient-index theorem using the full extension class rather than finitely many terminal periods.

---

## 13. The surviving extension-valued Cartier target

The homogeneous Faber form loses exactly the datum we need.  The product with the canonical Eichler extension does not:

\[
\mathcal E\Omega_m\frac{dq}{q}
=K(t)t^{-m}\frac{dt}{t}.
\]

For `0<=m<p`, Cartier therefore gives the exact coefficient section

\[
\boxed{
\mathcal C_p\!\left(\mathcal E\Omega_m\frac{dq}{q}\right)
=\sum_{j\ge0}\kappa_{m+pj}t^j\frac{dt}{t}.
}
\tag{13.1}
\]

At a common row,

\[
b_m=0,\qquad \kappa_m=0,
\]

so the homogeneous Cartier eigenvalue and the **constant term** of (13.1) vanish.  But there is no reason for the higher coefficients `kappa_{m+pj}` to vanish.

This is precisely why (13.1) escapes the earlier no-go results:

- it is not a homogeneous adjoint;
- it is not determined by finitely many terminal periods;
- it keeps the actual level-six source and fixed Eichler normalization;
- it is indexed by the coefficient label `m`, not by a geometric fiber value.

The next load-bearing theorem should therefore be a Dwork/Frobenius relation for the rank-two extension `(E,K)` which controls the family

\[
\{\mathcal C_p(Kt^{-m}dt/t):0\le m<p\}
\]

by a bounded-conductor extension class.  A sufficient form would imply that on the Hasse locus, the constant term of this extension section can vanish for at most one label.

I do not have that Frobenius-extension law.  This is a much narrower missing object than “find another determinant.”

---

## 14. Final status

No uniform multiplicity theorem is claimed.  The strongest theorem-level progress is the chain

\[
\boxed{
\begin{array}{l}
\text{actual }g_m\text{ has exact level-six residue and Franel-Cartier formulas};\\[1mm]
\mathcal C_p(\Omega_m\,dq/q)=b_m\Psi\,dq/q;\\[1mm]
D_q^3Y_{r,s}=r^3b_{r-1}\Omega_r-s^3b_{s-1}\Omega_s;\\[1mm]
\text{its fixed-Eichler pairing is exactly Green};\\[1mm]
\mathcal L(J\kappa-\kappa)=5(g_{p-r}+g_{r+1})\text{ but adds no return equation};\\[1mm]
S_{r,s}=[T^{p-1}]R_p(T)C_{r,s}(T),\\
C_{r,s}=T^{p-s-1}\sum_{j=0}^{s-r-2}b_{s-1-j}T^j
\text{ with unit endpoints.}
\end{array}}
\]

Thus EIS-2RET is reduced to the explicit nonvanishing `(RB-EIS)`.  The clean homogeneous Cartier/Hecke route, the canonical linear Eichler cocycle, and the natural inhomogeneous reflection cocycle are all proved insufficient.  The remaining credible global mechanism is the **inhomogeneous extension-valued Cartier/Frobenius law for `K`**, or an equivalent supersingular/intersection theorem that genuinely uses the distinguished condition `Xi_r=0`.
