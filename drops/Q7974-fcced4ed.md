ANSWER Q7974 fcced4ed

# HM3 semiprime \(TT^*\): exact primitive-\(p q\ell\) reduction and the remaining five-prime correlation

## Verdict

The ordered-pair operator can be reconstructed exactly, including the factor of two coming from \(R_{q,\ell}=R_{\ell,q}\). After conjugating away the harmless palindromic phases, it has the following concrete form:

```text
primitive Fourier tensor R_q \otimes R_l
        --inverse CRT Fourier transform-->
centered row product g_q(m) g_l(m) on 0 <= m < X^2
        --nonzero p-Fourier transform-->
output block B_p.
```

More precisely, with

\[
g_r(m)=1_{\{m\bmod r\in\mathcal Z_r\}}-z_r,
\qquad z_r=\frac{Z(r)}r,
\qquad M=X^2,
\]

one obtains the exact identity

\[
\boxed{
 B_p(c)=\frac1p\sum_{m<M}G_p(m)e_p(-cm),
 \qquad
 G_p(m)=\sum_{\substack{q<\ell\\q,\ell\ne p}}g_q(m)g_\ell(m),
 \quad c\in\mathbf F_p^\times .
}
\tag{1}
\]

Consequently

\[
\boxed{
 \sum_p\|B_p\|_2^2
 =\sum_p\frac1p\sum_{m,n<M}
 \left(1_{m\equiv n\pmod p}-\frac1p\right)
 G_p(m)G_p(n).
}
\tag{2}
\]

This is already a sharp actual-Apéry reduction: it is a centered congruence correlation of one output prime and four input-prime slots.

The full primitive-frequency calculation gives more structure. For one block \(p,q,\ell\), the pair \((c,a)\), with \(c\in\mathbf F_p^\times\) and \(a\in(\mathbf Z/q\ell\mathbf Z)^\times\), is in bijection with

\[
k=cq\ell-ap\in(\mathbf Z/pq\ell\mathbf Z)^\times.
\]

Thus the block really uses all primitive frequencies \(k/(pq\ell)\), not a fixed shift \(h\). Its exact Hilbert--Schmidt mass is

\[
\boxed{
 \|2T_{p,q\ell}\|_{\rm HS}^2
 =\frac1{p^2q^2\ell^2}
   \sum_{k\in(\mathbf Z/pq\ell\mathbf Z)^\times}
   \left|\mathcal D_M^\circ\!\left(\frac{k}{pq\ell}\right)\right|^2
 =\frac{\Xi(p,q,\ell)}{pq\ell},
}
\tag{3}
\]

where \(\Xi(p,q,\ell)\) is evaluated exactly below and satisfies

\[
\Xi(p,q,\ell)=M\bigl(1+O(X^{-1})\bigr).
\]

Hence one block has Hilbert--Schmidt square \(\asymp M/(pq\ell)\asymp X^{-1}\). Restricting all three primes to active Apéry columns gives

\[
\operatorname{Tr}(\mathsf S^*\mathsf S)
 \ll M\left(\sum_{Z(p)>0}\frac1p\right)^3
 \le M\lambda^3.
\tag{4}
\]

This is the full tensor orthogonality available for free. It does **not** imply the desired \(M X^{o(1)}\lambda^5\) bound. The actual input tensor has squared norm \(O(X^4\lambda^2)\), so the trace inequality gives only

\[
\sum_p\|B_p\|_2^2\ll X^6\lambda^5,
\tag{5}
\]

losing \(X^4\). The direct semiprime additive large sieve gives, for active outputs, only

\[
\sum_p\|B_p\|_2^2\ll X^5\lambda^3,
\tag{6}
\]

and has the wrong density exponent as well as three extra powers of \(X\).

The exact \(\mathsf S^*\mathsf S\) expansion isolates three classes:

* the same-semiprime block, with congruence \(k\equiv k'\pmod p\);
* a shared-prime block \(q\ell,qs\), with congruence \(ks\equiv k'\ell\pmod p\);
* a disjoint block \(q\ell,rs\), with congruence \(krs\equiv k'q\ell\pmod p\).

The last is the genuine five-distinct-prime term. The coefficients are reciprocal evaluations of the actual Apéry Fourier rows, for example

\[
\widehat g_q(-k\overline{p\ell})
\widehat g_\ell(-k\overline{pq})
\overline{
 \widehat g_r(-k'\overline{ps})
 \widehat g_s(-k'\overline{pr})}.
\]

No currently proved input in the repository controls this moving reciprocal correlation. Detecting its congruence by additive characters does not create a classical Kloosterman sum with fixed coefficients; Cauchy and Parseval return exactly the already-recorded \(\lambda^2\) input mass.

Therefore the full-\(pq\ell\) route gives a rigorous reduction but no unconditional saving to \(X^{2+o(1)}\lambda^5\). The missing theorem is the **Apéry five-prime centered-congruence estimate** stated in Section 10 below, equivalently its explicit reciprocal-frequency form in (24).

---

## 1. Prime sets, active outputs, and the exact row norm

Put

\[
\mathcal P_X=\{r\text{ prime}:X<r\le2X\},
\qquad
\mathcal A_X=\{r\in\mathcal P_X:Z(r)>0\}.
\]

All input tensors with an inactive prime vanish. For the HM3 output-energy Cauchy step one only needs output primes in \(\mathcal A_X\), because \(R_p=0\) when \(Z(p)=0\). To cover both conventions, let

\[
\mathcal O_X\subseteq\mathcal P_X
\]

be the output set. The literal all-prime version is \(\mathcal O_X=\mathcal P_X\), while the relevant active-output version is \(\mathcal O_X=\mathcal A_X\).

For \(r\in\mathcal P_X\), define the complete-period centered row

\[
g_r(x)=1_{\mathcal Z_r}(x)-z_r,
\qquad x\in\mathbf F_r,
\qquad z_r=\frac{Z(r)}r.
\]

Its unnormalized additive Fourier transform is

\[
\widehat g_r(a)=\sum_{x\bmod r}g_r(x)e_r(-ax).
\]

Since \(\sum_xg_r(x)=0\), one has \(\widehat g_r(0)=0\), while for \(a\ne0\), \(\widehat g_r(a)\) is the ordinary Fourier transform of the Apéry zero set. Parseval gives exactly

\[
\boxed{
 \eta_r:=\sum_{a\in\mathbf F_r^\times}|\widehat g_r(a)|^2
 =rZ(r)-Z(r)^2
 =r^2z_r(1-z_r).
}
\tag{7}
\]

Write

\[
\rho_r=z_r(1-z_r)\le z_r,
\qquad
\lambda=\sum_{r\in\mathcal P_X}z_r.
\]

The palindromic rotation used in the repository writes

\[
R_r(a)=\omega_r(a)\widehat g_r(a),
\qquad |\omega_r(a)|=1,
\]

with \(R_r(a)\) real after pairing conjugate frequencies. Let \(D_r\) be the corresponding diagonal unitary. Under CRT, \(D_{q\ell}=D_q\otimes D_\ell\). If \(T_{p,q\ell}\) is written in the paired real basis, then

\[
\widetilde T_{p,q\ell}=D_p^*T_{p,q\ell}D_{q\ell}
\]

is the complex Fourier-basis block below, and

\[
\widetilde B_p=D_p^*B_p.
\]

Therefore

\[
\|\widetilde B_p\|_2=\|B_p\|_2.
\]

All calculations may consequently be made with \(\widehat g_r\) and adjoints instead of \(R_r\) and transposes, with no change in the energy.

For \(v=q\ell\), CRT gives

\[
\widehat g_{q,\ell}(a)
 =\widehat g_q(a\overline\ell\bmod q)
  \widehat g_\ell(a\overline q\bmod\ell),
\qquad a\in(\mathbf Z/v\mathbf Z)^\times,
\tag{8}
\]

and hence

\[
\boxed{
 \|R_{q,\ell}\|_2^2
 =\eta_q\eta_\ell.
}
\tag{9}
\]

The normalized total input mass used by dual large sieve is

\[
\begin{aligned}
\mathcal E_2
&:=\sum_{q<\ell}\frac{\|2R_{q,\ell}\|_2^2}{q^2\ell^2}\\
&=4\sum_{q<\ell}\rho_q\rho_\ell\\
&=2\left[
 \left(\sum_r\rho_r\right)^2-\sum_r\rho_r^2
 \right]
\le2\lambda^2.
\end{aligned}
\tag{10}
\]

Two elementary active-prime bounds will be used repeatedly:

\[
\boxed{
 |\mathcal A_X|\le2X\lambda,
 \qquad
 \sigma_X:=\sum_{r\in\mathcal A_X}\frac1r\le\lambda.
}
\tag{11}
\]

Indeed, for \(r\in\mathcal A_X\), \(z_r\ge1/r\ge1/(2X)\).

---

## 2. The exact ordered symmetric operator

Let \(I=\{0,1,\ldots,M-1\}\). For squarefree \(d\), let

\[
\mathcal H_d=\ell^2((\mathbf Z/d\mathbf Z)^\times)
\]

and define the unit-row matrix

\[
U_d(a,m)=d^{-1/2}e_d(-am),
\qquad a\in(\mathbf Z/d\mathbf Z)^\times,\quad m\in I.
\tag{12}
\]

In the complex basis the block is

\[
T_{p,v}=\frac1{2\sqrt{pv}}U_pU_v^*,
\qquad v=q\ell,\quad p\nmid v.
\tag{13}
\]

The question uses ordered pairs and symmetric inputs:

\[
B_p=\sum_{\substack{q\ne\ell\\p\notin\{q,\ell\}}}
T_{p,q\ell}R_{q,\ell},
\qquad R_{q,\ell}=R_{\ell,q}.
\]

There is no normalization ambiguity. Pairing the two orders gives

\[
\boxed{
 B_p
 =2\sum_{\substack{q<\ell\\p\notin\{q,\ell\}}}
   T_{p,q\ell}R_{q,\ell}.
}
\tag{14}
\]

Equivalently, one may use unordered inputs \(x_{q\ell}=2R_{q,\ell}\) and the original block \(T_{p,q\ell}\), or unordered inputs \(R_{q,\ell}\) and the doubled block

\[
S_{p,q\ell}:=2T_{p,q\ell}
 =\frac1{\sqrt{pq\ell}}U_pU_{q\ell}^*.
\tag{15}
\]

I use the latter. Define

\[
\mathsf S:\bigoplus_{q<\ell}\mathcal H_{q\ell}
 \longrightarrow
 \bigoplus_{p\in\mathcal O_X}\mathcal H_p,
\]

by

\[
(\mathsf S Y)_p
 =\sum_{\substack{q<\ell\\p\notin\{q,\ell\}}}
 S_{p,q\ell}Y_{q\ell}.
\tag{16}
\]

Then the actual input is \(Y_{q\ell}=R_{q,\ell}\), and

\[
\boxed{
 \sum_{p\in\mathcal O_X}\|B_p\|_2^2
 =\|\mathsf S R\|_2^2
 =\langle R,\mathsf S^*\mathsf S R\rangle.
}
\tag{17}
\]

---

## 3. Exact inverse CRT transform: the operator becomes a pair field

The crucial point is that the tensor \(R_{q,\ell}\) is not an arbitrary vector. It is the primitive Fourier transform of the product \(g_qg_\ell\).

Because both factors have zero Fourier coefficient at frequency zero, the Fourier transform of \(g_qg_\ell\), viewed modulo \(q\ell\), is supported exactly on the unit frequencies. Fourier inversion therefore gives

\[
\begin{aligned}
(U_{q\ell}^*\widehat g_{q,\ell})(m)
&=\frac1{\sqrt{q\ell}}
  \sum_{a\in(\mathbf Z/q\ell\mathbf Z)^\times}
  \widehat g_{q,\ell}(a)e_{q\ell}(am)\\
&=\sqrt{q\ell}\,g_q(m)g_\ell(m).
\end{aligned}
\tag{18}
\]

Applying (15),

\[
\begin{aligned}
(S_{p,q\ell}\widehat g_{q,\ell})(c)
&=\frac1{\sqrt{pq\ell}}
 U_p\bigl(\sqrt{q\ell}\,g_qg_\ell\bigr)(c)\\
&=\frac1p\sum_{m<M}g_q(m)g_\ell(m)e_p(-cm).
\end{aligned}
\tag{19}
\]

Summing over \(q<\ell\) proves (1). This identity is exact on the finite interval; there is no completion of the interval and no CRT boundary estimate hidden in it.

For later use, put

\[
A_{p;q\ell}(u)
 =\sum_{\substack{m<M\\m\equiv u\pmod p}}
 g_q(m)g_\ell(m),
\qquad
A_{q\ell}=\sum_{m<M}g_q(m)g_\ell(m).
\tag{20}
\]

Nonzero-frequency Parseval modulo \(p\) gives

\[
\begin{aligned}
\langle S_{p,q\ell}R_{q,\ell},
        S_{p,rs}R_{r,s}\rangle
&=\frac1p\sum_{u\bmod p}
 A_{p;q\ell}(u)\overline{A_{p;rs}(u)}\\
&\quad-\frac1{p^2}A_{q\ell}\overline{A_{rs}}.
\end{aligned}
\tag{21}
\]

Equivalently, with

\[
\Delta_p(m,n)=1_{m\equiv n\pmod p}-\frac1p,
\]

one has

\[
\boxed{
\mathcal C_p(q,\ell;r,s)
:=\langle S_{p,q\ell}R_{q,\ell},
        S_{p,rs}R_{r,s}\rangle
=\frac1p\sum_{m,n<M}
\Delta_p(m,n)
 g_q(m)g_\ell(m)g_r(n)g_s(n).
}
\tag{22}
\]

Thus

\[
\boxed{
\sum_{p\in\mathcal O_X}\|B_p\|_2^2
=\sum_{p\in\mathcal O_X}
 \sum_{\substack{q<\ell,\ r<s\\p\notin\{q,\ell,r,s\}}}
 \mathcal C_p(q,\ell;r,s).
}
\tag{23}
\]

This is the exact centered five-slot correlation. It is not a fixed-gap statement: \(\Delta_p\) simultaneously sums every difference divisible by \(p\), with the zero Fourier mode removed exactly.

---

## 4. Exact \(\mathsf S^*\mathsf S\) kernel with the centered Dirichlet kernel

Let

\[
\mathcal D_M^\circ(t)
 =\sum_{m=0}^{M-1}e\!\left((m-(M-1)/2)t\right)
 =\frac{\sin(\pi Mt)}{\sin(\pi t)}.
\]

The centering changes only row phases. In this basis

\[
S_{p,v}(c,a)
 =\frac1{pv}\mathcal D_M^\circ
   \!\left(\frac cp-\frac av\right).
\]

For semiprimes \(v,w\), the block of \(\mathsf S^*\mathsf S\) is

\[
\boxed{
\mathcal K_{v,w}
=\sum_{\substack{p\in\mathcal O_X\\p\nmid vw}}
 S_{p,v}^*S_{p,w}
=\frac1{\sqrt{vw}}U_v
 \left(
  \sum_{\substack{p\in\mathcal O_X\\p\nmid vw}}
  \frac1pU_p^*U_p
 \right)U_w^*.
}
\tag{24}
\]

Entrywise,

\[
\boxed{
\mathcal K_{v,w}(a,b)
=\sum_{\substack{p\in\mathcal O_X\\p\nmid vw}}
 \frac1{p^2vw}
 \sum_{c\in\mathbf F_p^\times}
 \mathcal D_M^\circ\!\left(\frac cp-\frac av\right)
 \overline{
 \mathcal D_M^\circ\!\left(\frac cp-\frac bw\right)}.
}
\tag{25}
\]

The time-side output kernel is equally exact:

\[
(U_p^*U_p)(m,n)
=\frac1p\sum_{c\in\mathbf F_p^\times}e_p(c(m-n))
=1_{m\equiv n\pmod p}-\frac1p.
\tag{26}
\]

Equations (22), (24), and (25) are three exactly equivalent forms of the same \(\mathsf S^*\mathsf S\) operator.

---

## 5. Full primitive-frequency modulus \(pq\ell\)

Fix distinct \(p,q,\ell\), set \(v=q\ell\) and \(P=pv=pq\ell\). For

\[
c\in\mathbf F_p^\times,
\qquad a\in(\mathbf Z/v\mathbf Z)^\times,
\]

define

\[
k=cv-ap\pmod P.
\]

Then

\[
k\bmod p=cv\ne0,
\qquad
k\bmod q=-ap\ne0,
\qquad
k\bmod\ell=-ap\ne0.
\]

Hence \(k\in(\mathbf Z/P\mathbf Z)^\times\). Conversely,

\[
c\equiv k\overline v\pmod p,
\qquad
 a\equiv-k\overline p\pmod v,
\]

so \((c,a)\mapsto k\) is a bijection

\[
\mathbf F_p^\times\times(\mathbf Z/v\mathbf Z)^\times
 \simeq(\mathbf Z/P\mathbf Z)^\times.
\]

The two CRT component frequencies of \(a\) are

\[
a_q\equiv-k\overline{p\ell}\pmod q,
\qquad
a_\ell\equiv-k\overline{pq}\pmod\ell.
\]

Therefore one actual pair column has the exact primitive representation

\[
\boxed{
(S_{p,q\ell}R_{q,\ell})(c)
=\frac1{pq\ell}
 \sum_{\substack{k\in(\mathbf Z/pq\ell\mathbf Z)^\times\\
                  k\equiv cq\ell\pmod p}}
 \mathcal D_M^\circ\!\left(\frac{k}{pq\ell}\right)
 R_q(-k\overline{p\ell})
 R_\ell(-k\overline{pq}).
}
\tag{27}
\]

Here and below each inverse is taken in the modulus of the \(R\)-factor. In the unrotated basis replace \(R\) by \(\widehat g\) and insert the corresponding harmless unit phases.

Now let \(v=q\ell\) and \(w=rs\). Opening the square and using equality of the output frequency \(c\) gives

\[
 k\overline v\equiv k'\overline w\pmod p,
\]

or equivalently

\[
\boxed{
 kw\equiv k'v\pmod p.
}
\tag{28}
\]

Thus the complete reciprocal-frequency form of the energy is

\[
\boxed{
\begin{aligned}
\sum_{p\in\mathcal O_X}\|B_p\|_2^2
&=\sum_{p\in\mathcal O_X}
 \sum_{\substack{v=q\ell,\ w=rs\\q<\ell,\ r<s\\p\nmid vw}}
 \frac1{p^2vw}\\
&\quad\times
 \sum_{\substack{k\in(\mathbf Z/pv\mathbf Z)^\times\\
                   k'\in(\mathbf Z/pw\mathbf Z)^\times\\
                   kw\equiv k'v\ (p)}}
 \mathcal D_M^\circ\!\left(\frac{k}{pv}\right)
 \overline{\mathcal D_M^\circ\!\left(\frac{k'}{pw}\right)}\\
&\quad\times
 R_q(-k\overline{p\ell})R_\ell(-k\overline{pq})
 \overline{R_r(-k'\overline{ps})R_s(-k'\overline{pr})}.
\end{aligned}
}
\tag{29}
\]

This is the requested full primitive-modulus formula.

### 5.1 Same semiprime

If \(v=w=q\ell\), condition (28) is

\[
\boxed{k\equiv k'\pmod p.}
\tag{30}
\]

Equivalently \(k'=k+ph\) modulo \(pq\ell\), with \(h\bmod q\ell\). The frequency diagonal \(k=k'\) is only one part of this block; the other \(h\)'s are genuine self-correlations of the reciprocal Apéry Fourier tensor.

### 5.2 One shared input prime

If \(v=q\ell\) and \(w=qs\), with \(q,\ell,s\) distinct, (28) becomes

\[
kqs\equiv k'q\ell\pmod p.
\]

Since \(p\ne q\), this reduces exactly to

\[
\boxed{ks\equiv k'\ell\pmod p.}
\tag{31}
\]

The \(q\)-row appears twice, at two moving reciprocal arguments.

### 5.3 Disjoint semiprimes

If \(v=q\ell\) and \(w=rs\) are coprime, the condition is

\[
\boxed{krs\equiv k'q\ell\pmod p.}
\tag{32}
\]

Here \(p,q,\ell,r,s\) are five distinct primes. This is the genuine five-prime block.

---

## 6. Exact Hilbert--Schmidt trace

For squarefree \(d\),

\[
(U_d^*U_d)(m,n)
=\frac{c_d(m-n)}d.
\]

For a prime \(r\),

\[
\frac{c_r(h)}r=1_{r\mid h}-\frac1r,
\]

and for \(q\ell\), multiplicativity gives

\[
\frac{c_{q\ell}(h)}{q\ell}
=\left(1_{q\mid h}-\frac1q\right)
 \left(1_{\ell\mid h}-\frac1\ell\right).
\]

Therefore

\[
\boxed{
\|S_{p,q\ell}\|_{\rm HS}^2
=\frac1{pq\ell}\sum_{m,n<M}
 \prod_{r\in\{p,q,\ell\}}
 \left(1_{r\mid m-n}-\frac1r\right).
}
\tag{33}
\]

Define

\[
N_I(d)=\#\{(m,n)\in I^2:m\equiv n\pmod d\}.
\]

Expanding the product in (33), set

\[
\begin{aligned}
\Xi(p,q,\ell)
&=N_I(pq\ell)
 -\frac1pN_I(q\ell)
 -\frac1qN_I(p\ell)
 -\frac1\ell N_I(pq)\\
&\quad+\frac1{pq}N_I(\ell)
 +\frac1{p\ell}N_I(q)
 +\frac1{q\ell}N_I(p)
 -\frac{M^2}{pq\ell}.
\end{aligned}
\tag{34}
\]

Then

\[
\|S_{p,q\ell}\|_{\rm HS}^2
=\frac{\Xi(p,q,\ell)}{pq\ell}.
\]

Because \(p,q,\ell>X\) and \(M=X^2\), every pair product and the triple product exceed \(M\). Thus

\[
N_I(pq)=N_I(p\ell)=N_I(q\ell)=N_I(pq\ell)=M.
\]

For a single prime \(r\), write

\[
M=a_rr+b_r,
\qquad 0\le b_r<r.
\]

Every residue class occurs \(a_r\) or \(a_r+1\) times, so

\[
N_I(r)=b_r(a_r+1)^2+(r-b_r)a_r^2
=\frac{M^2}{r}+\varepsilon_r,
\]

where

\[
\varepsilon_r=b_r\left(1-\frac{b_r}{r}\right),
\qquad 0\le\varepsilon_r\le\frac r4.
\]

Substitution into (34) gives the exact closed form

\[
\boxed{
\begin{aligned}
\Xi(p,q,\ell)
&=M\left(1-\frac1p-\frac1q-\frac1\ell\right)
 +\frac{2M^2}{pq\ell}\\
&\quad+\frac{\varepsilon_p}{q\ell}
 +\frac{\varepsilon_q}{p\ell}
 +\frac{\varepsilon_\ell}{pq}.
\end{aligned}
}
\tag{35}
\]

In particular, uniformly for large \(X\),

\[
\frac{M}{2pq\ell}
\le\|S_{p,q\ell}\|_{\rm HS}^2
\le\frac{2M}{pq\ell}.
\tag{36}
\]

The primitive-frequency identity (3) follows directly from the bijection in Section 5:

\[
\begin{aligned}
\|S_{p,q\ell}\|_{\rm HS}^2
&=\frac1{p^2q^2\ell^2}
 \sum_{k\in(\mathbf Z/pq\ell\mathbf Z)^\times}
 \left|\mathcal D_M^\circ\!\left(\frac{k}{pq\ell}\right)\right|^2,\\
\sum_{k\in(\mathbf Z/pq\ell\mathbf Z)^\times}
 \left|\mathcal D_M^\circ\!\left(\frac{k}{pq\ell}\right)\right|^2
&=pq\ell\,\Xi(p,q,\ell).
\end{aligned}
\tag{37}
\]

Thus the full primitive-frequency trace has been computed without estimating a fixed gap.

Summing (36),

\[
\boxed{
\operatorname{Tr}(\mathsf S^*\mathsf S)
\asymp
M\sum_{p\in\mathcal O_X}\frac1p
 \sum_{\substack{q<\ell\\q,\ell\in\mathcal A_X\\p\notin\{q,\ell\}}}
 \frac1{q\ell}.
}
\tag{38}
\]

Hence

\[
\operatorname{Tr}(\mathsf S^*\mathsf S)
\ll M\sigma_{\mathcal O}\sigma_X^2,
\qquad
\sigma_{\mathcal O}=\sum_{p\in\mathcal O_X}\frac1p.
\tag{39}
\]

For active outputs,

\[
\boxed{
\operatorname{Tr}(\mathsf S^*\mathsf S)
\ll M\lambda^3.
}
\tag{40}
\]

For all outputs, Mertens on the dyadic block gives instead

\[
\operatorname{Tr}(\mathsf S^*\mathsf S)
\ll \frac{M\lambda^2}{\log X}.
\tag{41}
\]

---

## 7. Why the exact trace does not prove the target

The actual tensor input has

\[
\begin{aligned}
\|R\|_2^2
&=\sum_{q<\ell}\eta_q\eta_\ell\\
&\le\frac12\left(\sum_q\eta_q\right)^2\\
&\le\frac12(4X^2\lambda)^2
=8X^4\lambda^2.
\end{aligned}
\tag{42}
\]

Since

\[
\|\mathsf S\|_{\rm op}^2
\le\operatorname{Tr}(\mathsf S^*\mathsf S),
\]

(40) and (42) imply only

\[
\boxed{
\sum_{p\in\mathcal A_X}\|B_p\|_2^2
\le\operatorname{Tr}(\mathsf S^*\mathsf S)\|R\|_2^2
\ll X^6\lambda^5.
}
\tag{43}
\]

The desired estimate is

\[
X^{2+o(1)}\lambda^5=M X^{o(1)}\lambda^5.
\]

Thus exact Hilbert--Schmidt orthogonality has the correct fifth power only because the input norm contributes \(\lambda^2\), but it loses \(X^4\). Rescaling the input by \(q\ell\) merely transfers those \(X^4\) factors from the vector norm into the operator trace; it cannot improve (43).

---

## 8. What the additive large sieve gives after all parameters are inserted

The normalized semiprime coefficient mass is (10). All reduced fractions \(a/(q\ell)\), with \(q\ell\le4X^2\), are separated by at least

\[
\frac1{16X^4}.
\]

The dual additive large sieve on an interval of length \(M=X^2\) therefore has constant

\[
M-1+16X^4\ll X^4.
\]

The output operator satisfies \(\|(4p)^{-1}U_p^*U_p\|\le1/2\), exactly as in the prior semiprime calculation. Consequently, for any output set \(\mathcal O_X\),

\[
\boxed{
\sum_{p\in\mathcal O_X}\|B_p\|_2^2
\le |\mathcal O_X|
 (M-1+16X^4)\lambda^2.
}
\tag{44}
\]

For all output primes,

\[
\sum_{p\in\mathcal P_X}\|B_p\|_2^2
\ll\frac{X^5}{\log X}\lambda^2.
\tag{45}
\]

For active outputs, (11) gives

\[
\boxed{
\sum_{p\in\mathcal A_X}\|B_p\|_2^2
\ll X^5\lambda^3.
}
\tag{46}
\]

Relative to \(X^2\lambda^5\), this loses

\[
X^3\lambda^{-2}.
\]

The shared-prime subfamily has spacing \(\gg X^{-3}\), hence constant \(O(X^3)\), already one power above the interval scale. The disjoint-semiprime family has spacing \(\gg X^{-4}\), hence constant \(O(X^4)\), two powers above it. These are the same sharp row-geometric losses found in the preceding audit.

Using the full primitive fractions \(k/(pq\ell)\) does not turn this into a better generic large sieve:

* for one fixed \(P=pq\ell\asymp X^3\), primitive fractions have spacing \(P^{-1}\asymp X^{-3}\), so the exact complete-frequency constant is \(P\asymp X^3\), one power above \(M\);
* across distinct triple products \(P,P'\le8X^3\), distinct reduced fractions may be only \(1/(PP')\gg X^{-6}\) apart, so a direct all-triples large sieve has constant \(M+O(X^6)=O(X^6)\).

The tensor trace calculation in Section 6 is precisely what avoids this \(X^6\) catastrophe. It still stops at (43).

---

## 9. Diagonal, shared-prime, and disjoint density ledger

Write the column vector

\[
\beta_{p;q\ell}=S_{p,q\ell}R_{q,\ell}.
\]

Expanding the square gives

\[
\sum_p\|B_p\|_2^2
=\mathcal D_3+2\operatorname{Re}\mathcal S_4
 +2\operatorname{Re}\mathcal D_5,
\tag{47}
\]

where

\[
\mathcal D_3
=\sum_p\sum_{q<\ell}\|\beta_{p;q\ell}\|_2^2
\]

is the same-semiprime term, \(\mathcal S_4\) sums distinct semiprimes sharing one input prime, and \(\mathcal D_5\) sums disjoint semiprimes. The subscripts record the total number of distinct primes after including the output prime \(p\).

There is a rigorous inexpensive bound for \(\mathcal D_3\). For a fixed pair \(q,\ell\), the output fractions \(c/p\), \(p\in(X,2X]\), are \(1/(4X^2)\)-spaced. Hence

\[
\sum_{p\in\mathcal O_X}\|\beta_{p;q\ell}\|_2^2
\le5\sum_{m<M}g_q(m)^2g_\ell(m)^2.
\tag{48}
\]

To sum the right side, let

\[
f_r(m)=1_{\{m\bmod r\in\mathcal Z_r\}},
\qquad
A_r=\sum_{m<M}f_r(m).
\]

Since

\[
g_r(m)^2\le f_r(m)+z_r^2,
\]

and the exact one-column boundary estimate gives

\[
A_r\le(1+2/X)Mz_r,
\]

we obtain, using the proved second factorial moment,

\[
\begin{aligned}
H_2
&:=\sum_{q<\ell}\sum_{m<M}g_q(m)^2g_\ell(m)^2\\
&\le\frac12\sum_{m<M}(K_X(m))_2
 +(\sum_rA_r)(\sum_rz_r^2)
 +\frac M2(\sum_rz_r^2)^2\\
&\le\frac52M\lambda^2
 +(1+2/X)M\lambda^3
 +\frac12M\lambda^4.
\end{aligned}
\tag{49}
\]

Therefore

\[
\boxed{
\mathcal D_3
\ll M(\lambda^2+\lambda^4).
}
\tag{50}
\]

In the sparse range \(\lambda\le1\), this is only \(O(M\lambda^2)\), three density powers short of \(M\lambda^5\).

The raw \(m=n\) part of the same-semiprime block has one additional output reciprocal density:

\[
\begin{aligned}
\mathcal D_{3,0}
&=\sum_{p\in\mathcal A_X}\frac1p
  \sum_{q<\ell}\sum_{m<M}g_q(m)^2g_\ell(m)^2\\
&\le\sigma_XH_2\\
&\ll M(\lambda^3+\lambda^5).
\end{aligned}
\tag{51}
\]

Thus the natural three-modulus diagonal scale supplied by exact orthogonality is \(M\lambda^3\), not \(M\lambda^5\). The remaining same-semiprime shifts and the finite-mean subtraction are signed and could in principle cancel (51), but no available theorem forces such cancellation.

For comparison, a completely elementary column Cauchy bound uses

\[
\#\{q<\ell:q,\ell\in\mathcal A_X\}
\le2X^2\lambda^2
\]

and (50) to give

\[
\sum_p\|B_p\|_2^2
\ll X^4(\lambda^4+\lambda^6).
\tag{52}
\]

This is sometimes stronger than (46), but it still exceeds \(X^2\lambda^5\) by \(X^2/\lambda\) when \(\lambda\le1\).

The density bookkeeping therefore has a precise interpretation:

* the same-semiprime block visibly has only the three active prime slots \(p,q,\ell\);
* a shared-prime cross block has four distinct prime slots;
* only the disjoint block has five distinct prime slots.

Obtaining a fifth power of \(\lambda\) for the whole norm requires either two extra density gains in the same-semiprime sector and one in the shared sector, or a signed cancellation of those lower-overlap sectors against the disjoint sector. Neither follows from Parseval, the second moment, reflection, the gap-polynomial bounds, or the current one-prime Fourier estimates.

---

## 10. The exact missing correlation estimate

Define the pair field

\[
G_p(m)=\sum_{\substack{q<\ell\\q,\ell\ne p}}g_q(m)g_\ell(m).
\]

The desired operator estimate is exactly the following statement.

### Apéry five-prime centered-congruence estimate (AP5-CC)

Uniformly in dyadic \(X\),

\[
\boxed{
\sum_{p\in\mathcal A_X}
 \frac1p\sum_{u\bmod p}
 \left|
   \sum_{\substack{m<M\\m\equiv u\pmod p}}G_p(m)
   -\frac1p\sum_{m<M}G_p(m)
 \right|^2
\ll X^{2+o(1)}\lambda^5.
}
\tag{53}
\]

By (1)--(2), (53) is exactly equivalent to

\[
\sum_{p\in\mathcal A_X}\|B_p\|_2^2
\ll X^{2+o(1)}\lambda^5.
\]

Expanding \(G_p\), the same statement is

\[
\boxed{
\sum_{p\in\mathcal A_X}
 \sum_{\substack{q<\ell,\ r<s\\p\notin\{q,\ell,r,s\}}}
 \frac1p\sum_{m,n<M}
 \left(1_{m\equiv n\pmod p}-\frac1p\right)
 g_q(m)g_\ell(m)g_r(n)g_s(n)
\ll M X^{o(1)}\lambda^5.
}
\tag{54}
\]

Its full primitive reciprocal form is exactly (29), with the bound \(M X^{o(1)}\lambda^5\). Formula (29) names every modulus, inverse, frequency range, and congruence; no fixed \(h\) has been substituted.

If one insists on proving the three overlap sectors separately by absolute estimates, a sufficient package is

\[
\mathcal D_3\ll M X^{o(1)}\lambda^5,
\qquad
|\mathcal S_4|\ll M X^{o(1)}\lambda^5,
\qquad
|\mathcal D_5|\ll M X^{o(1)}\lambda^5.
\tag{55}
\]

This package is stronger than necessary. The weakest signed requirement, using (47), is

\[
\boxed{
2\operatorname{Re}(\mathcal S_4+\mathcal D_5)
\le-\mathcal D_3+M X^{o(1)}\lambda^5.
}
\tag{56}
\]

In a polynomially sparse block, (56) makes explicit what any proof of the unweighted \(B\)-norm target must accomplish: it must cancel the natural lower-overlap energy down to the five-density scale.

---

## 11. Why additive-character detection is not yet a Kloosterman theorem

The congruence in (28) can of course be detected by

\[
1_{kw\equiv k'v\ (p)}
=\frac1p\sum_{h\bmod p}e_p(h(kw-k'v)).
\tag{57}
\]

This identity alone gives no saving. The \(h=0\) term is nonoscillatory. For \(h\ne0\), the coefficients in (29) are

\[
R_q(-k\overline{p\ell})R_\ell(-k\overline{pq})
\overline{R_r(-k'\overline{ps})R_s(-k'\overline{pr})},
\]

so they depend simultaneously on all of \(p,q,\ell,r,s,k,k'\).

The parameter ranges are explicit:

\[
p,q,\ell,r,s\asymp X,
\qquad pq\ell,prs\asymp X^3,
\qquad M=X^2.
\]

Moreover

\[
\left|\mathcal D_M^\circ(k/P)\right|
\ll\min\left(M,\frac{P}{|k|_P}\right),
\]

so the principal Dirichlet lobe has

\[
|k|_P\lesssim P/M\asymp X,
\qquad
|k'|_{P'}\lesssim P'/M\asymp X.
\]

Thus a putative bilinear/Kuznetsov step would have modulus \(p\asymp X\), two short variables of length \(X\), and four additional prime variables of length \(X\), with modulus-dependent Fourier coefficients evaluated at moving inverses. The repository currently supplies only

\[
\sum_{a\ne0}|R_r(a)|^2=r^2\rho_r
\]

and one-prime Fourier nonconcentration in ordinary frequency intervals. It does not supply a horizontal trace-function family controlling

\[
R_q(k\overline{p\ell})
\]

as \(p,\ell,q\) vary. Cauchy followed by Parseval removes the reciprocal arguments but returns precisely the \(\lambda^2\) mass in (10), leading back to (44)--(46).

Accordingly, no DI/DFI/Kuznetsov conclusion is available from the present inputs. Invoking one would require first proving a theorem tailored to the coefficient family in (29), with the above lengths and moduli.

A usable new input could be stated as a reciprocal tensor large sieve: after dyadically decomposing the Dirichlet weights in \(k,k'\), the signed sum in (29) should be bounded by \(M X^{o(1)}\lambda^5\), uniformly for \(K,K'\le X^{1+o(1)}\). This is precisely a new cross-prime Apéry correlation theorem.

---

## Final conclusion

The full primitive-frequency modulus does not reveal a hidden standard large-sieve proof of the HM3 semiprime output norm. It gives three rigorous and useful facts:

1. the ordered symmetric operator is exactly the nonzero-\(p\) Fourier transform of the actual pair field \(G_p=\sum_{q<\ell}g_qg_\ell\);
2. each \(p,q,\ell\) block uses every primitive frequency modulo \(pq\ell\), and its Hilbert--Schmidt square is exactly \(\Xi(p,q,\ell)/(pq\ell)\asymp M/(pq\ell)\);
3. \(\mathsf S^*\mathsf S\) has the explicit reciprocal congruences (30)--(32), culminating in the five-prime condition \(krs\equiv k'q\ell\pmod p\).

The exact active-triple trace is \(O(M\lambda^3)\), but applying it to the actual Fourier tensors gives only \(O(X^6\lambda^5)\). Generic semiprime duality gives \(O(X^5\lambda^3)\). Neither reaches \(M X^{o(1)}\lambda^5\).

The unresolved input is therefore not another normalization or another application of Parseval. It is AP5-CC, equivalently the reciprocal tensor correlation (29): aggregate cancellation of the centered Apéry pair field across the moving output congruences, including the same-semiprime and shared-prime overlap sectors. Without such a new correlation theorem, the unweighted \(B\)-norm route remains open and, in polynomially sparse blocks, substantially stronger than what the currently proved row and moment estimates imply.
