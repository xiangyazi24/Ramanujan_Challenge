ANSWER Q6372 ce335333

# Self-twists of the Apéry rank-two system

## Verdict

The required high-order exclusion is true, in a form stronger than the proposed “order divides \(4\)” statement.

There are three objects that must not be conflated.

1. The **raw square-root system** \(\mathcal V_	au\), whose distinguished holomorphic solution is \(	au=\sqrt F\). Its determinant is a nontrivial quadratic Kummer system \(\delta\).
2. The other raw branch \(\mathcal V_\sigma\), whose distinguished solution is \(\sigma=\sqrt{F/q}\). It is the quadratic twist
   \[
   \mathcal V_\sigma\simeq\mathcal V_	au\otimes\delta\simeq\mathcal V_	au^ee.
   \]
3. The genuinely determinant-one system
   \[
   \widehat{\mathcal V}:=\mathcal V_	au\otimes\eta^{-1},
   \qquad \eta^2=\delta,
   \]
   where the canonical choice is the order-four Kummer system \(\eta=q(t)^{1/4}\).

For a **fixed normalization**, every geometric finite-order self-twist is trivial. A coefficient/Frobenius conjugate of \(\widehat{\mathcal V}\) differs from \(\widehat{\mathcal V}\) by at most the quadratic character \(\delta\). Thus the actual conjugate-twist list is
\[
oxed{\{1,\delta\}},
\]
not a high-order family. The four powers
\[
\{1,\eta,\eta^2=\delta,\eta^{-1}\}
\]
appear only when one changes between the raw and determinant-one normalizations; \(\eta\) itself is not a self-twist.

Moreover:

- \(\operatorname{Sym}^2\mathcal V_	au\) has **no nontrivial geometric self-twist or conjugate-twist**.
- The literal external product \(\mathcal V_1oxtimes\mathcal V_2\) on \(U	imes U\) is irreducible and has only factorwise twists from \(\{1,\delta\}\).
- After diagonal pullback, the rank-four tensor product decomposes only in the projectively equivalent case, and then only as the Clebsch–Gordan sum \(3+1\), with the same quadratic list. There is no \(2+2\) or other decomposition.

This proves that no high-order Mellin Kummer character can be a self-twist or cause the exceptional tensor decomposition used in the counting lemma.

A terminology correction is essential: \(oxtimes\) on \(U	imes U\) does **not** decompose; the familiar \(3+1\) decomposition occurs only after pullback to the diagonal, where \(oxtimes\) becomes an ordinary tensor product.

---

# 1. The rank-two differential operator

Let
\[
F(t)=\sum_{n\ge0}b_nt^n,\qquad q(t)=1-34t+t^2,\qquad 	au(t)=\sqrt{F(t)}=\sum_{n\ge0}	au_nt^n.
\]
The given half-index recurrence is
\[
4(j+2)^2	au_{j+2}
=2(68j^2+170j+107)	au_{j+1}-(2j+1)^2	au_j.
	ag{1.1}
\]
Put \(n=j+1\). Then
\[
(n+1)^2	au_{n+1}
-\left(34n^2+17n+rac52ight)	au_n
+\left(n-rac12ight)^2	au_{n-1}=0.
	ag{1.2}
\]
With \(	heta=t\,d/dt\), the coefficient of \(t^{n+1}\) in
\[
\mathcal D_	au
:=	heta^2-t\left(34	heta^2+17	heta+rac52ight)
+t^2\left(	heta+rac12ight)^2
	ag{1.3}
\]
is exactly the left side of (1.2). Hence \(	au\) is annihilated by \(\mathcal D_	au\).

Multiplying (1.3) by \(4\) gives the useful form
\[
oxed{
4q(t)	heta^2+4t(t-17)	heta+t(t-10).
}
	ag{1.4}
\]
In ordinary derivatives this is
\[
4t^2q(t)y''+4t(1-51t+2t^2)y'+t(t-10)y=0.
	ag{1.5}
\]
Thus the only singular points are
\[
0,\qquad \alpha=17+12\sqrt2,\qquad eta=17-12\sqrt2,\qquad \infty,
\]
where \(q(t)=(t-\alpha)(t-eta)\) and \(\alphaeta=1\).

For comparison, the \(\sigma=\sqrt{F/q}\) coefficients satisfy
\[
4(j+2)^2\sigma_{j+2}
=2(68j^2+238j+209)\sigma_{j+1}-(2j+3)^2\sigma_j.
	ag{1.6}
\]
Equivalently,
\[
\mathcal D_\sigma
=	heta^2-t\left(34	heta^2+51	heta+rac{39}{2}ight)
+t^2\left(	heta+rac32ight)^2.
	ag{1.7}
\]
The two operators are related by the gauge transformation
\[
	au=\sqrt q\,\sigma.
	ag{1.8}
\]

---

# 2. Exact local exponents

We compute the exponents of \(\mathcal D_	au\) directly.

## 2.1 At \(t=0\)

In (1.3), all terms containing \(t\) vanish in the indicial equation. Substituting \(y=t^\lambda\) gives
\[
\lambda^2=0.
\]
Therefore the exponents at \(0\) are
\[
oxed{0,0}.
	ag{2.1}
\]
The local monodromy eigenvalues are \(1,1\).

## 2.2 At a root \(a\in\{\alpha,eta\}\) of \(q\)

Write (1.5) as
\[
A_2(t)y''+A_1(t)y'+A_0(t)y=0,
\]
where
\[
A_2(t)=4t^2q(t),\qquad A_1(t)=4t(1-51t+2t^2).
\]
Since \(a\) is a simple root of \(q\), the indicial equation at \(a\) is
\[
\lambdaigl(A_2'(a)(\lambda-1)+A_1(a)igr)=0.
	ag{2.2}
\]
Using \(q(a)=0\), hence \(1=34a-a^2\), we obtain
\[
1-51a+2a^2=a(a-17),
\]
so
\[
A_1(a)=4a^2(a-17).
\]
Also
\[
A_2'(a)=4a^2q'(a)=8a^2(a-17).
\]
Thus
\[
rac{A_1(a)}{A_2'(a)}=rac12,
\]
and (2.2) gives
\[
oxed{0,rac12}.
	ag{2.3}
\]
The local monodromy is semisimple, with eigenvalues \(1,-1\), because the exponent difference is nonintegral.

## 2.3 At \(t=\infty\)

Put \(z=1/t\) and \(\delta_z=z\,d/dz=-	heta\). Equation (1.3), multiplied by \(z^2\), becomes
\[
z^2\delta_z^2-z\left(34\delta_z^2-17\delta_z+rac52ight)
+\left(-\delta_z+rac12ight)^2.
\]
At \(z=0\), the indicial equation is
\[
\left(-\lambda+rac12ight)^2=0.
\]
Hence the exponents at infinity are
\[
oxed{rac12,rac12}.
	ag{2.4}
\]
The two local monodromy eigenvalues are both \(-1\).

Consequently the Riemann scheme of the raw rank-two system is
\[
\mathcal P\left\{
egin{array}{cccc}
0&\alpha&eta&\infty\[2pt]
0&0&0&rac12\
0&rac12&rac12&rac12
\end{array}ight\}.
	ag{2.5}
\]
The exponent sum is \(0+rac12+rac12+1=2\), as required by the Fuchs relation for a second-order equation with four singular points.

For \(\mathcal D_\sigma\), the same calculation gives
\[
\mathcal P\left\{
egin{array}{cccc}
0&\alpha&eta&\infty\[2pt]
0&0&0&rac32\
0&-rac12&-rac12&rac32
\end{array}ight\}.
	ag{2.6}
\]
Modulo integral exponent shifts, the local monodromy eigenvalue multisets are the same as for \(\mathcal V_	au\).

Finally, taking symmetric squares of (2.5) gives the order-three Apéry exponents
\[
egin{array}{c|c}
0&0,0,0\
\alpha&0,rac12,1\
eta&0,rac12,1\
\infty&1,1,1.
\end{array}
	ag{2.7}
\]
These agree with the original operator
\[
	heta^3-t(34	heta^3+51	heta^2+27	heta+5)+t^2(	heta+1)^3.
\]

---

# 3. Determinant, the quadratic branch twist, and the \(SL_2\) normalization

Let
\[
U=\mathbf P^1-\{0,\alpha,eta,\infty\}.
\]
Work either with complex local systems or with tame \(\overline{\mathbf Q}_\ell\)-local systems, with residue characteristic and \(\ell\) different from \(2\).

Let \(\gamma_0,\gamma_\alpha,\gamma_eta,\gamma_\infty\) be standard positively oriented local loops, with
\[
\gamma_0\gamma_\alpha\gamma_eta\gamma_\infty=1.
\]
From (2.5), the determinant character of \(\mathcal V_	au\) has local values
\[
1,-1,-1,1.
	ag{3.1}
\]
Define
\[
\delta:=\mathcal L_{\sqrt{q(t)}}.
\]
This is the quadratic Kummer local system whose local monodromy tuple is exactly (3.1). Therefore
\[
oxed{\det\mathcal V_	au=\delta,\qquad \delta^2=1.}
	ag{3.2}
\]
The gauge identity \(	au=\sqrt q\,\sigma\) gives
\[
oxed{\mathcal V_\sigma\simeq\mathcal V_	au\otimes\delta.}
	ag{3.3}
\]
For any rank-two local system \(V\), there is a canonical identity
\[
V^ee\simeq V\otimes(\det V)^{-1}.
\]
Using (3.2),
\[
oxed{\mathcal V_\sigma\simeq\mathcal V_	au^ee.}
	ag{3.4}
\]

The raw system is therefore not literally \(SL_2\)-valued: its connected geometric monodromy is \(SL_2\), but its full geometric monodromy has the quadratic determinant component \(\delta\).

Set
\[
\eta:=\mathcal L_{q(t)^{1/4}}.
\]
Its local monodromy tuple is
\[
1,i,i,-1,
	ag{3.5}
\]
so
\[
\eta^2=\delta.
	ag{3.6}
\]
Define the determinant-one normalization
\[
oxed{\widehat{\mathcal V}:=\mathcal V_	au\otimes\eta^{-1}.}
	ag{3.7}
\]
Then
\[
\det\widehat{\mathcal V}=1.
\]
Its local exponents are
\[
egin{array}{c|c}
0&0,0\
\alpha&-rac14,rac14\
eta&-rac14,rac14\
\infty&0,0.
\end{array}
	ag{3.8}
\]
Thus its local eigenvalue multisets are
\[
\{1,1\},\quad\{i,-i\},\quad\{i,-i\},\quad\{1,1\}.
	ag{3.9}
\]
By the hypothesis in the question, the geometric monodromy of \(\widehat{\mathcal V}\) is Zariski-dense in \(SL_2\).

The cyclic order-four list
\[
\{1,\eta,\delta,\eta^{-1}\}
	ag{3.10}
\]
accounts for every change among the raw \(	au\)-branch, the raw \(\sigma\)-branch, and the determinant-one normalizations. It should not be called a self-twist group: the genuine self-twist group will be shown to be trivial.

---

# 4. A general finite-self-twist rigidity lemma

We use the following elementary lemma repeatedly.

## Lemma 4.1

Let \(\Gamma\) be a group, \(E\) an algebraically closed field of characteristic \(0\), and
\[
r:\Gamma\longrightarrow GL(V)
\]
an absolutely irreducible representation. Assume that the identity component of the Zariski closure of \(r(\Gamma)\) acts irreducibly on \(V\). Let \(\chi:\Gamma	o E^	imes\) be a finite-order character. If
\[
r\otimes\chi\simeq r,
\]
then \(\chi=1\).

### Proof

Let \(A
e0\) be an intertwiner:
\[
A\,r(g)\chi(g)=r(g)\,A\qquad(g\in\Gamma).
	ag{4.1}
\]
Put \(\Gamma_0=\ker\chi\). Since \(\Gamma_0\) has finite index in \(\Gamma\), the Zariski closures of \(r(\Gamma_0)\) and \(r(\Gamma)\) have the same identity component. For \(h\in\Gamma_0\), (4.1) says
\[
Ar(h)=r(h)A.
\]
By irreducibility of the connected monodromy action and Schur's lemma, \(A=cI\) for some \(c
e0\). Substituting into (4.1) gives
\[
\chi(g)r(g)=r(g)\qquad(g\in\Gamma),
\]
so \(\chi(g)=1\) for every \(g\). Hence \(\chi=1\). \(\square\)

For \(\widehat{\mathcal V}\), the connected group is \(SL_2\) acting through its standard representation, so the lemma applies. It also applies to \(\operatorname{Sym}^2\widehat{\mathcal V}\), because \(\operatorname{Sym}^2\) of the standard \(SL_2\)-representation is irreducible.

---

# 5. Complete rank-two twist classification

Let \(\mathcal L\) be a finite-order rank-one local system on \(U\). Write its local monodromy scalars as
\[
\lambda_0,\lambda_\alpha,\lambda_eta,\lambda_\infty,
\qquad
\lambda_0\lambda_\alpha\lambda_eta\lambda_\infty=1.
	ag{5.1}
\]
Suppose
\[
\mathcal V_	au\otimes\mathcal L\simeq\mathcal V_	au'
	ag{5.2}
\]
where \(\mathcal V_	au'\) is either \(\mathcal V_	au\), \(\mathcal V_\sigma\), or a coefficient/Frobenius conjugate having the same local exponent multisets.

At \(0\), the eigenvalue multiset of the left side is
\[
\{\lambda_0,\lambda_0\}.
\]
The target multiset is \(\{1,1\}\), so
\[
\lambda_0=1.
	ag{5.3}
\]
At infinity, the source multiset is
\[
\{-\lambda_\infty,-\lambda_\infty\},
\]
and the target multiset is \(\{-1,-1\}\), so
\[
\lambda_\infty=1.
	ag{5.4}
\]
At either root of \(q\), the source multiset is
\[
\{\lambda,-\lambda\},
\]
and the target multiset is \(\{1,-1\}\). Hence
\[
\lambda\in\{1,-1\}.
	ag{5.5}
\]
Using (5.1),
\[
\lambda_\alpha\lambda_eta=1,
\]
so
\[
\lambda_\alpha=\lambda_eta.
	ag{5.6}
\]
Therefore the only two local systems allowed by local monodromy are
\[
oxed{\mathcal L=1\quad	ext{or}\quad\mathcal L=\delta.}
	ag{5.7}
\]
This also follows from determinants:
\[
\det(\mathcal V_	au\otimes\mathcal L)
=\det(\mathcal V_	au)\otimes\mathcal L^2,
\]
and every conjugate under consideration has determinant \(\delta\), so \(\mathcal L^2=1\).

We now identify which candidate actually occurs.

## 5.1 Exact self-twists

If
\[
\mathcal V_	au\otimes\mathcal L\simeq\mathcal V_	au,
\]
Lemma 4.1 gives
\[
oxed{\mathcal L=1.}
	ag{5.8}
\]
In particular, although \(\delta\) passes the local spectral test, it is not a self-twist. Indeed,
\[
\mathcal V_	au\otimes\delta\simeq\mathcal V_\sigma\simeq\mathcal V_	au^ee
\]
is a different rank-two lift of the same symmetric square.

The same argument gives
\[
\operatorname{SelfTw}(\mathcal V_\sigma)=\{1\},\qquad
\operatorname{SelfTw}(\widehat{\mathcal V})=\{1\}.
	ag{5.9}
\]

## 5.2 The two branch systems

By (3.3),
\[
\mathcal V_	au\otimes\mathcal L\simeq\mathcal V_\sigma
\]
holds if and only if
\[
oxed{\mathcal L=\delta.}
	ag{5.10}
\]
The isomorphism space is one-dimensional, by absolute irreducibility.

Thus, up to a nonzero scalar multiple of the intertwiner, the complete branch-pair list is
\[
egin{array}{c|c|c}
	ext{source}&	ext{target}&\mathcal L\ \hline
\mathcal V_	au&\mathcal V_	au&1\
\mathcal V_	au&\mathcal V_\sigma&\delta\
\mathcal V_\sigma&\mathcal V_	au&\delta\
\mathcal V_\sigma&\mathcal V_\sigma&1.
\end{array}
	ag{5.11}
\]

## 5.3 Frobenius conjugates of the determinant-one normalization

The raw operator \(\mathcal D_	au\) has rational coefficients, so its coefficient conjugates are isomorphic to itself. The only nontrivial coefficient dependence in \(\widehat{\mathcal V}\) comes from the fourth root \(\eta\).

Let \(arphi\) be a coefficient automorphism. Since \(arphi(i)=i\) or \(-i\),
\[
arphi(\eta)=\eta\quad	ext{or}\quad\eta^{-1}.
\]
Define \(arepsilon(arphi)\in\{0,1\}\) by
\[
arphi(\eta)=\eta^{(-1)^{arepsilon(arphi)}}.
\]
Then
\[
oxed{\widehat{\mathcal V}^{\,arphi}\simeq
\widehat{\mathcal V}\otimes\delta^{arepsilon(arphi)}.}
	ag{5.12}
\]
Indeed, if \(arphi(\eta)=\eta^{-1}\), then
\[
\widehat{\mathcal V}^{\,arphi}\simeq\mathcal V_	au\otimes\eta
=\mathcal V_	au\otimes\eta^{-1}\otimes\eta^2
=\widehat{\mathcal V}\otimes\delta.
\]
For arithmetic Frobenius at an odd prime \(p\),
\[
\eta^{(p)}=\eta^p,
\]
so
\[
oxed{
\widehat{\mathcal V}^{(p)}\simeq
egin{cases}
\widehat{\mathcal V},&p\equiv1\pmod4,\
\widehat{\mathcal V}\otimes\delta,&p\equiv3\pmod4.
\end{cases}}
	ag{5.13}
\]
If
\[
\widehat{\mathcal V}\otimes\mathcal L\simeq\widehat{\mathcal V}^{\,arphi},
\]
then (5.12) and Lemma 4.1 imply
\[
oxed{\mathcal L=\delta^{arepsilon(arphi)}.}
	ag{5.14}
\]
Thus every conjugate-twist has order at most \(2\).

---

# 6. Symmetric-square self-twists

Put
\[
\mathcal A:=\operatorname{Sym}^2\mathcal V_	au.
\]
This is the rank-three Apéry local system.

At \(0\), its local eigenvalues are
\[
\{1,1,1\}.
\]
At either root of \(q\), the eigenvalues of \(\mathcal V_	au\) are \(1,-1\), so those of \(\mathcal A\) are
\[
\{1,-1,1\}=\{1,1,-1\}.
	ag{6.1}
\]
At infinity, the eigenvalues are
\[
\{1,1,1\}.
\]
Suppose
\[
\mathcal A\otimes\mathcal L\simeq\mathcal A^arphi.
	ag{6.2}
\]
All coefficient/Frobenius conjugates have the same local eigenvalue multisets, because the only eigenvalues involved are \(\pm1\).

At \(0\) and infinity, (6.2) forces
\[
\lambda_0=\lambda_\infty=1.
\]
At \(\alpha\), the twisted multiset is
\[
\{\lambda_\alpha,\lambda_\alpha,-\lambda_\alpha\}.
\]
For this to equal \(\{1,1,-1\}\), the repeated eigenvalue must be \(1\), hence
\[
\lambda_\alpha=1.
\]
Similarly,
\[
\lambda_eta=1.
\]
A tame geometric rank-one local system on a punctured projective line is determined by its local monodromy tuple, so
\[
oxed{\mathcal L=1.}
	ag{6.3}
\]

There is also a representation-theoretic proof: the connected monodromy acts on \(\mathcal A\) through the irreducible highest-weight-two representation of \(SL_2\), so Lemma 4.1 applies.

Notice the distinction
\[
\operatorname{Sym}^2(\mathcal V_	au\otimes\delta)\simeq
\operatorname{Sym}^2\mathcal V_	au,
	ag{6.4}
\]
because \(\delta^2=1\). Equation (6.4) says that \(\mathcal V_	au\) and \(\mathcal V_\sigma\) are two rank-two lifts of the same rank-three system. It does **not** give a nontrivial self-twist of the rank-three system itself.

For the determinant-one normalization,
\[
\operatorname{Sym}^2(\widehat{\mathcal V}^{\,arphi})\simeq
\operatorname{Sym}^2(\widehat{\mathcal V}\otimes\delta^{arepsilon(arphi)})\simeq
\operatorname{Sym}^2\widehat{\mathcal V}.
	ag{6.5}
\]
Hence every Frobenius conjugate of the symmetric square is already isomorphic to it, and its only twist is \(1\).

---

# 7. The rank-four external product

Let
\[
\mathcal R_{arphi,\psi}
:=\widehat{\mathcal V}^{\,arphi}oxtimes\widehat{\mathcal V}^{\,\psi}
\]
be the external product on \(U	imes U\).

In the Betti category,
\[
\pi_1(U	imes U)=\pi_1(U)	imes\pi_1(U).
\]
The same statement holds for the tame prime-to-characteristic category used by Kummer sheaves. Therefore every finite-order rank-one local system on the product splits as
\[
\mathcal L=\mathcal L_1oxtimes\mathcal L_2.
	ag{7.1}
\]

Because each factor is absolutely irreducible, the external tensor product is absolutely irreducible. Indeed,
\[
\operatorname{End}_{\pi_1(U)	imes\pi_1(U)}
(V_1\otimes V_2)
=\operatorname{End}_{\pi_1(U)}(V_1)\otimes
\operatorname{End}_{\pi_1(U)}(V_2)=E.
	ag{7.2}
\]
Thus the literal object denoted by \(oxtimes\) has no nontrivial direct-sum decomposition.

Now suppose
\[
\mathcal R_{arphi,\psi}\otimes(\mathcal L_1oxtimes\mathcal L_2)
\simeq
\mathcal R_{arphi',\psi'}.
	ag{7.3}
\]
The Hom space factors:
\[
egin{aligned}
&\operatorname{Hom}_{\pi_1(U)	imes\pi_1(U)}
igl((V_arphi\otimes L_1)oxtimes(V_\psi\otimes L_2),
V_{arphi'}oxtimes V_{\psi'}igr)\
&\qquad\simeq
\operatorname{Hom}_{\pi_1(U)}(V_arphi\otimes L_1,V_{arphi'})
\otimes
\operatorname{Hom}_{\pi_1(U)}(V_\psi\otimes L_2,V_{\psi'}).
\end{aligned}
	ag{7.4}
\]
Consequently (7.3) holds if and only if both factorwise twist equations hold. By (5.14),
\[
oxed{
\mathcal L_1=\delta^{arepsilon(arphi')-arepsilon(arphi)},\qquad
\mathcal L_2=\delta^{arepsilon(\psi')-arepsilon(\psi)},
}
	ag{7.5}
\]
where exponents are taken modulo \(2\).

Thus the complete external-product conjugate-twist list is
\[
oxed{
1oxtimes1,\quad
\deltaoxtimes1,\quad
1oxtimes\delta,\quad
\deltaoxtimes\delta.
}
	ag{7.6}
\]
Every member has order at most \(2\). For an exact self-twist, the source and target conjugates are the same, so (7.5) gives
\[
oxed{\operatorname{SelfTw}(\mathcal R_{arphi,\psi})=\{1oxtimes1\}.}
	ag{7.7}
\]

---

# 8. Diagonal tensor products and the only possible decomposition

The rank-four object decomposes only after pullback by the diagonal
\[
\Delta:U\longrightarrow U	imes U.
\]
Then
\[
\Delta^*(\mathcal V_1oxtimes\mathcal V_2)=\mathcal V_1\otimes\mathcal V_2.
\]

For the explicit Frobenius conjugates of \(\widehat{\mathcal V}\), write
\[
\widehat{\mathcal V}^{\,arphi}\simeq
\widehat{\mathcal V}\otimes\delta^{e_arphi},\qquad e_arphi\in\{0,1\}.
\]
Then
\[
egin{aligned}
\widehat{\mathcal V}^{\,arphi}\otimes\widehat{\mathcal V}^{\,\psi}
&\simeq
(\widehat{\mathcal V}\otimes\widehat{\mathcal V})\otimes
\delta^{e_arphi+e_\psi}\
&\simeq
\left(\operatorname{Sym}^2\widehat{\mathcal V}\oplus
igwedge^2\widehat{\mathcal V}ight)\otimes
\delta^{e_arphi+e_\psi}.
\end{aligned}
\]
Since \(\det\widehat{\mathcal V}=1\),
\[
oxed{
\widehat{\mathcal V}^{\,arphi}\otimes\widehat{\mathcal V}^{\,\psi}
\simeq
\left(\operatorname{Sym}^2\widehat{\mathcal V}\oplus1ight)\otimes
\delta^{e_arphi+e_\psi}.
}
	ag{8.1}
\]
This is the unique \(3+1\) decomposition.

For arithmetic Frobenius at a prime \(p\), (5.13) gives
\[
oxed{
\widehat{\mathcal V}\otimes\widehat{\mathcal V}^{(p)}\simeq
egin{cases}
\operatorname{Sym}^2\widehat{\mathcal V}\oplus1,&p\equiv1\pmod4,\[2pt]
(\operatorname{Sym}^2\widehat{\mathcal V}\otimes\delta)\oplus\delta,&p\equiv3\pmod4.
\end{cases}}
	ag{8.2}
\]

For the raw branches,
\[
oxed{
\mathcal V_	au\otimes\mathcal V_	au
\simeq\operatorname{Sym}^2\mathcal V_	au\oplus\delta,
}
	ag{8.3}
\]
and
\[
oxed{
\mathcal V_	au\otimes\mathcal V_\sigma
\simeq(\operatorname{Sym}^2\mathcal V_	au\otimes\delta)\oplus1.
}
	ag{8.4}
\]

## 8.1 Why there are no other decompositions

The following standard Goursat argument makes the classification intrinsic.

Let \(V_1,V_2\) be two rank-two local systems whose connected projective monodromy groups are \(PGL_2\). Let \(H\) be the connected Zariski closure of the joint projective representation
\[
\pi_1(U)\longrightarrow PGL_2	imes PGL_2.
\]
Both projections \(H	o PGL_2\) are surjective. Since \(PGL_2\) is simple, Goursat's lemma gives exactly two possibilities.

1. \(H=PGL_2	imes PGL_2\). Then the representation \(\mathrm{Std}\otimes\mathrm{Std}\) is irreducible, so \(V_1\otimes V_2\) is irreducible of rank \(4\).
2. \(H\) is the graph of an automorphism of \(PGL_2\). Every automorphism of \(PGL_2\) over an algebraically closed characteristic-zero field is inner. Therefore \(V_1\) and \(V_2\) are projectively isomorphic, so
   \[
   V_2\simeq V_1\otimes L
   \]
   for a rank-one local system \(L\). The rank-two classification above forces \(L\in\{1,\delta\}\) in the Apéry family. Then
   \[
   V_1\otimes V_2\simeq
   \left(\operatorname{Sym}^2V_1\oplus\det V_1ight)\otimes L,
   \]
   which is exactly the \(3+1\) decomposition already listed.

Thus there is no \(2+2\) decomposition and no decomposition caused by a high-order Kummer character.

## 8.2 Self-twists of the diagonal rank-four object

Let
\[
\mathcal T_e:=\left(\operatorname{Sym}^2\widehat{\mathcal V}\oplus1ight)\otimes\delta^e,\qquad e\in\{0,1\}.
\]
The rank-one constituent is unique, because the other constituent is irreducible of rank \(3\). If
\[
\mathcal T_e\otimes\mathcal L\simeq\mathcal T_e,
\]
the isomorphism must preserve the rank-one constituent. Hence
\[
\delta^e\otimes\mathcal L\simeq\delta^e,
\]
so
\[
oxed{\mathcal L=1.}
	ag{8.5}
\]
Between the two Frobenius-conjugate forms,
\[
\mathcal T_e\otimes\mathcal L\simeq\mathcal T_{e'}
\]
if and only if
\[
oxed{\mathcal L=\delta^{e'-e}.}
	ag{8.6}
\]
Again the order is at most \(2\).

---

# 9. Consequence for high-order Mellin characters

A finite-order Kummer character capable of producing a rank-two conjugate twist must belong to
\[
\{1,\delta\}.
\]
If one compares different determinant normalizations as well, the entire list is contained in
\[
\langle\etaangle=\{1,\eta,\delta,\eta^{-1}\},
\]
which has order \(4\).

The symmetric square admits no nontrivial geometric twist at all. The external two-factor object admits only the four factorwise quadratic possibilities (7.6), and its exact self-twist group is trivial. The diagonal rank-four tensor has only the \(3+1\) decomposition (8.1), caused by projective equivalence of the two factors, and the relevant rank-one twist is again \(1\) or \(\delta\).

Therefore:
\[
oxed{	ext{No Kummer character of order \(>4\) can be a self-twist or decomposition character.}}
	ag{9.1}
\]
In fact, after fixing one normalization, “\(>4\)” can be strengthened to “\(>2\).”

This is exactly the geometric prerequisite needed by the character-order counting lemma: all self-twist and tensor-decomposition exceptions lie in a fixed bounded-order set, while the surviving Mellin characters have growing order.

---

# 10. Category and normalization caveat

The proof above classifies **geometric finite-order Kummer local systems**, which is the category relevant to Mellin-character order. Two common broader phrases require care.

1. A rank-one **Weil** local system can be multiplied by a constant character of the arithmetic Galois group. Local exponents do not see that constant factor. Once a Frobenius normalization is fixed, Lemma 4.1 also rules out a nontrivial constant exact self-twist; however, changing the Weil normalization can artificially insert a constant character into an isomorphism with a conjugate. Such a change has no relation to high-order geometric Mellin characters.
2. The notation \(oxtimes\) means an external product on \(U	imes U\), which is irreducible. Any statement that it “decomposes as symmetric square plus determinant” is shorthand for its **diagonal pullback**, not for the external product itself.

With these conventions fixed, the theorem above is complete and has no residual geometric hypothesis beyond the stated Zariski-density of the connected monodromy.

# Least-confident step

The least-confident point is terminological rather than mathematical: the project notes use “Frobenius-twisted two-prime object” without always distinguishing the external product on \(U	imes U\) from its diagonal pullback, and without always specifying whether Frobenius conjugation includes a change of Weil normalization. Under the explicit conventions in Sections 5.3, 7, and 8, every proof step is formal and the twist list is exact.