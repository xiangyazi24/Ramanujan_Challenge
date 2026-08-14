ANSWER Q8215 abddd8b9

# P3.2 five-prime estimate: affine Hermite/Newton-shell audit

## Verdict

The special Apéry recurrence does **not** currently supply the missing affine closest-point theorem, and the proposed nonzero-mismatch route has a sharper obstruction than merely “Minkowski is homogeneous.”  For a fixed selected row \(m\), reflection reduces every even Hermite interpolant to one polynomial in \(Y=T(T+1)\).  After removing at most two central primes, **all derivative conditions collapse to one scalar Chinese-remainder class modulo the product of the hit primes**.  The corresponding affine coefficient lattice has exact index equal to that product.  No homogeneous Minkowski argument, Newton/Pascal basis change, or standard lattice transference theorem can make an arbitrary affine class short; a simple counting argument already forbids it at the requested exponents.

More decisively, if the proposed integer mismatch \(\Delta_m\neq0\) is divisible by every hit prime, then
\[
 \log |\Delta_m|>K(m)\log X.
\]
Hence for a selected row
\[
 \log |\Delta_m|\ge L\log X
 =X^{14/15+\eta}\log X.
\]
So an \(H\)-scale construction with \(\log |\Delta_m|=O(H\log X)\) is impossible, since
\[
 \frac{L}{H}=X^{7/15+\eta}(\log X)^{1/3}\to\infty.
\]
Moreover
\[
 (K)_4\log(2+|\Delta_m|)
 \ge (K)_4K\log X
 =\bigl((K)_5+4(K)_4\bigr)\log X.
\]
Thus the requested weighted \(\Delta_m\)-estimate is already at least as strong as the five-prime estimate it is supposed to prove.  It is a valid **conditional implication**, but not a height shortcut.

The special Apéry/Gessel structure does provide exact first-jet identities: the Gessel derivative satisfies an inhomogeneous differentiated Apéry recurrence, a first-order Casoratian identity, reflection oddness, and Gessel's \(p^2\) Lucas congruence.  Those identities are written out below.  They determine the local affine target exactly, but none collapses the cross-prime CRT class to an \(H\)-short integer.  In fact the natural global quantities produced by the \(p^2\) congruence have exponential height in \(r\) or \(m\), far above the root-lattice budget.

The route is therefore:

* **proved:** exact scalar CRT reduction, exact lattice index, exact product-modulus lower bound, and exact differentiated-recurrence/Casoratian identities;
* **conditional:** an Apéry-specific exceptional-coset theorem could still exist, but it would have to be a genuinely new cross-prime theorem, not pointed simplicity or a homogeneous geometry-of-numbers consequence;
* **scoped no-go:** no argument using only the existing homogeneous shortness, reflection, Newton/Pascal identities, or a uniform affine transference bound can produce the requested nonzero \(H\)-short mismatch;
* **different sufficient theorem:** a five-prime Fourier-dispersion estimate for the actual coefficient-zero sets \(I_p\) implies the desired fifth factorial moment directly.  I give its exact formula and implication below.  It is stronger than the target, not a renaming of it.

I audited the repository at main commit
`734a5a84c1e4fd8703a811aadaa2b4c7f532b20e`.  In particular I used the exact Newton/Smith/curvature scripts listed in §10.  Finite computations below are used only to **falsify universal identities or reproduce finite affine classes**; no asymptotic theorem is inferred from a scan.

---

## 1. Exponent ledger

The parameters in the question are
\[
 H=X^{7/15}(\log X)^{-1/3},
 \qquad
 L=X^{14/15+\eta}.
\]
Therefore
\[
 H^2=X^{14/15}(\log X)^{-2/3},
 \qquad
 \boxed{\frac{L}{H^2}=X^\eta(\log X)^{2/3}\to\infty.}
 \tag{1.1}
\]
This is the affine-lattice threshold: the selected multiplicity is already larger than the square of the isolation/Newton scale by the factor in (1.1).

Also
\[
 \frac{L}{H}
 =X^{7/15+\eta}(\log X)^{1/3},
 \tag{1.2}
\]
and hence
\[
 \frac{L\log X}{H}
 =X^{7/15+\eta}(\log X)^{4/3}.
 \tag{1.3}
\]

There is only a hard regime when
\[
 0<\eta<\frac1{15}.
\]
Indeed, if \(\eta\ge1/15\), then \(L\ge X\), whereas \(K(m)\) is at most the number of primes in \((X,2X]\), which is \(<X\) for large \(X\).  Thus the selected set is then empty without any Apéry input.

In what follows assume \(0<\eta<1/15\), which is the genuinely nontrivial regime.

---

# I. PROVED: the even Hermite problem is one scalar CRT problem

## 2. Exact reflection quotient

Let
\[
 \sigma(T)=-T-1,
 \qquad
 Y=T(T+1).
\]
Then
\[
 \boxed{\mathbf Z[T]^\sigma=\mathbf Z[Y].}
 \tag{2.1}
\]
A quick integral proof is useful here.  Since \(T^2=Y-T\), every polynomial has a unique decomposition
\[
 F(T)=A(Y)+T B(Y),\qquad A,B\in\mathbf Z[Y].
\]
Applying \(\sigma\) gives
\[
 F(-T-1)=A(Y)-(T+1)B(Y).
\]
Equality with \(F(T)\) forces \(B=0\).

For a fixed row \(m\), put
\[
 y_m=m(m+1).
\]
The reflected root factor is exactly
\[
 \boxed{
 G_m(T)=(T-m)(T+m+1)=Y-y_m.}
 \tag{2.2}
\]
If an integral polynomial \(P_m\) is reflection-even and divisible by \(G_m\), then the quotient is reflection-even as well, so there is a unique \(Q_m\in\mathbf Z[Y]\) with
\[
 \boxed{
 P_m(T)=(Y-y_m)Q_m(Y).}
 \tag{2.3}
\]
Differentiating gives
\[
 P_m'(T)=(2T+1)\left(Q_m(Y)+(Y-y_m)Q_m'(Y)\right),
\]
therefore at either root
\[
 \boxed{
 P_m'(m)=(2m+1)Q_m(y_m),
 \qquad
 P_m'(-m-1)=-(2m+1)Q_m(y_m).}
 \tag{2.4}
\]
This is the key exact reduction.  Regardless of the degree of \(Q_m\), all Hermite information at the two reflected roots is carried by **one integer** \(Q_m(y_m)\).

## 3. The Gessel jet has the required reflection parity

For the standard Gessel first jet write
\[
 D(r)=A'(r)
 =2\sum_{k=0}^r\binom rk^2\binom{r+k}{k}^2
      \bigl(H_{r+k}-H_{r-k}\bigr).
 \tag{3.1}
\]
For \(r<p\), its denominator is a \(p\)-unit.  Gessel's congruence, in the integer-\(n\) form proved by Rowland--Yassawi--Krattenthaler, is
\[
 \boxed{
 A(d+pn)\equiv\bigl(A(d)+pnD(d)\bigr)A(n)\pmod{p^2}}
 \tag{3.2}
\]
for \(0\le d<p\) and \(n\in\mathbf Z\).

The entire Apéry function satisfies
\[
 A(-1-z)=A(z).
 \tag{3.3}
\]
Let \(r'=p-1-r\).  Applying (3.2) with \(d=r,n=-1\), then using (3.3), gives
\[
 A(r')=A(r-p)\equiv A(r)-pD(r)\pmod{p^2}.
\]
Interchanging \(r,r'\) gives
\[
 A(r)\equiv A(r')-pD(r')\pmod{p^2}.
\]
Hence
\[
 \boxed{D(p-1-r)\equiv-D(r)\pmod p.}
 \tag{3.4}
\]
At the central residue \(r=(p-1)/2\), this gives
\[
 \boxed{D(r)\equiv0\pmod p.}
 \tag{3.5}
\]
Thus the Gessel derivative is exactly reflection-odd, as an even Hermite interpolant requires.

The same scalar reduction below applies to any other canonical project jet \(J_p(r)\) already proved to have this reflection-odd transformation law; simply replace \(D(r)\) by \(J_p(r)\).  I do **not** identify the repository's separate Eichler/inhomogeneous coordinates with \(D(r)\) unless their normalization says so.

## 4. At most two central hit primes

For a hit prime \(p\in(X,2X]\), the residue \(m\bmod p\) is central iff
\[
 2m+1\equiv0\pmod p.
\]
Since \(m<X^2\),
\[
 0<2m+1<2X^2+1.
\]
For \(X\ge3\), three distinct primes \(>X\) cannot all divide \(2m+1\), because their product exceeds \(X^3>2X^2+1\).  Therefore
\[
 \boxed{\#\{p\text{ hit}:p\mid2m+1\}\le2.}
 \tag{4.1}
\]
At such a central prime, both sides of the Gessel Hermite condition vanish automatically: (2.4) is \(0\pmod p\), and (3.5) says the Gessel jet is \(0\pmod p\).

Let \(S_m^*\) be the noncentral hit primes and
\[
 K^*(m)=|S_m^*|,\qquad
 M_m^*=\prod_{p\in S_m^*}p.
\]
Then
\[
 K^*(m)\ge K(m)-2.
 \tag{4.2}
\]

## 5. Exact scalar CRT class

For \(p\in S_m^*\), let \(r_p=m\bmod p\).  From (2.4), matching the Gessel derivative is equivalent to
\[
 (2m+1)Q_m(y_m)\equiv D(r_p)\pmod p.
\]
Since \(2m+1\not\equiv0\pmod p\), define
\[
 a_{p,m}:=(2r_p+1)^{-1}D(r_p)\in\mathbf F_p.
 \tag{5.1}
\]
Then every noncentral Hermite condition is simply
\[
 \boxed{Q_m(y_m)\equiv a_{p,m}\pmod p.}
 \tag{5.2}
\]
Chinese remaindering gives one and only one class
\[
 \boxed{c_m\pmod{M_m^*}}
 \tag{5.3}
\]
such that
\[
 Q_m(y_m)\equiv c_m\pmod{M_m^*}
 \tag{5.4}
\]
is equivalent to all noncentral derivative conditions.

This is the exact affine target.  It is **one scalar CRT class**, not \(K(m)\) independent real directions.  But its modulus is the full product of the hit primes.

Pointed simplicity says, at most, that some or all \(a_{p,m}\) are nonzero.  It gives no control on the balanced representative of \(c_m\).  The model example \(a_{p,m}=1\) for every \(p\) has perfectly nonzero local jets yet \(c_m=1\).  Conversely other nonzero residues can make the balanced CRT representative of order \(M_m^*\).  Therefore pointed simplicity is logically orthogonal to the needed affine closest-point statement.

---

# II. PROVED: exact affine lattice index and the product-modulus barrier

## 6. The affine coefficient lattice has determinant \(M_m^*\)

Fix a degree bound \(d\) and write
\[
 Q(Y)=q_0+q_1Y+\cdots+q_dY^d,
 \qquad q\in\mathbf Z^{d+1}.
\]
Let
\[
 v_m=(1,y_m,y_m^2,\ldots,y_m^d).
\]
The homogeneous evaluation lattice is
\[
 \mathcal L_{m,d}
 =\{q\in\mathbf Z^{d+1}:q\cdot v_m\equiv0\pmod{M_m^*}\}.
 \tag{6.1}
\]
Because the first coordinate of \(v_m\) is \(1\), the homomorphism
\[
 \mathbf Z^{d+1}\longrightarrow\mathbf Z/M_m^*\mathbf Z,
 \qquad q\longmapsto q\cdot v_m
\]
is surjective.  Therefore
\[
 \boxed{
 [\mathbf Z^{d+1}:\mathcal L_{m,d}]
 =\det\mathcal L_{m,d}=M_m^*.}
 \tag{6.2}
\]
The admissible Hermite coefficients form the single affine coset
\[
 \mathcal C_{m,d}
 =\{q:q\cdot v_m\equiv c_m\pmod{M_m^*}\}
 =q^{(0)}+\mathcal L_{m,d}.
 \tag{6.3}
\]

Nothing changes by moving to an integral falling-factorial/Newton basis in \(Y\): this is a unimodular integral basis change.  Even if one uses the integer-valued basis \(\binom{Y}{j}\), the evaluation map still contains the constant coordinate \(1\), hence is still surjective modulo \(M_m^*\).  The product-modulus index is not a monomial-basis artifact.

This is the precise form of the statement “homogeneous Minkowski does not solve the affine problem.”  Minkowski controls short nonzero vectors in (6.1); it does not select a short point in the one particular translate (6.3).

## 7. A uniform affine covering theorem is impossible at the requested scale

Let \(n=d+1\).  There are at most
\[
 (2B+1)^n
\]
integer coefficient vectors with \(\|q\|_\infty\le B\).  Each vector hits only one residue class modulo \(M_m^*\).  Therefore, in order that **every** affine class modulo \(M_m^*\) have such a representative, it is necessary that
\[
 (2B+1)^n\ge M_m^*.
 \tag{7.1}
\]
Equivalently, some affine class has
\[
 \boxed{
 B\ge\frac{(M_m^*)^{1/n}-1}{2}.}
 \tag{7.2}
\]
This is just pigeonhole counting; no lattice theorem can beat it.

Since every \(p\in S_m^*\) exceeds \(X\),
\[
 \log M_m^*>K^*(m)\log X\ge(K(m)-2)\log X.
 \tag{7.3}
\]
For a selected row \(K(m)\ge L\), and for an \(O(H)\)-dimensional Newton window \(n\le C H\), (7.2)--(7.3) give the worst-coset lower bound
\[
 \log B
 \ge \frac{L+o(L)}{CH}\log X-O(1)
 =\boxed{
 \frac1C X^{7/15+\eta}(\log X)^{4/3}(1+o(1)).}
 \tag{7.4}
\]
By contrast, an \(H\)-scale homogeneous coefficient budget of the usual form
\[
 B_{\rm hom}=X^{O(H)}
\]
has
\[
 \log B_{\rm hom}=O(H\log X)
 =O\!\left(X^{7/15}(\log X)^{2/3}\right).
 \tag{7.5}
\]
The ratio of the two logarithmic scales is exactly, up to constants,
\[
 \boxed{
 \frac{L}{H^2}=X^\eta(\log X)^{2/3}\to\infty.}
 \tag{7.6}
\]

This proves a **uniform/worst-coset no-go**: no affine transference theorem that depends only on the determinant/index and the homogeneous short-vector data can put every relevant coset into the existing short box.

It does **not** prove that the special Apéry class \(c_m\) is a worst coset.  A special-coset theorem remains logically possible.  But such a theorem must use new arithmetic information about the residues \(a_{p,m}\); it cannot be a consequence of Minkowski, reflection, root counts, or basis choice alone.

### Transference does not remove this distinction

Banaszczyk-type transference bounds relate covering radius to a minimum of the dual lattice.  Here
\[
 \mathcal L_{m,d}^*
 =\mathbf Z^n+\frac1{M_m^*}v_m\mathbf Z
 \tag{7.7}
\]
(up to the obvious identification of \(v_m\) with a column vector).  A useful covering-radius theorem would therefore require a strong dual Diophantine statement about multiples of
\[
 \frac1{M_m^*}(1,y_m,\ldots,y_m^d),
\]
or, even more specifically, about the displacement of the **particular** Apéry class \(c_m\).  No such statement follows from the current Newton-shell recurrence.  In any case, a theorem claiming a uniformly smaller box than (7.2) would contradict the counting argument before transference enters.

---

## 8. The nonzero mismatch has unavoidable product height

This is stronger than the covering-radius warning and applies directly to the proposed \(\Delta_m\).

Assume, as in the question, that for a selected row \(m\) one has an integer
\[
 \Delta_m\ne0
\]
such that **every hit prime** divides \(\Delta_m\).  The hit primes are distinct, hence their product divides \(\Delta_m\):
\[
 \prod_{p:\,m\bmod p\in I_p}p\mid\Delta_m.
\]
Since every factor is \(>X\),
\[
 \boxed{
 \log |\Delta_m|>K(m)\log X.}
 \tag{8.1}
\]
For selected rows,
\[
 \boxed{
 \log |\Delta_m|\ge
 X^{14/15+\eta}\log X.}
 \tag{8.2}
\]
Therefore a nonzero mismatch satisfying an \(H\)-scale pointwise height bound
\[
 \log(2+|\Delta_m|)=O(H\log X)
\]
**cannot exist** on a selected row.  The ratio between the mandatory lower bound and the \(H\)-scale upper bound is
\[
 \frac{L}{H}
 =X^{7/15+\eta}(\log X)^{1/3}\to\infty.
 \tag{8.3}
\]

For example, if \(P_m\) has degree \(D=O(H)\), coefficient sup-height \(B=X^{O(H)}\), and the integer jet representative against which one forms the mismatch has comparable height, then for \(m<X^2\)
\[
 \log |P_m'(m)|
 \le \log D+\log B+(D-1)\log(1+m)+O(1)
 =O(H\log X),
\]
so (8.2) gives an immediate contradiction.  If one instead uses the natural huge integer/rational lift of the Gessel jet, then the contradiction disappears--but so does the hoped-for short-height argument.

This is the exact reason the affine step is not a routine extension of the root-lattice theorem.

---

# III. CONDITIONAL IMPLICATION: the proposed weighted estimate really does imply the fifth moment

Let \((K)_j=K(K-1)\cdots(K-j+1)\).  From (8.1),
\[
 \log(2+|\Delta_m|)\ge K(m)\log X
\]
for every selected \(m\).  Hence
\[
 \begin{aligned}
 (K(m))_4\log(2+|\Delta_m|)
 &\ge (K(m))_4K(m)\log X\\
 &=\bigl((K(m))_5+4(K(m))_4\bigr)\log X.
 \end{aligned}
 \tag{9.1}
\]
Therefore the requested bound
\[
 \sum_{m\,\mathrm{selected}}
 (K(m))_4\log(2+|\Delta_m|)
 \le X^{o(1)}X^2\Lambda^5\log X
 \tag{9.2}
\]
would immediately imply
\[
 \boxed{
 \sum_{m\,\mathrm{selected}}(K(m))_5
 \le X^{o(1)}X^2\Lambda^5.}
 \tag{9.3}
\]
This implication is exact.

But (9.1) also shows that (9.2) is **not a weaker intermediate estimate**.  Its left side automatically contains the desired fifth factorial moment times \(\log X\), plus the positive term \(4(K)_4\log X\).  So the weighted mismatch formulation is a legitimate sufficient theorem, but proving it requires essentially the same five-prime control unless the \(\Delta_m\)'s possess additional global structure that is not present in the current homogeneous lattice argument.

If one normalizes away the at-most-two central primes and only knows divisibility by the \(K-2\) noncentral primes, the exact identity becomes
\[
 (K)_4(K-2)=(K)_5+2(K)_4,
\]
so the conclusion is unchanged at the required scale.

---

# IV. PROVED: what the special Apéry/Gessel recurrence actually gives

The preceding obstruction is generic.  The important remaining question is whether the Apéry recurrence forces the **special** class \(c_m\) to be exceptionally close.  The exact first-jet algebra shows what information is available.

## 10. Differentiated Apéry recurrence

Write
\[
 P(n)=34n^3+51n^2+27n+5
\]
and \(b_n=A(n)\).  The ordinary recurrence is
\[
 (n+1)^3b_{n+1}-P(n)b_n+n^3b_{n-1}=0.
 \tag{10.1}
\]
The entire Apéry interpolation satisfies an inhomogeneous recurrence whose right side and derivative both vanish at integer arguments.  Therefore \(d_n:=A'(n)\) obeys
\[
 \boxed{
 (n+1)^3d_{n+1}-P(n)d_n+n^3d_{n-1}=f_n,}
 \tag{10.2}
\]
where
\[
 \begin{aligned}
 f_n
 &=P'(n)b_n-3(n+1)^2b_{n+1}-3n^2b_{n-1}\\
 &=\boxed{
 \frac{3}{n+1}
 \bigl((17n^2+16n+4)b_n-n^2b_{n-1}\bigr).}
 \end{aligned}
 \tag{10.3}
\]
The initial values are
\[
 d_0=0,\qquad d_1=12.
 \tag{10.4}
\]

This is already enough to compute the exact Gessel jet with rational arithmetic and no harmonic-number summation.

## 11. Exact jet Casoratian identity

Define
\[
 W_n=b_nd_{n+1}-b_{n+1}d_n.
\]
Multiply (10.2) by \(b_n\), multiply (10.1) by \(d_n\), and subtract.  One gets
\[
 \boxed{
 (n+1)^3W_n-n^3W_{n-1}=b_nf_n.}
 \tag{11.1}
\]
Since \(W_0=12\), summing gives
\[
 \boxed{
 r^3W_{r-1}
 =12+\sum_{j=1}^{r-1}b_jf_j.}
 \tag{11.2}
\]
At a coefficient zero \(b_r\equiv0\pmod p\), consecutive zeros are impossible in the interior because (10.1) would propagate backwards to \(b_0\equiv0\).  Hence \(b_{r-1}\) is a \(p\)-unit, and (11.2) becomes
\[
 \boxed{
 d_r\equiv
 b_{r-1}^{-1}r^{-3}
 \left(12+\sum_{j=1}^{r-1}b_jf_j\right)
 \pmod p.}
 \tag{11.3}
\]
This is a concrete determinant identity for the canonical Gessel jet at a zero.

What it does **not** do is make the right side a bounded-degree or \(H\)-height function of the root location.  It is an accumulated prefix quantity.  The natural characteristic-zero height of the prefix is exponential in \(r\), not exponential in \(H\).

## 12. Gessel's \(p^2\) Lucas law is the strongest obvious cross-prime bridge--and it still leaves local quotient digits

For \(m=qp+r\), \(0\le r<p\), Gessel's congruence (3.2) gives
\[
 \boxed{
 b_m\equiv\bigl(b_r+pq\,d_r\bigr)b_q\pmod{p^2}.}
 \tag{12.1}
\]
If \(p\mid b_r\), Lucas modulo \(p\) also gives \(p\mid b_m\).  Dividing (12.1) by \(p\) yields the exact first-digit relation
\[
 \boxed{
 \frac{b_m}{p}
 \equiv b_q\left(\frac{b_r}{p}+q d_r\right)
 \pmod p.}
 \tag{12.2}
\]
Here the quotients are ordinary integers before reduction modulo \(p\).  On the subcase \(p\nmid q b_q\), one may solve for the jet:
\[
 \boxed{
 d_r\equiv q^{-1}
 \left(b_q^{-1}\frac{b_m}{p}-\frac{b_r}{p}\right)
 \pmod p.}
 \tag{12.3}
\]

This is important because it shows exactly where the affine information lives: in the **first \(p\)-adic quotient digits** \(b_m/p\) and \(b_r/p\).  Those change with \(p\).  Equation (12.3) is not a single small integer whose reductions give all the jets.

The obvious global lifts are also on the wrong height scale:

* from the harmonic formula, after clearing denominators up to \(2r\), the numerator/denominator height of \(d_r\) is \(\exp(O(r))\), hence \(\exp(O(X))\) in this prime range;
* \(b_m\) has \(\log b_m=O(m)\), hence up to \(O(X^2)\).

In the hard range \(0<\eta<1/15\), the mandatory product logarithm
\[
 L\log X=X^{14/15+\eta}\log X
\]
is \(o(X)\).  Thus the natural Gessel jet lift is already much too large to yield the required contradiction, while the \(b_m\) lift is vastly larger still.  The affine closest-point problem is precisely the demand to replace these natural lifts by a much smaller representative of the same simultaneous residue data.  Neither (10.2), (11.3), nor (12.2) performs that replacement.

---

# V. SCOPED NO-GO FROM THE CURRENT NEWTON-SHELL MACHINERY

The repository contains several exact audits that sharply delimit what the Newton identities themselves can do.  These are useful because they are not heuristic scans; they contain symbolic identities and/or explicit finite counterexamples to universal claims.

## 13. The recurrence residual has only a quartic quotient, but a hit does not force an extra factor

`problems/3.2/research/scripts/q32_newton_recurrence_residual.py` forms an interpolation polynomial \(F\) through consecutive actual Apéry values and its recurrence residual
\[
 \mathcal R_F(x)
 =(x+1)^3F(x+1)-P(x)F(x)+x^3F(x-1).
 \tag{13.1}
\]
Because \(F\) agrees with the Apéry sequence at every interior node, the full interior-node product divides \(\mathcal R_F\), leaving a quotient of degree at most \(4\).  That part is exact.

The same script then checks the most tempting extra-factor claim on the three reflected hits at \(n=321\):
\[
 (r,p)=(36,179),(64,193),(100,211).
\]
For all three, the quartic quotient evaluated at the target is a \(p\)-unit.  Thus **the statement “a coefficient-zero hit forces one more recurrence-residual factor” is false**.  This finite computation is used only as a counterexample to that universal identity; it is not evidence about asymptotic frequencies.

## 14. Exact Newton curvature prevents endpoint collapse from adjacent local data

`q32_cross_curvature_audit.py` proves the universal beta--Padé/Newton difference law
\[
 A_d(F)-A_{d+1}(F)=P_dJ_d(F)
\]
and the exact endpoint-curvature identity
\[
 \det(a_0,a_r)
 =\sum_sP_sE_s
 +\sum_{t<s}P_tP_s\det(j_s,j_t).
 \tag{14.1}
\]
It then constructs explicit integral arrays for which every adjacent normalized minor has a prescribed prime factor while the endpoint determinant does not.  Consequently an endpoint Hermite collapse cannot follow from the universal Newton/Pascal kernel alone.  An Apéry-specific identity would have to kill the curvature term in (14.1) for arithmetic reasons not present in the formal interpolation algebra.

## 15. Shared-node factors can be presentation content, not arithmetic information

`q32_beta_pade_packet_audit.py` proves the exact incomplete-beta packet identity and shows that, when two Newton carriers share a prime node, certain coefficientwise differences are divisible by that prime for **every input sequence**.  Such factors are presentation content.  They cannot be counted as a new Apéry/Gessel derivative condition.

## 16. The second-layer scalar is genuinely local

`q32_newton_second_layer_audit.py` proves for the shortest top-half carriers
\[
 G_{p-2,1}-G_{p-1,1}=p\,\Delta^2Y_{p-2},
\]
and, at a target,
\[
 \boxed{
 \frac{G_{p-2,1}-G_{p-1,1}}p
 \equiv Y_{p-2}+Y_p\pmod p,
 \qquad Y_p\equiv40b_{r-1}\pmod p.}
 \tag{16.1}
\]
A proof that this scalar never vanishes would be a useful local transversality/pointed-simplicity statement.  But it would still only say that each local affine residue is nonzero.  Section 5 shows why that gives no control on the simultaneous CRT representative \(c_m\).

## 17. Smith calculations confirm that the prime index survives the Newton family

`q32_newton_gcd_audit.py`, `q32_margin_lattice_audit.sage`, `q32_fixed_minor_jet_audit.sage`, and `q32_ghost_quotient_audit.sage` compute exact local Smith/determinantal divisors for several Newton families.  The common pattern is the same as (6.2): a target node contributes a genuine prime index after the presentation-content factors are removed.  Adding translated carriers or ghost coordinates can change the available lattice, but it does not turn the inhomogeneous target into the zero coset for free.

For example, `q32_newton_gcd_audit.py` has an exact hostile case at \(n=321\) where the gcd of the relevant carriers is exactly
\[
 179\cdot193\cdot211,
\]
while the carriers themselves have 528--529 decimal digits.  Again, that is a finite exact regression, not an asymptotic theorem; its relevance is that the formal Newton machinery does not manufacture a hidden extra common factor even in the canonical multi-hit example.

### Scoped no-go statement

The following is what is actually ruled out:

> **No-go.**  From the currently banked homogeneous root-lattice shortness, reflection \(T\mapsto-T-1\), the universal Newton/Pascal/beta--Padé identities, and standard lattice transference alone, one cannot deduce an \(H\)-short representative of the Apéry Hermite affine class.  Nor can one obtain a nonzero \(H\)-short integer mismatch divisible by all selected hit primes.  The latter is arithmetically impossible by (8.1).

What is **not** ruled out is a new theorem asserting that the one special Apéry class \(c_m\) satisfies an exceptional congruence across all its defining primes.  Such a theorem would be new mixed-characteristic arithmetic, not geometry-of-numbers bookkeeping.

---

# VI. What an Apéry-specific affine rescue theorem would have to say

The exact missing special-coset statement can now be written without ambiguity.

For each selected \(m\), let \(S_m^*\), \(M_m^*\), \(y_m\), and \(c_m\) be as in §§4--5.  A genuine affine rescue would need something of the form
\[
 \boxed{
 \exists\,Q_m\in\mathbf Z[Y],\quad
 \deg Q_m\le C H,\quad
 H_{\rm coeff}(Q_m)\le X^{C H},\quad
 Q_m(y_m)\equiv c_m\pmod{M_m^*}.}
 \tag{AR}
\]
This is **not** a consequence of the homogeneous lattice theorem.  By §7 it can only hold because the specific Apéry class is exceptional.

There is an even sharper issue.  If (AR) plus the rest of the construction made a nonzero integer mismatch \(\Delta_m\) with \(\log|\Delta_m|=O(H\log X)\) and all hit primes dividing it, §8 would force a contradiction.  Thus in the selected range an \(H\)-short affine lift cannot lead to a *small nonzero* mismatch: it would force the mismatch to be exactly zero.

So a successful special identity would more naturally look like an **exact compatibility identity**
\[
 \Delta_m=0,
\]
followed by a separate argument showing that this exact identity is impossible for a selected row, or that it forces an algebraic degeneracy with controlled multiplicity.  Merely asking for a short nonzero \(\Delta_m\) is inconsistent with the product modulus once \(K(m)\ge L\).

Equivalently, one could seek a direct cross-prime formula for the normalized jets
\[
 a_{p,m}=(2r_p+1)^{-1}D(r_p)
\]
of the form
\[
 a_{p,m}\equiv Z_m\pmod p
 \qquad(p\in S_m^*)
 \tag{17.1}
\]
with one integer \(Z_m\) of logarithmic height \(O(H\log X)\).  Neither the differentiated recurrence nor the \(p^2\) Gessel law provides such a \(Z_m\).  Formula (12.2) instead exhibits \(p\)-dependent quotient digits.

This is the exact theorem to search for if the affine route is pursued further.  It is a cross-prime Gessel-jet congruence theorem, not pointed simplicity.

---

# VII. A different exact theorem that really implies the five-prime bound

If the affine route is abandoned, the cleanest exact replacement is a **five-prime Fourier-dispersion theorem for the coefficient-zero sets themselves**.  This does not use cover-variable Hasse roots.

Put
\[
 z_p=|I_p|,
 \qquad
 F_p(a)=\sum_{r\in I_p}e_p(-ar),
 \qquad e_q(x)=e^{2\pi i x/q}.
\]
For an ordered tuple of five distinct primes
\[
 \mathbf p=(p_1,\ldots,p_5),
 \qquad Q=p_1\cdots p_5,
\]
let
\[
 N(\mathbf p)
 =\#\{0\le m<X^2:m\bmod p_i\in I_{p_i}\ \forall i\}.
 \tag{18.1}
\]
Let
\[
 W_M(h;Q)=\sum_{0\le m<M}e_Q(hm),
 \qquad M=X^2.
\]
CRT Fourier inversion gives the **exact identity**
\[
 \boxed{
 N(\mathbf p)
 =\frac{M}{Q}\prod_{i=1}^5z_{p_i}
 +E(\mathbf p),}
 \tag{18.2}
\]
where
\[
 \boxed{
 E(\mathbf p)
 =\frac1Q\sum_{h=1}^{Q-1}W_M(h;Q)
   \prod_{i=1}^5
   F_{p_i}\!\left(
      h\,(Q/p_i)^{-1}\bmod p_i
   \right).}
 \tag{18.3}
\]
Changing the sign convention in \(F_p\) only conjugates the local factors.

Now consider the following genuinely stronger statement:

> **(FD5) Five-prime Fourier dispersion.**
> \[
> \boxed{
> \sum_{\substack{p_1,\ldots,p_5\in\mathcal P_X\\
>                  p_i\ \mathrm{pairwise\ distinct}}}
> |E(p_1,\ldots,p_5)|
> \le X^{o(1)}X^2\Lambda^5.}
> \tag{FD5}
> \]

This is not the desired factorial-moment estimate under a new name.  It controls an **absolute five-linear Fourier discrepancy tuple by tuple before the tuple sum**.  It is strictly stronger than what the nonnegative fifth moment asks for.

Nevertheless it implies the desired estimate immediately.  Indeed
\[
 \sum_{0\le m<M}(K(m))_5
 =\sum_{\mathbf p}^{*}N(\mathbf p),
 \tag{18.4}
\]
where \(^*\) means ordered distinct tuples.  Summing the main term in (18.2),
\[
 \sum_{\mathbf p}^{*}
 \frac{M}{Q}\prod_i z_{p_i}
 \le
 M\left(\sum_{p\in\mathcal P_X}\frac{z_p}{p}\right)^5
 =X^2\Lambda^5.
 \tag{18.5}
\]
Under (FD5), (18.2)--(18.5) give
\[
 \sum_{m<X^2}(K(m))_5
 \le X^{o(1)}X^2\Lambda^5.
 \tag{18.6}
\]
The selected sum is smaller, hence
\[
 \boxed{
 \sum_{m\,\mathrm{selected}}(K(m))_5
 \le X^{o(1)}X^2\Lambda^5.}
\]

This is a real alternate theorem: it moves the missing arithmetic from an affine closest-point problem into a mixed-prime five-linear correlation estimate.  The special Apéry recurrence can legitimately enter there through the local transforms \(F_p\), Casoratian gap identities, or higher local energies.  But the repository's earlier Fourier audit (`drops/Q7901-eed6cfc4.md`) already shows that reflection plus Parseval alone does not force the analogous mixed-prime cancellation even in lower order.  So (FD5) needs a genuinely new cross-prime decorrelation theorem; it is not supplied by one-prime energy bounds.

A less absolute five-linear operator estimate could suffice and might be closer to sharpness, but (FD5) has the advantage that its implication is one line and contains no disguised use of the desired fifth moment.

---

# VIII. Reproducible exact affine/jet calculation

The following Sage code decides the finite affine CRT class for any supplied list of **actual coefficient-zero hits** \((p,r)\) of a fixed \(m\).  It never substitutes cover-variable/Hasse roots.  It computes the Gessel derivative from the differentiated recurrence (10.2), checks the coefficient-zero condition, checks reflection oddness, forms the normalized scalar CRT class \(c_m\), and optionally solves the exact finite \(\ell_\infty\) closest-point problem for a chosen degree using Sage's integer linear programming backend.

The `hits` list should come from the already-defined \(I_p\) data; this script deliberately does not guess the repository's precise cyclic/noncyclic convention for \(H\)-isolation.

```sage
# q8215_affine_crt_decider.sage
# Exact arithmetic only.  Finite decision aid, not an asymptotic proof.

from sage.all import *


def apery_and_gessel_jet(N):
    """Return b_n=A(n), d_n=A'(n) in QQ for 0<=n<=N."""
    b = [QQ(1), QQ(5)]
    d = [QQ(0), QQ(12)]
    if N == 0:
        return b[:1], d[:1]
    for n in range(1, N):
        P = 34*n^3 + 51*n^2 + 27*n + 5
        bn1 = (P*b[n] - n^3*b[n-1]) / (n+1)^3
        assert bn1.denominator() == 1
        f = (
            (102*n^2 + 102*n + 27)*b[n]
            - 3*(n+1)^2*bn1
            - 3*n^2*b[n-1]
        )
        # Independent check of the simplified forcing term.
        f2 = QQ(3, n+1) * (
            (17*n^2 + 16*n + 4)*b[n] - n^2*b[n-1]
        )
        assert f == f2
        dn1 = (P*d[n] - n^3*d[n-1] + f) / (n+1)^3
        b.append(bn1)
        d.append(dn1)
    return b, d


def mod_QQ(x, p):
    x = QQ(x)
    num = ZZ(x.numerator()) % p
    den = ZZ(x.denominator()) % p
    assert den != 0
    return ZZ(num * inverse_mod(den, p) % p)


def crt_incremental(residues, moduli):
    c = ZZ(0)
    M = ZZ(1)
    for a, p in zip(residues, moduli):
        p = ZZ(p)
        a = ZZ(a) % p
        t = ((a - c) % p) * inverse_mod(M % p, p) % p
        c += M*t
        M *= p
        c %= M
    return c, M


def balanced(c, M):
    c = ZZ(c % M)
    return c if 2*c <= M else c-M


def analyze(m, hits, degree=None):
    """
    hits = [(p,r), ...] with p in (X,2X], r=m mod p,
    and r an ACTUAL coefficient-zero residue in I_p.
    """
    maxp = max(p for p, _ in hits)
    b, d = apery_and_gessel_jet(maxp-1)

    residues = []
    moduli = []
    central = []

    for p, r in hits:
        p, r = ZZ(p), ZZ(r)
        assert m % p == r
        assert mod_QQ(b[r], p) == 0      # actual Apéry coefficient zero

        rp = p-1-r
        assert mod_QQ(d[rp] + d[r], p) == 0  # reflection oddness

        if (2*r+1) % p == 0:
            assert mod_QQ(d[r], p) == 0
            central.append(p)
            continue

        a = mod_QQ(d[r], p) * inverse_mod((2*r+1) % p, p) % p
        residues.append(ZZ(a))
        moduli.append(p)

    c, M = crt_incremental(residues, moduli)
    cb = balanced(c, M)
    y = ZZ(m)*(m+1)

    print("m", m)
    print("K", len(hits), "noncentral", len(moduli), "central", central)
    print("M_noncentral", M)
    print("log_M", RR(log(M)) if M > 1 else RR(0))
    print("balanced_CRT", cb)
    print("log_abs_balanced_CRT",
          RR(log(max(1, abs(cb)))))

    # Exact product-modulus lower bound for any nonzero all-hit mismatch.
    M_all = prod(ZZ(p) for p, _ in hits)
    print("M_all", M_all)
    print("log_M_all", RR(log(M_all)))

    if degree is None:
        return

    # Exact finite closest-point problem:
    # minimize B subject to sum_j q_j*y^j = c (mod M), |q_j|<=B.
    # This checks the affine evaluation coset only.  Any additional
    # root-lattice constraints can be added as further linear equations.
    n = degree + 1
    mip = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    q = mip.new_variable(integer=True)
    t = mip.new_variable(integer=True)
    B = mip.new_variable(integer=True, nonnegative=True)

    lhs = sum(q[j] * y^j for j in range(n)) - M*t[0]
    mip.add_constraint(lhs, min=c, max=c)
    for j in range(n):
        mip.add_constraint(q[j] - B[0], max=0)
        mip.add_constraint(-q[j] - B[0], max=0)
    mip.set_objective(B[0])
    optimum = mip.solve()
    coeffs = [ZZ(round(mip.get_values(q[j]))) for j in range(n)]
    Bmin = ZZ(round(mip.get_values(B[0])))

    assert (sum(coeffs[j]*y^j for j in range(n)) - c) % M == 0
    assert max(abs(v) for v in coeffs) <= Bmin
    print("degree", degree, "exact_min_supnorm", Bmin)
    print("coeffs", coeffs)


# Example usage: replace by one selected row and its actual I_p hit list.
# analyze(m=..., hits=[(p1,r1), (p2,r2), ...], degree=...)
```

For a selected row with \(K\) hits, this finite calculation should be read against the theorem-level lower bounds
\[
 \log M_{\rm all}>K\log X,
 \qquad
 \log M_m^*>(K-2)\log X.
\]
A surprisingly small `balanced_CRT` or small exact `Bmin` for one row is **not** a theorem; the asymptotic affine route would require a proof of such smallness uniformly or with the exact weighted aggregate needed for all selected rows.

### Existing exact repository regressions

From the repository root, the following already-existing commands reproduce the relevant finite exact checks:

```bash
python3 problems/3.2/research/scripts/q32_newton_recurrence_residual.py
python3 problems/3.2/research/scripts/q32_newton_gcd_audit.py --extended --blocks
python3 problems/3.2/research/scripts/q32_newton_second_layer_audit.py --limit 5000
sage problems/3.2/research/scripts/q32_margin_lattice_audit.sage
```

Interpretation is deliberately narrow:

* the first command falsifies a specific universal “hit gives an extra quartic-residual factor” identity;
* the gcd/Smith commands verify exact finite Newton-lattice structure and hostile examples;
* the second-layer command tests the distinguished local divided digit;
* none of them is used to infer a density, moment, or asymptotic theorem.

---

# IX. Decision tree for the affine route

The issue can now be decided by one of three theorem-level outcomes.

### A. Exact special compatibility

Prove a formula forcing the normalized Gessel residues \(a_{p,m}\) into one \(H\)-height integer class, or prove (AR) directly from Apéry-specific arithmetic.  Because of §8, a selected row cannot then produce a small **nonzero** all-hit mismatch.  The next step must be to exploit the resulting exact zero identity.

### B. Exceptional-coset lower bound

Prove that for the actual Apéry class \(c_m\), not merely for a worst class,
\[
 \operatorname{dist}_\infty(\mathcal C_{m,d},0)>X^{C H}
\]
for every selected row (or for all but a weighted-negligible set).  This would rigorously kill the affine route.  The dual form would be a Diophantine approximation theorem for (7.7) correlated with the target \(c_m\).  No such theorem is currently in the repository.

### C. Bypass affine CVP

Prove a mixed-prime theorem such as (FD5).  This attacks exactly the five-prime factorial moment and does not require turning homogeneous root polynomials into affine Hermite interpolants.

Given the present code and identities, **C is the cleaner strategic target**.  The existing Newton machinery is very effective at local algebra, Smith structure, and exact target preservation, but every audited attempt to promote that local information to a cross-prime endpoint/affine collapse encounters either presentation content, curvature, or the product-modulus index.  The missing theorem is therefore mixed-prime distribution, not pointed simplicity.

---

# X. Source/provenance audit

Repository sources read at main commit `734a5a84c1e4fd8703a811aadaa2b4c7f532b20e` include:

* `problems/3.2/AVENUES.md` -- Apéry recurrence/reflection campaign state;
* `problems/3.2/RUN_LOG_P32.md` -- order-drop/Newton campaign audit trail;
* `problems/3.2/FABLE_SECTION_coinc_target.tex` -- exact gap-Casoratian and mixed-correlation target;
* `problems/3.2/pairpalm_result.tex` -- exact factorial/Palm identities and warning against equivalent renamings;
* `SUBMIT/3.2/FINAL_AUDIT.md` -- current theorem-status audit;
* `drops/Q7359-662fb1bd.md` -- interpolation/Padé primitivity obstruction;
* `drops/Q7901-eed6cfc4.md` -- reflection + Parseval does not supply mixed-prime dispersion;
* `problems/3.2/research/scripts/q32_newton.py`;
* `q32_newton_recurrence_residual.py`;
* `q32_newton_second_layer_audit.py`;
* `q32_newton_gcd_audit.py`;
* `q32_cross_curvature_audit.py`;
* `q32_beta_pade_packet_audit.py`;
* `q32_margin_lattice_audit.sage`;
* `q32_fixed_minor_jet_audit.sage`;
* `q32_ghost_quotient_audit.sage`.

Primary external facts used:

* Ira Gessel, **Some congruences for Apéry numbers**, *Journal of Number Theory* 14 (1982), 362--368: Lucas and \(p^2\) Gessel congruences.
* Eric Rowland, Reem Yassawi, Christian Krattenthaler, **Lucas congruences for the Apéry numbers modulo \(p^2\)**, *Integers* 21 (2021), A20, arXiv:2005.04801: the entire Apéry interpolation, reflection \(A(-1-z)=A(z)\), denominator \(p\)-integrality of \(A'(d)\), and the integer-\(n\) form of the \(p^2\) congruence.
* W. Banaszczyk, **New bounds in some transference theorems in the geometry of numbers**, *Mathematische Annalen* 296 (1993), 625--635: standard primal/dual transference context.  The no-go in §7 does not depend on its constants; the pigeonhole lower bound is elementary.

---

# Bottom line

The answer to the strategic question is **no, not from the currently available special Newton-shell recurrence**.  The recurrence gives a real first-jet structure, but it does not solve the affine closest-point problem.  Reflection actually makes the obstruction more transparent: for a fixed \(m\), all noncentral Hermite conditions are one scalar CRT class modulo a modulus of logarithm at least \((K(m)-2)\log X\).  At the selected threshold \(K\ge X^{14/15+\eta}\), that modulus is far beyond the \(H\)-scale homogeneous height.

Most importantly, any nonzero integer mismatch divisible by all hits already has logarithmic height at least \(K\log X\).  Consequently the proposed weighted mismatch estimate is a **valid but stronger reformulation of the five-prime burden**, not a route around it.  A successful affine attack now requires a genuinely new Apéry-specific cross-prime compatibility theorem for the special CRT class and, because a short all-hit mismatch would then have to vanish exactly, an additional rigidity theorem after that.

Absent such an identity, the mathematically clean alternative is the five-prime Fourier-dispersion theorem (FD5), whose exact CRT expansion proves the required factorial moment without cover-variable substitutions and without hiding the desired estimate behind a new name.