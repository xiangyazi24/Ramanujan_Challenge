# Terminal Turan and Hankel target carriers

This note records a scalar cumulative-family construction following the
selector theorem (68.205).  It uses only
\[
Y_d=C_M(d),
\]
not the two-coordinate \(Y/W\) shell.

Put
\[
n=M+1,\qquad d_0=\lfloor M/2\rfloor+1,\qquad L_0=M-d_0,
\]
\[
f_L=G_{M-L,L}(C_M),\qquad F_j=f_{L_0-j}.
\]
With
\[
c_L=(-1)^L\binom nL,\qquad
B_L=\Delta^LC_M(M-L),
\]
the exact terminal difference law is
\[
f_L-f_{L-1}=c_LB_L. \tag{1}
\]

## 1. Adjacent Turan carrier

For \(1\le j<L_0\), set \(L=L_0-j\) and define
\[
\mathcal E_j=
\frac{F_{j-1}F_{j+1}-F_j^2}
{\gcd\{\binom nL,\binom n{L+1}\}}. \tag{2}
\]
This is an integer.  Indeed,
\[
F_{j-1}F_{j+1}-F_j^2
=f_{L+1}f_{L-1}-f_L^2
\]
is an integer linear combination of \(c_L\) and \(c_{L+1}\) by
(1).

Let \(p=n-r\) be a top-half target.  If
\[
p>d_0+K+1, \tag{3}
\]
then \(p\) divides \(F_0,\ldots,F_{K+1}\) by the exact selector
theorem.  Equivalently \(r<L_0-K\).  Also
\[
v_p\binom nL=v_p\binom n{L+1}=1
\]
for every \(j\le K\).  The numerator in (2) has valuation at least
two and the denominator has valuation exactly one.  Therefore
\[
\boxed{\quad
p\mid\gcd(\mathcal E_1,\ldots,\mathcal E_K)
\quad(p>d_0+K+1).
\quad} \tag{4}
\]
Thus a fixed \(K\) misses only a fixed-width strip above \(n/2\).

For the four hostile rows the last two normalized minors have:
\[
\begin{array}{c|c}
n&\gcd(\mathcal E_1,\mathcal E_2)\\ \hline
200&2^2\cdot5\cdot139\cdot181\\
272&2\cdot191\cdot233\\
300&11^2\cdot191\cdot227\\
321&179\cdot193\cdot211.
\end{array} \tag{5}
\]
The individual minors have respectively about
\(2112,2885,3186,3411\) bits, while their gcds have
\(19,17,23,23\) bits.

## 2. Exact first boundary resultant

Put
\[
a=F_0-F_1=c_{L_0}B_{L_0},\quad
b=F_1-F_2=c_{L_0-1}B_{L_0-1},\quad
c=F_2-F_3=c_{L_0-2}B_{L_0-2},
\]
and let
\[
A_0=F_0F_2-F_1^2,\qquad A_1=F_1F_3-F_2^2.
\]
Writing \(S=F_0\) gives
\[
A_0=S(a-b)-a^2,\qquad
A_1=(S-a)(b-c)-b^2.
\]
Eliminating \(S\) yields
\[
\boxed{\quad
(b-c)A_0-(a-b)A_1=b(ac-b^2).
\quad} \tag{6}
\]

Let
\[
q_0=\gcd(c_{L_0},c_{L_0-1}),\qquad
q_1=\gcd(c_{L_0-1},c_{L_0-2}).
\]
Since \(A_0=q_0\mathcal E_1\), \(A_1=q_1\mathcal E_2\),
\(q_1\mid b-c\), and \(q_0\mid a-b\), (6) becomes
\[
\boxed{\quad
\frac{b-c}{q_1}\mathcal E_1
-\frac{a-b}{q_0}\mathcal E_2
=\frac{b(ac-b^2)}{q_0q_1}.
\quad} \tag{7}
\]

This exact Bezout identity also gives a sharp negative verdict for
the first boundary elimination.  Every candidate prime
\(p>d_0+3\), target or not, divides each of \(a,b,c\) and each of
\(q_0,q_1\) exactly once.  The right side of (7) therefore retains
one universal copy of \(p\).  Its candidate-prime part is the full
Pascal tail, not the target-selected part.

## 3. Full-family rank-one alias

Put
\[
a_j=F_{j-1}-F_j,\qquad
A_j=F_{j-1}F_{j+1}-F_j^2,\qquad
q_j=\gcd(c_{L_0-j},c_{L_0-j+1}).
\]
The exact identity
\[
A_j=F_j(a_j-a_{j+1})-a_ja_{j+1} \tag{8}
\]
controls every scale.  At a candidate prime in the common selector
tail, write
\[
a_j=p\alpha_j,\qquad q_j=pu_j
\]
with \(u_j\) a unit.  Since all \(F_j\)'s select the same scalar
\(\tau_p=C_M(p-1)\), division of (8) gives
\[
\boxed{\qquad
\mathcal E_j\equiv
\tau_pu_j^{-1}(\alpha_j-\alpha_{j+1})\pmod p.
\qquad} \tag{9}
\]

Hence the whole normalized Turan family remains a rank-one alias.
Any target-blind linear combination which eliminates \(\tau_p\)
vanishes for every candidate prime and reacquires the common
primorial.  A combination which does not eliminate it is merely
another lift of the original marked scalar.

The general adjacent resultant is
\[
(a_{j+1}-a_{j+2})A_j-(a_j-a_{j+1})A_{j+1}
=a_{j+1}(a_ja_{j+2}-a_{j+1}^2). \tag{10}
\]
After its two Pascal normalizations, the right side retains one
universal candidate copy.  This proves a sharp no-go for scalar
linear elimination across the full terminal family.

## 4. Multiscale Hankel normalization

Let
\[
D_k=\det(F_{i+j})_{0\le i,j<k}.
\]
After unitriangular finite differences in rows and columns, the
\((i,j)\)-entry is \(\Delta^{i+j}F_0\).  Put
\[
q_0=1,\qquad
q_s=\gcd(c_{L_0},c_{L_0-1},\ldots,c_{L_0-s+1})
\quad(s\ge1).
\]
Then \(q_s\mid\Delta^sF_0\).  Since \(q_{s+1}\mid q_s\), a
primewise assignment argument shows that the certified Pascal
divisor
\[
U_k=\prod_{s=k}^{2k-2}q_s \tag{11}
\]
divides \(D_k\).  This is a sufficient universal divisor; additional
coefficient content may occur.

If
\[
p>d_0+2k-2 \tag{12}
\]
is a target, all entries \(F_0,\ldots,F_{2k-2}\) are divisible by
\(p\), so \(p^k\mid D_k\).  On the other hand, every \(q_s\) in
(11) has \(p\)-valuation one, and hence
\[
v_p(U_k)=k-1,\qquad
p\mid T_k:=D_k/U_k. \tag{13}
\]
The hostile examples have \(v_p(T_k)=1\) for every target and
\(2\le k\le7\).

The gcd of \(T_2,T_3\) is again tiny, but this is not an independent
mechanism.  If
\[
A_i=F_iF_{i+2}-F_{i+1}^2,
\]
Desnanot--Jacobi gives
\[
F_2D_3=A_0A_2-A_1^2. \tag{14}
\]
Thus, away from the middle normalization singularity, a common prime
of \(T_2,T_3\) is forced into the adjacent Turan family.  At a target
the middle normalization is singular, which is exactly where the one
target copy survives.

## 5. Present verdict

Equations (2)--(13) give a rigorous \(J\)-dependent normalization:
the Pascal tail is removed to multiplicity \(k-1\), while each target
is retained at least once.  It does not yet give the required height
bound.  Direct calculation shows
\[
\log_2|T_k|\asymp k n
\]
on the tested range; increasing \(k\) makes the individual carrier
worse.  The small characteristic-zero gcd in (5) remains an
Apéry-specific gcd problem.

There is also a sharp formal limitation.  For an arbitrary terminal
array, after imposing the marked scalar \(F_j\equiv0\pmod p\), the
once-divided tail coordinates can be prescribed freely through the
unitriangular binomial transform \(B_L\mapsto C_M(M-L)\).  Hence
Pascal/Hankel algebra alone supplies no second target equation.
Any proof that the gcd in (5) has subexponential height must use the
distinguished Laurent coefficient coupling or an arithmetic theorem
for the actual terminal sequence.

The exact reproducible audit is
`problems/3.2/research/scripts/q32_terminal_turan_hankel_audit.py`.
