PARTIAL BREAKTHROUGH — THE GLOBAL LIMITING APERY CELL AND \(\cosh(2\pi\eta)=2\) ARE PROVED, BUT THE POINTWISE GCD THEOREM, THE FIXED-PRIME \(S_D=o(p)\) ESTIMATE, AND ALL-HEIGHT CRITICAL-VALUE SEPARATION REMAIN OPEN.

# Post-submission unconditional report

## 1. What changed

This report records the conclusions that survived the final audit.  It is not a
campaign ledger and it does not promote numerical fits to theorems.  There are
four distinct endpoints:

1. The original pointwise statement
   \[
   G_n:=\gcd(d_na_n,d_nb_n)=e^{o(n)},\qquad
   d_n=\operatorname{lcm}(1,\ldots,n)^3,
   \]
   remains open.
2. The local Apéry connection curvature is now exact:
   \[
   (\log\kappa)''(0)=-4\zeta(2)=-\frac{2\pi^2}{3}.
   \]
3. Restoring the companion term in the cell and combining the exact Apéry
   difference equation with the global Frobenius function of
   Golyshev--Zagier and Bloch--Vlasenko closes that bridge.  Unconditionally,
   \[
   J(z)=\pi^3\cot^3(\pi z)+\frac{\pi^3}{3}\cot(\pi z).
   \]
4. This global limiting formula does not supply the uniform finite-height
   remainder needed for all-height critical-value separation; that problem
   remains open beyond the certified range \(2\le h\le60\).

The verification program `RC_BREAKTHROUGH_verify.py` checks all finite algebra
used here.  Its final line is `PASS`.  Analytic continuation and convergence
inputs are cited to the primary sources, while the additional vertical-growth
bound is derived explicitly below; none is inferred from finite sampling.

## 2. The arithmetic wall did not move

For a fixed prime \(p\ge7\), put \(N=p-2\).  Let \(c_0=0,c_1=6\)
be the second solution of the same Apéry recurrence as \(b_n\), and set
\(\pi(r)=[b_r:c_r]\in\mathbf P^1(\mathbf F_p)\).  Write

\[
 R_h=\#\{1\le r\le N-h:\pi(r)=\pi(r+h)\},\qquad
 S_D=\sum_{h\le D}R_h.
\]

Equivalently, if

\[
 d_D(r)=\#\{1\le h\le D:r+h\le N,\ \pi(r)=\pi(r+h)\},
 \qquad Q_D=\sum_r\binom{d_D(r)}2,
\]

then \(S_D=\sum_r d_D(r)\) and
\(\sum_r d_D(r)^2=S_D+2Q_D\).

At the mesoscopic scale \(D=\sqrt N L\), where \(L\to\infty\),
\(L=N^{o(1)}\), and \(D=o(N)\), the current unconditional estimate is

\[
 S_D\ll N L\sqrt{\log N}.
\]

The target needed to improve the fixed-prime energy record is \(S_D=o(N)\).
Thus the remaining loss is subpower but unbounded.  The exact identity above
and Cauchy--Schwarz show why merely deleting the logarithm from the abstract
bound for \(Q_D\) is not enough: one must also gain against the factor \(L\).

Double additive-character orthogonality gives an exact decomposition of
\(Q_D\), but it retains two unsigned boundary-weight terms.

Reflection interchanges those terms and makes the centered bias real; it does
not make the bias negative or small.  Because no estimate for those terms is
proved, the decomposition is not promoted into the theorem list here.  None
of parity, renewal, squarefreeness, generic symplectic transport, or the
finite-height critical certificates supplies the missing estimate.

The proposed uniform cross-cell expansion also needs a correction.  A term
\(h^{-1}G(\alpha)\) has adjacent variation of order \(h^{-2}G'(\alpha)\) in
general.  Reflection cancels it only in the central symmetric comparison.
Consequently an assertion of a uniform
\(O(h^{-2-\delta})\) remainder over every nonmirror adjacent pair is false
without subtracting this branchwise term.

## 3. Exact connection curvature

Let

\[
 P(x)=34x^3+51x^2+27x+5
\]

and define the shifted Frobenius solution by

\[
 (n+1+z)^3F_{n+1}(z)-P(n+z)F_n(z)
 +(n+z)^3F_{n-1}(z)=0,
 \quad F_{-1}(z)=0,\ F_0(z)=1.
\]

The recurrence and initial conditions agree term for term with the regular
Frobenius deformation used by Golyshev and Zagier.  With

\[
 \lambda_+=17+12\sqrt2=(1+\sqrt2)^4,
 \qquad \lambda_-=\lambda_+^{-1}=17-12\sqrt2,
\]

Put

\[
 \phi_J(z)=\lim_{n\to\infty}\frac{F_n(z)}{F_n(0)},
 \qquad \kappa(z)=\lambda_-^z\phi_J(z).
\]

Their exact telescoping calculation gives the Taylor germ

\[
 \log\kappa(z)=-2\zeta(2)z^2+O(z^3).
\]

Therefore

\[
 \boxed{c:=(\log\kappa)''(0)=-4\zeta(2)=-\frac{2\pi^2}{3}.}
\]

This is not a decimal recognition.  The verifier independently checks the
rational telescoping certificate through \(z^2\), including the gauge
\(\lambda_+^{-z}\).  The normalization and coefficient are also stated for the Apéry
operator in Golyshev--Zagier and in the later global Frobenius-constant
framework of Bloch--Vlasenko.

The cited primary sources are:

- V. Golyshev and D. Zagier, [*Proof of the gamma conjecture for Fano
  3-folds of Picard rank one*, Izvestiya: Mathematics 80 (2016), 24--49,
  DOI 10.1070/IM8343](https://www.mathnet.ru/eng/im8343).
- V. Golyshev and D. Zagier, [*Interpolated Apéry numbers, quasiperiods
  of modular forms, and motivic gamma functions*](https://people.mpim-bonn.mpg.de/zagier/files/tex/DubrovinVolume/GZ_AperyInterpolation.pdf),
  Sections 6--7.
- S. Bloch and M. Vlasenko, [*Gamma functions, monodromy and Frobenius
  constants*](https://math.uchicago.edu/~bloch/paper-revision-final3.pdf),
  especially Propositions 8 and 15, Example 29, Theorem 30, Corollary 33,
  and Remark 32.

## 4. Global limiting cell: unconditional closure

The global bridge is available once the complete two-component cell is used.
Section 7 of Golyshev--Zagier's later interpolation paper defines the same
function by the asymptotic limit below, where \(L_A\) is the explicit
asymptotic factor defined in their Section 6:

\[
 \kappa(z)=\lim_{n\to\infty}
 \frac{F_n(z)}{2^{-9/4}\pi^{-3/2}L_A(n+\tfrac12+z)}
\]

and states that it is meromorphic on \(\mathbf C\) and holomorphic at zero.
Their Section 6 gives

\[
 b_n=F_n(0)\sim2^{-9/4}\pi^{-3/2}L_A(n+\tfrac12),
 \qquad
 \frac{L_A(x+z)}{L_A(x)}\longrightarrow \lambda_+^z.
\]

Thus this published asymptotic normalization is exactly
\(\kappa(z)=\lambda_-^z\phi_J(z)\), not merely a function sharing its first few
Taylor coefficients.  Bloch--Vlasenko Example 29 uses the same Apéry operator
and Frobenius normalization; Theorem 30, Corollary 33, and Remark 32 identify
its direct-positive-path Frobenius germ with a global meromorphic gamma
function.  Equality of the normalized germs and the meromorphic identity
theorem therefore identify the two global functions.  This also fixes the
otherwise possible factor \(\lambda e^{2\pi i m z}\): the chosen direct-path
branch fixes \(m=0\), and \(\kappa(0)=1\) fixes \(\lambda=1\).

This normalization does **not** use Bloch--Vlasenko Lemma 24.  Its boundary
divergence hypothesis fails here because \(b_n\lambda_-^n\asymp n^{-3/2}\).  The
identification instead comes directly from the Golyshev--Zagier asymptotic
definition above and the direct-path normalization in Theorem 30.

To avoid collision with the gcd notation below, denote the companion solution
by \(H_n\), with

\[
 H_0(z)=0,\qquad H_1(z)=\frac1{(z+1)^3},
\]

and put \(\chi_J(z)=\lim H_n(z)/b_n\).  The finite identity

\[
 F_{n-1}(z+1)=(z+1)^3H_n(z)
\]

and \(b_{n-1}/b_n\to\lambda_-\) give

\[
 \psi(z):=\frac{\kappa(z+1)}{(z+1)^3}
          =\lambda_-^z\chi_J(z).                    \tag{4.1}
\]

Here is the finite identity that connects this pair of Jost limits to the
two-block cell.  Set

\[
 N_0(X)=0,\quad N_1(X)=1,\quad
 N_{m+1}(X)=P(X+m)N_m(X)-(X+m)^6N_{m-1}(X),
\]

and \(q_h(X)=\prod_{a=1}^h(X+a)\),
\(\delta_h(X)=N_h(X)/q_h(X)^3\).  If \(1\le j\le h\),
\(r=j-1\), and \(s=h-j\), the continuant addition law is

\[
 N_h(X)=N_j(X)N_{s+1}(X+j-1)
 -(X+j)^6N_{j-1}(X)N_s(X+j).                        \tag{4.2a}
\]

For fixed \(X,j\), both terms on the right obey the defining recurrence in
\(s\); the cases \(s=0,1\) are exactly the two initial values.  Reversing the
tridiagonal continuant and using \(P(-X-1)=-P(X)\) gives, for every \(m\ge1\),

\[
 N_m(-m-1-X)=(-1)^{m-1}N_m(X).                      \tag{4.2b}
\]

After substituting \(X=-j+z\), splitting the denominator at its \(j\)-th
factor, and applying (4.2b), (4.2a) becomes the exact all-height formula

\[
 \delta_h(-j+z)
 =z^{-3}F_r(-z)F_s(z)+z^3H_r(-z)H_s(z).             \tag{4.2c}
\]

Divide by \(b_rb_s\) and let both \(r,s\to\infty\).  Equation (4.1) and
the bilateral cancellation of \(\lambda_-^{\pm z}\) give the normalized
two-block limit

\[
 \boxed{J(z)=z^{-3}\kappa(-z)\kappa(z)
              +z^3\psi(-z)\psi(z).}                 \tag{4.2d}
\]

The second term is essential.  An audit of only
\(z^{-3}\kappa(-z)\kappa(z)\) misses exactly the shift defect canceled by the
companion.

Bloch--Vlasenko Proposition 8 applied to

\[
 L=D^3-tP(D)+t^2(D+1)^3
\]

gives the meromorphic identity below.  The adjoint convention causes no sign
gap: formal adjunction sends \(D\) to \(-D\), and
\(P(-D-1)=-P(D)\), so \(L^\vee=-L\).  Hence the adjoint gamma module in
Theorem 30 has the same scalar equation:

\[
 \kappa(z)-\frac{P(z)}{(z+1)^3}\kappa(z+1)
 +\frac{(z+1)^3}{(z+2)^3}\kappa(z+2)=0.              \tag{4.3}
\]

Equivalently,

\[
 \kappa(z+1)=(z+1)^3\psi(z),\qquad
 (z+1)^3\psi(z+1)=P(z)\psi(z)-\kappa(z).             \tag{4.4}
\]

Since \(P(-z-1)=-P(z)\), equations (4.4) at \(z\) and \(-z-1\) give

\[
 \psi(-z-1)=-\frac{\kappa(-z)}{z^3},\qquad
 \kappa(-z-1)=z^3\psi(-z)-P(z)\psi(-z-1).
\]

Substitution in (4.2d) makes the two \(P(z)\)-terms cancel and proves

\[
 \boxed{J(z+1)=J(z),\qquad J(-z)=-J(z).}              \tag{4.5}
\]

This cancellation is checked symbolically by
`exact_global_cell_cocycle_gate`; the finite shift identity is checked
independently by `exact_shift_and_casoratian_gate`.  The addition,
reflection, and normalized cell formulas (4.2a)--(4.2c) are checked with
generic symbolic induction variables by `exact_universal_continuant_gate`;
small-height reconstruction is retained only as an independent regression.

There are no hidden nonintegral poles.  Corollary 33 writes

\[
 \kappa(z)=\frac{z^3\Gamma(z)}{(1-e^{-2\pi i z})^3}
\]

with \(\Gamma\) entire, so \(\kappa\), and hence \(\psi\), can have poles only
at integers, of order at most three.  Remark 32 gives the finite value
\(\kappa(1)=\zeta(3)/6\), so \(\psi\) is regular at zero.  The exact local
germ from Section 3 now yields

\[
 J(z)=z^{-3}-\frac{2\pi^2}{3}z^{-1}+O(z).             \tag{4.6}
\]

Periodicity transports this complete principal part to every integer and in
particular cancels the apparent sixth-order poles in the separate products.

The required vertical bound also follows from the cited construction rather
than being an extra conjecture.  Proposition 15, on the direct positive-real
path \(0<t<\lambda_-\), applies to the original Apéry cyclic vector without an
extra endpoint multiplier.  Indeed, in ordinary-derivative form the two top
coefficients are

\[
 a_3(t)=t^3(t^2-34t+1),\qquad
 a_2(t)=3t^2(2t^2-51t+1).
\]

At \(\lambda_-=17-12\sqrt2\),

\[
 a_3(\lambda_-)=0,\qquad
 \frac{a_2(\lambda_-)}{a_3'(\lambda_-)}=\frac32,
\]

so the local indicial polynomial is
\(\rho(\rho-1)(\rho-\tfrac12)\).  The one-dimensional monodromy-variation
line is the nonintegral \(\rho=\tfrac12\) branch: its local-monodromy
eigenvalue is \(-1\), whereas the two integral branches have eigenvalue
\(1\).  Hence its pairing is
\(O((\lambda_--t)^{1/2})\), which satisfies Proposition 15's endpoint hypothesis
with \(\alpha=3/2>0\).  Cancelling the factor
\(R(e^{-2\pi iz})=(1-e^{-2\pi iz})^3\) between Proposition 15 and
Corollary 33 gives, on a sufficiently far-right strip,

\[
 \kappa(z)=A z^3\int_0^{\lambda_-} f(t)t^{z-1}\,dt.  \tag{4.7}
\]

Here \(f(t)=\langle m,\delta(t)\rangle\) is the direct-path
monodromy-variation pairing in Proposition 15, and \(A\ne0\) is the scalar
fixed by the normalization \(\kappa(0)=1\).

Its endpoint hypotheses make the integral absolutely convergent there, and
\(|t^{i\operatorname{Im}z}|=1\) on this path.  Hence
\(\kappa(z)=O((1+|\operatorname{Im}z|)^3)\) on a right strip of width two.
Solving (4.3) for \(\kappa(z)\) and shifting left finitely many times propagates
the same polynomial bound to any fixed vertical strip; the rational
coefficients are uniformly bounded on its vertical tails.  Equation (4.1)
then gives \(\psi(z)=O(1)\), and therefore

\[
 J(z)=O((1+|\operatorname{Im}z|)^3)                  \tag{4.8}
\]

on a fundamental strip away from its pole.

Subtract the periodic polar model

\[
 E(z)=J(z)-\pi^3\cot(\pi z)\csc^2(\pi z)
          +\frac{2\pi^3}{3}\cot(\pi z).
\]

Equations (4.5)--(4.8) make \(E\) entire, odd, one-periodic, and polynomially
bounded in both vertical directions.  Shifting the Fourier integral downward
for positive frequencies and upward for negative frequencies kills every
nonzero Fourier coefficient; oddness kills the constant coefficient.  Thus
\(E=0\), and the global identity is unconditional:

\[
 \boxed{
 J(z)=\pi^3\cot^3(\pi z)+\frac{\pi^3}{3}\cot(\pi z)
 =\frac{\pi^3}{3}
   \frac{\cos(\pi z)(\cos(2\pi z)+2)}{\sin^3(\pi z)}.} \tag{4.9}
\]

Its nonreal zeros are

\[
 z=\frac12\pm i\eta\pmod{\mathbf Z},\qquad
 \boxed{\cosh(2\pi\eta)=2},
 \qquad
 \eta=\frac{\log(2+\sqrt3)}{2\pi}.                  \tag{4.10}
\]

The nonreal critical points of the global function \(J\) lie at the distinct
height

\[
 a=\frac{\log2}{2\pi}.                               \tag{4.11}
\]

The verifier checks the recurrence signs, complete two-term cocycle, gauge
cancellation, adjoint reordering, reflection-point indicial exponents,
Laurent coefficient, factorization, derivative, zero equation, and both
distinct heights.  The analytic inputs in this section are imported
with their exact normalizations from the cited primary sources; no finite
sampling is used to establish continuation or growth.

## 5. An independent exact hypergeometric anomaly

Define

\[
 B(x)={}_4F_3\!\left(\begin{matrix}-x,-x,x+1,x+1\\1,1,1\end{matrix};1\right)
 =\sum_{k\ge0}\frac{(-x)_k^2(x+1)_k^2}{(k!)^4}
\]

and

\[
 (\mathcal Ly)(x)=(x+1)^3y(x+1)-P(x)y(x)+x^3y(x-1).
\]

Termwise Gosper telescoping proves

\[
 \boxed{(\mathcal LB)(x)
 =8(2x+1)\left(\frac{\sin\pi x}{\pi}\right)^2.}
\]

This identity explains why the interpolation is homogeneous on the integer
lattice and inhomogeneous off it.  In particular, \(B\) is not the normalized
homogeneous Jost connection: \(B''(0)=\pi^2/3\), which has the wrong sign and
magnitude for \(c\).  The explicit rational certificate is checked identically
by `exact_hypergeometric_anomaly_gate`; three complex evaluations at 60-digit
precision are only an independent regression, not the proof.

## 6. Positive Apéry remainder

Let \(q_0=0,q_1=8\) and

\[
 (n+1)^3q_{n+1}-P(n)q_n+n^3q_{n-1}=8(2n+1).
\]

For integers \(m\ge1\), put

\[
 W_n(m)=\frac1{m^2}\prod_{k=0}^{n}
 \left(\frac{m-k}{m+k}\right)^2.
\]

The \(k=0\) factor is one.  A second exact Gosper identity gives

\[
 \boxed{b_n\zeta(2)-q_n=\sum_{m\ge1}W_n(m).}
\]

For \(n\ge1\), one has \(W_n(m)=0\) for \(m\le n\), while
\(0<W_n(m)<m^{-2}\) for \(m>n\).  Hence

\[
0<b_n\zeta(2)-q_n
 <\sum_{m>n}\frac1{m^2}<\frac1n,
\]

and therefore \(q_n/b_n\to\zeta(2)\) from below.  Monotonicity also follows
directly: the displayed product gives \(W_{n+1}(m)\le W_n(m)\) termwise,
while the binomial formula gives \(b_{n+1}>b_n\), so
\((b_n\zeta(2)-q_n)/b_n\) decreases.  The proof is
not an interchange of divergent series: sum the rational identity to a finite
cutoff, use its explicit boundary limit, and only then let the cutoff tend to
infinity.  The verifier checks the cleared rational identity, boundary sign,
initial values, and recurrence.

This analytic result is independent of the pointwise gcd theorem.  It is best
presented in a companion note rather than inserted into the arithmetic proof
chain.

## 7. All-height critical separation remains open

Arb interval arithmetic certifies the finite range \(2\le h\le60\).  The
global formula makes the limiting central critical height exact:

\[
 a=\frac{\log2}{2\pi}.
\]

This limiting constant does not prove the tail.  A complete theorem still
requires:

- a uniform \(C^2\) Jost/Riccati remainder with one discrete height
  difference;
- a rigorous inventory of every analytic critical branch;
- NEAR, FAR, EDGE, and branch-switch separation with explicit constants;
- comparison of the resulting cutoff with the finite Arb range.

There is also a noncompact obstruction that the global formula makes visible.
In the cylinder coordinate \(q=e^{2\pi iz}\), the upper end
\(\operatorname{Im}z\to+\infty\) has \(\cot(\pi z)\to-i\), and (4.9) gives

\[
 J(z)\longrightarrow \frac{2i\pi^3}{3},\qquad J'(z)\longrightarrow0.
\]

Thus convergence on every fixed compact set, even with a sharp central-cell
error, supplies no uniform derivative gap at the cylinder end and cannot by
itself inventory all finite-height critical branches.  The missing analytic
input is a moving cylinder-chart expansion, together with a uniform
transversality/noncollision estimate for those branches.

No finite scan, fitted remainder, or formal (1/h)-series can replace those
uniform gates.

The finite Arb certificate is produced separately by `CRON_kinf_branch.py` and
recorded in `CRON_kinf_results.json`.  `RC_BREAKTHROUGH_verify.py` audits every
archived certificate row but does not pretend that a manifest check replaces
the interval computation.

## 8. Relation to the main gcd problem

The corrected every-\(n\) master sum is

\[
 M(n)=\sum_{\sqrt n<p\le n}\log p\,
       \mathbf1_{p\mid b_{n\bmod p}},
\]

with

\[
 \log G_n\le6M(n)+O(n^{2/3}),\qquad
 M(n)\le\log G_n+O(1).
\]

The exact global cell calculation and the fixed-prime \(S_D\) problem do not control
the diagonal choice \(n\bmod p\) simultaneously across different primes.
The strongest proved pointwise-in-\(X\) substitute remains

\[
 \#\{n\le X:\log G_n>\varepsilon n\}
 =O_\varepsilon((\log X)^2).
\]

Thus neither the exact global connection formula nor the finite-height certificates
change the terminal verdict for Problem 3.2.

## 9. Claim status

**PROVED:** the exact connection curvature, the complete global limiting-cell
formula (4.9), the zero and critical heights (4.10)--(4.11), the
hypergeometric anomaly, and the positive Apéry remainder.

**CERTIFIED:** critical-value separation for each integer height
\(2\le h\le60\), by the archived Arb interval computation.

**OPEN:** the every-\(n\) gcd theorem, the fixed-prime estimate \(S_D=o(p)\),
and critical-value separation for all heights.

## 10. Reproducibility

Run

```bash
python3 RC_BREAKTHROUGH_verify.py
python3 CODEX_MAINTHM_verify.py
python3 -u CRON_kinf_branch.py --min-h 2 --max-h 60 --digits 120 \
  --cert-digits 110 --skip-kinf --report /tmp/RC_KINF_recheck_report.md \
  --json /tmp/RC_KINF_recheck_results.json
```

The first two programs terminate with `PASS`.  The first checks the analytic finite
algebra in this report, including the complete cell cocycle; the second checks
the corrected Wronskian, master sum,
determinant reduction, saturation examples, exceptional-set algebra, digit
identities, and exponent budgets.  The third command independently recomputes
the finite Arb certificate at every integer height \(2\le h\le60\), writing
only to temporary files.  Diagnostics printed by these programs are not used
as substitutes for asymptotic proofs.
