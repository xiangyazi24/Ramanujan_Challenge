# Q7699: p-adic horizontal Mellin attack

## Verdict

The defining-characteristic Mellin transform does admit an exact p-adic/Gross--Koblitz analysis, but it does **not** produce a valuation-only horizontal sieve.

The strongest exact statement I obtain is this.  For an interior top-half prime
\(n/2<p<n\), put
\[
r=n-p,\qquad N=p-1,\qquad
\rho=\min(r,N-r)=\min(n-p,2p-n-1).
\]
Q7690 gives
\[
M_p(n-1)\equiv-b_r\pmod{\mathfrak p_p}.
\]
The Gross--Koblitz/Morita expansion gives
\[
b_r=\sum_{k=0}^r p^{2\mathbf1_{r+k\ge p}}U_{p,r,k}^2,
\qquad
U_{p,r,k}=\frac{\Gamma_p(r+k+1)}{\Gamma_p(k+1)^2\Gamma_p(r-k+1)}\in\mathbf Z_p^\times.
\]
Hence exactly \(\rho+1\) terms have minimum valuation zero, every remaining term has valuation two, and
\[
\boxed{
v_p(M_p(n-1))>0
\iff
\Theta_{p,r}:=\sum_{k=0}^{\rho}\overline U_{p,r,k}^{\,2}=0\quad\text{in }\mathbf F_p.}
\tag{A}
\]
For every interior \(r\), \(\rho+1\ge2\).  Thus there is **never a unique minimum-valuation term**.  Cancellation among equal-slope units is not a residual nuisance; it is the entire defining-characteristic zero event.

Moreover the minimum-slope sum is not a new small-complexity object.  If \(r\le N/2\), it is just \(b_r\bmod p\).  If \(r>N/2\), the elementary reflection of the no-carry terms gives
\[
\Theta_{p,r}\equiv b_{N-r}=b_\rho\pmod p.
\]
Thus in all cases
\[
\boxed{
v_p(M_p(n-1))>0\iff p\mid b_{\rho_n(p)},\qquad
\rho_n(p)=\min(n-p,2p-n-1).}
\tag{B}
\]
This reduces the moving coefficient to the first half of the mod-\(p\) Apéry row, but \(\rho_n(p)\) is still linear in \(n\) (as large as about \(n/3\)).  It is an exact reflection/reindexing, not a power-saving horizontal theorem.

The literal finite-field/Jacobi transform has **exactly \(p-1\) character summands**.  Gross--Koblitz leaves \(\rho+1\) minimum-slope unit contributions.  There is no \(O(1)\) or \(O(\gcd(m,p-1))\) decomposition in this representation.  The Kummer twist order \((p-1)/\gcd(m,p-1)\) is not the number of Mellin/Jacobi terms.

A targeted Sage run at fixed \(n=321\) confirms the obstruction at genuine horizontal zeros.  The 29 top-half primes are all interior; the three bad primes are
\[
(p,r)=(179,142),(193,128),(211,110).
\]
They have respectively 37, 65, and 101 **valuation-zero** unit terms whose residues cancel to zero.  Nearby nonzeros have the same type of valuation profile.  For every prime \(p\ge223\) in this sample all \(r+1\) terms have valuation zero, so the Newton/carry polygon has literally no separation to exploit.

The p-adic route therefore fails at the slope/carry stage.  The exact remaining theorem is a **horizontal unit-cancellation sieve**, not a Newton-polygon theorem.

Artifacts:

- `problems/3.2/research/scripts/q7699_padic_horizontal_mellin.sage`
- `problems/3.2/research/scripts/q7699_padic_horizontal_mellin_n321.csv`

The Sage 10.6 container run ended with `Q7699_SAGE_VERIFY PASS`.

---

## 1. Normalization and endpoints from Q7690

Let
\[
\Lambda(x,y,z)=\frac{(1+x)(1+y)(1+z)}{xyz}
\bigl((1+y)(1+z)+xyz\bigr)
\]
be the fixed Laurent polynomial used by the toric/K3 realization, and let \(\theta_p(a)\) be the primitive rank-three K3 trace for the fiber \(\Lambda=a\), \(a\in\mathbf F_p^\times\).

Fix the Teichmuller lift
\[
\omega_p:\mathbf F_p^\times\longrightarrow\mu_{p-1}\subset\mathbf Z_p^\times
\]
with \(\omega_p(a)\equiv a\pmod p\), and define the **unnormalized** Mellin transform
\[
M_p(m)=\sum_{a\in\mathbf F_p^\times}\theta_p(a)\omega_p(a)^m.
\tag{1}
\]
There is no factor \(1/(p-1)\) in (1).  The exact Q7690 theorem, with this choice of Teichmuller generator and exponent direction, is
\[
\boxed{M_p(m)\equiv-b_r\pmod{\mathfrak p_p}}
\tag{2}
\]
when
\[
1\le r\le p-2,\qquad m\equiv r\pmod{p-1}.
\]
The minus sign and exponent direction are part of the normalization.  Replacing the character generator by another generator reindexes exponents by a unit modulo \(p-1\); one may not silently keep the same integer \(r\).

Most importantly, \(r\) is a **Mellin exponent / Apéry coefficient index**.  It is not the geometric fiber parameter \(a\).  The transform sums over every nonzero fiber \(a\), and only after multiplicative Fourier projection does the coefficient \(b_r\) appear.

### Fixed \(n\), top-half primes

Set \(m=n-1\).  If
\[
\frac n2<p<n
\]
and neither endpoint degeneration below occurs, then
\[
n-1=(p-1)+(n-p),
\]
so
\[
r=n-p\in[1,p-2],\qquad m\equiv r\pmod{p-1}.
\tag{3}
\]
Therefore Q7690 gives
\[
M_p(n-1)\equiv-b_{n-p}\pmod{\mathfrak p_p}.
\tag{4}
\]
Because \(\theta_p(a)\in\mathbf Z\) and each Teichmuller value is a \(p\)-adic unit, \(M_p(n-1)\) is \(p\)-integral.  Hence
\[
\boxed{v_p(M_p(n-1))>0\iff p\mid b_{n-p}.}
\tag{5}
\]

There are two boundary cases outside the interior theorem.

1. **Upper endpoint \(p=n\)** (possible only when \(n\) itself is prime): then \(n-p=0\), and \(b_0=1\), so this prime is never bad.
2. **Lower Mellin endpoint \(n=2p-1\)**: then \(n-p=p-1\) while \(m=n-1\equiv0\pmod{p-1}\).  Q7690's interior statement does not identify this with the \(r=p-1\) coefficient.  Directly,
   \[
   b_{p-1}\equiv1\pmod p,
   \]
   because the \(k=0\) Apéry summand is 1 and for every \(k\ge1\), \(\binom{p-1+k}{k}\) contains a factor \(p\).  Thus this endpoint is also never bad.

No endpoint is smuggled into (2).

---

## 2. Exact finite-field/Jacobi complexity: \(p-1\), not \(O(1)\)

The earlier finite-field hypergeometric reduction gives a second exact expression for \(b_r\bmod p\).  Let
\[
N=p-1,
\]
choose a generator \(T\) of the multiplicative character group compatibly with the Teichmuller prime, put
\[
\varphi=T^{N/2},\qquad \chi_r=T^r,
\]
and let \(\eta\) run over **all \(N\) multiplicative characters**.  Characters are extended by 0 at 0, and
\[
J(A,B)=\sum_{x\in\mathbf F_p}A(x)B(1-x).
\]
Then, for \(1\le r\le p-2\),
\[
\boxed{
b_r\equiv-rac1N\sum_{\eta}
J(\varphi\eta\chi_r,\varphi\eta^{-1})^2
J(\varphi\eta\chi_r^{-1},\varphi\eta^{-1})^2
\pmod{\mathfrak p_p}.}
\tag{6}
\]
The factor in (6) is \(-1/N\); it should not be confused with the unnormalized geometric Mellin sum (1).

Writing \(\eta=T^t\), \(0\le t<N\), makes the complexity literal:

\[
\boxed{\text{(6) has exactly }N=p-1\text{ character summands}.}
\tag{7}
\]

The dependence on \(m\) or \(r\) only shifts the two numerator characters.  It does not reduce the \(t\)-sum to the \(m\)-th roots of unity.  In particular, \(\gcd(m,p-1)\) controls the **order** of the Kummer character
\[
\operatorname{ord}(\omega_p^m)=\frac{p-1}{\gcd(m,p-1)},
\]
but it is not the number of summands in (6).

This distinction is visible numerically at \(n=321\).  For \(p=179\),
\[
\gcd(320,178)=2,
\]
yet the literal Jacobi formula has 178 terms and the minimum-slope packet below has 37 terms.  So an \(O(\gcd(m,p-1))\) term count is false for this representation.

The direct truncated hypergeometric representation
\[
\boxed{
b_r=\sum_{k=0}^r
\left(\binom rk\binom{r+k}{k}\right)^2}
\tag{8}
\]
has exactly \(r+1\) terms.  Thus the two exact explicit complexities are

- \(p-1\) terms in the all-character Jacobi/Gauss transform;
- \(r+1\) terms in the coefficient/Gross--Koblitz-Morita truncation.

Both are \(O(p)\) in the top-half regime in the worst case.  Neither is uniformly \(O(1)\) or \(O(\gcd(m,p-1))\).

---

## 3. Stickelberger/Gross--Koblitz valuation of the Jacobi terms

Choose the Gross--Koblitz exponent convention \(T=\bar\omega_p=\omega_p^{-1}\).  For \(0\le j\le N-1\),
\[
g(T^j)=-\pi^j\Gamma_p(j/N),
\qquad v_p(\pi)=1/N.
\]
Consequently, whenever \(T^{a+b}\) is nontrivial,
\[
v_p J(T^a,T^b)
=\frac{[a]_N+[b]_N-[a+b]_N}{N},
\tag{9}
\]
which is the carry indicator 0 or 1.

For the \(t\)-th term in (6), put
\[
a_+(t)=[N/2+t+r]_N,
\quad
a_-(t)=[N/2+t-r]_N,
\quad
c(t)=[N/2-t]_N.
\]
The two products of Jacobi characters are \(\chi_r\) and \(\chi_r^{-1}\), both nontrivial for \(1\le r\le p-2\).  Hence there is no trivial-product exception in (9).  Define
\[
\epsilon_+(t)=\frac{a_+(t)+c(t)-r}{N}\in\{0,1\},
\]
\[
\epsilon_-(t)=\frac{a_-(t)+c(t)-(N-r)}{N}\in\{0,1\}.
\]
Then the whole squared Jacobi product has exact valuation
\[
\boxed{v_p(\mathcal J_t(r))=2\bigl(\epsilon_+(t)+\epsilon_-(t)\bigr)\in\{0,2,4\}.}
\tag{10}
\]

Let
\[
\rho=\min(r,N-r).
\]
A direct carry count gives
\[
\boxed{\#\{t:v_p(\mathcal J_t(r))=0\}=\rho+1.}
\tag{11}
\]
For example, when \(r\le N/2\), the minimum terms are exactly
\[
t=N/2-r,\,N/2-r+1,\ldots,N/2.
\]
The case \(r>N/2\) follows by the symmetry \(r\leftrightarrow N-r\), which interchanges the two Jacobi factors.

Thus the all-character GK calculation already has an unavoidable equal-slope packet of size \(\rho+1\).  The prefactor \(-1/N\) is a \(p\)-adic unit and does not affect valuations.

---

## 4. Exact Morita form and the minimum-unit cancellation criterion

For \(0\le r\le p-1\) and \(0\le k\le r\), define Morita's integer gamma by
\[
\Gamma_p(a)=(-1)^a\prod_{\substack{1\le j<a\\p\nmid j}}j
\]
and set
\[
U_{p,r,k}
=\frac{\Gamma_p(r+k+1)}
{\Gamma_p(k+1)^2\Gamma_p(r-k+1)}\in\mathbf Z_p^\times.
\]
Since \(r+k\le2p-2\), there is at most one factorial carry, and the exact identity is
\[
\boxed{
\binom rk\binom{r+k}{k}
=p^{\mathbf1_{r+k\ge p}}U_{p,r,k}.}
\tag{12}
\]
All signs cancel.  Squaring and summing yields
\[
\boxed{
b_r=\sum_{k=0}^r
p^{2\mathbf1_{r+k\ge p}}U_{p,r,k}^2.}
\tag{13}
\]
Therefore
\[
\underbrace{0,\ldots,0}_{\rho+1\text{ terms}},
\quad
\underbrace{2,\ldots,2}_{\max(0,2r+1-p)\text{ terms}}
\tag{14}
\]
is the complete termwise valuation profile of (13), and
\[
\boxed{
p\mid b_r
\iff
\Theta_{p,r}:=\sum_{k=0}^{\rho}\overline U_{p,r,k}^{\,2}=0
\quad\text{in }\mathbf F_p.}
\tag{15}
\]
Combining with Q7690,
\[
\boxed{
v_p(M_p(n-1))>0\iff\Theta_{p,n-p}=0.}
\tag{16}
\]

### Cancellation is the only remaining issue

For every interior Mellin index
\[
1\le r\le p-2,
\]
one has
\[
\rho+1\ge2.
\]
So Stickelberger/Gross--Koblitz **never produces a unique lowest term** on the interior range.  No interior top-half prime is excluded by the Newton polygon alone.

For fixed \(n\), the counts become especially transparent.  With \(r=n-p\),
\[
u_p(n)=\rho+1
=\min(n-p+1,\,2p-n),
\tag{17}
\]
\[
h_p(n)=\max(0,2n-3p+1).
\tag{18}
\]
Thus:

- if \(p\ge(2n+1)/3\), then \(h_p(n)=0\): **all \(r+1=n-p+1\) terms are p-adic units**;
- if \(n/2<p<(2n+1)/3\), then exactly \(2p-n\) unit terms survive modulo \(p\), and the rest have valuation two;
- at the excluded lower endpoint \(n=2p-1\), the unit packet has size one and gives \(b_{p-1}\equiv1\), so the only unique-minimum case is harmless.

Hence the carry polygon cannot even exclude a positive-density subinterval of the interior top-half primes: on the upper third it has no slope separation at all, and on the lower part it still has at least two equal-minimum units.

---

## 5. The minimum-unit packet is just the reflected Apéry coefficient

There is one exact simplification, but it is a reindexing rather than a horizontal bound.

If \(r\le N/2\), then \(\rho=r\) and (15) is simply the original coefficient modulo \(p\).

If \(r>N/2\), write
\[
r=N-\rho=p-1-\rho.
\]
For \(0\le k\le\rho\), reduction modulo \(p\) gives
\[
\binom{p-1-\rho}{k}
\equiv(-1)^k\binom{\rho+k}{k},
\]
\[
\binom{p-1-\rho+k}{k}
\equiv(-1)^k\binom{\rho}{k}.
\]
After multiplying and squaring,
\[
\overline U_{p,r,k}^{\,2}
\equiv
\left(\binom\rho k\binom{\rho+k}{k}\right)^2.
\]
Consequently
\[
\boxed{\Theta_{p,r}\equiv b_\rho\pmod p.}
\tag{19}
\]
This is the basic reflection
\[
b_r\equiv b_{p-1-r}\pmod p.
\tag{20}
\]
No stronger reflected-depth law is being asserted.

For the fixed top-half problem,
\[
\boxed{
\rho_n(p)=\min(n-p,2p-n-1),
\qquad
v_p(M_p(n-1))>0\iff p\mid b_{\rho_n(p)}.}
\tag{21}
\]
The tent index \(\rho_n(p)\) is at most about \(n/3\), but still has linear scale.  Equation (21) therefore does not supply a power-saving prime count by itself.

It does explain the Sage profiles below: the number of slope-zero units is exactly \(\rho_n(p)+1\).

---

## 6. Explicit elliptic--Asai cover: exact mod-p reduction, still O(p)

The exact toric-fiber descent gives
\[
D_a=a^2-34a+1
\]
and the quadratic cover
\[
8x^2+(a-1)x+a=0,
\tag{22}
\]
with elliptic curve
\[
E_x:\quad Y^2+(1-2x)XY+x^2Y=X^3.
\tag{23}
\]
For the primitive trace \(\theta_p(a)\):

- if \(D_a\) is a nonzero square,
  \[
  \theta_p(a)=A_p(E_x)^2-p-\mathbf1_{a=1}p;
  \]
- if \(D_a\) is a nonsquare,
  \[
  \theta_p(a)=\chi_p(-3)A_{p^2}(E_x)-p-\mathbf1_{a=1}p;
  \]
- if \(D_a=0\),
  \[
  \theta_p(a)=A_p(E_x)^2-p-\chi_p(-6)p.
  \]

Therefore, modulo the Teichmuller prime, every displayed Tate correction vanishes and
\[
\theta_p(a)\equiv
\begin{cases}
A_p(E_x)^2,&D_a\text{ square or }0,\\
\chi_p(-3)A_{p^2}(E_x),&D_a\text{ nonsquare}
\end{cases}
\pmod p.
\tag{24}
\]
Substitution in (1) is an exact p-adic test for positivity of \(v_p(M_p(m))\), but it leaves a sum over the \(p-1\) nonzero parameters \(a\).  The degree-two cover (22) replaces each K3 trace by one elliptic/Asai trace; it does not diagonalize the multiplicative parameter.  Its Mellin complexity is still \(O(p)\).

The Jacobi formula (6) is precisely the more useful Fourier-diagonalized form, and even there the literal term count is \(p-1\).

Thus the elliptic--Asai cover gives an exact geometric interpretation of each summand but no \(O(1)\) or \(O(\gcd(m,p-1))\) p-adic decomposition.

---

## 7. Dwork/unit-root interpretation: equivalent obstruction

The rank-three middle-extension hypergeometric/K3 stalk has the usual first Hasse invariant equal, in this normalization, to the Apéry coefficient modulo \(p\).  Equivalently:

- \(p\nmid b_r\): a slope-zero/unit-root direction is present;
- \(p\mid b_r\): the first Hasse invariant vanishes and the slope-zero direction disappears.

This is a useful conceptual restatement of (15), but it does not strengthen it horizontally.  The Dwork unit-root criterion asks whether the same unit residue \(\Theta_{p,r}\) vanishes.  Since the Kummer character order varies with \(p\), there is no fixed p-adic family in which the Hasse invariant becomes a bounded-degree polynomial in one fixed parameter.

So the Dwork route reaches the same exact wall: control the horizontal zero set of a moving Hasse/unit residue, not merely its Newton polygon.

---

## 8. Targeted Sage experiment at fixed n = 321

The experiment was run in Sage 10.6 using

```bash
sage -python problems/3.2/research/scripts/q7699_padic_horizontal_mellin.sage --n 321
```

The script does **not** merely count Apéry zeros.  For every top-half prime it independently checks:

1. the exact integer hypergeometric summands;
2. the valuation vector \(0^{u_p}2^{h_p}\);
3. the Morita-gamma unit for every minimum term;
4. equality between the normalized gamma-unit square and the direct binomial residue;
5. \(\Theta_{p,r}=b_r\bmod p\);
6. the zero event.

It records the full per-prime profile in

`q7699_padic_horizontal_mellin_n321.csv`.

The run summary was

```text
Q7699_FIXED_N 321
top_half_primes 29
interior 29
bad_count 3
bad_pairs [(179, 142), (193, 128), (211, 110)]
complexity_min 5
complexity_max 159
min_slope_terms_min 5
min_slope_terms_max 101
Q7699_SAGE_VERIFY PASS
```

The three genuine cancellations are:

| p | r=n-p | rho | valuation-0 terms | valuation-2 terms | Theta |
|---:|---:|---:|---:|---:|---:|
| 179 | 142 | 36 | 37 | 106 | 0 |
| 193 | 128 | 64 | 65 | 64 | 0 |
| 211 | 110 | 100 | 101 | 10 | 0 |

The first and last normalized minimum residues are nonzero in every case, e.g.

- \(p=179\): first 1, last 29; 37 units cancel;
- \(p=193\): first 1, last 143; 65 units cancel;
- \(p=211\): first 1, last 122; 101 units cancel.

Nearby controls show no valuation distinction:

- \(p=181,r=140\): 41 minimum units, 100 valuation-two terms, \(\Theta=137\ne0\);
- \(p=191,r=130\): 61 minimum units, 70 valuation-two terms, \(\Theta=24\ne0\);
- \(p=197,r=124\): 73 minimum units, 52 valuation-two terms, \(\Theta=14\ne0\).

For \(p=223,227,\ldots,317\) in this sample the high-count column is zero: every hypergeometric summand is a p-adic unit.  The corresponding \(\Theta\) values are nonzero in this finite sample, but **the valuation profile gives no reason for that fact**.

This is the exact obstruction requested by the task: not just many terms, but many **equal minimum-valuation units** whose residual sum is uncontrolled by Stickelberger/Gross--Koblitz.

---

## 9. Can the carry/digit condition exclude a positive-density subset?

Not by valuation alone.

For every interior top-half prime,
\[
u_p(n)=\min(n-p+1,2p-n)\ge2.
\]
Therefore the unique-minimum mechanism excludes **no** interior prime.

On the entire subinterval
\[
p\ge\frac{2n+1}{3},
\]
one has \(h_p(n)=0\), so all terms lie on the same valuation level.  This is a positive-proportion region of the top-half interval where the carry polygon has no separating power whatsoever.

On the lower subinterval, the only information is that the high terms disappear modulo \(p\), leaving \(2p-n\) units.  Their sum is exactly the reflected coefficient \(b_{2p-n-1}\bmod p\).  This is still a moving-index horizontal zero problem.

The first nontrivial packet size illustrates the limitation.  If \(2p-n=2\), then \(\rho=1\) and the residue is \(b_1=5\), so all primes \(p>5\) in this one lattice slice are excluded.  More generally a bounded value \(2p-n=d\) reduces the zero test to the fixed integer \(b_{d-1}\).  But for fixed \(d\) there is at most one candidate prime \(p=(n+d)/2\); allowing \(d=O(1)\) removes only \(O(1)\) candidates, not a positive-density subset.  For a positive proportion of primes, \(d\asymp n\) and the cancellation problem returns at linear scale.

So the exact carry law yields only a thin endpoint exclusion, not a power-saving horizontal sieve.

---

## 10. Precise sufficient p-adic horizontal sieve theorem

Define the interior bad set
\[
\mathcal H(n)=\left\{
 p\text{ prime}:\frac n2<p<n,\ 1\le n-p\le p-2,\
 \Theta_{p,n-p}=0
\right\}.
\]
By Q7690 and (15), this is exactly
\[
\mathcal H(n)=\left\{
 p:\frac n2<p<n,\ v_p(M_p(n-1))>0
\right\}
\]
after removing the two harmless endpoints.

The **minimal weighted horizontal theorem sufficient for the top-half defining-characteristic factor** is
\[
\boxed{
\sum_{p\in\mathcal H(n)}\log p=o(n).}
\tag{PHS}
\]
Equivalently, using (21),
\[
\sum_{\substack{n/2<p<n\\p\mid b_{\rho_n(p)}}}\log p=o(n).
\]
A convenient stronger power-saving form would be: there exists \(\delta>0\) such that
\[
\boxed{
\#\mathcal H(n)\ll \frac{n^{1-\delta}}{\log n}.}
\tag{PHS}_\delta
\]
Then
\[
\sum_{p\in\mathcal H(n)}\log p\ll n^{1-\delta}=o(n).
\]

This is genuinely weaker than a theorem controlling the full Apéry radical: it only addresses

- primes in the top-half interval;
- the single moving Mellin exponent \(m=n-1\);
- support of the event \(v_p(M)>0\), not higher multiplicity;
- the defining-characteristic slice isolated by Q7690.

If the preceding reduction has already isolated this slice as the remaining radical contribution, (PHS) is exactly sufficient for that missing step.  It makes no claim about other prime ranges or the full factorization of any \(b_j\).

Crucially, (PHS) is **not** a Newton-polygon statement.  By (15), it is a horizontal anti-cancellation theorem for the unit traces \(\Theta_{p,n-p}\).

---

## 11. Exact obstruction and surviving route

The p-adic attack succeeds up to the following exact normal form:

\[
M_p(n-1)\bmod p
\quad\longleftrightarrow\quad
-b_{n-p}\bmod p
\quad\longleftrightarrow\quad
-\Theta_{p,n-p}
\quad\longleftrightarrow\quad
-b_{\rho_n(p)}\bmod p.
\]

The literal Jacobi transform has \(p-1\) terms.  Gross--Koblitz/Stickelberger removes all positive-slope terms but leaves \(\rho_n(p)+1\) equal-slope units, which can be \(\asymp p\).  At genuine horizontal zeros these units really do cancel, as the \(n=321\) Sage data demonstrate.

Therefore:

- **no valuation/digit/Newton-polygon condition obtained here gives a power-saving count;**
- **cancellation among minimum-valuation terms is the only remaining p-adic issue;**
- the surviving p-adic theorem must control the horizontal distribution of the unit residue \(\Theta_{p,n-p}\) itself (or an equivalent Dwork Hasse invariant / Jacobi-unit trace), not just its slopes;
- the elliptic--Asai cover and bounded geometric conductor do not reduce the character sum below \(O(p)\) because the Kummer twist varies with \(p\).

A successful next theorem would therefore have to be genuinely horizontal, e.g. a cross-characteristic equidistribution/large-sieve statement for these **unit residues**.  Gross--Koblitz has already extracted all the valuation information available from the termwise carry polygon.