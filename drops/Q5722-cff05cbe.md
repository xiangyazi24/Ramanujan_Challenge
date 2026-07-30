ANSWER Q5722 cff05cbe

# P3.2 route-selection audit: what a first-cell recurrence would and would not buy

## Verdict

A fixed-order P-recursive equation for

\[
F_M(r)=C_M(M-r)-b_M
\]

would be useful, but it would **not by itself** control the common gcd of the pair-preserving margin family.  The precise result is:

1. After the known core-prime Smith factors are retained and all outside formal Smith content is saturated away, the margin family becomes an observation matrix on the finite recurrence state.
2. If the recurrence has order \(R\), then:
   - \(S+T=R-1\) is the dimensionally sharp minimum for a **projective** exceptional integer.  This integer necessarily involves the actual initial state and can have linear height or vanish identically.
   - \(S+T=R\) is the sharp minimum for a clean **primitive-state** theorem: outside primes divide a recurrence singularity product, an augmented observability determinant, or the content of the augmented state \((b_M,F_M(r),\ldots,F_M(r+R-1))\).
3. Consequently a separate small Bézout identity

   \[
   U_{-1}(M,r)b_M+\sum_{i=0}^{R-1}U_i(M,r)F_M(r+i)=J_M(r)
   \]

   with \(J_M(r)\ne0\) and \(\log|J_M(r)|=\operatorname{polylog}M\) is load-bearing.  The recurrence alone does not imply it.
4. For one short core with \(L=O(\log M)\), order, degree, and coefficient height all polylogarithmic, the resulting determinant has polylogarithmic logarithmic height, hence \(o(M^{1/3})\).  But recurrence theory alone does **not** make that determinant reusable across an entire \(M^{1/3}\)-block.  Multiplying the translated determinants, or using one stencil of length \(M^{1/3}\), restores at least linear block height.
5. In the exact first-cell decomposition \(F_M=K_M+L_M+H_M\), the hard squared-binomial part of \(K_M\) already has the explicit degree-\(M+1\) Green/continuant denominator from Q5716.  A fixed-order recurrence can arise only through near-total cancellation with the two boundary pieces.  Even if that cancellation is proved, the augmented primitive-state theorem remains a separate arithmetic obligation.

Thus the route-selection answer is:

\[
\boxed{\text{A full recurrence is not a closure theorem.  A uniform primitive-state/Bézout theorem is load-bearing.}}
\]

---

## 1. Live repository audit

I used the current project state and the exact scripts/artifacts relevant to this question:

- [`Q32_SEPARATION_ANALYSIS.md`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/main/problems/3.2/research/working_notes/Q32_SEPARATION_ANALYSIS.md), current fetched `main` blob `76f820a89dd64deccfc53f78bb069baccb132a28`;
- [`q5711_first_cell_audit.py`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/chatgpt-drop/problems/3.2/research/scripts/q5711_first_cell_audit.py), blob `13ff6cf3fccf9c7d47239aa7ff8fbb54f68ba98e`;
- [`q5714_margin_audit.py`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/q5714-audit-d816e1d2/problems/3.2/research/scripts/q5714_margin_audit.py), blob `5dc5da4f3d5b07979222023ad4a7598d23726dde`;
- [`q5716_cd_audit.py`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/chatgpt-drop/problems/3.2/research/scripts/q5716_cd_audit.py), blob `b8f9ac02b7f9bde483f339f88e65cb36a59ea973`;
- the completed Q5710/Q5715/Q5716/Q5717 analyses.

The checked-in `main` note still ends at Section 52.  The exact decomposition referred to as Section 57 is nevertheless source-verifiable in `q5711_first_cell_audit.py`: with \(r=M-d\), it is precisely

\[
C_M(M-r)=b_M+K_M(r)+L_M(r)+H_M(r).
\]

The formulas are reproduced in Section 9 below.  I do not rely on an unobserved recurrence computation or on unexecuted numerical output.

---

## 2. Exact margin calculus and the correct saturation

Write

\[
Y_d=C_M(d),\qquad
G_{d,L}=\sum_{i=0}^L(-1)^i
 \binom{d+i}{i}\binom{d+L+1}{L-i}Y_{d+i}.
\]

For a fixed core \([d,d+L]\), put

\[
H_{s,t}=G_{d-s,L+s+t},\qquad 0\le s\le S,\quad0\le t\le T,
\]

and set

\[
m=S+T,\qquad \ell=d+L+1.
\]

The exact Pascal identities checked in `q5714_margin_audit.py` are, for \(n=L+s+t\),

\[
\boxed{
H_{s,t}-H_{s-1,t}
=(-1)^n\binom{\ell+t}{n}\Delta^nY_{d-s}}
\tag{2.1}
\]

when \(s>0\), and

\[
\boxed{
H_{s,t}-H_{s,t-1}
=(-1)^n\binom{\ell+t-1}{n}\Delta^nY_{d-s}}
\tag{2.2}
\]

when \(t>0\).  The grid has rational rank exactly

\[
\boxed{m+1=S+T+1.}
\tag{2.3}
\]

Let

\[
\mathcal Q=\prod_{\substack{q\ {m prime}\\q-1\in[d,d+L]}}q.
\]

All core primes are assumed larger than every stencil length.  For each such \(q\), every row is the same selector modulo \(q\):

\[
H_{s,t}\equiv Y_{q-1}\pmod q.
\tag{2.4}
\]

The local Smith form is

\[
\operatorname{diag}(1,q,\ldots,q)
\tag{2.5}
\]

with \(m\) copies of \(q\).  Hence

\[
q\mid\gcd_{s,t}H_{s,t}
\iff q\mid Y_{q-1},
\tag{2.6}
\]

so the core-prime part is exactly the target part.

The appropriate integral operation is therefore **saturation away from \(\mathcal Q\)**.  If \(\mathcal L\) is the row lattice of the margin coefficient matrix, define

\[
\mathcal L^{(\mathcal Q)}
=\{x\in\mathbf Z^N:\exists c,\ \gcd(c,\mathcal Q)=1,\ cx\in\mathcal L\}.
\tag{2.7}
\]

This removes every formal Smith factor at primes outside the core and leaves (2.5) untouched.  Let \(K^{\rm sat}_{d,L;S,T}\) be the positive generator of the evaluated ideal of this saturated lattice.  All the theorems below concern this post-Smith gcd.  The deleted factors are known target-blind presentation content.

For later use, a path through the grid gives a basis consisting of one base carrier and \(m\) high-difference rows.  Taking the one-sided path \(S=0\) is especially clean:

\[
H_{0,t}-H_{0,t-1}
=(-1)^{L+t}\binom{\ell+t-1}{L+t}\Delta^{L+t}Y_d,
\qquad1\le t\le T.
\tag{2.8}
\]

After stripping the known core-prime factor and the outside Smith content, these are primitive consecutive high differences.

---

## 3. Recurrence state and the augmented observation matrix

Assume that on the relevant first-cell interval

\[
\boxed{
\sum_{j=0}^{R}a_j(M,r)F_M(r+j)=0,}
\tag{3.1}
\]

where \(a_j\in\mathbf Z[M,r]\), \(a_0a_R\ne0\), and the order is \(R\).

The key point is that the margin family is formed from \(Y\), not from \(F\).  Every Newton row has coefficient sum one, so

\[
H_{s,t}(Y)=b_M+H_{s,t}(F).
\tag{3.2}
\]

Thus the relevant state has dimension \(R+1\), not \(R\):

\[
\boxed{
z_c=igl(b_M,F_M(c),F_M(c+1),\ldots,F_M(c+R-1)\bigr)^t.}
\tag{3.3}
\]

The union of all \(r\)-indices used by the margin rectangle is

\[
u=M-d-L-T,\qquad v=M-d+S.
\tag{3.4}
\]

Choose \(c\) with \(u\le c\le v-R+1\).  Let \(T_r\) be the usual companion transition matrix taking

\[
(F(r),\ldots,F(r+R-1))^t
\longmapsto
(F(r+1),\ldots,F(r+R))^t.
\]

Its denominators are supported on \(a_R(M,r)\), and those of \(T_r^{-1}\) are supported on \(a_0(M,r)\).  Define the explicit recurrence-singularity integer

\[
\boxed{
\Sigma_{u,c,v}(M)=
\prod_{r=u}^{c-1}a_0(M,r)
\prod_{r=c}^{v-R}a_R(M,r),}
\tag{3.5}
\]

with empty products equal to one.

For every prime \(p\nmid\Sigma_{u,c,v}\), every value used by the margins is an invertible linear image of the state at \(c\) over \(\mathbf F_p\).  Reduce a basis of \(\mathcal L^{(\mathcal Q)}\) through these transition matrices, clear the common denominator supported on \(\Sigma_{u,c,v}\), and remove all powers supported on \(\mathcal Q\Sigma_{u,c,v}\).  This gives an integral **outside observation matrix**

\[
W_{d,L;S,T}(M,c)\in
\mathbf Z^{(m+1)\times(R+1)}
\tag{3.6}
\]

such that, for every prime \(p\nmid\mathcal Q\Sigma_{u,c,v}\),

\[
p\mid K^{\rm sat}_{d,L;S,T}
\quad\Longrightarrow\quad
W_{d,L;S,T}(M,c)z_c\equiv0\pmod p.
\tag{3.7}
\]

This matrix is completely explicit from:

- the Newton/Pascal weights;
- the leading and trailing recurrence coefficients;
- the chosen initial state location \(c\).

Equivalently, if \(E\) denotes the shift Ore operator, reduce the saturated margin operators modulo the recurrence operator

\[
P_M(r,E)=\sum_{j=0}^{R}a_j(M,r)E^j.
\]

The coefficient rows of the remainders of degree below \(R\), augmented by the constant coordinate, are the rows of \(W\).  Its maximal minors are therefore genuine recurrence/margin subresultants, not a repackaging of the evaluated gcd.

---

## 4. Sharp abstract recurrence-to-gcd theorems

Let

\[
n=R+1
\]

be the augmented state dimension.

### 4.1 The projective theorem: \(S+T=R-1\)

Assume

\[
m=R-1,
\]

so \(W\) has \(R\) rows and \(R+1\) columns, and assume it has rational rank \(R\).  Let

\[
\kappa_j=(-1)^j\det W_{\widehat j},
\qquad0\le j\le R,
\tag{4.1}
\]

where \(W_{\widehat j}\) deletes column \(j\).  Then

\[
W\kappa^t=0.
\]

Define the explicit projective state resultant

\[
\boxed{
\Xi(W,z_c)=\gcd_{0\le i<j\le R}
(\kappa_i z_{c,j}-\kappa_j z_{c,i}).}
\tag{4.2}
\]

#### Theorem 4.1

Let \(R_{\rm core}\) be the product of the actual core targets.  If \(\Xi(W,z_c)\ne0\), then

\[
\boxed{
\operatorname{rad}
\frac{K^{\rm sat}_{d,L;S,T}}{\gcd(K^{\rm sat}_{d,L;S,T},\mathcal Q^\infty)}
\mid
\operatorname{rad}\bigl(\Sigma_{u,c,v}\,\Xi(W,z_c)\bigr).}
\tag{4.3}
\]

Together with (2.6), this gives

\[
\boxed{
\operatorname{rad}K^{\rm sat}_{d,L;S,T}
\mid
R_{\rm core}\,
\operatorname{rad}\bigl(\Sigma_{u,c,v}\,\Xi(W,z_c)\bigr).}
\tag{4.4}
\]

#### Proof

Let \(p\nmid\mathcal Q\Sigma_{u,c,v}\) divide the saturated gcd.  Then \(Wz_c=0\) modulo \(p\).  If \(W\) has rank below \(R\) modulo \(p\), all cofactors \(\kappa_j\) vanish modulo \(p\), so every integer in (4.2) is divisible by \(p\).  If the rank is \(R\), its kernel is the line generated by \(\kappa\), hence \(z_c\) is proportional to \(\kappa\) modulo \(p\), and again all wedges in (4.2) vanish.  Thus \(p\mid\Xi\).  The excluded primes divide \(\mathcal Q\Sigma\).  Core primes are handled exactly by (2.6).  ∎

This is the sharpest theorem with only \(R\) observations.  It is also usually quantitatively useless: \(\Xi\) contains the full actual initial state, whose entries can have logarithmic height \(\Theta(M)\), and \(\Xi\) can be zero if the state is rationally aligned with the kernel line.

### 4.2 The primitive-state theorem: \(S+T\ge R\)

Assume

\[
m\ge R.
\]

Let \(\Delta_n(W)\) be the gcd of all \(n\times n\) minors of \(W\), with \(\Delta_n(W)=0\) if the rational rank is below \(n\).  Suppose there is an explicit Bézout identity

\[
\boxed{
U_{-1}(M,c)b_M+
\sum_{i=0}^{R-1}U_i(M,c)F_M(c+i)=J_M(c),}
\tag{4.5}
\]

where all \(U_i,J_M(c)\) are integers and \(J_M(c)\ne0\).  Equivalently,

\[
\gcd z_c\mid J_M(c).
\tag{4.6}
\]

#### Theorem 4.2

If \(\Delta_n(W)\ne0\), then

\[
\boxed{
\operatorname{rad}K^{\rm sat}_{d,L;S,T}
\mid
R_{\rm core}\,
\operatorname{rad}\left(
\Sigma_{u,c,v}(M)\,\Delta_n(W)\,J_M(c)
\right).}
\tag{4.7}
\]

For the minimal clean choice \(m=R\), one may take a path of \(R+1\) independent margin rows; then \(W\) is square and \(\Delta_n(W)=|\det W|\).

#### Proof

Take a prime \(p\nmid\mathcal Q\Sigma\Delta_n(W)\) dividing the saturated gcd.  Some \(n\times n\) submatrix of \(W\) is invertible modulo \(p\).  Equation (3.7) therefore gives

\[
z_c\equiv0\pmod p.
\]

The Bézout identity (4.5) gives \(p\mid J_M(c)\), a contradiction unless \(p\) is in the right side of (4.7).  Core primes are again exactly the targets.  ∎

### 4.3 Why this is genuinely a Casoratian/resultant theorem

For the one-sided choice \(S=0,T=R\), (2.8) shows that, after the known Smith factors are removed, the homogeneous rows are

\[
\nabla^{L+1}F_M(r_0),\ldots,
\nabla^{L+R}F_M(r_0)
\tag{4.8}
\]

up to signs, where \(r_0=M-d\).  Reducing the corresponding Ore operators

\[
(1-E^{-1})^{L+1},\ldots,(1-E^{-1})^{L+R}
\]

modulo \(P_M\) produces the homogeneous block of \(W\).  Its determinant is the discrete observability Casoratian; its maximal-minor ideal is the corresponding Ore subresultant ideal.  Thus (4.7) is exactly the desired mechanism: an outside common prime is forced into a recurrence singularity, a Casoratian/resultant, or the primitive-state scalar.

---

## 5. Exactly how many margins are necessary

The answer has two levels.

### 5.1 Dimensionally sharp projective control

The augmented state dimension is \(R+1\), while the margin rank is \(S+T+1\).

- If \(S+T\le R-2\), the observation kernel has generic dimension at least two.  No state-independent nonzero scalar can control all primitive states in that kernel.
- If \(S+T=R-1\), the generic kernel is a line.  The projective resultant (4.2) is the sharp exceptional integer.

Hence

\[
\boxed{S+T=R-1}
\tag{5.1}
\]

is the minimum for any nontrivial projective elimination.

### 5.2 Clean primitive-state control

To force the entire augmented state to vanish modulo a prime outside the determinant, one needs \(R+1\) independent observations.  Therefore

\[
\boxed{S+T=R}
\tag{5.2}
\]

is the minimum for Theorem 4.2.

The split between left and right margins is irrelevant to the dimension count.  A one-sided family \((S,T)=(0,R)\) is algebraically simplest.  A balanced family may have a smaller numerical determinant, but that is additional arithmetic, not a rank improvement.

If one derives a recurrence directly for \(Y\), rather than for \(F=Y-b_M\), and its order is \(R_Y\), then the state dimension is only \(R_Y\), and the corresponding clean threshold is \(S+T=R_Y-1\).  Under the hypothesis stated in the question, however, the constant \(b_M\) is an extra state coordinate and (5.2) is the correct threshold.

---

## 6. A fixed-order recurrence alone is insufficient

### 6.1 Scalar-content counterexample

Let \(f(r)\) be any integer P-recursive sequence satisfying a primitive order-\(R\) recurrence.  For an arbitrary integer \(C\), put

\[
b=C b_0,\qquad F(r)=C f(r),\qquad Y(r)=C(b_0+f(r)).
\]

The recurrence operator, its order, its degree, and its coefficient height are unchanged, but every Newton extrapolant is multiplied by \(C\):

\[
H_{s,t}(Y)=C H_{s,t}(b_0+f).
\]

Hence the whole margin gcd contains \(|C|\).  Since \(C\) is arbitrary, no theorem depending only on the recurrence can bound the gcd.

The extreme order-one example is the constant shell

\[
Y_d=C,\qquad F(r)=0,\qquad F(r+1)-F(r)=0.
\]

Every Newton row has coefficient sum one, so every \(H_{s,t}=C\), for every number of margins.

### 6.2 Why primitive recurrence coefficients do not help

Making \(\gcd(a_0,\ldots,a_R)=1\) controls the operator, not the solution lattice.  Multiplying one solution by \(C\) still gives a solution.  The missing datum is the content of the actual augmented state.

### 6.3 Why \(R-1\) margins need more than state primitivity

With \(S+T=R-1\), let \(\kappa\) be the primitive kernel vector of \(W\).  There is an integer solution of the recurrence whose augmented initial state equals \(\kappa\).  That state can be primitive, yet all \(R\) margin observations vanish identically.  Thus the condition \(\gcd z_c=1\) does not rescue the projective theorem; one needs the nonalignment scalar \(\Xi(W,z_c)\ne0\) with a small height bound.

### 6.4 Minimal useful extra hypothesis

The minimal clean hypothesis is the **augmented primitive-state Bézout property**:

> There are explicit integer polynomials or integer-valued functions
> \(U_{-1},U_0,\ldots,U_{R-1},J\), uniform on every relevant block, such that
> \[
> U_{-1}b_M+\sum_{i=0}^{R-1}U_iF_M(r+i)=J_M(r),
> \]
> \(J_M(r)\ne0\), and both the coefficient heights and \(\log|J_M(r)|\) are polylogarithmic in \(M\).

The strongest form is \(J_M(r)=1\).  A radical version is enough for P3.2 if

\[
\log\operatorname{rad}J_M(r)=\operatorname{polylog}M.
\]

It is important that this identity be **local and uniform in \(r\)**.  Propagating a primitive state from \(r=0\) to a block at distance \(\Theta(M)\) introduces the product of \(\Theta(M)\) leading/trailing recurrence coefficients, which has linear or worse logarithmic height and destroys the block budget.

---

## 7. Complete height ledger

Let

\[
N=L+S+T,\qquad w=v-u+1=L+S+T+1.
\]

Suppose:

- recurrence order \(R\);
- coefficient degree at most \(\delta\);
- logarithmic coefficient height at most \(h\);
- \(|r|,M\le M+O(w)\).

Put

\[
A=h+\delta\log(M+w+2).
\tag{7.1}
\]

Then every evaluated recurrence coefficient has logarithmic size \(O(A)\).

### 7.1 Newton/Pascal rows

For every stencil of length at most \(N\),

\[
\log\left|
\binom{d+i}{i}\binom{d+N+1}{N-i}
\right|
\le 2N\log(M+N+2).
\tag{7.2}
\]

After exact outside Smith saturation, the difference rows are primitive finite-difference rows; their coefficient height is \(O(N)\).  The base row remains bounded by (7.2).  Thus one may take

\[
B_{\rm Newt}=O(N\log(M+N)).
\tag{7.3}
\]

### 7.2 Transfer matrices

A product of at most \(w\) companion or inverse-companion matrices has numerator-entry logarithmic height

\[
B_{\rm tr}=O\bigl(w(A+\log(R+1))\bigr),
\tag{7.4}
\]

and

\[
\log|\Sigma_{u,c,v}|=O(wA).
\tag{7.5}
\]

### 7.3 Observation determinant

After known core-prime powers and recurrence-singularity powers are removed, every entry of the outside observation matrix has height

\[
B_W=O\bigl(N\log(M+N)+w(A+\log(R+1))\bigr).
\tag{7.6}
\]

Hadamard's inequality gives, with \(n=R+1\),

\[
\boxed{
\log|\Delta_n(W)|
\le n\left(B_W+\tfrac12\log n\right).}
\tag{7.7}
\]

Consequently Theorem 4.2 gives the block-local exceptional height

\[
\boxed{
\begin{aligned}
\log|\mathcal E_{d,L}|
\ll{}&wA\\
&+(R+1)\left[
N\log(M+N)+w(A+\log(R+1))+\log(R+1)
\right]\\
&+\log|J_M(c)|.
\end{aligned}}
\tag{7.8}
\]

If

\[
L=O(\log M),
\quad R,\delta,h=\operatorname{polylog}M,
\quad \log|J_M(c)|=\operatorname{polylog}M,
\]

then

\[
\boxed{
\log|\mathcal E_{d,L}|=\operatorname{polylog}M=o(M^{1/3}).}
\tag{7.9}
\]

This is a genuine improvement over an exponential-height Newton carrier.

### 7.4 Why the projective threshold is not enough quantitatively

For \(S+T=R-1\), the cofactor vector \(\kappa\) has polylogarithmic determinant height, but (4.2) also contains the actual state.  Since

\[
\log|b_M|,\ \log|F_M(r)|=\Theta(M)
\]

in the relevant range, the generic bound is

\[
\log|\Xi(W,z_c)|=\Theta(M),
\]

not \(o(M^{1/3})\).  Therefore the extra margin in (5.2) is not cosmetic: it replaces full state height by state **content**, which a primitive Bézout theorem can make small.

---

## 8. Can one determinant be reused for an entire \(H=M^{1/3}\) block?

### Conditional positive statement

If one proves a single integer \(\mathcal E_{M,B}\), depending only on the \(H\)-block \(B\), such that for every translated short core inside \(B\)

\[
\Sigma_{d}\,\Delta_{R+1}(W_d)\,J_M(c_d)
\mid \mathcal E_{M,B}
\]

up to the known core-prime factors, and

\[
\log|\mathcal E_{M,B}|=\operatorname{polylog}M,
\]

then yes: one reusable determinant per block has height \(o(H)\), and the global product over \(O(M/H)\) blocks has logarithmic height \(o(M)\).

### What a recurrence alone actually gives

A polynomial-coefficient recurrence produces a different observation determinant

\[
\Delta_{R+1}(W_d)
\]

at each translation \(d\).  A nuisance prime for one core need divide only that one value.  There is no reason for it to divide the determinants of the neighboring cores.

There are three generic ways to combine the translations, and all lose the required saving:

1. **Multiply the short-core determinants.**  There are \(\asymp H/L\) target-blind cores.  The generic bound is
   \[
   \frac HL\operatorname{polylog}M,
   \]
   which is not automatically \(o(H)\); with the natural \(L=\Theta(\log M)\) and the Newton/transfer factor in (7.8), it is typically \(H\operatorname{polylog}M\).
2. **Use one core of length \(H\).**  Then \(N\asymp H\), and (7.8) gives at least an \(O(H\operatorname{polylog}M)\) upper ledger, not \(o(H)\).
3. **Take a symbolic norm/resultant over all translations.**  If \(\Omega_M(d)\) is the symbolic observability determinant, the natural interval norm is
   \[
   \prod_{d\in B}\Omega_M(d),
   \]
   or the equivalent resultant against \(\prod_{d\in B}(X-d)\).  Its logarithmic height is again proportional to \(H\) in general.

Moreover, taking the gcd of the translated determinants controls only primes common to **every** core, whereas a nuisance prime may occur in one core only.

Therefore

\[
\boxed{
\text{One reusable determinant per }M^{1/3}\text{-block is a separate uniform fixed-divisor theorem.}}
\tag{8.1}
\]

It is not a consequence of fixed recurrence order.  The required additional statement can be phrased as a uniform divisibility of all translated observability determinants and primitive-state scalars by one polylogarithmic-height block integer.

This also explains why the recurrence route is better than one \(O(\log^2M)\) certificate per **realized** target pair only if that extra reuse theorem is proved.  Without reuse, it merely replaces realized-pair dependence by a target-blind but still linear-size covering family.

---

## 9. Application to the exact first-cell decomposition

The source-verified decomposition is as follows.  Put

\[
A_t=\binom Mt,\quad N_t=2M-t,\quad
B_t=\binom{N_t}{M},\quad P_t(r)=\binom{N_t}{r}.
\]

Then

\[
\boxed{
K_M(r)=\sum_{t=0}^{M}A_t^2
\left(2B_tP_t(r)+P_t(r)^2\right).}
\tag{9.1}
\]

For \(t\le r\), put

\[
Q_t(r)=\binom{N_t}{r-t},\qquad
U_t(r)=\binom M{r-t}.
\]

The low boundary is

\[
\boxed{
\begin{aligned}
L_M(r)=\sum_{t=0}^{r}A_t\Bigl[
&A_t\bigl(2(B_t+P_t)Q_t+Q_t^2\bigr)\\
&+U_t(B_t+P_t+Q_t)^2
\Bigr].
\end{aligned}}
\tag{9.2}
\]

The high boundary is

\[
\boxed{
H_M(r)=\sum_{k=0}^{r}
\binom Mk\binom M{r-k}
\left(\binom{M+k}{k}+\binom{M+k}{r}\right)^2.}
\tag{9.3}
\]

Thus

\[
\boxed{F_M(r)=K_M(r)+L_M(r)+H_M(r).}
\tag{9.4}
\]

### 9.1 The long core is not presently fixed-order

Split (9.1) into

\[
K_M=2K_M^{\rm lin}+K_M^{\rm sq}.
\]

After substituting \(k=M-t\),

\[
K_M^{\rm sq}(r)
=\sum_{k=0}^{M}\binom Mk^2\binom{M+k}{r}^2.
\tag{9.5}
\]

Its generating polynomial is exactly

\[
\mathcal S_M(z)=
\sum_{k=0}^{M}\binom Mk^2J_{M+k}(z),
\qquad
J_n(z)=\sum_r\binom nr^2z^r.
\tag{9.6}
\]

Q5716 proves the exact Green/Christoffel--Darboux identity

\[
E_M(z)\mathcal S_M(z)
=-(2M+1)\eta_M(z)J_{2M+1}(z)
-M(1-z)^2\zeta_M(z)J_{M-1}(z),
\tag{9.7}
\]

where \(E_M\) is an explicit tridiagonal continuant of degree at most \(M+1\).  At \(z=1\),

\[
\boxed{
E_M(1)=(-1)^{M+1}2^{M+1}
\prod_{k=0}^{M}(2M+2k+1).}
\tag{9.8}
\]

Hence

\[
\log H(E_M)=\Theta(M\log M).
\tag{9.9}
\]

The linear part \(K_M^{\rm lin}\) is much smaller and has the explicit transformed terminating \({}_3F_2\) operator of order three.  The growing obstruction is (9.5), not the linear term.

### 9.2 What the boundary terms would have to do

The sums \(L_M,H_M\) are triangular boundary forcings.  They do not alter the tridiagonal bulk operator behind (9.7).  A fixed-order recurrence for the complete \(F_M\) can therefore occur only if their forcing numerators cancel almost all of the bulk continuant.

The exact scalar that must become small is a full-forcing analogue of

\[
\bar E_M(z)=
\frac{E_M(z)}{\gcd(E_M(z),\eta_M(z),\zeta_M(z),	ext{boundary forcing numerators})}.
\tag{9.10}
\]

A fixed-order, polylogarithmic-height recurrence would imply a very strong version of this cancellation.  Generic holonomic closure does not prove it; the explicit descriptor construction in Q5716 has rank growing linearly with \(M\).

Thus deriving the recurrence is a meaningful test:

- if the minimal order or primitive denominator still grows like \(M\), the route is quantitatively closed;
- if a genuinely fixed-order primitive operator appears, it records a new, highly nontrivial cancellation between \(K_M,L_M,H_M\).

### 9.3 Even a miraculous recurrence leaves the primitive-state gap

At \(r=0\), the exact endpoint formula gives

\[
\boxed{
F_M(0)=
2\sum_{k=0}^{M}\binom Mk^2\binom{M+k}{k}
+\binom{2M}{M}^2
+7\binom{2M}{M}+11.}
\tag{9.11}
\]

This has exponential size.  No current identity proves that

\[
\gcd\bigl(b_M,F_M(r),\ldots,F_M(r+R-1)\bigr)
\]

has polylogarithmic radical or equals one uniformly in \(r\).  The constant \(b_M\) cannot be omitted: by (3.2), it survives in every \(H_{s,t}\).

Therefore the complete recurrence route still requires a theorem of the form

\[
\boxed{
U_{-1}(M,r)b_M+
\sum_{i=0}^{R-1}U_i(M,r)F_M(r+i)=J_M(r),
\qquad
\log\operatorname{rad}J_M(r)=\operatorname{polylog}M.}
\tag{9.12}
\]

This identity must come from the distinguished Apéry state, not from P-recursiveness.

---

## 10. Final route selection

The recurrence program has three logically separate stages.

### Stage A: primitive recurrence

Prove a recurrence for the **complete** \(F_M=K_M+L_M+H_M\) with

\[
R,\ \deg a_j,\ \log H(a_j)=\operatorname{polylog}M.
\]

This already requires the near-total cancellation of the Q5716 continuant obstruction.

### Stage B: local observability and primitive state

For \(S+T=R\), prove uniformly that

\[
\Delta_{R+1}(W_{d,L;S,T})\ne0
\]

and has polylogarithmic height after the known core and recurrence-singularity factors are removed.  Prove the augmented Bézout identity (9.12).  Then Theorem 4.2 gives a genuine \(o(M^{1/3})\) nuisance bound for one short target-blind core.

### Stage C: block reuse

Prove that all translated short-core exceptional integers in one \(M^{1/3}\)-block divide one polylogarithmic-height block integer, or prove an equivalent fixed-divisor/resultant statement.  Without this stage, summing local certificates does not automatically give \(o(M^{1/3})\).

The exact decision is therefore:

\[
\boxed{
\begin{array}{c}
\text{A fixed-order recurrence would be a major structural breakthrough,}\
\text{but it would not close P3.2 by itself.}\[1mm]
\text{The augmented primitive-state theorem is load-bearing,}\
\text{and block-level reuse is a further independent globalization step.}
\end{array}}
\]

In particular, “derive the recurrence” is a plausible research route only if it is pursued together with the two arithmetic outputs it must expose: a small primitive observability determinant and a uniform small Bézout scalar for \((b_M,F_M(r),\ldots,F_M(r+R-1))\).  Otherwise the recurrence merely repackages the same common-gcd problem in a finite state space.