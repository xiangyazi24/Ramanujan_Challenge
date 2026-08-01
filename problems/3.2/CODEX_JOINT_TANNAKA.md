# Joint Tannakian moment computation: stall report

## Verdict

The finite-field Mellin inversion gate passes at all eight requested primes,
but it does **not** determine the characteristic-zero traces of the two deck
eigensheaves.  The eigen-trace equation in the specification conflicts with
the supplied source arithmetic at split fibres, and the requested Mellin
normalization has the wrong Weil weight.  Following the specification's stall
protocol, no moments were computed.

## Exact convention checked

For each rational source fibre I used the rank-three symmetric-square trace

\[
 f_p(u)=T_F(u)^2-p.
\]

Thus the exact symmetric-square pushforward trace on
\(\mathbf F_p^\times\) is

\[
 P_p(t)=\sum_{\phi(u)=t}f_p(u)
       =T_G(t)-pN(t).
\]

The subtraction is once per rational source fibre.  It removes the
determinant/Tate summand from the literal tensor-square trace and agrees
modulo \(p\) with `T_G`.  The reused source arithmetic proves

\[
 P_p(t)\equiv(1+\chi_2(q(t)))A_p(t)\pmod p. \tag{1}
\]

The corrected pullback \(A_p(t)\) supplied by the source script is an element
of \(\mathbf F_p\).  On nonsplit fibres it is recovered in
\(\mathbf F_{p^2}\) using

\[
 A_p(t)=\frac{H_p(x)^2}{(1+x)^{p-1}},\qquad
 x=\frac{1-t+\sqrt{q(t)}}{16}.
\]

This is an exact residual value, not an integral or complex Frobenius trace.

## Mellin gate

Formally imposing the specification's difference residue

\[
 D_p(t)\equiv\chi_2(q(t))A_p(t)\pmod p
\]

gives

\[
 -P_p(t)+D_p(t)\equiv-A_p(t)\pmod p.
\]

Consequently its Mellin inversion is exactly the already-verified Apéry
identity.  For \(N=p-1\), the raw transform has the two endpoint aliases

\[
 R_p(0)=R_p(N)=b_0+b_N,
\]

and subtracting \(b_N\) at \(r=0\), respectively \(b_0\) at \(r=N\),
reproduces every \(b_r\).  The gate-only script verifies all
\(0\le r\le p-1\) exactly.

| p | exponents verified | raw failures before endpoint correction |
|---:|---:|:---|
| 29 | 29 | `0, 28` |
| 37 | 37 | `0, 36` |
| 41 | 41 | `0, 40` |
| 53 | 53 | `0, 52` |
| 61 | 61 | `0, 60` |
| 73 | 73 | `0, 72` |
| 89 | 89 | `0, 88` |
| 101 | 101 | `0, 100` |

This gate is residual and tautological after (1): it cannot select an
integral lift of \(A_p(t)\), hence cannot validate complex moment data.

## Exact conflict in the eigen-trace formula

Section Z.1 states the projection-formula correction: the rank-six
pushforward consists of the two rank-three deck descents.  Equation (1) then
has the pointwise residual decomposition

\[
 A_p(t)+\chi_2(q(t))A_p(t), \tag{2}
\]

so the quadratic companion is one summand, not the difference of the two
summands.  At every split fibre, \(\chi_2(q(t))=1\), the two exact source
`Sym²` traces agree.  Their deck-descended traces therefore agree and have
difference zero.

The specification instead requires

\[
 T_+(t)-T_-(t)=\chi_2(q(t))A_p(t). \tag{3}
\]

This is already false at the first small witness:

\[
 p=29,\quad t=2,\quad \phi^{-1}(2)=\{8,10\},\quad
 f_{29}(8)=f_{29}(10)=7.
\]

Hence the exact deck traces are equal at this split fibre, their sum is 14,
and their difference is 0.  Formula (3) instead asks for difference 7.  With
the only evident lift \(A_{29}(2)=7\), solving the two equations in the spec
would give \(T_+=21/2\), \(T_-=7/2\), not integral Frobenius traces.

The same contradiction occurs at every requested prime:

| p | split witness t | A_p(t) mod p | each exact source `Sym²` trace | centered-lift parity failures |
|---:|---:|---:|---:|---:|
| 29 | 2 | 7 | 7 | 16 |
| 37 | 2 | 4 | -33 | 15 |
| 41 | 1 | 36 | -5 | 23 |
| 53 | 2 | 36 | -17 | 30 |
| 61 | 4 | 3 | 3 | 17 |
| 73 | 1 | 4 | -69 | 29 |
| 89 | 1 | 57 | 235 | 18 |
| 101 | 3 | 21 | 223 | 46 |

Adding multiples of \(p\) to the residual companion can repair parity at an
individual point, but the supplied sources give no rule selecting those
lifts.  Such arbitrary choices change every complex Mellin moment while
leaving the mod-\(p\) gate unchanged.

The source-compatible exact identities should instead have the form

\[
 P_p=T_++T_-,\qquad T_-=\chi_2(q)T_+
\]

away from branch stalk corrections, with
\(T_+\equiv A_p\pmod p\).  Then the Apéry virtual trace is
\(-P_p+T_-=-T_+\), not \(-P_p+(T_+-T_-)\).  Computing \(T_+\) on nonsplit
fibres requires genuine integral Frobenius data; the Hasse--Witt residue in
the supplied script is insufficient.

## Independent normalization conflict

A pure local trace of weight 2 has eigenvalues of size \(p\).  Its
multiplicative Mellin sum is a Frobenius trace of weight 3, hence has natural
size \(p^{3/2}\).  Therefore the bounded normalized trace is

\[
 s(\chi)=\frac{S(\chi)}{p^{3/2}},
\]

not \(S(\chi)/p\).  This is also forced by exact Parseval:

\[
 \frac1{p-1}\sum_\chi |S(\chi)|^2
 =\sum_{t\in\mathbf F_p^\times}|T(t)|^2\asymp p^3.
\]

Thus \(\operatorname{avg}|S/p|^2\asymp p\), rather than tending to 1, and
the asserted bound \(|S/p|\le3\) cannot hold uniformly.  The repository's
existing Mellin ledger likewise records bounds proportional to
\(p^{3/2}\).

## What is needed to resume

Two corrections are required before numerical Goursat moments are meaningful:

1. replace the difference equation (3) by an exact definition of the two
   deck-descended integral Frobenius traces, including nonsplit and branch
   fibres; and
2. replace the normalization \(S/p\) by \(S/p^{3/2}\), or explicitly explain
   an additional prior division by \(\sqrt p\) in the local traces.

Once those are fixed, the relevant theoretical fourth moment is 3 for both
the standard representations of \(O_3\) and \(SO_3\): it is
\(\dim(V^{\otimes4})^G\), generated by the three pair contractions.  For a
product group, \(M_{2,2}=1\); for the diagonal identity-component graph it is
3.  Those comparisons were not applied to unverified numerical traces.

## Limitations

Even corrected moment data could only fingerprint the compact trace
distribution.  It could not prove equality of arithmetic and geometric
Tannakian groups, rule out all small-index or arithmetically exceptional
subgroups, determine integral models, or promote finitely many prime samples
to a compatible-system theorem.  Here the obstruction occurs earlier: the
available finite-field residues do not determine the complex traces whose
moments the specification asks to average.
