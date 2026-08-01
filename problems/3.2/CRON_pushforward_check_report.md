# Rank-corrected Franel pushforward: numerical verification and self-twist audit

## Scope and verdict

This report executes `CODEX_SPEC_CRON_pushforward_check.md` on

\[
P=\{29,37,41,53,61,73,89,101\}.
\]

No requested prime degenerates: all are prime to `6`, the degree-two map is
separable, and `q(t)=t^2-34t+1` has no repeated root.  The computation uses
only Python's standard library and exact integer/finite-field arithmetic.

The two sanity gates pass for all eight primes.  The small-order Kummer test
finds only the trivial character for both the pushforward trace and its
`q(t)`-quadratic partner.  The unique inversion shift is

\[
\boxed{c=0}.
\]

The inherited **rank-18** claim is refuted, not merely left unverified.  The
ledger's `[FALSE-3F2]` correction says that the Franel local system has rank
two.  Therefore the literal tensor-square pushforward has generic rank

\[
  2\cdot \operatorname{rank}(\mathcal F\otimes\mathcal F)
  =2\cdot 4=8,
\]

while the minimal symmetric-square pushforward has rank `2 * 3 = 6`.
The rank-18 arithmetic in Q6375 used the already-refuted rank-three Franel
premise.

## Exact conventions

The cover and its branch polynomial are

\[
 \phi(u)=\frac{u(1-8u)}{1+u},\qquad q(t)=t^2-34t+1.
\]

The exact model located in the life-side scripts is

\[
 E_u:\quad y^2+(1-2u)xy+u^2y=x^3.
\]

For `u != -1`, the script defines

\[
 T_F(u)=p+1-\#E_u(\mathbf F_p)
\]

by direct point count and then

\[
 T_G(t)=\sum_{\substack{u\in\mathbf F_p\\ \phi(u)=t}}T_F(u)^2.
\]

For every smooth source fibre, the script checks both

\[
 T_F(u)\equiv H_p(u)\pmod p,
 \qquad T_F(u)^2<4p,
\]

and verifies that `T_F(u)` is the unique centered lift of `H_p(u)`.
The nodal source fibres `u=0,1/8` both lie over `t=0`; direct point count
still agrees with the centered `H_p`-lift there.  Every calculation on
`G_m` uses only smooth source fibres.

For each `t`, the actual list of rational preimages is constructed and
checked against

\[
 N(t)=1+\chi_2(q(t)).
\]

This is checked for every `t`, not just sampled values.

## Task 1: pushforward trace substrate

The following columns give the numbers of `t`'s with `N(t)=0,1,2`, the
number of smooth elliptic point counts checked, and summary invariants of
the vector `(T_G(t))`.  `TG hash` is the first 16 hexadecimal digits
of SHA-256 applied to the comma-separated full integer trace vector in
increasing `t`-order.  The full vectors can be reproduced as CSV with
`--dump-traces`.

| p | N=0/1/2 | smooth checks | support on G_m | sum T_G | sum T_G^2 | TG hash |
|---:|---:|---:|---:|---:|---:|:---|
| 29 | 15/0/14 | 26 | 10 | 720 | 51,840 | `e86e77e94e42ab98` |
| 37 | 19/0/18 | 34 | 17 | 1,000 | 155,456 | `1251d28ea2fdfe0e` |
| 41 | 20/2/19 | 38 | 15 | 1,512 | 233,280 | `baac7ed46778e3b2` |
| 53 | 27/0/26 | 50 | 22 | 2,448 | 425,088 | `0146201f75cf5146` |
| 61 | 31/0/30 | 58 | 29 | 4,072 | 961,856 | `f9f93dbf87af0c21` |
| 73 | 36/2/35 | 70 | 36 | 4,888 | 1,533,600 | `54d87bb022a98180` |
| 89 | 44/2/43 | 86 | 37 | 6,840 | 2,731,968 | `3d5edac05be17a08` |
| 101 | 51/0/50 | 98 | 42 | 10,512 | 4,987,008 | `ec7463d2d71c5d90` |

For every `p,t`, the integer trace data also satisfy

\[
 T_G(t)\equiv(1+\chi_2(q(t)))A_p(t)\pmod p.
\]

**Task 1: PASS for every `t` and every requested prime.**

## Task 2: Mellin inversion and the normalization correction

For nonsplit fibres, the script works in
\(\mathbf F_p[\sqrt{q(t)}]\), takes

\[
 x=\frac{1-t+\sqrt{q(t)}}{16},
\]

checks `8x^2+(t-1)x+t=0`, and verifies the corrected pullback identity

\[
 \boxed{
 A_p(t)=\frac{H_p(x)^2}{(1+x)^{p-1}}.}
\]

The quotient is fixed by Frobenius, lies in \(\mathbf F_p\), and equals direct
evaluation of the Apéry truncation.  The normalization factor is genuinely
needed: the table records `nonsplit / naive failures`, where a naive failure
means `H_p(x)^2 != A_p(t)` before division.

| p | nonsplit / naive failures | Mellin-vector hash |
|---:|---:|:---|
| 29 | 15 / 12 | `2c0ec3a5b4710a77` |
| 37 | 19 / 16 | `817c82f76eeef70a` |
| 41 | 20 / 18 | `484af910bd19a70d` |
| 53 | 27 / 24 | `a530e9763e89d6c4` |
| 61 | 31 / 28 | `687339f10397943a` |
| 73 | 36 / 34 | `5ea71ea79de8bd93` |
| 89 | 44 / 40 | `e9d51a260d235968` |
| 101 | 51 / 48 | `1573e87e558e8eec` |

Put, modulo `p`,

\[
 V_p(t)=-T_G(t)+\chi_2(q(t))A_p(t).
\]

The pointwise test gives `V_p(t)=-A_p(t)` for every `t`.  Thus, for
`N=p-1`,

\[
 R_p(r)=\sum_{t\in\mathbf F_p^*}V_p(t)t^{-r}
\]

satisfies

\[
 R_p(r)=b_r\quad(1\le r\le N-1).
\]

There is one necessary endpoint convention that was absent from the literal
specification.  On \(\mathbf F_p^*\), the exponents `0` and `N` are the
same character, so orthogonality sees both coefficients `b_0` and `b_N`:

\[
 R_p(0)=R_p(N)=b_0+b_N.
\]

The all-(r) corrected formula is therefore

\[
 \boxed{
 b_r=R_p(r)-\mathbf1_{r=0}b_N-\mathbf1_{r=N}b_0,
 \qquad 0\le r\le N.}
\]

The raw formula fails exactly at `{0,p-1}` for every requested prime;
the displayed endpoint correction passes for every `0 <= r < p`.  This is
character aliasing, not another `(1+x)^(p-1)` error.

**Task 2: PASS for every requested `p,r`, with the explicit endpoint
correction above.**

## Task 3: small Kummer self-twists

For each prime, the least primitive root `g` is used.  A character is
represented exactly by `k mod (p-1)`, with
\(\chi_k(g)=\zeta_{p-1}^k\).  Its order is
\((p-1)/\gcd(k,p-1)\).  Every character of order at most 30 is tested
(this includes every character whose order divides 24).  No complex
rounding is involved: the equality
\(T(t)\chi(t)=T(t)\) is checked as
\(k\log_g(t)=0\pmod{p-1}\) on the exact nonzero support.

Let

\[
 T_{G_q}(t)=\chi_2(q(t))T_G(t).
\]

This pointwise identity is checked for the full arrays.  It is important that
\(\chi_2(q(t))\) is a character of the argument \(q(t)\), not a
multiplicative character of \(t\): exhaustive multiplication tables show
that \(t\mapsto\chi_2(q(t))\) is nonmultiplicative for every prime in \(P\).

| p | least g | characters tested | passing for G | passing for G_q | chi_2(q(t)) multiplicative? |
|---:|---:|---:|:---|:---|:---:|
| 29 | 2 | 28 | trivial | trivial | no |
| 37 | 2 | 24 | trivial | trivial | no |
| 41 | 6 | 24 | trivial | trivial | no |
| 53 | 2 | 28 | trivial | trivial | no |
| 61 | 2 | 44 | trivial | trivial | no |
| 73 | 5 | 36 | trivial | trivial | no |
| 89 | 3 | 28 | trivial | trivial | no |
| 101 | 2 | 40 | trivial | trivial | no |

Hence the only candidate that passes at every prime is the trivial character.
There are no nontrivial candidates to report in the tested range.

**Task 3: PASS; no nontrivial `t`-Kummer self-twist of order at most 30
was found for either trace function.**

### Q6375's direct pair formula is refuted

The finite-field correction in Task 2 is
`chi_2(q(t)) A_p(t)`, not
`chi_2(q(t)) T_G(t)`.  These cannot be interchanged.  Indeed, away from
the branch locus, `T_G` is supported only where `chi_2(q(t))=1`, so
`T_G-T_{G_q}` is zero away from rational branch points.

The literal Q6375 formula

\[
 b_r\stackrel?=M(r;G)-M(r;G\otimes\mathcal L_q)
\]

was tested for all \(0\le r<p\).  The table gives accidental matches and the
exact support size of `T_G-T_{G_q}`.

| p | matches / p | support of T_G-T_Gq |
|---:|---:|---:|
| 29 | 0/29 | 0 |
| 37 | 2/37 | 0 |
| 41 | 2/41 | 0 |
| 53 | 0/53 | 0 |
| 61 | 2/61 | 0 |
| 73 | 2/73 | 2 |
| 89 | 0/89 | 0 |
| 101 | 0/101 | 0 |

Thus the Q6375 direct-difference packaging is **REFUTED** at the trace level.
This does not refute the corrected K.2 virtual trace used in Task 2.

## Task 4: inversion symmetry

The script checks the stronger pointwise integer identities

\[
 T_G(t)=T_G(t^{-1}),\qquad
 T_{G_q}(t)=T_{G_q}(t^{-1})
\]

for every \(t\in\mathbf F_p^*\), as well as the corresponding identity for
the corrected virtual trace modulo `p`.  The twist factor is compatible
because

\[
 t^2q(t^{-1})=q(t),
\]

and \(t^2\) is a square.  The mod-\(p\) Mellin vectors were then searched over
all shifts `c mod (p-1)`.  Pointwise integer inversion proves the `c=0`
identity for exact character values, while the mod-\(p\) search excludes
every other shift.  For `G`, `G_q`, and the virtual trace, every prime
has the unique shift set

\[
 \{c:M(r)=M(-r-c)\ \forall r\}=\{0\}.
\]

This is precisely the character-space form of
`b_{p-1-r}=b_r`: modulo `p-1`, the reflected exponent is `-r`.

**Task 4: PASS with the unique fixed shift `c=0`.**

## Task 5: rank probe

The `m=1` integer trace substrate and all mod-`p` Mellin vectors are
computed; their hashes are listed above.  Determining the generic
L-function degree from Frobenius power sums requires genuine
\(T_G^{(m)}\) over extension fields.  The \(\mathbf F_{p^2}\) arithmetic used
in Task 2 only evaluates the Hasse polynomial and is not an extension-field
Frobenius trace.  A full point-count implementation over all
\(\mathbf F_{p^m}\)-parameters would be substantially heavier and was skipped
as explicitly allowed by the specification.

**Task 5: SKIPPED beyond `m=1`; the generic Mellin L-function degree remains
`[待验]`.**  This does not rescue rank 18: generic sheaf rank is already fixed
by the corrected rank-two input and finite pushforward degree.

## P.3 claim audit

### VERIFIED by this computation

- The degree-two cover has the exact fibre counts
  \(1+\chi_2(t^2-34t+1)\) for every tested `t`.
- The literal tensor-square integer trace `T_G` is reproducibly constructed
  from the Beauville-IV elliptic fibres.
- The corrected Franel/Apéry pushforward and Mellin identities hold for all
  requested primes and all exponents, with the stated endpoint convention.
- No nontrivial multiplicative self-twist of `t`-character order at most 30
  occurs for `G` or its `q(t)`-quadratic partner.
- The partner really is the pointwise `q(t)`-twist, and the twisting function
  is not a multiplicative character of `t`.
- The numerical duality shift is uniquely `c=0`.

### REFUTED

- `rank(G)=18`: the corrected literal tensor-square rank is 8, and the
  symmetric-square reduction has rank 6.
- Consequently the stated
  \(G_{\mathrm{geom}}\simeq SL_3\rtimes\mathbf Z/2\subset GL_{18}\) package is
  false as written; it is built on the false rank-three Franel premise.
- The direct trace-level formula
  \(b_r=M(G)-M(G\otimes\mathcal L_q)\) is false.  The Task-2 correction uses
  the Apéry companion `chi_2(q(t)) A_p(t)`, not
  `chi_2(q(t)) T_G(t)`.

### Still `[待验]`

- The corrected full geometric/arithmetic Tannakian group for the rank-eight
  tensor-square pushforward or rank-six symmetric-square pushforward.
- Local monodromy and determinant/autoduality sign at the sheaf level.  The
  branch locus is verified, but a trace table does not prove the local
  monodromy representation.
- Exclusion of Kummer self-twists of every possible order.  The present result
  is exhaustive only through order 30 at the requested primes.
- Any claimed exact size-4 full symmetry group: inversion, the
  `q(t)`-twist, and `t`-Kummer self-twists are different kinds of
  operations and must not be conflated.
- Generic Mellin L-function degree and extension-field power sums.
- The characteristic-zero compatible-system and middle-extension comparison
  needed to promote these finite-field trace identities to the full
  Tannakian conclusion `[GAP-2]`.

## Reproduction

Run from the repository root:

```text
python3 -B problems/3.2/CRON_pushforward_check.py
python3 -B problems/3.2/CRON_pushforward_check.py --dump-traces
```

The first command runs every assertion and prints the summary.  The second
prints 485 CSV lines (one header plus one row for every `p,t`) containing
`p,t,N,T_G,T_G_chi_q,A_from_Franel,virtual_mod`.
