# Publication Architecture for Problem 3.2

## Publication thesis

The submitted archive should not be converted into a paper by incremental
editing. Its strongest coherent publication is a short article about an exact
large-prime reduction and a quantitative density-one theorem for common
divisors in Apéry's recurrence. The pointwise conjecture remains open.

### Proposed title

**Large-Prime Reductions and a Quantitative Density-One Theorem for Apéry
GCDs**

### Proposed abstract

Let \(a_n/b_n\) be Apéry's rational approximants to \(\zeta(3)\), let
\(d_n=\operatorname{lcm}(1,\ldots,n)^3\), and put
\(G_n=\gcd(d_na_n,d_nb_n)\).  We prove that for every fixed
\(\varepsilon>0\),

\[
\#\{n\le X:\log G_n>\varepsilon n\}
=O_\varepsilon((\log X)^2).
\]

The proof begins with an exact large-prime reduction of \(\log G_n\) to the
sparse master sum

\[
M(n)=\sum_{\sqrt n<p\le n}\log p\,
      \mathbf1_{p\mid b_{,n\bmod p}},
\]

including a uniform large-prime valuation cap and an \(O(n^{2/3})\) remainder.  Gap
continuants then provide nonzero integer carriers of controlled height for
the primes shared by two nearby indices.  A localized codegree argument turns
this control into the stated polylogarithmic exceptional set.  We also derive
fixed-Apéry-value radical corollaries and an exact digit description of the
top-prime window.  The conjectural pointwise estimate \(\log G_n=o(n)\) for
every \(n\) remains open.

Let

\[
d_n=\operatorname{lcm}(1,\ldots,n)^3,
\qquad
G_n=\gcd(d_na_n,d_nb_n).
\]

The conjectural pointwise estimate

\[
\log G_n=o(n)
\]

is not proved. The primary article should instead establish the following
quantitative substitute as its main theorem: for every fixed
\(\varepsilon>0\),

\[
\#\{n\le X:\log G_n>\varepsilon n\}
=O_\varepsilon((\log X)^2).
\]

The conceptual thesis is that a corrected large-prime decomposition converts
the gcd problem into a sparse prime-indexed master sum, while the continuant
geometry of two nearby indices supplies enough codegree control to rule out
all but polylogarithmically many exceptional integers. This is a genuine
density-one theorem, but it does not control any prescribed exceptional
integer.

The abstract, introduction, theorem statements, and conclusion must all say
explicitly that the every-\(n\) conjecture remains open.

## Contribution hierarchy

The paper should distinguish one principal theorem, three supporting theorem
packages, and one independent arithmetic corollary.

### I. Principal theorem

**Quantitative exceptional set.** For every fixed \(\varepsilon>0\),

\[
\#\{n\le X:\log G_n>\varepsilon n\}
=O_\varepsilon((\log X)^2).
\]

This is the result to feature in the title, abstract, and first paragraph of
the introduction. It is the strongest unconditional progress in the archive
toward the original pointwise problem.

### II. Foundational arithmetic reduction

The exact Wronskian is

\[
a_nb_{n-1}-a_{n-1}b_n=\frac{6}{n^3}.
\]

With \(A_n=d_na_n\) and \(B_n=d_nb_n\), the safe adjacent determinant is

\[
A_{n+1}B_n-A_nB_{n+1}
=\frac{6d_nd_{n+1}}{(n+1)^3},
\]

so that

\[
G_nG_{n+1}\mid\frac{6d_nd_{n+1}}{(n+1)^3},
\qquad
v_p(G_n)\le 6
\quad(\sqrt n<p\le n,\ p\ge5).
\]

The associated master sum is exactly

\[
M(n)=\sum_{\sqrt n<p\le n}
\log p\;\mathbf 1_{p\mid b_{\,n\bmod p}},
\]

and the corrected comparison is

\[
\log G_n\le6M(n)+O(n^{2/3}),
\qquad
M(n)\le\log G_n+O(1).
\]

Thus \(M(n)=o(n)\) is equivalent to the original pointwise conjecture. This is
an exact reduction, not a solution of that conjecture.

The same section may state the fixed-integer corollaries

\[
\mathcal S(n)=\log\operatorname{rad}\gcd\!\left(b_n,
                    \prod_{p\le n}p\right),
\qquad
0\le\mathcal S(n)-M(n)\le12n^{2/3},
\]

and

\[
C_n^{\mathrm{core}}
=\gcd\!\left(\operatorname{rad}(b_n),
             \prod_{n^{2/3}<p\le n}p\right),
\qquad
\left|\log C_n^{\mathrm{core}}-M(n)\right|\le8n^{2/3}.
\]

They turn the moving-residue sum into one fixed Apéry value, but do not
reduce the difficulty: proving either radical has logarithm \(o(n)\) is
equivalent to the open pointwise conjecture.

They nevertheless inherit a genuine unconditional consequence of the
principal theorem.  For every fixed \(\varepsilon>0\),

\[
\#\{n\le X:\mathcal S(n)>\varepsilon n\}
=O_\varepsilon((\log X)^2),
\]

and the same bound holds with \(\mathcal S(n)\) replaced by
\(\log C_n^{\mathrm{core}}\).  This fixed-Apéry-value radical corollary belongs
in the theorem hierarchy; it must not be advertised as pointwise control.

### III. Gap-carrier and codegree theorem

For two indices \(m<n\) in one dyadic block, put \(h=n-m\). The paper should
develop only the continuant facts required to prove that the number of primes
\(p>\sqrt N\) appearing in both large-prime counts is \(O(h)\):

- in the no-wrap case, a common prime divides the nonzero gap carrier
  \(N_h(m)\);
- in the wrap case, it divides the explicit carrier
  \(\prod_{j=1}^{h}(m+j)\);
- both carriers have logarithmic height \(O(h\log N)\).

For the no-wrap carrier, the positivity argument rests on the exact identity

\[
P(k)-k^3-(k+1)^3=4(2k+1)^3>0.
\]

The reflection, renewal, pole-value, root-strip, and characteristic-zero
squarefreeness results are valuable structural strengthenings. They should be
stated together only if their proofs can be included without obscuring the
codegree argument; otherwise the primary article should prove the
dependency-critical carrier statements and move the stronger root theory to a
companion note or supplement.

### IV. Localized packing mechanism

For \(n\in(N,2N]\), let

\[
t(n)=\#\{\sqrt n<p\le n:p\mid b_{\,n\bmod p}\}.
\]

If \(\log G_n>\varepsilon n\), then for sufficiently large \(N\),

\[
t(n)\gg\varepsilon\frac{N}{\log N}.
\]

On an interval of length

\[
Y=c_1\varepsilon^2\frac{N}{\log N},
\]

Cauchy's inequality gives a lower bound for the total prime-pair load, while
the codegree theorem bounds the same load by the sum of pairwise distances.
For sufficiently small absolute \(c_1\), only \(O_\varepsilon(1)\) exceptional
integers can lie in one such interval. Summation over the intervals and dyadic
blocks gives the principal theorem.

The manuscript must call this a **localized codegree-controlled second
factorial moment**. It must not attribute the theorem to the separate global
Chinese-remainder second moment.

### V. Digit and top-window theorem

For

\[
Z_p=\{0\le r<p:p\mid b_r\},
\]

Gessel's congruence gives

\[
b_{mp+r}\equiv b_mb_r\pmod p,
\]

and hence the base-\(p\) digit criterion: \(p\mid b_n\) if and only if at least
one base-\(p\) digit of \(n\) belongs to \(Z_p\). In the top window
\(n/2<p\le n\), this reduces, apart from the explicit finite small-prime
correction, to the digit \(n-p\).

The resulting exact averaging identity identifies

\[
\sum_{p\le X}|Z_p|=o(X)
\]

with the unnormalized Cesàro \(o(X)\) statement for the quotient-one/top-window
contribution. It does not see the quotient ranges \(q\ge2\), is not equivalent
to the pointwise conjecture, and gives no anti-concentration for a fixed
integer. The current unconditional bound remains

\[
\sum_{p\le X}|Z_p|\ll\frac{X^{5/3}}{\log X}.
\]

## Exact dependency chain

The proof architecture is intentionally linear.

1. The Apéry recurrence and integrality give the Wronskian and the integer
   determinant divisibility.
2. The determinant gives the large-prime valuation cap \(6\).
3. The block congruence
   \[
   p^3a_{qp+r}\equiv a_qb_r\pmod p
   \]
   together with the companion numerator-height estimate removes the
   companion channel at total cost \(O(n^{2/3})\).
4. These facts give the exact master sum and its two-sided comparison with
   \(\log G_n\).
5. Gap continuants provide nonzero no-wrap carriers; the explicit wrap carrier
   handles the complementary case. Their height bounds give the linear
   cross-level codegree estimate.
6. The codegree estimate and localized prime-load packing give the
   \(O_\varepsilon((\log X)^2)\) exceptional-set theorem.

Two branches are parallel, not dependencies of this proof:

- Gessel's digit criterion gives the exact quotient-one/top-window averaging
  theorem.
- The global Chinese-remainder identity
  \[
  \sum_{0\le m\le X^2}(K_X(m))_2\le4X^2\lambda_X^2
  \]
  is an independent second-factorial-moment theorem. A hypothetical
  independence-scale moment of one fixed order \(k>6\) would prove the
  pointwise conjecture, but only \(k=2\) is currently proved.

The determinant/symplectic, critical-value, limiting-cell, and regular-gcd
triangle programs have no edge into the unconditional chain above.

## Primary article: section and page budget

The target is a self-contained article of approximately 27 pages, including
references. The budget is a constraint on exposition: additional campaign
material should be removed rather than compressed into footnotes.

| Section | Content | Pages |
|---|---|---:|
| 1. Introduction and theorem statements | Open pointwise problem; principal theorem; three supporting packages and one arithmetic corollary; status box; proof overview | 3 |
| 2. Apéry arithmetic and the Wronskian | Recurrence, integrality, corrected determinant, and valuation cap \(6\) | 3 |
| 3. The exact large-prime reduction | Block congruence, companion height, small/medium-prime disposal, master-sum and smooth-radical equivalences | 4 |
| 4. Gap carriers and cross-level codegree | Continuants, no-wrap and wrap carriers, nonvanishing, height bounds, \(O(h)\) codegree theorem | 5 |
| 5. Exceptional integers | Prime loads, local interval decomposition, Cauchy lower bound, codegree upper bound, dyadic summation | 5 |
| 6. Digits and the top window | Gessel/Lucas criterion, finite small-prime correction, exact Cesàro identity, current \(X^{5/3}/\log X\) bound | 3 |
| 7. Boundary of the method | Global CRT moment, the precise higher-moment gateway, and a short statement of the remaining pointwise obstruction | 2 |
| References and acknowledgements | Only sources used in proofs or historical framing | 2 |
| **Total** |  | **27** |

The introduction should contain one theorem-status box and one dependency
paragraph. Section 7 should be no more than a concise research outlook; it is
not a catalogue of failed campaigns.

## Companion-note split

Material not used in the primary dependency chain should be separated by
mathematical maturity rather than collected in one miscellaneous appendix.

### Companion note A: certified critical-value separation

Proposed title:

**Certified Critical-Value Separation for Apéry Gap Rational Functions Through
Height 60**

This note may contain:

- the exact reconstruction of the critical polynomials;
- the Arb interval-Newton isolation and nonvanishing checks;
- the certified theorem for every \(2\le h\le60\);
- the full reproducibility package, including software versions, source hashes,
  exact inputs, interval outputs, and a deterministic checker;
- the reflection, root-strip, and all-height characteristic-zero
  squarefreeness theory when it improves the exposition.

It must not claim an all-height theorem, a cutoff \(h_0=60\), or a validated
two-regime NEAR/FAR extrapolation. It may cite the proved global limiting-cell
formula from the analytic companion, but every passage from that limit to a
uniform finite-height tail must remain explicitly open.

### Companion note B: the global cell and exact analytic identities

Proposed title:

**The Global Apéry Cell, Its Exact Zero Height, and Positive Remainders**

A separate analytic note should lead with the unconditional connection theorem

\[
J(z)=z^{-3}\kappa(-z)\kappa(z)+z^3\psi(-z)\psi(z)
=\pi^3\cot^3(\pi z)+\frac{\pi^3}{3}\cot(\pi z),
\]

where \(\psi(z)=\kappa(z+1)/(z+1)^3\). Its source-normalization proof must
identify the Golyshev--Zagier asymptotic kappa function with the
direct-positive-path Bloch--Vlasenko Frobenius function, retain the companion
term, derive the exact shift cocycle, and prove polynomial vertical growth from
the Mellin representation. The resulting nonreal zero height is

\[
\cosh(2\pi\eta)=2,
\qquad \eta=\frac{\log(2+\sqrt3)}{2\pi}.
\]

The same note may present the exact positive Apéry remainder,

\[
0<\zeta(2)-\frac{q_n}{b_n}<\frac{1}{nb_n}\qquad(n\ge1),
\]

together with its Gosper/WZ certificate and the exact normalization of the
inhomogeneous companion.  This material is mathematically clean but
independent of the gcd exceptional-set theorem.

The Taylor-germ scalar

\[
(\log\kappa)''(0)=-4\zeta(2)=-\frac{2\pi^2}{3}
\]

is the residue input for the global theorem, not a substitute for it. The
manuscript must show the full two-component cocycle: the one-sided product
\(z^{-3}\kappa(-z)\kappa(z)\) is not periodic by itself. This analytic theorem
is exact but remains logically independent of the gcd exceptional-set theorem.

### Material retained only as research notes

The determinant/sheaf program and the regular-gcd triangle program should
remain research notes until they yield a new theorem at the target scale. A
repaired, saturated radical bound of order \(D^{8/3}\) may be recorded in a
specialized algebraic note with the hypotheses \(p\ge7\), \(D<p\), the radical
and algebraic-closure conventions, and the cut/infinity saturations stated in
full. It must not be described as a subquadratic energy saving or as progress
past the current endpoint wall.

## Status vocabulary

Every theorem, computation, conjecture, and route in the publication package
must carry one of the following labels in the internal claim ledger.

| Label | Meaning |
|---|---|
| **PROVED** | A complete mathematical proof is present, with every imported result and normalization identified. |
| **CERTIFIED** | A finite-range exact or interval computation has a reproducible checker and archived proof objects. |
| **CONDITIONAL** | The implication is proved, but at least one explicitly named hypothesis is open. |
| **HEURISTIC** | Numerical evidence, an asymptotic model, a fitted constant, or an unproved identification. |
| **FAILED** | The proposed route or statement is refuted by an explicit counterexample or complete exponent bookkeeping. |

**OPEN** may be used for a target problem, but not as a substitute for labelling
the status of a claimed result. A statement should never move to **PROVED** or
**CERTIFIED** merely because it appears in a campaign ledger or passed a
non-independent consistency check.

## Claims excluded from the primary article

The following claims are excluded from the abstract, theorem list, and
conclusion unless their named gaps are independently closed before revision:

1. \(\log G_n=o(n)\) for every \(n\).
2. The assertion that \(\sum_{p\le X}|Z_p|=o(X)\) is proved, or that it is
   equivalent to the pointwise conjecture.
3. The assertion that the global Chinese-remainder second moment proves either
   the exceptional-set theorem or the pointwise conjecture.
4. Any unconditional \(O(p\log p)\), \(O(p)\), or fixed-power-saving
   determinant completion; any claim that the required bounded-conductor sheaf
   has been constructed.
5. All-height critical-value separation, a certified cutoff \(h_0=60\), or an
   inference from a finite scan to an infinite tail.
6. Any inference from the exact global limiting-cell function to exact
   finite-height root lattices, an all-height critical-value theorem, or fitted
   parity constants presented as certified finite-height values.
7. A regular-gcd triangle theorem claimed to yield a saving at the target
   scale; pairwise scalar resultants claimed to preserve common root-component
   labels.
8. Fixed-prime coprimality inferred from characteristic-zero coprimality
   without saturation.
9. Restart or renewal described as amplifying one collision into many; it is a
   detector and factorization identity, not a multiplicity pump.
10. All-height absolute irreducibility of the critical-value curves.
11. A blanket claim that the project is formally verified in Lean.
12. Empirical Poisson laws, bounded-looking maxima, fitted constants, or finite
    root tables stated as theorems.

The repaired \(D^{8/3}\) saturated radical bound, the finite
\(2\le h\le60\) certificate, and the exact global cell theorem are valid
results. They are excluded from the primary article because they are independent of its proof
and are easy to overstate when presented beside the open pointwise problem.

## Verification and publication package

The public package should make the textual proof primary and computation
auditable.

1. **Source hierarchy.** The current `CODEX_MAINTHM_report.md`, overridden by
   `ERRATA.md`, governs all formulas and statuses in the primary article.
   `RC_BREAKTHROUGH_report.md` governs the analytic companion.  This
   architecture governs editorial scope. Historical campaign ledgers and
   external-model answers are not authoritative sources.
2. **Proof-to-gate table.** Maintain a table mapping each displayed exact
   identity or finite claim to its textual proof and, where applicable, to an
   independent verification gate.
3. **Main verifier.** Run `python3 CODEX_MAINTHM_verify.py` from a clean
   checkout. Archive the complete output ending in `PASS`, the Python version,
   platform information, source commit, and a hash of the script.
4. **Analytic companion verifier.** Run
   `python3 RC_BREAKTHROUGH_verify.py` and archive its complete output, runtime
   versions, source commit, and script hash.  Its exact symbolic gates support
   the continuant, cocycle, residue, telescoping, and WZ identities; cited
   analytic continuation and growth theorems remain textual source inputs.
5. **Independent ground truths.** Gates must compare genuinely independent
   constructions, use several primes and scales where finite tests are
   relevant, and contain no placeholder, skipped, or default-success branch.
   A script supports a proof; it does not replace the written argument.
6. **Regression cases.** Preserve explicit tests for the corrected factor
   \(6\), valuation cap \(6\), the false \(|R_p(n)|\) incidence formula, the
   \(p=5\) top-window correction, wrap/no-wrap carriers, and saturation at a
   fixed working prime.
7. **Computer-assisted companion.** Archive exact input reconstruction, the
   Arb version and precision, interval-Newton outputs, branch counts, source
   and generated-file hashes, and a small deterministic checker. The theorem
   statement must say exactly \(2\le h\le60\).
8. **Formalization ledger.** List Lean-checked lemmas individually, together
   with every imported recurrence or private-helper dependency. Do not use
   “formally verified” for the end-to-end gcd theorem.
9. **Immutable submission record.** Leave the submitted 136-page ZIP unchanged
   at `SUBMIT/dist/ramanujan-3.2-huang.zip`.  Its SHA-256 checksum is
   `9923f47c4614ab15c2c1f7320bfceed13a66576bbdb82e5d819c550ea00e7270`.
   Build the revised article and companions from new source files so that the
   historical submission and the corrected publication package cannot be
   confused.

Before circulation, a final claim audit should search the manuscript for every
use of “prove,” “equivalent,” “all heights,” “formally verified,” and “main
theorem,” and reconcile each occurrence with the status ledger.

## Immediate post-submission correction note

The following note should be sent promptly if the submitted title, abstract,
or cover letter can be read as claiming a solution of Problem 3.2.

> **Correction concerning the scope of the submitted archive.** A
> post-submission audit has shown that the 136-page archive does not prove the
> conjectural pointwise estimate
> \(\log\gcd(d_na_n,d_nb_n)=o(n)\) for every \(n\). That problem remains open.
> The archive should therefore be regarded as a partial-results dossier that
> records several exact reductions, structural results, computational
> certificates, conditional gateways, and unsuccessful approaches; it is not
> the manuscript we intend to publish in its present form. The strongest
> unconditional result relevant to the gcd problem is the quantitative
> exceptional-set theorem
> \[
> \#\{n\le X:\log G_n>\varepsilon n\}
> =O_\varepsilon((\log X)^2)
> \]
> for every fixed \(\varepsilon>0\). We are preparing a new, substantially
> shorter article centered on that theorem, the corrected master-sum
> reduction, and the continuant codegree argument. Independent critical-value
> and analytic results will be separated into companion notes and labelled by
> their precise proof status. We ask that no claim in the submitted archive be
> interpreted as an announced solution of the original every-\(n\) problem.

This correction changes the advertised scope, not the historical submission:
the submitted ZIP remains frozen as the record of what was sent. The future
paper should be written anew according to the architecture above.
