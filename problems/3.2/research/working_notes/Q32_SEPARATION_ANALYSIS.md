# Q3.2 quotient-separation audit: Casoratian and defect-dimension routes

**Date:** 2026-07-29

**Source under audit:** `/tmp/q_separation_sol.txt` and the scripts in the same
temporary scratch directory.

**Goal:** assess whether the proposed CRT quotient-separation mechanism can prove
the fully pointwise estimate `log G_n=o(n)`, with emphasis on Route A
(Casoratian) and Route B (defect-dimension counting).

## 1. Executive verdict

The proposed separation mechanism does **not** currently give a route to the
pointwise theorem.

1. The Casoratian in the source note has the wrong sign:

   \[
   a_n b_{n+1}-a_{n+1}b_n=-\frac{6}{(n+1)^3},
   \]

   with no alternating factor. The generating script prints failures at every
   odd index but continues and later tests only even outer indices, thereby
   masking the error.

2. After the sign is repaired and denominators are treated correctly, the
   Casoratian gives an exact formula for the local quotient `b_n/p`. It does
   not give an independent condition on that quotient. Replacing `b_n` by
   `R_nX` produces a defect which factors **exactly** as a fixed integer times
   `R_n(A_n-X)`. If `X≡A_n (mod R_n^k)`, the cleared defect contains
   `R_n^{k+1}` for the tautological reason that `A_n-X` already contains
   `R_n^k`. There is no compression.

3. The cited Liu 2024 theorem is materially different from the claimed
   endpoint tower. Ji-Cai Liu, arXiv:2404.16636, Theorem 1.1 proves one extra
   Gauss digit (here `\mathbf A` is the generalized Apéry sequence, not the
   quotient defined below):

   \[
   \mathbf A_{np^m}^{(r,s,t)}
   \equiv \mathbf A_{np^{m-1}}^{(r,s,t)}
      +p^{3m}B_{p-3}\,\mathcal A_n^{(r,s,t)}
      \pmod {p^{3m+1}}.
   \]

   It concerns indices `np^m`, adds only `B_{p-3}`, and for `m=1` is a
   congruence modulo `p^4`, not `p^5`. It does not state an arbitrary-order
   shifted endpoint expansion for `n=p+r`, and it does not supply the asserted
   sequence `B_{p-3},B_{p-5},...`.

4. Sun's relevant statements do not supply the missing tower either.
   Conjecture 5.1 of arXiv:2005.02081 is only modulo `p^4` for a few fixed
   multiples. Conjecture 2.4 of arXiv:2409.06544 predicts a modulo-`p^5`
   formula at multiplicative indices, but it is unproved and contains one
   fixed combination of `B_{2p-4}` and `B_{p-3}`, not the claimed sequence
   `B_{p-3},B_{p-5},...`. Thus those papers do not define the asserted
   dimensions at `p^5,p^7,p^9`. The temporary dimension script is a
   speculative free-polynomial model, not a computation of an Apéry defect
   module, and it is internally inconsistent even as a speculative model.
   A later independent block calculation, synchronized and audited in
   Section 11 below, does prove concrete target carriers through `p^7`;
   this is new input, not support for the proposed dimension table.

5. The inverse-limit formulation is backwards. If each finite congruence
   really computes the residue of the fixed integer `A_n` modulo `R_n^k`,
   then the compatible inverse-limit element is **automatically** the diagonal
   image of `A_n`, and its least nonnegative representatives eventually
   stabilize at `A_n`. Proving non-diagonality would contradict the congruence
   input, not prove the desired estimate.

6. A fixed level such as `k=7` could at best improve a linear constant. To
   prove `log R_n=o(n)` by reconstruction, one needs a separation theorem
   uniformly at arbitrarily large fixed levels `k` (or a controlled growing
   level), plus a genuinely independent nonzero low-height certificate.

7. The later small-zero-fiber computation does not change this verdict.
   Even \(|Z_p|\le2\) for every prime would not bound one fixed integer
   column: reflection-invariant two-point rows can all be aligned at the
   same \(n\). The exact scan through two million instead supports a
   Poisson law with slowly unbounded row and column maxima. The claimed
   evenness also has the central exceptions \(p=11,3137\).

8. The clean surviving target is horizontal. It is enough to prove, for
   every fixed quotient \(q\), a power saving for the shell pair energy of
   the affine zero sets \(qp+Z_p\). This is a cross-prime small-CRT theorem,
   not a consequence of a local Casoratian or of defect-space dimension.

9. There is nevertheless one unconditional pointwise pruning. A
   divisor-sensitive Kummer-order sieve shows that all middle primes with
   \[
   \gcd(p-1,n-q)>
   \tau(n-q)\log n\,F(n),\qquad q=\lfloor n/p\rfloor,
   \]
   have total logarithmic weight \(o(n)\) for any fixed subpolynomial choice
   \(F(n)\to\infty\), even before the Apéry zero condition is imposed. Thus
   the unresolved primes select characters of order \(p^{1-o(1)}\). This
   sharpens the target, but does not supply the missing cross-prime
   nonalignment.

10. The safe large-prime depth bound is \(v_p(G_n)\le6\), not 3. The
    claimed bound 3 in `proof.tex` uses \(a_{n-1}\) as an integer Bézout
    coefficient although it is generally rational. Clearing both adjacent
    denominators gives the valid constant 6. This changes no support-to-gcd
    little-oh implication.

11. The entire individual-gap factorial-moment variant of Route A can be
    optimized explicitly. In the direct prefix window it gives
    \[
    \sum_{z_p(H)>B}\log p\ll H^3\log H/B^3.
    \]
    At the available \(B\asymp H^{2/3}\) scale this is only
    \(O(H\log H)\). Higher moments are worse; the missing improvement is a
    joint short-gap compression identity.

12. A sharper top-half reduction is available. A prime-gap collision is
    **degenerate** if either characteristic sees zeros at both endpoint
    indices; the continuant then proves the pointwise bound
    \(M_h^{\rm deg}(n)\ll h\). On every remaining **pure-cross** collision,
    the same continuant is a unit modulo both primes. Consequently the full
    pair-energy theorem can be weakened to one logarithmic-gap statement:
    adjacent pure-cross target pairs of gap \(O(\log n)\) must be
    \(o(n/\log n)\). Replacing the index gap \(h\) by \(qh\) extends the
    argument to every fixed quotient \(q\); after the Kummer pruning,
    (PC.21) is one minimal all-middle obstruction.

13. The CFVZ square-root structure does not yield a hidden low-order
    Casoratian. The reversed convolution \(b_r=\sum_i s_i s_{r-i}\) admits
    no uniformly bounded rational boundary telescoper: the natural
    two-term certificate fails by a translation-by-two divisor imbalance,
    while a general bounded-width collapse would contradict the
    \({\rm SL}_2\) differential Galois group. Its minimal outer telescoper
    is exactly the original symmetric-square Apéry operator.

14. The divisor-sensitive order parametrization has an exact \(L^2\)
    consequence, but it does not close the horizontal gap. Low--low
    collision energy is negligible through order \(N^{1/2-o(1)}\), while
    collisions with one unrestricted partner are controlled only through
    \(N^{1/3-o(1)}\). The pointwise sieve should instead be applied first:
    it reduces the sufficient pair-energy theorem to two nearly primitive
    characters of order \(N^{1-o(1)}\). No Cauchy--Schwarz interpolation
    bridges these scales.

15. This nearly primitive restriction does not make the known local
    structure horizontally selective. Taking centers with \(n-q\) prime
    gives an exact reflected two-point adversary with
    \(\gcd(p-1,n-q)=1\) for every selected \(q\)-arc prime. It also obeys
    nonconsecutivity and the universal reflected-continuant factor. Thus
    any proof of the residual pair energy must use arithmetic of the
    distinguished Apéry initial state, not these local axioms alone.

16. In the exact scan through two million, all 2,764 unordered top-half
    target pairs are pure cross and none is degenerate. Only 16 of the
    2,711 adjacent edges have gap at most \(10\log n\). This is finite
    evidence, but it shows that the pure-cross logarithmic-gap criterion
    isolates the observed residual without leakage.

17. Gessel--Lucas collapses every top-half event to \(p\mid b_n\), and
    the central binomial carrier gives an exact holonomic-gcd formulation.
    More generally, it suffices to prove
    \(\log\gcd(b_n,\binom n{\lfloor n/d\rfloor})=o_d(n)\) for every
    fixed prime \(d\): a fixed quotient \(q\) is captured by any
    \(d\mid q+1\). This is a useful Path D interface, not yet a proof: checked
    moving-target gcd theorems cover \(S\)-unit or constant-coefficient
    recurrence orbits, not an Apéry \(P\)-recursive coefficient paired
    with a factorial ratio.

18. A forbidden spike has a second deterministic amplification. After
    kernel partitioning it forces \(\Omega(n)\) pure-cross target pairs
    with a common kernel and gap
    \(O(\tau(n-q)\log^2n)=n^{o(1)}\). This yields the structured
    same-kernel criterion (SK.5), complementary to the adjacent-gap
    criterion (PC.21).

19. There is now a genuine positive local extension of the quotient tower.
    Endpoint block identities give the common cofactor modulo
    `R,R^2,...,R^7`; at the seventh grade the apparent new scalar satisfies
    an unconditional quadratic H6 relation and can be removed by a fixed
    nonlinear direct/reflected combination. This repairs the former
    order-seven barrier. It still does not prove the pointwise theorem:
    one fixed level changes only a linear reconstruction threshold, while
    even an arbitrary compatible tower needs an independent small-lift or
    non-stabilization theorem to constrain the diagonal integer cofactor.

20. The direct/reflected Casoratian extends unconditionally through the
    eighth quotient digit.  With the raw endpoint residual
    \[
      {\cal F}_p=b_{2p-1}-5-8\Delta_p-\frac{336}{5}H_p,
    \]
    the fixed law is
    \[
      (166144+33296H_p+{\cal F}_p)D_8+5Z_8
      \equiv166149\,b_{p+r}/p\pmod {p^8}.
    \]
    It is valid at \(p=769\); the former \(769\) exception came only from
    splitting \({\cal F}_p\) into two normalized coordinates.  All \(163\)
    target rows through \(p\le1000\), including the two rows at \(769\),
    pass the exact audit.  The sole fixed inversion exception is \(18461\).
    This gains one local digit but remains CRT-saturated.

21. The proposed four-coordinate effective-weight-seven finite-MHS module
    collapses rigorously to one line:
    \[
      3\eta=2\xi,\qquad3A=14\xi,\qquad2B=-7\xi.
    \]
    Thus the next local obstruction is not defect dimension.  A precise
    all-\(m\) endpoint rank-one formula at precision eight passes \(6520\)
    exact divisibility checks, but its symbolic proof still needs the
    lifted lower-coordinate change and a Gosper certificate.  Even that
    proof would not supply horizontal separation, because all fixed-order
    local quotient rows reconstruct the same diagonal integer.

22. Path D can be stated exactly as a lower-denominator problem for an
    explicit Hadamard \(G\)-function.  Every fixed residue section already
    has asymptotically maximal prefix denominator lcm:
    \[
      \log\operatorname {lcm}_{m\le M}q_{dm+a,d}=dM+o_d(M).
    \]
    This does not control one coefficient.  In the central case a target is
    precisely an isolated denominator hole with valuation pattern
    \((1,0,1)\), and the integer
    \[
      \mathfrak d_n=
      \frac{\gcd(q_{n-1},q_{n+1})}
           {\gcd(q_{n-1},q_n,q_{n+1})}
    \]
    has exact top-half radical \(\prod_{p\in T_n}p\).

23. The hole carrier is selective but not yet a height compression.  It
    satisfies
    \[
      \mathfrak d_n\mid2\gcd\!\left(
        b_n,\binom n{\lfloor n/2\rfloor}\right),
    \]
    with no factor \(2\) for even \(n\).  Moreover the exact normalized
    primitive normalized recurrence has unit polar coefficients at every
    top-half candidate and admits a rational solution with an arbitrary
    prescribed squarefree product of simultaneous \((1,0,1)\) holes.  Its
    exact primitive state-dependent carrier is divisible by
    \(\mathfrak d_n\), but has logarithmic height at least
    \[
      n\log\!\left(\frac{17+12\sqrt2}{2}\right)-O(\log n).
    \]
    Consequently no initial-state-uniform theorem based only on the
    recurrence, generic \(G\)-function axioms, or local denominator
    propagation can prove \(\log\mathfrak d_n=o(n)\).  A theorem for the one
    distinguished Apéry initial state is not excluded.

24. Higher universal tensors satisfy an exact conservation law.  On a
    pure-cross pair, if \({\cal C}_{p,\ell}\) is the formal cross-divisibility
    ideal and \(\Delta\) is the endpoint determinant, then
    \[
      ({\cal C}_{p,\ell}:\Delta^m)={\cal C}_{p,\ell}.
    \]
    Thus primitive selectivity leaves an exponential slot, while complete
    contraction leaves coefficient content \(p\ell\).  The only fixed-order
    escape is an initial-state-specific leading cancellation.  Its leading
    symbol would give a nonzero polynomial over
    \(\mathbb Q(\sqrt2)\) vanishing at \(\zeta(3)\); growing templates require
    a corresponding quantitative algebraic-approximation measure.  Neither
    input is currently known.

25. There is one unconditional pointwise gain on the horizontal side.
    For each fixed even gap \(h\le A\log n\), the classical
    dimension-two Selberg sieve gives
    \[
      M_h^\times(n)
      \ll {\mathfrak S}(h)\frac n{\log^2n},
      \qquad
      {\mathfrak S}(h)\ll\log\log(3h).
    \]
    Hence the product of all ambient prime pairs \(p(p+h)\) at that exact
    gap is an \(o(n)\)-height integer containing every pure-cross target
    pair of gap \(h\).  The mean law
    \(\sum_{h\le H}{\mathfrak S}(h)\ll H\) shows the precise stopping
    point: multiplying these carriers over \(h\le A\log n\) costs
    \(O_A(n)\), not \(o(n)\).  Any Apéry-specific saving factor tending to
    infinity over this ambient sieve scale would complete (PC.21).

26. The later \(S(p)\)/Katz handoff does not close the problem.  Its
    multiplier set describes the companion \(a_q\) channel, which is
    already pointwise \(O(n^{2/3})\), while omitting the lower \(b_r\)
    channel.  The exact counterexample is \(n=16,p=11,r=5\).  Moreover a
    bounded row set cannot control a fixed column, and the variable-length
    recurrence prefix product is not a Frobenius statistic of the
    Beukers--Peters K3 sheaf.  Ordinary square-root character sums would
    give \(1+O(\sqrt p)\) zeros, not \(O(1)\).

27. Universal direct carriers are now saturated in every degree.  The
    intersection of the moving target ideals evaluates exactly to the
    unknown radical \(R_n\mathbb Z\), even after imposing the full
    recurrence.  Signed Newton carriers cannot cancel the central \(4^k\)
    alias at subexponential cost: exact cancellation forces full candidate
    primorial content, and all-degree alias saturation changes \(R_n\) by
    at most a fixed factor \(24\).

28. Neither the standard irrationality proof nor a factorial-gcd
    reformulation supplies the missing estimate.  Van der Poorten,
    Beukers, and Zudilin use unreduced lcm-cleared linear forms, integrals,
    or hypergeometric denominator estimates and do not prove a
    subexponential common divisor.  The full \(\gcd(b_n,n!)\) statement is
    stronger because Lucas congruences do not control small-prime
    valuation depth.

29. The first \(p\)-adic lift gives no second target equation.  With the
    harmonic derivative \(D_r\),
    \[
      b_{p+r}\equiv5(b_r+pD_r),\qquad
      p^3a_{p+r}\equiv6(b_r+pD_r)\pmod {p^2}.
    \]
    On the target locus the two divided coordinates remain proportional to
    \((5,6)\).  The example \(p=73,r=2\) has valuation exactly one, so the
    residual scalar is not forced to vanish.

30. The two sharp positive interfaces are now explicit.  One is the
    shell-uniform fixed-\(q\) pair energy
    \(o_q(N^2/\log^2N)\); the other is a pointwise
    two-characteristic dispersion estimate over prime gaps
    \(h\le A\log N\).  Standard \(L^2\) large sieves are critical even for
    bounded rows, and reflected two-point masks have optimal individual
    Fourier bounds while still violating the pointwise estimate.  A proof
    must therefore introduce actual-Apéry cross-prime cancellation, not a
    stronger one-characteristic estimate.

31. The two exceptional targets detected by the recurrence factor
    \(2n+1\) have an exact modular explanation.  For \(n=p+r\),
    \[
      p\mid2n+1\quad\Longleftrightarrow\quad
      r=(p-1)/2,\qquad n=(3p-1)/2.
    \]
    Beukers' congruence identifies this central zero with
    \(p\mid a_p(\eta(2z)^4\eta(4z)^4)\).  This classifies the factor-covered
    targets, but it does **not** obstruct a pointwise bound
    \(z_p=O(1)\): such a bound may allow the one optional central zero.
    The modularity observation rules out eventual central nonvanishing,
    not every pointwise zero-fibre estimate.

32. Lifting the full \(2\times2\) transfer frame to \(p^2\) does not by
    itself add a selective equation.  The integral crossing block at
    \(p-1\) has rank one modulo \(p\), the two divided endpoint columns are
    proportional to \((5,6)\), and every nonsingular neighboring row
    propagates the same free scalar.  A local frame deformation preserves
    the recurrence, determinant, and zeroth layer while varying that
    scalar.  Since it need not preserve the distinguished initial frame,
    this is a no-go for the proposed **local** matrix argument, not a
    theorem excluding new global connection data.

33. The latest uisai2 `dm` and `family` computations sharpen both sides of
    the saturation ledger.  A half-prefix has a universal third layer
    \(p^3\) on \(n/2<p\le(3n+1)/4\), but the divided carrier still has
    positive linear height and an actual target has no fourth layer.
    Two Laurent models have cross-model Smith form
    \({\rm diag}(1,U_n)\); their primitive mutation commutator is not
    target-divisible.  The second selected Kummer-trace digit is a free
    Witt coordinate.  These are genuine structural refinements, but none
    supplies horizontal separation.

34. There is a new experimental Route A reformulation.  Writing the
    rational companion in lowest terms as \(a_n=A_n/C_n\), exact
    computation through \(n=10000\) found
    \(\gcd(A_n,b_n)=1\), which would imply
    \[
      G_n=\frac{d_n^3}{C_n}.
    \]
    The two adjacent reduced Wronskians only give
    \[
      \gcd(A_n,b_n)\mid
      \gcd\!\left(
        b_n,\frac{6C_{n-1}C_n}{n^3},
        \frac{6C_nC_{n+1}}{(n+1)^3}\right).
    \]
    This carrier is usually nontrivial: it exceeds one for 1952 of
    \(1\le n<2000\), including a residual factor \(17\) at \(n=20\).
    Thus the coprimality pattern is promising finite evidence, not a
    consequence of two neighboring Casoratians.

35. The zero-fibre parity is now an unconditional two-line theorem:
    \[
      z_p\equiv{\bf1}_{p\mid a_p(f)}\pmod2.
    \]
    Reflection pairs every noncentral zero, and (21.2) identifies the
    fixed point.  The scan through \(p\le200000\) found mean \(z_p=1.00801\),
    maximum \(12\), and only the two odd fibres \(11,3137\).  These are
    strong data for (20.1), not a proof.

36. The strongest present unconditional prime average obtained from the
    continuant theorem is
    \[
      \sum_{p\le x}z_p\ll x^{5/3}/\log x.
    \]
    Averaging root counts of each fixed gap polynomial can improve this
    exponent only with uniformity for a growing family; even \(O(1)\)
    average roots for every gap through \(p^{1/2}\) yields only the
    \(x^{3/2+o(1)}/\log x\) scale.  The conjectural bounded mean requires
    compression of the whole fibre, not a refinement of the one-gap
    ledger.

37. The all-\(n\) ``L2 reduction'' is valid only after strengthening L2
    from the top-half assertion \(K_1(n)=o(n/\log n)\) to the family
    \[
      K_q(n)=o_q(n/\log n)
      \quad\text{for every fixed }q.
    \]
    The valuation cap makes the subsequent \(\eta\)-truncation uniform,
    but it does not derive the \(q\ge2\) hypotheses from \(q=1\).
    Consequently the quotientwise L2 family is a clean sufficient
    interface, while Fable's original L2 alone is insufficient.

38. The reduced denominator defect is itself an exact form of the hard
    channel.  If \(D_n=d_n^3\), \(P_n=D_na_n\), and \(a_n=A_n/C_n\) is
    reduced, then
    \[
      \frac{D_n}{C_n}=\gcd(D_n,P_n).
    \]
    For \(p>\sqrt n\), \(n=qp+r\),
    \[
      p\mid D_n/C_n\Longleftrightarrow p\mid a_qb_r,
    \]
    and in the top half this is equivalent to \(p\mid b_n\).  Thus even a
    universal proof of \(\gcd(A_n,b_n)=1\) would only convert P3.2 into the
    same unresolved pointwise denominator-defect estimate; it would not
    remove the cross-prime obstruction.

39. The neighboring reduced companion denominators give a sharper exact
    local dictionary.  For \(p\ge7\), \(n=p+r\), and \(p\mid b_n\), put
    \(t=v_p(b_n)\).  Then
    \[
      (v_p(C_{n-1}),v_p(C_n),v_p(C_{n+1}))
      =(3,\,3-\min(t,3),\,3).
    \]
    Hence the positive denominator curvature
    \[
      \kappa_n=
      \frac{\gcd(C_{n-1},C_{n+1})}
           {\gcd(C_{n-1},C_n,C_{n+1})}
    \]
    satisfies \(v_p(\kappa_n)=\min(t,3)\) throughout the top half.  This
    proves the proposed implication
    \(p\mid A_n,b_n\Rightarrow v_p(b_n)=3\), but also shows why it does not
    solve P3.2: \(\kappa_n\) exactly re-encodes the full top-half
    denominator defect.  In fact \(\kappa_n\mid b_n\) for every \(n\), and
    the `uisai2` denominator-valley argument proves
    \[
      \log\kappa_n=o(n)
      \quad\Longleftrightarrow\quad
      \log\operatorname {rad}\gcd(b_n,d_n)=o(n),
    \]
    hence equivalently P3.2 after the already controlled channels.  This
    is the cleanest target-free reformulation, not yet a height bound.

40. There is a useful global common-carrier equivalence.  Let
    \[
      B(n)=\{p>\sqrt n:p\mid b_{n\bmod p}\},\qquad
      C(n)=\{p>\sqrt n:p\mid b_n\}.
    \]
    Lucas gives \(B(n)\subset C(n)\).  The logarithmic weight of
    \(C(n)\setminus B(n)\) is
    \[
      O(Q^2+n/Q)=O(n^{2/3})
      \quad(Q=n^{1/3}),
    \]
    because a false positive with \(\lfloor n/p\rfloor=q\le Q\) divides
    \(b_q\), while one with \(q>Q\) has \(p\le n/Q\).  Thus the remaining
    channel is equivalent, up to an unconditional \(o(n)\) error, to
    \[
      \log\operatorname {rad}_{\sqrt n<p\le n,\ p\mid b_n}b_n=o(n).
    \]
    This removes the moving remainder from the statement, but current
    prime-factor theorems for Apéry or \(G\)-function coefficients do not
    prove this all-\(n\) large-prime-radical estimate.

41. Delaygue's valuation theorem does not control the new multiplicity
    issue.  It proves
    \(v_p(b_n)\ge\alpha_p(b,n)\), where \(\alpha_p\) counts base-\(p\)
    digits in \(Z_p\).  At \(n=p+r\), \(p\ge7\), this gives only
    \(v_p(b_{p+r})\ge1\) when \(p\mid b_r\), exactly the Lucas
    information.  It supplies no upper bound and cannot exclude cubic
    divisibility.  An independent exact scan through \(p\le50000\) found
    5165 shifted targets of valuation one, five of valuation two, and
    none of valuation at least three.  Even a theorem excluding every
    cubic event would control only \(\gcd(A_n,b_n)\), not the simple
    targets which dominate the denominator defect.

42. The Hasse--Witt Mellin normalization is now exact.  If
    \(\Phi=(1+x)(1+y)(1+z)(1+1/(xyz))\), then
    \(b_m=\operatorname {CT}\Phi^m\) and
    \(F_{<p}(u)=\sum_{m<p}b_mu^m\) is the first Hasse--Witt scalar for
    \(1-u\Phi=0\).  In the coordinate \(\Phi=t\), one has \(u=t^{-1}\),
    up to the monomial unit determined by the Hodge-line
    trivialization.  Exact zero detection of a Mellin coefficient opens
    into a product over all \(p-1\) fibers, or equivalently a moment of
    order \(p-1\); it does not reduce to an ordinary two-fiber
    Katz--Deligne correlation.  This closes that naive fixed-order
    expansion, not every possible geometric argument.

43. Plain binomial-window carriers have an exact top-band obstruction.
    For \(p>n/2\),
    \[
      p\mid\binom nk\quad\Longleftrightarrow\quad
      p>\max(k,n-k).
    \]
    Therefore every prime \(p>n/2+H\) divides every coefficient in
    \(|k-n/2|\le H\), and
    \[
      \log\gcd_{|k-n/2|\le H}\binom nk
      \ge\vartheta(n)-\vartheta(n/2+H)
      =n/2-H+o(n).
    \]
    Making this carrier subexponential forces
    \(H=n/2-o(n)\), at which point it retains only the vanishing tail
    \(p>n-o(n)\).  Thus the whole one-parameter window family has a rigid
    height-versus-coverage tradeoff.

44. The proved row bound already gives an unconditional density-one
    theorem.  On \(N<n\le2N\), the total logarithmic channel-\(B\) mass is
    \(O(N^{5/3})\).  Hence for every \(2/3<\alpha<1\), all but
    \(O(N^{5/3-\alpha})\) integers in the shell have mass at most
    \(N^\alpha=o(N)\); together with the pointwise small-prime and
    companion bounds this proves P3.2 outside a density-zero set.  It
    does not remove the sparse exceptional columns.  Direct computation
    of the top-half falling second moment through \(N=128000\) is close
    to the reflected-Poisson prediction \(N\lambda_N^2\), many orders
    below the sufficient threshold
    \(o(N^2/\log^2N)\), but remains finite evidence.

Thus Routes A and B are closed **as presently formulated**. This is not a
disproof of `log G_n=o(n)`. It identifies the exact missing ingredient:
cross-prime, Apéry-specific Archimedean/p-adic coupling not already contained
in the local residue statement.

## 2. Scope and notation

Let

\[
b_n=\sum_{j=0}^n\binom nj^2\binom{n+j}j^2,
\qquad (b_0,b_1)=(1,5),
\]

and let the rational companion solution have `(a_0,a_1)=(0,6)`. Both satisfy

\[
(n+1)^3u_{n+1}=P(n)u_n-n^3u_{n-1},
\quad
P(n)=34n^3+51n^2+27n+5.
\]

Following the source note, let

\[
R_n=\prod_{\substack{\sqrt n<p\le n\\p\mid b_n}}p,
\qquad A_n=b_n/R_n.
\]

Let \(d_n=\operatorname{lcm}(1,\ldots,n)^3\) and
\(G_n=\gcd(d_na_n,d_nb_n)\). Since
\(d_na_{n-1},d_nb_{n-1}\in\mathbf Z\), the valid integral combination is

\[
 (d_nb_{n-1})(d_na_n)-(d_na_{n-1})(d_nb_n)
 =\frac{6d_n^2}{n^3}.
\]

Consequently, for \(p\ge5\),

\[
 v_p(G_n)\le6\lfloor\log_p n\rfloor-3v_p(n)
 \le6\lfloor\log_p n\rfloor.
\]

In particular \(v_p(G_n)\le6\) for \(p>\sqrt n\); the total contribution of
\(p\le\sqrt n\) is \(O(\sqrt n)\). Hence an \(o(n)\) logarithmic bound for
the relevant large-prime radical is enough for \(\log G_n=o(n)\). The issue
is identifying which factors of the stronger proxy \(R_n\) are actually
relevant to \(G_n\).

There is a scope caveat. For the full interval `sqrt(n)<p<=n`, the block
criterion proved in the repository is

\[
p\mid G_n\quad\Longleftrightarrow\quad D_qb_r\equiv0\pmod p,
\]

where `D_q` is the companion-block constant. The repository's stronger
identification `D_q≡a_q (mod p)` gives the two-channel form

\[
p\mid G_n\quad\Longleftrightarrow\quad a_qb_r\equiv0\pmod p.
\]

The resume note records that the Wronskian argument alone does not prove
this identification at Hasse-zero indices; it points instead to the audited
explicit companion-sum/Kummer proof of the full block congruence. This
foundation issue does not affect the top-half calculation below, where
`q=1` and `D_1=a_1=6`.

The one-channel equivalence with `p|b_r` is valid for `p>n/2`, where `q=1`
and `a_1=6`. The companion `a_q` channel is already controllable with
`O(n^{2/3})` logarithmic weight, so a sufficiently uniform direct-channel
bound would still close the main problem. Nevertheless, `R_n` as defined
above is not literally the middle-prime part of `G_n`. The source note
silently switches between:

- all middle primes `sqrt(n)<p<=n`;
- only top-half primes `n/2<p<=n`;
- the lower-digit support relevant to `G_n`.

The accompanying `crt_tower_digits.py` claims the first definition in its
docstring but implements only the second.

## 3. Reconstruction threshold: the correct calculation

The dominant root of the Apéry recurrence is

\[
\lambda=17+12\sqrt2=(1+\sqrt2)^4,
\qquad c=\log\lambda=3.52549\ldots,
\]

and

\[
\log b_n=cn-\frac32\log n+O(1).
\]

Put `rho_n=(log R_n)/n`. Since `A_n=b_n/R_n`,

\[
\log A_n=(c-\rho_n)n+o(n),
\]

not `cn+o(n)` as stated in `/tmp/q_separation_sol.txt`.

If a level-`k` congruence determines `A_n mod R_n^k`, its least representative
is forced to equal `A_n` once

\[
A_n<R_n^k.
\]

The correct threshold is therefore

\[
cn+o(n)<(k+1)\log R_n,
\qquad\text{or}\qquad
\rho_n>\frac{c}{k+1}+o(1).
\]

The extra `+1` comes from the factor of `R_n` already removed from `b_n`.
This explains the numerical comment in the source note:

\[
c-\frac{7}{2}=0.02549\ldots.
\]

For a top-half radical, the prime number theorem gives the ceiling
`log R_n<=theta(n)-theta(n/2)=n/2+o(n)`. Thus level `k=6` misses reconstruction by
`0.02549n+o(n)`, while level `k=7` is the first level that could improve
the `1/2` constant. The source note's displayed threshold
`cn<k log R_n` is inconsistent with its own statement that `k=7` is first.

More importantly, ruling out

\[
\rho_n>\frac{c}{k+1}+o(1)
\]

for one fixed `k` gives only another positive linear constant. A
little-oh theorem requires the argument for arbitrarily large fixed `k`:
given every `epsilon>0`, choose `k>c/epsilon`, prove the relevant separation
uniformly in the moving `n,p,r` at that level, and only then let `n` grow.
No rate uniform as `k` itself tends to infinity is needed.

## 4. Route A: exact Casoratian audit

### 4.1 Correct Casoratian

Set

\[
W_n=a_nb_{n+1}-a_{n+1}b_n.
\]

Substitution of the recurrence gives

\[
\begin{aligned}
W_n
&=\frac{n^3}{(n+1)^3}
  (a_{n-1}b_n-a_nb_{n-1})\\
&=\frac{n^3}{(n+1)^3}W_{n-1}.
\end{aligned}
\]

Since `W_0=-6`,

\[
\boxed{W_n=-\frac6{(n+1)^3}.}
\]

Equivalently, in the orientation used in the repository,

\[
a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}.
\]

The temporary script `apery_casoratian.py` makes the sign error exactly
between its lines 47 and 48: line 47 has the expression equal to
`W_{n-1}`, but line 48 inserts an extra minus sign. When executed, it prints

```text
n=1: W_n/W_{n-1} = 1/8, expected -1/8 [FAIL]
n=2: W_n/W_{n-1} = 8/27, expected -8/27 [FAIL]
...
```

and its proposed closed form fails at every odd index. Its later assertions
use only `n=200,300,500`, all even, so the error happens not to trigger an
exception.

### 4.2 Exact local `p`-adic content

Fix a prime

\[
p\in(n/2,n],\qquad p\ge7,\qquad p\mid b_n,
\]

and put `r=n-p`. Then `0<=r<=p-2`: if `r=p-1`, Lucas and reflection give
`b_n≡5b_{p-1}≡5 (mod p)`, contrary to `p|b_n`.

Lucas gives

\[
b_n\equiv5b_r\pmod p,
\]

so `p|b_r`. Two consecutive Apéry values below `p` cannot both vanish,
hence

\[
p\nmid b_{r+1},
\qquad
p\nmid b_{n+1}.
\]

Let

\[
R_n=pS,\qquad
B=b_n=pSA_n,\qquad
C=b_{n+1},
\]

and define the `p`-integral companion values

\[
U=p^3a_n,\qquad V=p^3a_{n+1}.
\]

The block congruence in the audited repository gives

\[
U\equiv6b_r\equiv0\pmod p,
\qquad
V\equiv6b_{r+1}\not\equiv0\pmod p.
\]

Thus `u:=U/p` is `p`-integral and `V,C,S` are `p`-adic units. Multiplying
the correct Casoratian by `p^3` yields

\[
UC-VpSA_n=-\frac{6p^3}{(n+1)^3}.
\]

After division by `p`,

\[
\boxed{
VS A_n
=uC+\frac{6p^2}{(n+1)^3}
\quad\text{in }\mathbf Z_p.
}
\tag{A.1}
\]

Consequently,

\[
\boxed{
A_n
=(VS)^{-1}
\left(uC+\frac{6p^2}{(n+1)^3}\right).
}
\tag{A.2}
\]

Modulo `p`, this reduces to

\[
A_n\equiv (VS)^{-1}uC\pmod p.
\]

This is the exact local content. It expresses the quotient digit in terms of
the normalized companion digit `u=(p^3a_n)/p`; it does not eliminate that
new digit or constrain `A_n` independently.

The formula remains valid when `v_p(b_n)>=2`. In that case `A_n` itself
may still be divisible by `p`, because `R_n` removes only one copy of each
prime.

### 4.3 Exact numerical checks

Independent exact-rational recurrence arithmetic gives:

| `n` | `p` | `S=R_n/p` | `u mod p` | `V mod p` | `C mod p` | `A_n mod p` |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 11 | 5 | 9 | 6 | 5 | 7 |
| 20 | 17 | 1 | 6 | 7 | 3 | 5 |
| 27 | 19 | 11 | 8 | 17 | 11 | 15 |

In each case `(VS)^{-1}uC` gives the last column. For example, at
`n=16,p=11`, the full middle radical is `R_16=5*11`, so

\[
A_{16}\equiv(6\cdot5)^{-1}\cdot9\cdot5\equiv7\pmod {11}.
\]

These checks also give `v_p(U)>=1` and `v_p(V)=0`, exactly as the derivation
requires.

### 4.4 What happens to a fake quotient

Suppose `X` is any integer with

\[
X\equiv A_n\pmod {R_n^k}.
\]

Equation (A.1) gives

\[
uC-VSX+\frac{6p^2}{(n+1)^3}
=VS(A_n-X),
\]

so the left side has `p`-adic valuation at least `k`. In the un-divided
scaled Casoratian, the valuation is at least `k+1`. This is not a new
constraint: it is exactly the original congruence multiplied by a unit.

There is an equally sharp global formulation. Let

\[
L_{n+1}=\operatorname{lcm}(1,\ldots,n+1)^3
\]

and replace only `b_n` by `R_nX` in the Casoratian. The cleared defect is

\[
\begin{aligned}
E_n(X)
&=L_{n+1}\left(
 a_nb_{n+1}-a_{n+1}R_nX+\frac6{(n+1)^3}
\right)\\
&=(L_{n+1}a_{n+1})R_n(A_n-X).
\end{aligned}
\tag{A.3}
\]

All quantities in the last line are integers. If
`A_n-X=R_n^kt`, then

\[
\boxed{
E_n(X)=(L_{n+1}a_{n+1})R_n^{k+1}t.
}
\tag{A.4}
\]

Therefore:

- if `X=A_n`, the defect is exactly zero;
- if `X!=A_n`, then `|E_n(X)|>=R_n^{k+1}`.

The divisibility and the Archimedean lower bound saturate at the same power.
Indeed, in the nonzero case,

\[
\log|E_n(X)|
=(k+1)\log R_n+\log|t|+\log|L_{n+1}a_{n+1}|
\ge (k+1)\log R_n.
\]

The standard denominator and recurrence-growth bounds give
`\log|L_{n+1}a_{n+1}|=O(n)`; they cannot reverse this inequality. The strict
inequality needed for a contradiction,

\[
0<|E_n(X)|<R_n^{k+1},
\]

is impossible from this construction.

### 4.5 Why recurrence rigidity and Padé uniqueness do not repair this

The integer `X` is a proposed value at one index, not a second solution of
the Apéry recurrence.

- Keeping `b_{n-1}` and `b_{n+1}` fixed while replacing `b_n` by `R_nX`
  normally breaks the recurrence.
- Requiring the modified value to satisfy the recurrence forces
  `R_nX=b_n`, hence `X=A_n`; this simply assumes the desired equality.
- Padé uniqueness concerns a whole solution or a polynomial approximant
  with prescribed degree and normalization. A single CRT representative is
  neither.
- Interlacing or real positivity does not constrain the `p`-adic quotient
  digit `(b_n/p) mod p^j`.

The repository already contains the structural versions of this no-go:

- `Q32_CODEX_RESUME_2026-07-23.md`, Section 13: consecutive Casoratian
  certificates are an invertible change of coordinates, with linear
  height;
- Section 15: every multi-shift Casoratian family has rank two;
- Section 143: a comparison section synchronized with a star of bad primes
  carries the full star radical in its determinant;
- Section 149: the rational horizontal endomorphism ring of the canonical
  recurrence is scalar.

Thus “more Casoratians” do not add independent rows. They either remain in
the same rank-two module or encode the target radical in their normalization.

There is one exact many-prime identity, but it does not by itself change
this verdict.  If the target primes at a fixed \(n\) are
\(p_1,\ldots,p_K\), put

\[
 R=\prod_i p_i,\qquad A=b_n/R,\qquad
 P(T)=\prod_i(T-p_i).
\]

Then

\[
 \frac{b_n}{p_i}
 =A\prod_{j\ne i}p_j
 \equiv(-1)^{K-1}A\,P'(p_i)\pmod {p_i}.
\tag{A.5}
\]

All factors \(P'(p_i)\) are \(p_i\)-units in the top-half interval.
Thus the normalized divided residues recover the same common coordinate:

\[
 (-1)^{K-1}\frac{b_n/p_i}{P'(p_i)}
 \equiv A\pmod {p_i}.
\tag{A.6}
\]

This coherence is exact and could be a useful interface for a genuinely
new many-prime jet identity.  It is not yet quotient separation.  For any
chosen distinct primes and any chosen integer \(A\), the integer
\(A\prod_i p_i\) satisfies (A.5) identically.  Moreover,
\(P'(p_i)\) depends on the entire target set, and its unreduced product
has the same linear-height scale one is trying to beat.  CRT reconstruction
of the normalized residues is simply reconstruction of \(A\bmod R\), the
original quotient problem.  Transferring \(b_n/p_i\bmod p_i\) to a lower
Apéry quotient through a modulo-\(p_i^2\) Lucas jet does not remove this
free common \(A\).  A new relation among those transferred jets, with a
nonzero sublinear-height elimination, is still required.

### 4.6 Route A verdict and reopen condition

**Verdict:** the ordinary Casoratian route is exact but non-compressive.
It reconstructs the same local quotient coordinate and gives no
non-coincidence theorem.

A genuine reopening would require a second global observable `C_n` with all
of the following properties:

1. every target prime contributes a power `p^k` to `C_n`;
2. that divisibility is not an algebraic multiple of
   `X-A_n (mod p^k)` and not a rank-two change of coordinates;
3. `C_n` is proved nonzero;
4. `log|C_n|=o(k n)` or, more precisely, its Archimedean height is strictly
   below the accumulated target-prime valuation;
5. the construction is uniform for arbitrarily large fixed `k`.

The exact cubic `P(n)`, the canonical low-height gauge, or the
hypergeometric/Dwork realization would have to supply the independence.
The determinant line alone cannot.

## 5. Route B: literature and dimension audit

### 5.1 What Liu 2024 actually proves

There is a notation collision to avoid: these papers write `A_j` for the
ordinary Apéry number that this note calls `b_j`. In this section,
`\mathbf A_j` denotes the papers' Apéry number; the quotient
`A_n=b_n/R_n` retains its earlier meaning.

The relevant primary source is:

- Ji-Cai Liu, *An extension of Gauss congruences for Apéry numbers*,
  [arXiv:2404.16636](https://arxiv.org/abs/2404.16636), Theorem 1.1.

For `p>=5`, positive `N,ell`, and generalized-Apéry parameters
`r>=2`, `s,t>=0`, it proves

\[
\mathbf A_{Np^\ell}^{(r,s,t)}
\equiv \mathbf A_{Np^{\ell-1}}^{(r,s,t)}
+p^{3\ell}B_{p-3}\mathcal A_N^{(r,s,t)}
\pmod {p^{3\ell+1}}.
\tag{B.1}
\]

The earlier conjectures cited by Liu do not say more. Zhi-Hong Sun,
*Congruences for two types of Apéry-like sequences*,
[arXiv:2005.02081](https://arxiv.org/abs/2005.02081), Conjecture 5.1, gives
only special cases `\mathbf A_p,\mathbf A_{2p},\mathbf A_{3p} (mod p^4)`,
each involving
`B_{p-3}`. Liu's theorem is the uniform generalized-Apéry extension of this
kind of statement.

A later paper of Sun,
[arXiv:2409.06544](https://arxiv.org/abs/2409.06544), both quotes (B.1) and
states the stronger Conjecture 2.4:

\[
\mathbf A_{Mp^\ell}-\mathbf A_{Mp^{\ell-1}}
\equiv
2C_M p^{3\ell}
\left(
\frac{B_{2p-4}}{2p-4}
-2\frac{B_{p-3}}{p-3}
\right)
\pmod {p^{3\ell+2}}.
\tag{B.2}
\]

Here `p>5` and `M,ell` are positive. At `ell=1`, (B.2) is indeed a
conjectural congruence modulo `p^5`. It still concerns multiplicative
indices, not `p+r`; all `M` share the same one displayed prime-local scalar
combination; and it gives neither an arbitrary-order tower nor terms
`B_{p-5},B_{p-7},...`.

The differences from the source note are decisive:

- (B.1) is at multiplicative indices `Np^ell`, not shifted indices `p+r`;
- it is one digit deeper than the order-three Gauss congruence;
- for `ell=1`, its modulus is `p^4`;
- its only Bernoulli coordinate is `B_{p-3}`;
- the correction is additive with an explicit coefficient
  `mathcal A_n`, not the asserted general multiplicative formula
  `b_{mp}=b_m(1+C_m Delta_p) (mod p^5)`;
- it provides no formulas involving `B_{p-5},B_{p-7},...`.

No source for the claimed all-order tower or for “Q2672” was present in the
canonical repository or the temporary research directory.

### 5.2 Honest status of the requested level table

From the cited theorem, the only defensible table is:

| level | unconditional information | conjectural information relevant here |
|---|---|---|
| `p^3` | Gauss congruence `\mathbf A_{Np}≡\mathbf A_N`; 0 new coordinates | none needed |
| `p^4` | Liu (B.1); at most 1 coordinate, `B_{p-3}` | Sun 2020 special cases |
| `p^5` | no formula from Liu; dimension unknown | Sun (B.2), at most one displayed scalar combination |
| `p^7` | no cited formula; dimension unknown | no tower supplied |
| `p^9` | no cited formula; dimension unknown | no tower supplied |

Therefore the requested unconditional dimensions at `p^5,p^7,p^9` have not
been computed. They first require actual congruence theorems defining a
defect module. Even assuming (B.2) would specify only one
`\mathbf Z/p^2\mathbf Z`-valued scalar combination—at most one new
`p`-adic digit beyond (B.1)—not the free Bernoulli algebra used by the
temporary script. Turning that combination into an `F_p`-dimension statement
would itself require its Kummer and higher-digit relations.

**Later update.** Section 11 records unconditional in-house endpoint and
target identities through \(p^7\), proved after this source audit.  They do
not retroactively validate the table proposed in the source note: the new
coordinates satisfy rank-one and nonlinear relations rather than forming
the free Bernoulli polynomial algebra assumed by
`defect_dimension_count.py`.

### 5.3 Internal errors in `defect_dimension_count.py`

The temporary script assumes, without derivation, that the defect algebra is
the free polynomial algebra on symbols

\[
X_j=B_{p-(2j+1)},\qquad \operatorname{wt}(X_j)=j.
\]

Even under this unproved model, its table is inconsistent.

The number of monomials of **exact** weight `w` is the partition number
`p(w)`. If all weights up to `j` are retained, the nonconstant cumulative
dimension is

\[
D(j)=\sum_{w=1}^j p(w).
\]

Thus the speculative counts would be

| putative level | maximum weight | exact new monomials | cumulative nonconstant monomials |
|---|---:|---:|---:|
| `p^3` | 0 | 0 | 0 |
| `p^5` | 1 | 1 | 1 |
| `p^7` | 2 | 2 | 3 |
| `p^9` | 3 | 3 | 6 |

The script instead:

- says in prose that the `p^3` dimension is zero but prints one;
- calls `p(w)` a cumulative dimension although it counts exact weight;
- obtains “new dimension” by subtracting consecutive partition numbers,
  which has no meaning for either grading;
- counts `m=0` as an equation although Liu's theorem assumes positive `n`,
  and the `n=0` identity would be tautological;
- doubles every row by a “reflected endpoint” without proving that the
  reflected row is independent;
- invokes “Bernoulli independence at different primes”, for which no such
  theorem is supplied or known in the required form.

The output is consequently not evidence for rank growth.

### 5.4 Why equations from varying `m` do not overdetermine a fixed `A_n`

There are four separate quantifier failures.

#### (i) Every interpretation of “varying `m`” changes the sequence index

Varying the exponent `ell` in (B.1) concerns

\[
\mathbf A_{Np},\mathbf A_{Np^2},\mathbf A_{Np^3},\ldots.
\]

The temporary dimension script instead varies the base multiplier
`M=0,1,2,...`. That produces values
`\mathbf A_{Mp^\ell}` at different outer indices (and `M=0` is
tautological). Neither variation supplies several equations for one fixed
Apéry value.

For the moving-prime problem, the outer index is fixed and a top-half bad
prime has

\[
n=p+r.
\]

Changing either the exponent or the base multiplier does not give additional
equations for the same integer `A_n=b_n/R_n`.

#### (ii) Each equation has new left-side data

If the values `\mathbf A_{Np^\ell}` or `\mathbf A_{Mp^\ell}` are treated as
unknowns, every new index adds a new unknown along with the equation. If they
are treated as known Apéry values, the equations are already consistent
identities. More rows can determine or verify `B_{p-3}`; “more equations
than Bernoulli symbols” does not imply a contradiction.

An overdetermined linear system is contradictory only after one proves that
its augmented right-hand side is outside the coefficient column space.
Dimension counting alone says nothing of the kind.

#### (iii) Different primes have different local coordinates

For each `p`, the residues

\[
B_{p-3},B_{p-5},\ldots\pmod p
\]

live in a different residue field and form different local data. For a set
of bad primes `p|R_n`, the natural defect space is a direct product of the
`p`-local spaces. The number of local unknowns therefore grows with the
number of primes at the same rate as the number of local congruence rows.

To turn them into shared variables one needs a proved cross-prime
interpolation, motive, or low-height rational relation. That is precisely
the missing horizontal theorem; it cannot be supplied by counting.

#### (iv) Reflection does not automatically double rank

Modulo `p`, the direct and reflected Apéry values satisfy

\[
b_{p-1-r}\equiv b_r\pmod p,
\]

so they are the same condition. At the next digit the exact repository
calculation is

\[
b_{p-1-r}=b_r-2pW_r\pmod {p^2},
\]

where `W_r` is a new harmonic coordinate. A target condition `p|b_r` does
not force `p|W_r`. Thus reflection either duplicates the first row or
introduces a new local unknown; it does not provide a free independent
equation.

### 5.5 No Bernoulli-independence contradiction

Even a rigorously derived polynomial relation among
`B_{p-3},B_{p-5},... (mod p)` would require a theorem showing that the
relation is nonzero for the relevant primes. Statements such as “false for
almost all primes” are insufficient for a pointwise theorem: a hypothetical
large `R_n` could be supported on the exceptional primes unless their total
logarithmic weight along the moving diagonal is itself proved to be `o(n)`.
That is another form of the original problem.

There is also no general algebraic-independence theorem for these varying
Bernoulli residues across varying finite fields that could be inserted at
this point. Kummer congruences and finite harmonic identities create
dependencies rather than generic independent coordinates.

### 5.6 The inverse-limit no-go

The logical obstruction can be stated independently of Apéry numbers.

**Lemma (diagonal tautology).** Let `R>=2` and `A` be a nonnegative integer.
For every `k>=1`, let `x_k` be the least nonnegative representative of
`A mod R^k`. Then:

1. `x_{k+1}≡x_k (mod R^k)`;
2. `(x_k)_k` is the diagonal image of `A` in
   `lim_k Z/R^kZ`;
3. `x_k=A` for every `k` with `R^k>A`.

The proof is immediate from the definition.

Apply this with `R=R_n` and `A=A_n`. If a proposed endpoint tower is proved
to compute `A_n mod R_n^k` at every level, its compatible inverse-limit point
is diagonal and its canonical lifts stabilize. It is impossible to prove
non-diagonality from the same congruences.

Writing

\[
d_k=\frac{x_{k+1}-x_k}{R^k}
\]

only produces the base-`R` digits of `A`. Showing one digit `d_k` is nonzero
does not contradict eventual stabilization. One would need infinitely many
nonzero digits for the same fixed `A_n`, but an ordinary integer has only
finitely many.

If a symbolic Bernoulli tower instead defines a non-diagonal profinite
integer, the conclusion is that the symbolic tower is not the residue tower
of the actual `A_n`; it is not a contradiction to the size of `R_n`.

### 5.7 Route B verdict and reopen condition

**Verdict:** Route B has no defined defect module at the claimed levels, and
the proposed equation count has the wrong variables and quantifiers. Even a
correct all-order local expansion would reconstruct the actual integer rather
than prove its non-diagonality.

A genuine reopening needs all of the following.

1. A proved shifted expansion for `b_{qp+r}` or `b_{qp+r}/p` to arbitrary
   fixed order, uniformly in the moving residue `r`, not only a Gauss
   congruence at `np^m`.
2. A rigorously defined defect module after all harmonic, Kummer, reflection,
   and integrality relations are imposed.
3. An elimination identity for a **fixed outer `n`** whose rows do not add
   new Apéry-value unknowns as fast as they add equations.
4. A proof that the eliminated integer is nonzero.
5. A cross-prime height bound strong enough that its target-prime divisibility
   exceeds its Archimedean size.

Items 4 and 5 are the actual separation theorem. Rank growth without them is
only local reconstruction.

## 6. The V1 defect identity does not create a growing-rank tower

Let `mathcal P_p(z)` be the truncated Apéry polynomial. The proved V1 identity
over `F_p[z]` is

\[
(z+1)^3\mathcal P_p(z+1)-P(z)\mathcal P_p(z)
 +z^3\mathcal P_p(z-1)
=-16(2z+1)(z^p-z)^2.
\tag{V1}
\]

At a residue `a in F_p`, put `z=a+t`. Then

\[
z^p-z=t^p-t
\]

and the defect is

\[
-16(2a+1+2t)(t^p-t)^2.
\tag{V2}
\]

For Hasse orders below `p+1`, (V2) has only:

\[
[t^2]=-16(2a+1),
\qquad
[t^3]=-32.
\]

Thus:

- the value and first-derivative defects vanish;
- the next two local jets are explicitly forced;
- there is no sequence of new low-order independent Bernoulli coordinates;
- at the central residue `a=-1/2`, the vanishing order rises from two to
  three, but still by only one fixed order.

This explains why V1 is powerful for close-pair transport and the
`|Z_p|<<p^{2/3}` theorem: it gives homogeneous recurrence transport through
first order. It does not give higher `p`-adic divisibility of the normalized
quotient `b_n/p`.

To use V1 modulo `p^k`, one would first need an integral lift of its boundary
term. The factorization by `(z^p-z)^2` is a characteristic-`p` identity; its
higher lifts introduce Wilson quotients, harmonic sums, and normalized
Frobenius jets. Existing exact audits in
`Q32_CODEX_RESUME_2026-07-23.md` show:

- Section 14: the natural `q=1` Lucas defect normally has valuation exactly
  one;
- Section 20/Q376: after the forced reflected factor, the next quotient is
  Gessel's derivative defect, a new unconstrained coordinate;
- Section 118: mod-`p^2` reflection introduces the independent harmonic
  value `W_r`;
- Section 150/Q824: higher Frobenius jets saturate at homogeneous degree;
  one bad equation gives at most one selective prime factor per row degree.

Therefore repeated differentiation or exterior powers increase selective
valuation and algebraic/Archimedean degree together. A new algebraic relation
among the normalized singular jets would be required to break saturation.

## 7. Audit of the remaining temporary scripts

### `crt_tower_digits.py`

This script:

- computes `A_n` exactly first;
- prints the tautological residues `A_n mod R_n^k`;
- never derives those residues from Liu/Sun endpoint formulas;
- uses only top-half primes despite claiming all middle primes;
- correctly observes that no distinct congruent integer can lie within
  distance `<R_n^k` of `A_n`.

It is useful as a base-`R_n` digit display but supplies no separation
evidence.

### `verify_quotient_digit.py`

This script fails when it tries to invert `(r+1)^3` modulo `p^2` at the
singular step `r=p-1`:

```text
ValueError: base is not invertible for the given modulus
```

This is exactly why the scaled companion values and the block theorem are
needed.

### `verify_quotient_digit_v2.py`

The repaired script stays in exact integer arithmetic through the singular
step and then reduces modulo `p^2`. It successfully reproduces
`(b_n/p) mod p`, but this is not an independent endpoint formula: it has
already computed the exact Apéry sequence through the problematic index. It
verifies arithmetic implementation, not the claimed Bernoulli tower.

## 8. What the finite-level countermodel teaches

The source note asks what genuinely global property of the Apéry numbers is
absent from an artificial sequence that matches all finite local data while
prescribing the quotient.

It is **not** merely “more `p`-adic levels”. Any fixed integer already has a
compatible tower at every level, and CRT can impose arbitrarily many finite
local digits. The audited countermodels in the repository show that the
following data are also insufficient by themselves:

- a rank-two recurrence module;
- the determinant/Casoratian law;
- a common initial vector;
- finitely many comparison solutions;
- generic interpolation, Smith, Fitting, or resultant data;
- correct first-order exponential growth.

The feature not preserved cheaply by the strongest countermodels is the
simultaneous combination of:

1. the exact canonical cubic `P(n)`;
2. the canonical bounded-height integral gauge and initial state `(1,5)`;
3. the hypergeometric/G-function or Dwork realization across all indices;
4. a global Archimedean height bound coupled to the same nonzero integer that
   accumulates the local valuations.

The fourth item is essential. The first three are structural sources from
which such an integer might be built, but they do not themselves constitute
a separation theorem. In the fixed-state lattice calculation of repository
Section 143, the unknown radical is exactly the second successive minimum;
changing gauge can realize the prescribed local star only by writing that
radical into the coefficient height.

Thus the true missing theorem is not “the CRT lift is non-diagonal”. It is
closer to:

> Construct an Apéry-specific nonzero integer `C_{n,k}` whose target-prime
> valuation is at least `k log R_n`, while
> `log|C_{n,k}|=o(k n)`, uniformly for arbitrarily large fixed `k`, and prove
> that its divisibility is not a tautological multiple of `A_n-x_{n,k}`.

No such integer is produced by Routes A or B.

## 9. Final assessment

### Route A

**Closed in its present form.** The exact local formula is (A.1), and the
global fake-value defect is (A.3). Both are saturated reformulations of
`X≡A_n`; neither yields non-coincidence, additional `p`-adic depth, or a
sublinear-height cross-prime certificate.

### Route B

**Not yet mathematically instantiated.** The cited unconditional result gives
only one extra Gauss digit involving `B_{p-3}` at multiplicative indices;
Sun's next digit remains conjectural. The claimed high-order endpoint tower
and its partition dimensions are unsupported; the temporary count is
internally incorrect; varying `m` does not preserve the fixed outer index;
and correct residue formulas force a diagonal inverse limit rather than
exclude it.

### Route C

**Impossible as stated.** The true residue tower of the integer `A_n` is
diagonal and eventually stabilizes by definition. A useful replacement would
need a second, independently constructed lift with a strict Archimedean
separation, not the canonical residue lift itself.

### Research implication

These verdicts leave the conjecture open. They narrow the viable frontier:
any further quotient-separation attempt must start with an explicit
Apéry-specific global nonzero certificate and its height ledger. More local
digits, more Casoratians, or symbolic Bernoulli dimension counts do not
address the missing cross-prime Archimedean inequality.

## 10. Follow-up audit: sparse zero fibers do not give pointwise separation

**Date:** 2026-07-29.

This section audits the subsequent computational observation that the
single-prime zero fibers

\[
 Z_p=\{0\le r<p:b_r\equiv0\pmod p\}
\]

are very small, while the top-half target set

\[
 T_n=\{p:n/2<p\le n,\ p\text{ prime},\ n-p\in Z_p\}
\tag{10.1}
\]

has at most three elements in the tested range. The proposed inference was

\[
 |Z_p|\le D\text{ uniformly}\quad\Longrightarrow\quad |T_n|=O(1),
\tag{10.2}
\]

which would make the top-half logarithmic weight \(O(\log n)\).

The data are real and striking, but (10.2) is false. It confuses a bound on
the degree of every **prime row** with a bound on the degree of an
**integer column**. Reflection, nonconsecutivity, the Casoratian, and the
currently known defect identities do not repair that change of quantifiers.
Moreover, the observed Poisson law predicts that both \(\max |Z_p|\) and
\(\max |T_n|\) grow slowly rather than remain absolutely bounded.

The useful outcome is a corrected target: a horizontal, cross-prime
collision estimate for the translated sets \(p+Z_p\). This is exactly the
kind of datum missing from Routes A and B.

### 10.1 Exact data audit

The audit used the independently validated binary scan

`problems/3.2/data_zp_pairs.bin`

and the scanner report

`problems/3.2/bn_bigscan_report.md`.

The file contains every pair \((p,r)\), for \(7\le p\le2,000,000\), with
\(p\) prime, \(0\le r<p\), and \(p\mid b_r\). Its SHA-256, recorded by the
scanner, is

```text
8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d
```

An independent parser reproduced the following.

For \(p\le10,000\), among the 1,226 primes \(p\ge7\):

| \(|Z_p|\) | number of primes |
|---:|---:|
| 0 | 756 |
| 1 | 2 |
| 2 | 369 |
| 4 | 83 |
| 6 | 14 |
| 8 | 2 |

Thus the exact mean is \(1172/1226=0.955954\ldots\). The two omitted rows in
the originally quoted even-valued histogram are

\[
 Z_{11}=\{5\},\qquad Z_{3137}=\{1568\}.
\tag{10.3}
\]

They are both central fixed-point zeros. Consequently, “\(|Z_p|\) is always
even” is already false below \(10^4\); the correct parity statement is
proved in Section 10.2.

Under the scanner's convention \(p\ge7\), the exact top-half histogram for
\(1\le n\le10,000\) is:

| \(|T_n|\) | number of \(n\) |
|---:|---:|
| 0 | 9,204 |
| 1 | 769 |
| 2 | 26 |
| 3 | 1 |

For the literal all-prime definition (10.1), \(Z_5=\{1,3\}\) adds singleton
hits at \(n=6,8\), changing only the first two counts to 9,202 and 771.
Primes 2 and 3 add nothing. All maxima, double hits, triple hits, and
asymptotic statements below are unchanged.

The phrase “27 examples with two primes” is correct only if it means “at
least two”: there are 26 exact double hits and the one triple hit

\[
 n=321,\qquad T_{321}=\{179,193,211\}.
\]

At this \(n\),

\[
 \frac{\log\prod_{p\in T_n}p}{\log n}
 =2.737952929\ldots .
\]

The full scan through two million gives:

| \(|Z_p|\) | number of primes |
|---:|---:|
| 0 | 90,375 |
| 1 | 2 |
| 2 | 45,020 |
| 4 | 11,375 |
| 6 | 1,875 |
| 8 | 257 |
| 10 | 24 |
| 12 | 2 |

There are 148,930 primes and 149,112 zero pairs, so the mean is

\[
 1.001222050628\ldots .
\]

The maximum is now 12, at

\[
 p=159,977,\qquad p=1,823,963.
\]

Thus the earlier maximum 8 was a finite-range value, not a stable ceiling.
The assertion in `Q3.2_density_theorem.tex` that the maximum zero count is
16 through \(80,000\) is also false: the exact full scan gives maximum 8 in
that range. The value 16 occurs in a separate maximum-**fiber** computation
\(\max_a N_p(a)\); it is not the zero-fiber maximum.

For \(1\le n\le2,000,000\), again with the report's \(p\ge7\) convention,
the top-half histogram is:

| \(|T_n|\) | number of \(n\) |
|---:|---:|
| 0 | 1,896,672 |
| 1 | 100,670 |
| 2 | 2,605 |
| 3 | 53 |

Including \(p=5\) changes only the first two counts to 1,896,670 and
100,672.

There is no fourfold hit in this range. The largest logarithmic-radical
exponent is

\[
 2.960803420\ldots
\]

at

\[
 n=1,689,988,\qquad
 T_n=\{1,253,839,\ 1,328,449,\ 1,651,781\}.
\]

This confirms the quoted finite-range phenomenon while extending it by a
factor 200.

### 10.2 The parity is elementary reflection, not a new K3 signal

For \(m<p\), write

\[
 \binom zm\binom{z+m}{m}
 =\frac{1}{(m!)^2}
   \prod_{j=0}^{m-1}(z-j)(z+j+1)
 =\frac{1}{(m!)^2}
   \prod_{j=0}^{m-1}\bigl(z(z+1)-j(j+1)\bigr).
\tag{10.4}
\]

Define

\[
 H_p(X)=
 \sum_{m=0}^{(p-1)/2}
 \frac{\prod_{j=0}^{m-1}(X-j(j+1))^2}{(m!)^4}
 \in\mathbf F_p[X].
\tag{10.5}
\]

For every \(r\in\{0,\ldots,p-1\}\), terms with

\[
 m>\min(r,p-1-r)
\]

vanish at \(X=r(r+1)\). The remaining terms are exactly the defining
Apéry sum, so

\[
 b_r=H_p(r(r+1))\quad\text{in }\mathbf F_p.
\tag{10.6}
\]

Since

\[
 (-1-r)(-r)=r(r+1),
\]

(10.6) immediately gives

\[
 b_{p-1-r}\equiv b_r\pmod p.
\tag{10.7}
\]

Notice that (10.5) has degree \(p-1\) in \(X\), not bounded degree: the
top summand has degree \(2\cdot(p-1)/2=p-1\) and nonzero leading
coefficient.

The involution \(r\mapsto p-1-r\) has the unique fixed point

\[
 c_p=(p-1)/2.
\]

Therefore the exact parity law is

\[
 |Z_p|\equiv
 \mathbf 1_{\{b_{(p-1)/2}\equiv0\pmod p\}}
 \pmod2.
\tag{10.8}
\]

Ahlgren--Ono identify the central value with the \(p\)-th coefficient
\(\gamma(p)\) of

\[
 \eta(2\tau)^4\eta(4\tau)^4
\]

(in fact via a stronger supercongruence). Thus the fixed weight-four form
controls precisely the possible unpaired root. It says nothing about the
off-center reflection pairs which account for essentially all of the data.

There is an additional pointwise reason that the central slice is harmless.
If \(n=qp+r\) and \(r=(p-1)/2\), then

\[
 2n+1=(2q+1)p.
\tag{10.9}
\]

Hence all central bad primes at a fixed \(n\) divide the single integer
\(2n+1\), and their total logarithmic weight is at most

\[
 \log(2n+1)=O(\log n).
\tag{10.10}
\]

For the top-half arc \(q=1\), (10.9) is \(3p=2n+1\), so there is at most one
central candidate. Consequently, a theorem about the nonordinary primes of
the fixed eta-product cannot be the missing pointwise theorem: that entire
fixed-form slice is already negligible. The residual problem consists of
the moving off-center roots.

### 10.3 The row-versus-column counterexample

The logical failure of (10.2) is absolute and does not depend on a weakness
of the known bound for \(Z_p\).

Fix a large integer \(n_0\). For every prime

\[
 n_0/2<p\le n_0,
\]

put

\[
 r_p=n_0-p
\]

and define an artificial reflected zero set

\[
 S_p=\{r_p,\ p-1-r_p\}.
\tag{10.11}
\]

Then:

1. \(|S_p|\le2\) for every \(p\);
2. \(S_p\) is invariant under \(r\mapsto p-1-r\);
3. \(n_0-p=r_p\in S_p\) for every prime in the interval.

Thus the column over \(n_0\) contains every prime in
\((n_0/2,n_0]\), even though every row has size at most two. Its logarithmic
weight is

\[
 \sum_{n_0/2<p\le n_0}\log p
 \sim n_0/2.
\tag{10.12}
\]

The possible coincidences in which the two entries of (10.11) collapse, or
become consecutive, occur for at most \(O(1)\) primes. Omitting those primes
preserves (10.12). Thus the model can also obey the central exception rule
and nonconsecutivity, up to an immaterial bounded deletion.

This countermodel does not claim that the actual Apéry zero sets are
adversarial. It proves the precise logical statement:

> No theorem whose hypotheses see only \(|Z_p|\), reflection, central
> parity, and within-row spacing can bound \(|T_n|\) pointwise.

The reflection partner of an actual edge \((p,n=p+r)\) lies at

\[
 n'=p+(p-1-r)=3p-1-n.
\tag{10.13}
\]

It is a second edge in a generally different integer column. Formula
(10.13) explains why reflection doubles row counts but does not constrain
the number of distinct primes meeting one fixed \(n\).

In bipartite-graph language, a uniform degree bound on the prime side gives
no maximum-degree bound on the integer side. The proposed deduction was the
forbidden transposition of those two degrees.

### 10.4 The Poisson model predicts slow unbounded growth

The finite data do not merely fail to prove an absolute constant. Their
best-fitting model predicts that no such constant exists.

After removing the two central exceptions, put \(Y_p=|Z_p|/2\). The scan is
extremely close to

\[
 Y_p\sim\operatorname{Poisson}(1/2).
\tag{10.14}
\]

For 148,930 primes, (10.14) predicts:

| threshold | observed number | Poisson prediction |
|---:|---:|---:|
| \(|Z_p|\ge8\) | 283 | 260.87 |
| \(|Z_p|\ge10\) | 26 | 25.63 |
| \(|Z_p|\ge12\) | 2 | 2.11 |
| \(|Z_p|\ge14\) | 0 | 0.15 |

The agreement at the extreme tail is especially relevant: the two rows of
size 12 are almost exactly the predicted number. For independent
Poisson variables of fixed mean, the sample maximum grows like

\[
 \frac{\log M}{\log\log M}
\]

over \(M\) samples. In this model, \(\max_{p\le P}|Z_p|\) is therefore
unbounded (with an extra factor two from reflection), although it grows very
slowly.

For the top-half target count, a random sparse-row model gives

\[
 \lambda_n:=\mathbf E|T_n|
 \simeq\sum_{n/2<p\le n}\frac1p
 \simeq\frac{\log2}{\log n}.
\tag{10.15}
\]

Using the literal varying-mean Poisson mixture in (10.15), the expected
numbers of exceptional columns are:

| range | expected \(|T_n|\ge3\) | observed | expected \(|T_n|\ge4\) | observed |
|---:|---:|---:|---:|---:|
| \(n\le10,000\) | 1.251 | 1 | 0.055 | 0 |
| \(n\le2,000,000\) | 45.606 | 53 | 0.636 | 0 |

Under an independence heuristic, seeing no fourfold hit through two million
has probability approximately

\[
 e^{-0.636}=0.529\ldots,
\]

so it is not evidence for a hard ceiling at three. Numerical integration of
the same model makes the expected cumulative number of fourfold hits cross
one near \(N=4,000,000\). This is a heuristic prediction, not a theorem, but
it shows that a modest extension of the finite range could naturally
produce a fourfold hit.

More generally, (10.15) predicts a slowly unbounded top-half maximum,
roughly of order

\[
 \frac{\log N}{2\log\log N},
\]

not \(O(1)\). This growth is still vastly smaller than what Problem 3.2
requires. At the extreme-value scale it predicts

\[
 \max_{n\le N}\sum_{p\in T_n}\log p
 \asymp \frac{(\log N)^2}{2\log\log N},
\]

up to lower-order factors. Thus the finite-range law
\(\log R_n=O(\log n)\) is also not the natural asymptotic prediction, even
though the predicted polylogarithmic growth is still \(o(n)\).

### 10.5 A new full-middle scatter: maximum seven and Poisson mean \(\log2\)

The top-half statistic sees only the quotient \(q=1\). For the actual
middle-prime \(b\)-channel, define

\[
 K_b(n)=
 \#\{\sqrt n<p\le n:p\text{ prime},\ b_{n\bmod p}\equiv0\pmod p\}.
\tag{10.16}
\]

The saved pair file suffices to compute (10.16) exactly through two million.
For each recorded \((p,r)\), scatter it to

\[
 n=qp+r,\qquad q\ge1,\qquad n\le2,000,000,\qquad n<p^2.
\tag{10.17}
\]

The last inequality is exactly \(p>\sqrt n\).

The result is

\[
 \max_{n\le2,000,000}K_b(n)=7,
\tag{10.18}
\]

attained at 14 integers, beginning with

```text
208253, 501092, 513596, 855688, 939041, 985143, 1097287, ...
```

For example, the seven hits at \(n=208,253\) are

```text
(p,r,q) =
(499,170,417), (1381,1103,150), (6949,6732,29),
(22697,3980,9), (47017,20185,4), (54629,44366,3),
(195319,12934,1).
```

The histogram on the final million integers is:

| \(K_b(n)\) | count for \(1,000,001\le n\le2,000,000\) |
|---:|---:|
| 0 | 510,950 |
| 1 | 343,332 |
| 2 | 114,928 |
| 3 | 25,808 |
| 4 | 4,341 |
| 5 | 572 |
| 6 | 61 |
| 7 | 8 |

Its mean, variance, and variance-to-mean ratio are

\[
 0.671258,\qquad
 0.671072697436,\qquad
 0.999723947328.
\tag{10.19}
\]

The falling second moment on the same million-point interval is exactly

\[
 \sum_n(K_b(n))_2=450,402,
\qquad
 \frac1{10^6}\sum_n(K_b(n))_2=0.450402.
\]

The Poisson prediction is
\(0.671258^2=0.450587\ldots\). Thus the specific cross-prime collision
statistic isolated later in (10.23), not just the one-point histogram,
already matches the independence model numerically.

This is the predicted scale:

\[
 \sum_{\sqrt n<p\le n}\frac1p
 =\log\log n-\log\log\sqrt n+o(1)
 =\log2+o(1).
\tag{10.20}
\]

A Poisson variable of mean \(0.671258\), sampled one million times,
predicts 6.79 values at least seven and 0.56 values at least eight. The
observed counts are eight and zero. Thus the full-middle data agree with the
same Poisson picture down to the extreme tail.

This new scatter is a useful warning. The maximum has already grown from
five below \(10^4\), to six below \(10^5\), to seven below \(10^6\). The
natural prediction is

\[
 \max_{n\le N}K_b(n)\asymp
 \frac{\log N}{\log\log N},
\tag{10.21}
\]

not a uniform constant. Of course (10.21), if proved even with a much weaker
subpower upper bound, would solve the \(b\)-channel with enormous room.

The separate companion \(a_q\)-channel in the exact support law is not
included in (10.16). Existing divisor counting already gives that channel
sublinear logarithmic weight. The moving \(b_{n\bmod p}\)-channel remains
the pointwise obstruction.

### 10.6 Why K3/Frobenius does not supply the missing implication

There are three distinct arithmetic objects which must not be conflated.

1. **The fixed central modular form.** Ahlgren--Ono control
   \(b_{(p-1)/2}\). By (10.9), this slice is already pointwise harmless.
2. **The geometric Apéry/K3 family parameter.** Ordinary-reduction theorems
   for a fixed K3 surface vary the prime while holding the surface fixed.
3. **The coefficient-index interpolation (10.5).** Here \(p\) is fixed and
   the character/index \(r\) varies. The polynomial itself depends on \(p\)
   and has degree \(p-1\).

The distinction between the last two objects can be made exact. Put

\[
 \mathcal A_p(t)=\sum_{r=0}^{p-1}b_rt^r\in\mathbf F_p[t].
\]

This truncated generating series is the genuine first Hasse--Witt scalar
for the geometric Apéry K3 pencil (in the standard local trivialization).
For \(1\le r\le p-2\), finite Mellin inversion gives

\[
 b_r=-\sum_{t\in\mathbf F_p^\times}\mathcal A_p(t)t^{-r}.
\]

Thus \(b_r=0\) is the vanishing of one multiplicative Fourier coefficient
of the whole Hasse section. It is not the geometric nonordinary-fiber
condition \(\mathcal A_p(t_0)=0\).

Caruso--Fürnsinn--Vargas-Montoya--Zudilin, Theorem 2, prove the precise
factorization

\[
 \mathcal A_p(t)=
 \begin{cases}
  Q_p(t)^2,&p\equiv1,5,7,11\pmod {24},\\
  (t^2-34t+1)Q_p(t)^2,&p\equiv13,17,19,23\pmod {24}.
 \end{cases}
\]

This is strong structure in the **geometric generating variable** \(t\).
It gives quadratic convolution identities among the coefficients, but no
coefficient-support bound by itself: even a square such as
\((1+t^m)^2\) can have arbitrarily many zero coefficients. In particular,
the theorem neither bounds \(Z_p\) nor identifies (10.5) with a bounded
degree K3 Hasse polynomial in the spectral variable \(X\).

There is a stronger exact coefficient package, but it leads to a
convolution problem rather than a zero bound. Let

\[
 F(t)=\sum_{m\ge0}b_mt^m,\qquad D(t)=t^2-34t+1,
\]

and normalize \(Q_p(0)=1\). Uniqueness of a formal square root in odd
characteristic shows that \(Q_p\) is the appropriate truncation modulo
\(p\) of

\[
 S_+(t)=\sqrt{F(t)}
 \quad\text{or}\quad
 S_-(t)=\sqrt{F(t)/D(t)},
\]

according to the two CFVZ branches. Self-reciprocity of
\(\mathcal A_p\) and \(D\) also gives \(Q_p^*=\pm Q_p\).

Direct symmetric-square reduction of the Apéry differential operator gives
the independently checked equations

\[
\begin{aligned}
 4tD S_+''+4(2t^2-51t+1)S_+'+(t-10)S_+&=0,\\
 4tD S_-''+4(4t^2-85t+1)S_-'
   +3(3t-26)S_-&=0.
\end{aligned}
\tag{K3.1}
\]

Writing \(S_\pm=\sum s_m^\pm t^m\), \(s_{-1}^\pm=0\), these are

\[
\begin{aligned}
 4(m+1)^2s_{m+1}^+
 &-(136m^2+68m+10)s_m^+
 +(2m-1)^2s_{m-1}^+=0,\\
 4(m+1)^2s_{m+1}^-
 &-(136m^2+204m+78)s_m^-
 +(2m+1)^2s_{m-1}^-=0.
\end{aligned}
\tag{K3.2}
\]

Thus an Apéry coefficient is exactly a Cauchy-convolution coefficient

\[
 b_r=\sum_{i=0}^r s_i^+s_{r-i}^+,
\]

or, in the other branch, the fixed three-term \(D\)-filter of the analogous
convolution. This rank-two square-root period is genuine extra structure.
It is not yet transverse structure: the convolution map is triangular,
with \(b_r=2s_r^++P_r(s_1^+,\ldots,s_{r-1}^+)\), and the third-order Apéry
operator is precisely the symmetric square of (K3.1). Substitution therefore
recovers the original recurrence rather than a second hit-specific equation.
A useful escape would require a new anti-cancellation or bounded-boundary
telescoping theorem for this canonical convolution; factorization and
reciprocity alone permit reciprocal squares with linearly many zero
coefficients.

The most natural bounded-boundary escape can in fact be ruled out.  After
the self-adjoint gauges

\[
 s_m^\varepsilon=
 (2m+1)^\varepsilon\frac{\binom{2m}{m}}{4^m}y_m^\varepsilon,
 \qquad \varepsilon\in\{0,1\},
\]

the two square-root recurrences have a genuine same-index
Christoffel--Darboux identity.  It sums
\(y_i(\lambda)y_i(\mu)\) at two spectral parameters.  The Apéry
convolution instead pairs the reversed indices \(s_i s_{r-i}\), so this
identity does not apply.

The obstruction is exact.  Put \(u_i=s_i\), \(v_i=s_{r-i}\) and try the
natural bilinear concomitant

\[
 {\cal W}_i=X_i u_{i+1}v_i-Y_i u_iv_{i+1}.
\]

Cancellation of its two cross terms forces

\[
 \frac{X_{i+1}}{X_{i-1}}
 =
 \left\{
 \frac{(i+2)(2r-2i+2\varepsilon-1)}
 {(2i+2\varepsilon+1)(r-i+1)}
 \right\}^{\!2}.
\tag{K3.3}
\]

No nonzero \(X(i,r)\in\overline{\mathbf Q(r)}(i)\) satisfies (K3.3).
For a rational shift quotient \(X(i+1)/X(i-1)\), the zero-minus-pole
multiplicity must balance on every translation-by-two orbit.  The right
side of (K3.3) has integer zero/pole locations \(-2,r+1\) and
half-integer locations \(r+\varepsilon-\tfrac12,
-\varepsilon-\tfrac12\); exactly one of the two pairs lies on the same
parity orbit.  The other orbit is unbalanced.  Thus even a rationally
weighted reversed Green identity already fails before one asks for unit
weight.

Solving (K3.3) separately on the even and odd subsequences does produce a
hypergeometric weight, but its Pochhammer products have length
\(\asymp r\).  Clearing its moving zeros and poles has degree
\(\asymp r\) and logarithmic height \(O(r\log r)\).  It merely encodes the
whole convolution and is not a bounded boundary certificate.

The obstruction extends beyond this two-term ansatz.  A uniformly
bounded-width rational adjacent-state telescoper would sum to an identity

\[
 Q(r)[t^r]S_\varepsilon(t)^2
 =P_0(r)s_r^\varepsilon+P_1(r)s_{r-1}^\varepsilon
\tag{K3.4}
\]

up to finitely many initial terms (and with the fixed \(D\)-filter in the
minus branch).  Generating functions would turn (K3.4) into a nonzero
quadratic--linear differential relation between
\(S_\varepsilon\) and \(\theta S_\varepsilon\).
The CFVZ algebraic pullback identifies this rank-two equation with the
Gauss equation for
\({}_2F_1(\tfrac13,\tfrac23;1;y)\).  Its irreducible monodromy contains
two nontrivial unipotents with distinct fixed lines, so its connected
differential Galois group is \({\rm SL}_2\).  Hence
\(S_\varepsilon,\theta S_\varepsilon\) are algebraically independent over
\(\overline{\mathbf Q(t)}\), contradicting such a relation.

Accordingly the minimal bounded outer telescoper for the observable
\(F=D^\varepsilon S_\varepsilon^2\) is its irreducible
three-dimensional symmetric-square operator.  In coefficient form it is
exactly the original Apéry recurrence.  The square-root factorization
therefore supplies neither a second zero condition nor a lower-degree
gap certificate.  What remains open is a batch anti-cancellation theorem
for actual convolution returns; it cannot be obtained by collapsing each
coefficient to \(O(1)\) boundary states.

Bogomolov--Zarhin prove density-one ordinary reduction, after a finite
extension, for a fixed K3 surface over a number field. That theorem has
neither of the quantifiers needed here:

- it does not bound the number of \(r\)-values in one field
  \(\mathbf F_p\);
- it does not control the moving specialization \(r=n\bmod p\) while both
  \(p\) and the character vary.

Likewise, generic ordinarity of a family does not make its nonordinary
divisor bounded-degree uniformly in \(p\). Classical Hasse invariants
already provide counterexamples to that inference. In the present problem
the mismatch is sharper: \(H_p\) is an interpolation in the coefficient
index, not the geometric Hasse invariant in the K3 pencil parameter.

Archimedean Frobenius equidistribution and Sato--Tate control fixed
continuous test functions of normalized traces. The exact event
\(b_r\equiv0\pmod p\) is a same-characteristic, \(p\)-adically shrinking
target. Passing from equidistribution to an exact residue-zero count would
require a local limit theorem uniform in the moving character. No cited
Katz/Deligne, ordinary-reduction, or Chebotarev result supplies such a
theorem.

Therefore K3/Frobenius currently explains why modular and hypergeometric
structures are present, but it proves neither:

\[
 |Z_p|=O(1)
\]

nor the genuinely needed horizontal statement

\[
 \#\{p:n\bmod p\in Z_p\}=o(n/\log n)
\quad\text{uniformly in }n.
\]

The reflection parity is already fully explained by (10.4), with no
algebraic geometry.

#### An unconditional divisor-sensitive Kummer-order pruning

The fixed-\(n\) relation does give one useful cross-prime reduction that is
independent of the zero-fiber bound. Write

\[
 n=qp+r,\qquad q=\lfloor n/p\rfloor,\qquad0\le r<p,
\]

and let

\[
 d_{p,n}=\operatorname{ord}(\omega_p^r)
 =\frac{p-1}{\gcd(p-1,r)}
 =\frac{p-1}{\gcd(p-1,n-q)}.
\tag{KO.1}
\]

The last equality follows from \(r\equiv n-q\pmod {p-1}\). For fixed
\(n,q,D\), one has the purely divisor-theoretic bound

\[
 \boxed{
 \#\left\{p:\left\lfloor\frac np\right\rfloor=q,\ d_{p,n}\le D\right\}
 \le
 \tau(n-q)\left(\frac Dq+2\right).
 }
\tag{KO.2}
\]

Indeed, put \(g=\gcd(p-1,r)\), \(p-1=dg\), and \(r=jg\). Then
\(0\le j\le d\), \(\gcd(d,j)=1\), and

\[
 n-q=(qd+j)g.
\]

Writing \(h=qd+j\), one has

\[
 h\mid n-q,\qquad qd\le h\le(q+1)d,\qquad
 p=1+d\frac{n-q}{h}.
\]

For each divisor \(h\mid n-q\), the possible integers \(d\le D\) lie in
\([h/(q+1),h/q]\), an interval containing at most \(D/q+2\) integers.
This proves (KO.2) without using either primality or the Apéry condition.
It sharpens the earlier \(O(QD\,n^{o(1)})\) order split by retaining the
factor \(1/q\).

Let

\[
 M(n)=\max_{m\le n}\tau(m),\qquad
 L_n=4M(n)(\log n)^2=n^{o(1)}.
\tag{KO.3}
\]

For \(p>\sqrt n\), if
\(\gcd(p-1,n-q)>L_n\), then

\[
 d_{p,n}<\frac{n}{qL_n}.
\]

Applying (KO.2) with \(D=n/(qL_n)\), summing over
\(q<\sqrt n\), and bounding every prime by \(\log n\), gives

\[
\begin{aligned}
 &\sum_{\substack{p>\sqrt n\\
   \gcd(p-1,n-\lfloor n/p\rfloor)>L_n}}\log p\\
 &\quad\le
 \sum_{q<\sqrt n}
 \left(\frac{n}{q^2L_n}+2\right)\tau(n-q)\log n\\
 &\quad\ll
 \frac{nM(n)\log n}{L_n}
 +\sqrt n\,M(n)\log n
 =O\!\left(\frac n{\log n}\right)+n^{1/2+o(1)}
 =o(n).
\end{aligned}
\tag{KO.4}
\]

This is an unconditional, all-\(n\) improvement over the older
bounded-order pruning. The full unresolved lower-digit channel may therefore
be restricted to

\[
 \gcd\!\left(p-1,n-\lfloor n/p\rfloor\right)\le L_n,
\qquad
 d_{p,n}\ge\frac{p-1}{L_n}=p^{1-o(1)}.
\]

The uniform threshold (KO.3) is convenient but not optimal. The same proof
allows a quotient-adaptive cutoff. Let \(F(n)\to\infty\) with
\(F(n)=n^{o(1)}\), and put

\[
 L_{n,q}=\tau(n-q)\log n\,F(n).
\tag{KO.5}
\]

For primes with \(\gcd(p-1,n-q)>L_{n,q}\), (KO.2) with
\(D=n/(qL_{n,q})\) gives total logarithmic weight

\[
\begin{aligned}
 &\sum_{q<\sqrt n}
 \left(\frac{n}{q^2L_{n,q}}+2\right)
 \tau(n-q)\log n\\
 &\qquad\le
 \frac n{F(n)}\sum_{q\ge1}\frac1{q^2}
 +2\log n\sum_{q<\sqrt n}\tau(n-q)\\
 &\qquad=O\!\left(\frac n{F(n)}\right)+n^{1/2+o(1)}
 =o(n).
\end{aligned}
\tag{KO.6}
\]

Thus the sharper canonical residual is

\[
 \gcd(p-1,n-q)\le\tau(n-q)\log n\,F(n),
\qquad
 d_{p,n}\ge
 \frac{p-1}{\tau(n-q)\log n\,F(n)}
 =p^{1-o(1)}.
\tag{KO.7}
\]

For the top half this becomes
\(\gcd(p-1,n-1)\le\tau(n-1)\log n\,F(n)\). The reduction is sharp in the
sense relevant here: small kernels appear frequently in the computed target
set, and the sieve gives no estimate for their Apéry-zero condition. It
removes bounded- and moderately growing-order motives; it does not replace
the horizontal collision theorem.

#### What the order sieve does and does not give in pair energy

There is a deterministic \(L^2\) version of the divisor parametrization,
but it stops far below the nearly primitive scale.  This is worth recording
because otherwise (KO.6) can easily be misused through an invalid
\(L^1\)-to-\(L^2\) inference.

Fix \(q\), let \(N<n\le2N\), put \(m=n-q\), and let
\(K_{q,\le D}(n)\) count target primes in the \(q\)-arc whose selected
character order \(d_{p,n}\) is at most \(D\).  Set

\[
 H=(q+1)D,\qquad
 a_{q,D}(h)=
 \#\left\{d\le D:
       \frac{h}{q+1}\le d\le\frac{h}{q},\
       (d,h)=1\right\}.
\tag{KO.8}
\]

The parametrization used in (KO.2) injects every such prime into one pair
\((h,d)\), with \(h\mid m\).  Consequently

\[
 K_{q,\le D}(n)
 \le A_{q,D}(m):=
 \sum_{\substack{h\mid m\\h\le H}}a_{q,D}(h),
 \qquad
 a_{q,D}(h)\le\frac{h}{q(q+1)}+1.
\tag{KO.9}
\]

Expanding the square and counting common multiples gives

\[
\begin{aligned}
 \sum_{N<n\le2N}\bigl(K_{q,\le D}(n)\bigr)_2
 &\le \sum_{N-q<m\le2N-q}A_{q,D}(m)^2\\
 &\ll
 N\left\{
  \frac{D^2}{q^2}\log(2H)
  +\frac Dq\log^2(2H)
  +\log^3(2H)
 \right\}
 +(D^2+qD)^2.
\end{aligned}
\tag{KO.10}
\]

Indeed,

\[
 \sum_{h,k\le H}
 \frac{a_{q,D}(h)a_{q,D}(k)(h,k)}{hk}
 =
 \sum_{e\le H}\varphi(e)
 \left(
   \sum_{\substack{h\le H\\e\mid h}}
   \frac{a_{q,D}(h)}h
 \right)^2,
\]

and the inner sum is

\[
 \ll \frac1e
 \left\{\frac Dq+\log\frac{2H}{e}\right\}.
\]

The last term in (KO.10) is the accumulated \(+1\) in the count of
multiples.  If \(H\le\sqrt N\), it is absorbed into the main
least-common-multiple term.  In particular, for every fixed \(q\),

\[
 D\le\frac{\sqrt N}{\log^2N}
 \quad\Longrightarrow\quad
 \sum_{N<n\le2N}\bigl(K_{q,\le D}(n)\bigr)_2
 =o_q\!\left(\frac{N^2}{\log^2N}\right).
\tag{KO.11}
\]

Thus a pair in which **both** characters are low-order can be removed
almost to the square-root scale.

The mixed estimate is weaker.  Write

\[
 M_N=\max_{m\le2N}\tau(m).
\]

The vertical bound \(|Z_p|\ll p^{2/3}\) gives the first moment

\[
 \sum_{N<n\le2N}K_q(n)
 \ll
 \frac{N^{5/3}}{q^{5/3}\log(N/q)}.
\tag{KO.12}
\]

Since (KO.2) gives
\(K_{q,\le D}(n)\le M_N(D/q+2)\), the ordered-pair energy with
at least one low-order member satisfies

\[
\begin{aligned}
 \mathcal E_q^{L*}(N;D)
 &\le
 2\max_n K_{q,\le D}(n)
   \sum_{N<n\le2N}K_q(n)\\
 &\ll
 \frac{M_N(D+q)N^{5/3}}
 {q^{8/3}\log(N/q)}.
\end{aligned}
\tag{KO.13}
\]

For fixed \(q\), this is negligible at the sufficient pair-energy scale
when \(D\le N^{1/3-\epsilon}\).  More exactly, it is
\(o(N^2/\log^2N)\) whenever

\[
 D+q=
 o\left(
 \frac{q^{8/3}N^{1/3}\log(N/q)}
 {M_N\log^2N}
 \right).
\tag{KO.14}
\]

The exponent \(1/3\) is \(1-2/3\), and therefore records precisely the
present vertical zero-fiber exponent.  No optimization of the divisor
cutoff changes it.

There is also a uniform harmless range in \(q\).  Brun--Titchmarsh bounds
the number of candidate primes in one \(q\)-arc, while (KO.12) bounds the
first moment, giving, uniformly for
\(q\le N^{1/2-\delta}\),

\[
 \mathcal E_q(N)
 \ll_\delta
 \frac{N^{8/3}}
 {q^{8/3}(q+1)
  \log(N/q)\log(N/q^2)}.
\tag{KO.15}
\]

Hence every \(q\ge N^{2/11+\epsilon}\) in this range is already negligible.
For the pointwise theorem this is secondary, because a hypothetical
linear-weight spike can first be localized to finitely many fixed
quotients.

The correct way to combine these statements with (KO.6) is not
interpolation.  Define the nearly primitive residual count

\[
\begin{aligned}
 K_q^{\rm np}(n;F)=
 \#\bigl\{p:\;& n/(q+1)<p\le n/q,\quad
                 n-qp\in Z_p,\\
              &\gcd(p-1,n-q)
                 \le\tau(n-q)\log n\,F(n)\bigr\}
\end{aligned}
\]

and its energy

\[
 \mathcal E_q^{\rm np}(N;F)=
 \sum_{N<n\le2N}\bigl(K_q^{\rm np}(n;F)\bigr)_2.
\]

Then the strictly weaker sufficient theorem is

\[
 \boxed{
 \text{for every fixed }q,\qquad
 \mathcal E_q^{\rm np}(N;F)
 =o_q\!\left(\frac{N^2}{\log^2N}\right).
 }
\tag{KO.16}
\]

To prove this implication, apply (KO.6) pointwise first.  Its discarded
primes have \(o(n)\) total logarithmic weight, so any forbidden spike
leaves \(\gg n/\log n\) nearly primitive targets.  The same finite-\(q\)
pigeonhole argument used for (10.35) then contradicts (KO.16).

Conversely, (KO.6) does not imply that the discarded incidences have small
shell pair energy: they could be concentrated on a sparse set of rows and
meet many fresh high-order zeros there.  Nor do (KO.10)--(KO.13) touch
the residual in (KO.16).  For fixed \(q\), its orders are
\(d_{p,n}=N^{1-o(1)}\), whereas the mixed deterministic estimate ends at
\(N^{1/3-o(1)}\).  The exact remaining problem is therefore correlation
between two nearly primitive, different-characteristic Apéry
coefficients.

This residual still admits the same logical adversary as the earlier
reflection star.  Fix \(q\), choose primes \(m_j\to\infty\) so rapidly
that the intervals

\[
 I_j=\left(\frac{m_j+q}{q+1},\frac{m_j+q}{q}\right]
\]

are pairwise disjoint, and put \(n_j=m_j+q\).  For every prime
\(p\in I_j\), apart from the \(O(1)\) central and boundary coincidences,
set

\[
 r_p=n_j-qp,\qquad
 S_p=\{r_p,p-1-r_p\}.
\tag{KO.17}
\]

Here \(p-1<m_j=n_j-q\), so

\[
 \gcd(p-1,n_j-q)=1,\qquad d_{p,n_j}=p-1.
\tag{KO.18}
\]

The sets (KO.17) are reflected, have no adjacent distinct elements, and
put \(\asymp_q n_j/\log n_j\) primes on one column.  They also obey the
known unsaturated continuant divisibility: for the reflected interval
with lower endpoint \(a\) and gap \(h\),

\[
 2a+h+1=p
\]

is the universal reflection factor of \(N_h(a)\), independently of any
zero condition.  Saturating the continuant removes exactly this factor.
Thus even full character order, reflection, nonconsecutivity, and the
forced continuant factor do not logically imply horizontal dispersion.
This is an abstract zero-set adversary, not a construction of actual
Apéry zeros; its purpose is to identify the additional arithmetic input
that (KO.16) must use.

The exact scan confirms that the hard sector is already populated. Among
106,039 top-half hits with \(n\le2,000,000\), 39,804 have
\(\gcd(p-1,n-1)=1\). Of the 53 triple-hit columns, 12 have kernel 1 at all
three primes, and 21 have the same kernel at all three primes. For example,

\[
 n=11576,\qquad p=8893,9319,11437
\]

is a triple with \(\gcd(p-1,n-1)=1\) for every \(p\). Thus the
nearly-primitive restriction is a genuine theorem-level pruning, not an
explanation of the observed low column multiplicity.

### 10.7 Consequences for Routes A and B

The new zero-fiber data do not reopen either audited route.

There is a stronger exact saturation statement. Let
\(\Omega_p=\{z(z+1):z\in Z_p\}\), so
\(q_p=|\Omega_p|\) is the number of reflection orbits. For a fixed target
column \(n\), let \(\mathcal L_d(n)\) be the lattice of integral
polynomials \(F(Y)\) of degree at most \(d\) which vanish modulo each target
prime \(p\) at every point of \(\Omega_p\). Since
\((n-p)(n-p+1)\equiv n(n+1)\pmod p\), every such polynomial satisfies

\[
 \prod_{p\in T_n}p\mid F(n(n+1)).
\]

For \(d\ge\max_p(q_p-1)\), interpolation and CRT give the exact local/global
indices

\[
 [\mathbf Z^{d+1}:\mathcal L_d(n)]
 =\prod_{p\in T_n}p^{q_p},
\qquad
 [(Y-N)\mathbf Z[Y]_{\le d-1}:
   \mathcal L_d(n)\cap(Y-N)\mathbf Z[Y]_{\le d-1}]
 =\prod_{p\in T_n}p^{q_p-1},
\tag{S.1}
\]

where \(N=n(n+1)\). Their transverse quotient is exactly
\(\prod_{p\in T_n}p\). Thus the entire zero fiber contributes \(q_p\)
conditions, but the fake-evaluation subspace \(F(N)=0\) consumes
\(q_p-1\) of them. The observed values \(q_p\le6\) through two million
cancel completely from the only quotient capable of producing a nonzero
global evaluation.

#### Route A: Casoratian

The Casoratian is a same-prime, same-recurrence determinant. It can prove:

- nonvanishing of the local solution pair;
- absence of simultaneous local zeros;
- transport and close-pair certificates inside one \(Z_p\);
- saturation of a locally imposed fake quotient.

More precisely, use integral recurrence rows \(L_m\) which evaluate an
arbitrary initial state at \(m\). If \(p>b\) and \(a,b\in Z_p\), then
\(L_a\) and \(L_b\) are nonzero rows with the same kernel, namely the line
spanned by the Apéry initial state. Hence all rows indexed by \(Z_p\) have
rank exactly one over \(\mathbf F_p\), regardless of whether the fiber has
2, 6, or 12 elements. Every nonzero exterior construction therefore
reduces to a \(2\times2\) row minor; all higher homogeneous exterior powers
vanish.

The forced reflected minor is also saturated. If the two indices are
\(t,t+h\) with \(2t+h+1=p\), its gap continuant has the universal factor
\(2t+h+1\), whether or not \(b_t\) vanishes. Dividing by this factor exposes
the first reflection/derivative defect, not a second zero condition.
Consequently a one-orbit fiber supplies no selective Casoratian divisor
after saturation.

Even a hypothetical proof \(|Z_p|\le2\) obtained from the Casoratian would
not compare the selected root at \(p\) with the selected root at \(q\).
The adversarial section (10.11) obeys the same row bound and reflection.

For a top-half target prime, Lucas gives

\[
 b_n\equiv b_1b_{n-p}=5b_{n-p}\equiv0\pmod p.
\tag{10.22}
\]

Thus the target radical divides the fixed integer \(b_n\). This is real
global coupling, but \(\log b_n=\Theta(n)\), exactly the saturated height
scale. The Casoratian formulas audited in Section 4 reconstruct the same
local quotient once its residue is imposed; they do not construct a
second nonzero integer of \(o(n)\) logarithmic height.

There is one exact narrow reopen condition. A non-reflected simultaneous
zero pair \(a,b\in Z_p\) forces \(p\) into the corresponding gap continuant
after all universal reflection factors have been divided out. Hence a
prime-independent collection of such pairs covering every target prime
would give one positive integer divisible by the target radical, with
logarithmic height

\[
 O\!\left(\log n\sum_{\{a,b\}}|a-b|\right).
\]

To make this height \(o(n)\), however, the total gap budget must already be
\(o(n/\log n)\). Minimal one-orbit fibers provide no non-reflected pair,
and bounded fiber size gives neither a cover nor reuse of one continuant by
many primes. This is a precise collision-cover reformulation, not a height
gain from the Casoratian.

The capacity ledger is sharp. A saturated continuant with gap \(h\) has
degree \(3(h-1)\), less only the one universal reflection factor when
applicable, and

\[
 |\Sigma_n(a,h)|\log(n/2)
 \le \log N_h^\sharp(a)
 \le(3(h-1)+O(1))\log n.
\]

It can therefore cover at most \((3+o(1))(h-1)\) target primes. Any cover
of \(K\) targets already has total gap \(\Omega(K)\); the determinant
product is not an amplification mechanism.

For completeness, the same continuant gives a sharp aggregate vertical
repulsion statement.  For \(1\le h<p\), let

\[
 C_p^{\rm full}(h)=
 \#\{0\le r<p-h:r,r+h\in Z_p\}.
\]

Conditional on \(r\in Z_p\), the return equation is equivalent to
\(N_h(r)\equiv0\pmod p\).  The reduced continuant is nonzero and has
degree at most \(3(h-1)\), so

\[
 C_p^{\rm full}(h)\le3(h-1),\qquad
 \sum_{h\le H}C_p^{\rm full}(h)\le\frac32H(H-1).
\tag{VR.1}
\]

This bounds the number of short differences, not the minimum gap.  A
bounded number of gap-two or other fixed-gap pairs is allowed, and two
individually well-spaced rows can still be exact translates at the one
cross-prime shift being tested.  Thus (VR.1) strengthens the vertical
\(p^{2/3}\) theory but supplies no horizontal decorrelation by itself.

#### Quantitative closure of gap-continuant factorial moments

The continuant obstruction can be made quantitative.  This also records the
strongest new result recovered from the uncommitted `uisai2` notes Q910 and
Q913.

Fix the direct prefix window

\[
 n=3H+1,\qquad
 2H<p\le3H+1,\qquad
 Z_p(H)=Z_p\cap[0,H],\qquad z_p(H)=|Z_p(H)|.
\]

For \(h\ge1\), let

\[
 C_p(h)=\#\{0\le s\le H-h:s,s+h\in Z_p(H)\}.
\]

The cleared recurrence gives a positive gap continuant \(N_h(s)\), and two
zero endpoints imply \(p\mid N_h(s)\).  With

\[
 \Lambda_H=\log118+6\log(H+1),
\]

the elementary tridiagonal recurrence gives

\[
 0<N_h(s)\le
 \exp\bigl((h-1)\Lambda_H\bigr).
\tag{CM.1}
\]

Consequently, for every \(L\le H+1\),

\[
 \begin{aligned}
 W_2(H,L)
 &:=
 \sum_{2H<p\le3H+1}
 \sum_{1\le h<L}C_p(h)\log p\\
 &\le
 \sum_{h=2}^{L-1}\sum_{s=0}^{H-h}\log N_h(s)
 \le \frac12HL^2\Lambda_H .
 \end{aligned}
\tag{CM.2}
\]

This is a weighted factorial-moment inequality, not a heuristic.  If
\(8\le B\le H\), choose

\[
 L=\left\lceil\frac{4(H+1)}B\right\rceil.
\]

Partition \([0,H]\) into at most \(B/2\) intervals of length at most \(L\).
Cauchy--Schwarz on the numbers of zeros in these intervals shows that
\(z_p(H)>B\) creates at least \(B/2\) pairs of zeros at distance \(<L\).
Thus (CM.2) gives the explicit tail bound

\[
 \boxed{
 \sum_{\substack{2H<p\le3H+1\\z_p(H)>B}}\log p
 \le
 \frac{25H(H+1)^2\{\log118+6\log(H+1)\}}{B^3}.
 }
\tag{CM.3}
\]

The same inequality remains valid after imposing the moving target condition
\(p\mid b_{3H+1-p}\), since that only restricts the prime set.  At the
available vertical scale \(B=C H^{2/3}\), however, (CM.3) is only

\[
 O_C(H\log H),
\tag{CM.4}
\]

whereas the desired conclusion is \(o(H)\).  It becomes sublinear only when

\[
 \frac{B^3}{H^2\log H}\longrightarrow\infty,
\tag{CM.5}
\]

which lies beyond the proved \(O(H^{2/3})\) row bound.

Higher moments do not remove the logarithm.  Multiplying the \(r-1\)
consecutive gap certificates over an \(r\)-tuple yields, for fixed \(r\ge2\),

\[
 \sum_{2H<p\le3H+1}(z_p(H))_r\log p
 \ll_r H^{r+1}\log H.
\tag{CM.6}
\]

At \(B=H^{2/3}\), ordinary Markov bounds worsen as \(r\) grows.  Localizing
all \(r\) zeros to short blocks is also minimized by the pair case: the
height pays for each possible gap location at exactly the rate at which
the forced multiplicity grows.  Hence (CM.3) closes the entire method that
multiplies individual gap continuants and estimates them by absolute
height.

The precise reopening statement is a joint compression theorem.  For
\(L_H\asymp H^{1/3}\), let \(c_p(H)\) count simultaneous zero pairs with
gap \(<L_H\).  It would suffice to construct one nonzero integer \(C_H\),
independent of \(p\), such that

\[
 v_p(C_H)\ge c_p(H)
 \quad(2H<p\le3H+1),\qquad
 \log|C_H|=o(H^{5/3}).
\tag{CM.7}
\]

Indeed, \(z_p(H)>H^{2/3}\) forces \(c_p(H)\gg H^{2/3}\), so (CM.7) would
turn the high-fiber tail into \(o(H)\).  The naive product of the individual
continuants has logarithmic height \(O(H^{5/3}\log H)\).  Thus the missing
object is a genuinely Apéry-specific common compression identity, not
another Casoratian or another factorial moment.

#### Cross-prime collisions: the degenerate part is now controlled

There is a sharper pointwise decomposition for the top-half arc. Order two
target primes as \(p<\ell\), and write

\[
 h=\ell-p,\qquad s=n-\ell.
\]

The collision conditions are exactly

\[
 \ell=p+h,\qquad
 p\mid b_{s+h},\qquad
 \ell\mid b_s,\qquad
 0\le s<s+h<p<\ell.
\tag{PC.1}
\]

Call the pair **degenerate** if

\[
 p\mid b_s\quad\hbox{or}\quad \ell\mid b_{s+h};
\tag{PC.2}
\]

otherwise call it **pure cross**. The terminology records the essential
field distinction. In a degenerate pair, one of the two characteristics
sees two zero endpoints. In a pure-cross pair, each characteristic sees
only the endpoint belonging to the other modulus.

Let \(M_h^{\rm deg}(n)\) count the degenerate target pairs with prime gap
\(h\). Then the same-prime continuant gives the new unconditional pointwise
estimate

\[
 \boxed{M_h^{\rm deg}(n)\ll h.}
\tag{PC.3}
\]

Indeed, if \(p\mid b_s,b_{s+h}\), then \(p\mid N_h(s)\). Since
\(s=n-h-p\), this implies

\[
 p\mid N_h(n-h).
\]

If instead \(\ell\mid b_s,b_{s+h}\), then \(s=n-\ell\) gives

\[
 \ell\mid N_h(n).
\]

Both continuants are positive nonzero integers and

\[
 \log\{N_h(n-h)N_h(n)\}\ll h\log n.
\tag{PC.4}
\]

Every target prime exceeds \(n/2\), so (PC.3) follows by counting their
distinct large prime divisors. In particular,

\[
 \sum_{h\le H}M_h^{\rm deg}(n)\ll H^2.
\tag{PC.5}
\]

On the pure-cross locus the same continuant is forced to be a unit in both
fields:

\[
 p\nmid N_h(s),\qquad \ell\nmid N_h(s).
\tag{PC.6}
\]

For \(\ell\), divisibility of \(N_h(s)\), together with
\(\ell\mid b_s\), would propagate the zero to \(b_{s+h}\). For \(p\),
divisibility of \(N_h(s)\), together with \(p\mid b_{s+h}\), would force
\(p\mid N_{h-1}(s+1)\). The continuant Dodgson identity has right-hand
side \(-\prod_{j=2}^{h-1}(s+j)^6\), a \(p\)-unit because \(s+h<p\),
and gives the contradiction. Thus the ordinary
Casoratian detector does not merely fail to prove divisibility in the hard
case: its relevant factor is provably nonzero modulo both moving primes.

Writing

\[
 b_{s+h}=pu,\qquad b_s=\ell v
\]

turns divisor switching into the exact equation

\[
 u b_s-vb_{s+h}=huv.
\tag{PC.7}
\]

The new co-divisors \(u,v\) can have exponential height. Equivalently, the
transfer identity leaves the two opposite quotient residues

\[
 b_{s+h}/p\pmod\ell,\qquad b_s/\ell\pmod p
\]

uncontrolled. This is the precise different-characteristic obstruction.
A crossed companion/Casoratian expression can be made divisible by
\(p\ell\), but its low-gap term is accompanied by a multiple of
\(b_sb_{s+h}\); its logarithmic height is \(\Theta(s+h)\), not
\(O(h\operatorname{polylog}n)\).

This decomposition yields a sufficient pointwise theorem strictly weaker
than the shell pair-energy hypothesis. Let \(A_H^\times(n)\) be the number
of **adjacent**, in the ordered target-prime set \(T_n\), pure-cross pairs
with gap at most \(H\). It is enough to prove that, for every fixed \(A>0\),

\[
 \boxed{
  \sup_{N<n\le2N}A_{A\log N}^{\times}(n)
  =o_A(N/\log N).
 }
\tag{PC.8}
\]

To see this, a target set of size \(\gg_\varepsilon n/\log n\) has, after
choosing \(A>2/\varepsilon\), all but \(O(n/(A\log n))\) of its adjacent
gaps at most \(A\log n\). Equation (PC.5) makes the degenerate short gaps
only \(O_A(\log^2n)\). The remaining
\(\gg_\varepsilon n/\log n\) adjacent gaps would be pure cross,
contradicting (PC.8).

An algebraic statement far stronger than necessary, but exactly targeted,
would be a pure-cross gap certificate: for every \(h\le A\log n\),
construct one nonzero integer \(D_{n,h}\), independent of the pair location
\(s\), such that every pure-cross gap-\(h\) pair satisfies

\[
 p(p+h)\mid D_{n,h},\qquad
 \log|D_{n,h}|\ll h(\log n)^C.
\tag{PC.9}
\]

The degenerate certificate \(N_h(n-h)N_h(n)\) has the required height but
misses every pure-cross pair by (PC.6). Conversely, all currently known
certificates which see (PC.1) either depend on \(s\), retain the full
exponential Apéry height, or insert the moving product \(p(p+h)\)
universally. Thus (PC.8), rather than the full second moment, is the
narrowest clean top-half horizontal target exposed so far.

The complete scan through two million shows that this decomposition isolates
the empirical obstruction with no leakage.  The 2,605 exact double-hit
columns and 53 triple-hit columns contain

\[
 2605+3\cdot53=2764
\]

unordered target pairs.  Direct lookup of both cross endpoints in the
binary \(Z_p\) bank found **zero** degenerate pairs: all 2,764 are pure
cross.  Likewise all

\[
 2605+2\cdot53=2711
\]

adjacent target edges are pure cross.  Among those adjacent edges, only one
has \(h\le\log n\), two have \(h\le2\log n\), six have
\(h\le5\log n\), and sixteen have \(h\le10\log n\).  The smallest normalized
gap is

\[
 n=587250,\qquad (p,\ell)=(296507,296519),\qquad
 h/\log n=0.903396\ldots .
\tag{PC.10}
\]

This is finite evidence only, but it confirms that (PC.8) targets exactly
the observed rare event.  The unconditional bound (PC.3) is not being
mistaken for an explanation of the existing double and triple columns:
every one of them lies in the different-characteristic residual.

There is also an unconditional pointwise theorem for each **fixed exact
gap**.  Let
\[
 {\mathfrak S}(h)=
 \prod_{\substack{q\mid h\\q>2}}\frac{q-1}{q-2}
\tag{PC.11}
\]
with the empty product equal to one; an absolute twin-prime constant is
absorbed in the implied constants.  The Selberg upper-bound sieve applied
to \(m(m+h)\) gives, uniformly for even
\(1\le h\le A\log n\),
\[
 M_h^\times(n)
 \le
 \#\{n/2<p\le n-h:p,\ p+h\ {\rm prime}\}
 \ll_A {\mathfrak S}(h)\frac n{\log^2n}.
\tag{PC.12}
\]
For odd \(h\), the set is empty once both primes exceed \(2\).  The
standard maximal-order and mean estimates are
\[
 {\mathfrak S}(h)\ll\log\log(3h),
 \qquad
 \sum_{\substack{h\le H\\2\mid h}}{\mathfrak S}(h)\ll H.
\tag{PC.13}
\]
For completeness, the local sieve density at an odd prime \(q\) is
\(\nu_q=1\) when \(q\mid h\) and \(\nu_q=2\) otherwise.  Thus
\[
 \prod_{q\le z}\left(1-\frac{\nu_q}{q}\right)
 \asymp
 \frac{{\mathfrak S}(h)}{\log^2z}.
\]
The interval remainder for a squarefree sieve modulus is
\(O(\prod_{q\mid d}\nu_q)\), so the standard Selberg weights with
\(z=n^{1/4}(\log n)^{-B}\) give (PC.12), uniformly in the stated
logarithmic range.  The maximal-order estimate in (PC.13) follows by
putting the smallest possible primes into
\[
 {\mathfrak S}(h)=
 \prod_{\substack{q\mid h\\q>2}}
 \left(1+\frac1{q-2}\right).
\]
For the mean estimate, expand this product over squarefree divisors of
\(h\); after summing multiples of each divisor, the remaining Euler
product converges because its \(q\)-factor differs from one by
\(O(q^{-2})\).

This estimate has an exact integer-certificate interpretation.  Define
\[
 D^{\rm pp}_{n,h}=
 \prod_{\substack{n/2<p\le n-h\\p,\ p+h\ {\rm prime}}}p(p+h).
\tag{PC.14}
\]
Every pure-cross target pair of exact gap \(h\) has
\(p(p+h)\mid D^{\rm pp}_{n,h}\), and
\[
 \log D^{\rm pp}_{n,h}
 \ll_A {\mathfrak S}(h)\frac n{\log n}
 \ll_A\frac{n\log\log\log n}{\log n}
 =o_A(n).
\tag{PC.15}
\]
Thus a location-independent sublinear-height carrier already exists for
each exact logarithmic gap.  It is an ambient prime-pair carrier, not an
Apéry compression identity.

The obstruction occurs only when the gaps are combined.  From the mean
bound in (PC.13),
\[
 \sum_{\substack{h\le A\log n\\2\mid h}}
       \log D^{\rm pp}_{n,h}
 \ll_A n.
\tag{PC.16}
\]
This exactly recovers the scale of the adjacent-pair criterion and gives
no little-oh.  It would suffice to improve (PC.12), after summing the
relevant gaps, by any factor \(L(n)\to\infty\) coming from the Apéry
conditions.  The missing horizontal theorem is therefore not a
fixed-gap prime-pair estimate; it is an Apéry-specific saving over the
ambient sieve **across the full logarithmic gap range**.

The same reduction extends without loss to every fixed quotient \(q\).
Define

\[
 T_{n,q}=
 \left\{p:\frac{n}{q+1}<p\le\frac nq,\quad
              p\mid b_{n-qp}\right\}.
 \tag{PC.17}
\]

For \(p<\ell=p+h\) in \(T_{n,q}\), put \(s=n-q\ell\).  The two
conditions become

\[
 p\mid b_{s+qh},\qquad \ell\mid b_s.
 \tag{PC.18}
\]

Declare the pair degenerate if \(p\mid b_s\) or
\(\ell\mid b_{s+qh}\).  The same transfer argument, now across the index
gap \(qh\), shows that every degenerate pair contributes one of

\[
 p\mid N_{qh}(n-qh),\qquad
 \ell\mid N_{qh}(n).
\]

For fixed \(q\),

\[
 \log\{N_{qh}(n-qh)N_{qh}(n)\}\ll_q h\log n.
\]

Both primes are \(\asymp_q n\), and a prime determines its mate once
\(h\) is fixed.  Therefore

\[
 M_{q,h}^{\rm deg}(n)\ll_q h,\qquad
 \sum_{h\le H}M_{q,h}^{\rm deg}(n)\ll_q H^2.
 \tag{PC.19}
\]

On the complementary pure-cross locus, \(N_{qh}(s)\) is a unit modulo
both \(p\) and \(\ell\), by the same transfer and Dodgson proof as for
(PC.6).

Let \(A_{q,H}^{\times}(n)\) count adjacent primes in the ordered set
\(T_{n,q}\) which are pure cross and have prime gap at most \(H\).
Then the following fixed-\(q\) family already implies the full middle-prime
little-oh:

\[
 \boxed{
 \text{for every fixed }q,A>0,\qquad
 \sup_{N<n\le2N}A_{q,A\log N}^{\times}(n)
 =o_{q,A}(N/\log N).
 }
 \tag{PC.20}
\]

Indeed, a failure of the desired logarithmic-weight estimate first
localizes, after discarding \(p\le\eta n\), to one fixed \(q\) with
\(|T_{n,q}|\gg n/\log n\).  The sum of all adjacent prime gaps in that
arc is at most

\[
 \frac nq-\frac{n}{q+1}=\frac{n}{q(q+1)}.
\]

Choosing \(A\) sufficiently large leaves
\(\gg n/\log n\) adjacent gaps at most \(A\log n\); (PC.19) makes only
\(O_{q,A}(\log^2n)\) of them degenerate, contradicting (PC.20).
The order pruning can be incorporated before taking adjacency. Fix any
subpolynomial \(F(n)\to\infty\), let

\[
 T_{n,q}^{\rm np}(F)=
 \left\{p\in T_{n,q}:
 \gcd(p-1,n-q)\le\tau(n-q)\log n\,F(n)\right\},
\]

and let \(A_{q,H}^{{\rm np},\times}(n;F)\) count adjacent pure-cross
pairs in this smaller ordered set. By (KO.6), the still weaker family

\[
 \boxed{
 \text{for every fixed }q,A>0,\qquad
 \sup_{N<n\le2N}A_{q,A\log N}^{{\rm np},\times}(n;F)
 =o_{q,A}(N/\log N)
 }
 \tag{PC.21}
\]

also suffices. A forbidden spike retains linearly many vertices after the
pointwise pruning, and the identical span and degeneracy argument applies
inside the residual set. Thus (PC.21) is one of the two most focused
all-middle horizontal theorems isolated in this audit. It is strictly weaker than
the shell pair-energy condition (10.35): it asks only for short adjacent
pairs, after all non-nearly-primitive vertices have already been discarded.

There is a complementary amplification which trades adjacency for a shared
kernel.  It gives more pairs at a longer, but still subpolynomial, gap.
Fix \(q\), write \(m=n-q\), and for a target prime put

\[
 g_p=\gcd(p-1,m),\qquad p-1=g_pd_p.
\tag{SK.1}
\]

Suppose that a residual \(q\)-arc contains

\[
 K\ge c\,n/\log n
\tag{SK.2}
\]

targets.  Partition the prime interval separately for every divisor
\(g\mid m\) into cells of length \(\Delta\).  There are at most

\[
 B\le\tau(m)
 \left\{\frac{n}{q(q+1)\Delta}+2\right\}
\]

cells.  If their occupancies are \(k_C\), Cauchy--Schwarz gives

\[
 \#\{p<\ell:g_p=g_\ell,\ \ell-p\le\Delta\}
 \ge
 \sum_C\binom{k_C}{2}
 \ge\frac{K^2}{2B}-\frac K2.
\tag{SK.3}
\]

Taking, for a sufficiently large constant depending only on \(c,q\),

\[
 \Delta_{n,q}=C_{c,q}\tau(n-q)\log^2n=n^{o(1)}
\]

forces \(\Omega_{c,q}(n)\) same-kernel close pairs.  The degenerate
pairs among all gaps \(h\le\Delta_{n,q}\) number only

\[
 \sum_{h\le\Delta_{n,q}}M_{q,h}^{\rm deg}(n)
 \ll_q\Delta_{n,q}^2=n^{o(1)}
\]

by (PC.19).  Hence a forbidden spike actually forces
\(\Omega_{c,q}(n)\) **pure-cross, same-kernel** pairs.

Their arithmetic normal form is explicit.  If
\(\ell-p=gt\) and the common kernel is \(g\), then for some \(d,e\)

\[
\begin{aligned}
 m&=g(qd+e),\\
 p&=gd+1,& n-qp&=ge,\\
 \ell&=g(d+t)+1,& n-q\ell&=g(e-qt).
\end{aligned}
\tag{SK.4}
\]

Thus both zero indices and the prime gap carry the same divisor \(g\).
The following is consequently another sufficient all-middle theorem:

\[
 \boxed{
 \begin{aligned}
 \#\{p<\ell:\;&p,\ell\in T_{n,q}^{\rm np}(F),\
 g_p=g_\ell,\ \ell-p\le\Delta_{n,q},\\
 &\text{the pair is pure cross}\}=o_{c,q}(n)
 \end{aligned}
 }
\tag{SK.5}
\]

uniformly for \(N<n\le2N\), for each fixed \(q,c>0\).
One may omit \(c\) by requiring the statement for every fixed constant
in the definition of \(\Delta_{n,q}\).  Conditions (PC.21) and (SK.5)
are not ordered by logical strength.  The former counts only adjacent
\(O(\log n)\)-gap pairs and asks for \(o(n/\log n)\); the latter counts
all same-kernel \(n^{o(1)}\)-gap pairs and asks for \(o(n)\).  The common
kernel in (SK.4) is extra arithmetic structure which may make (SK.5) the
more accessible target.  That possible advantage must not be overstated:
the generic class is \(g=1\), where “same kernel” adds no nontrivial common
modulus or bounded-order motive.  Of the 2,764 observed top-half target
pairs, 1,606 have equal kernels and 826 of those have common kernel one.
Among the 53 triple columns, 21 have one common kernel and 12 have kernel
one at all three primes.  Thus (SK.5) remains genuinely
cross-characteristic even in its arithmetically simplest class.

There is also a hard uncovered component. If

\[
 Z_p=\{r_p,p-1-r_p\},
\]

then \(p\) has no non-reflected zero pair at all, so no saturated
collision-cover factor is guaranteed to contain it. In the scan through
two million, 45,020 of the 58,555 nonempty rows are such exact doublets,
about 76.9 percent of the active vertical population. This does not prove
the same proportion in a fixed target column, but it shows that the
primitive obstruction is structurally dominant, not exceptional.

Consequently a collision-cover proof must split into two genuinely new
horizontal assertions:

1. the total logarithmic weight of target primes with exact-doublet fibers
   is \(o(n)\);
2. the remaining rich-fiber targets admit a non-reflection cover with
   total gap \(o(n/\log n)\).

Neither assertion follows from the current zero-fiber histogram.

#### Route B: defect dimension

A bound \(|Z_p|\le D\) reduces the number of choices at each prime, but it
does not reduce the number of prime coordinates. The local defect spaces
still form a direct product over \(p\), and a section may choose one of the
\(D\) roots in every row. Varying the local level introduces new local
coordinates just as in Section 5.

Dimension growth becomes relevant only after producing one of:

1. a common rational coordinate in which roots for different primes are
   reductions of bounded-height global points;
2. a nonzero eliminated integer whose target-prime valuation accumulates;
3. a horizontal equidistribution theorem showing that the CRT section has
   no exceptional low-height lift.

None follows from the row dimension. Hence the new computation reinforces,
rather than changes, the earlier verdict: local rank is not global
separation.

The V1 jet tower has the same exact defect. At derivative order \(j\), once
all lower jets are fixed, the highest jet satisfies an inhomogeneous copy of
the original second-order Apéry operator. Its solution space is therefore
an affine torsor under a fresh two-dimensional homogeneous state. Each new
equation arrives with new local coordinates; reflection leaves one free
scalar at every order. Thus higher jets identify successive quotient and
reflection-defect coordinates but do not create growing codimension. This
is the structural reason the speculative dimension table in Route B cannot
be repaired merely by adding more formal defect levels.

#### An exact fixed-\(n\) carrier, but no new height saving

There is a useful deterministic reformulation of the top-half branch. Set

\[
 C_n=\binom n{\lfloor n/2\rfloor}.
\]

For a prime \(p\in(n/2,n]\), the numerator \(n!\) contains exactly one
factor \(p\). Both denominator arguments are below \(p\), except when
\(n=2p-1\), where \(\lceil n/2\rceil=p\) cancels it. Consequently

\[
 v_p(C_n)=
 \begin{cases}
  1,&n\ne2p-1,\\
  0,&n=2p-1.
 \end{cases}
\]

The boundary prime is never a target for \(p\ge7\), since its residue is
\(p-1\) and

\[
 b_{p-1}\equiv1\pmod {p^2}.
\]

The stronger congruence is elementary: for \(1\le k\le p-1\), the factor
\(\binom{p-1+k}{k}\) contains exactly one \(p\), so every nonzero-\(k\)
summand in \(b_{p-1}\) contains \(p^2\), while the \(k=0\) summand is 1.

For \(n>10\), Gessel--Lucas also gives

\[
 b_n\equiv b_1b_{n-p}=5b_{n-p}\pmod p.
\]

It follows in the stronger valuation-by-valuation form that, for every
prime \(p\in(n/2,n]\),

\[
 v_p\!\left(\gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right)\right)
 =
 \mathbf1_{\{p\mid b_{n-p}\}}.
\]

Thus the interval part of the gcd is already squarefree. With
\(\operatorname{rad}_{(n/2,n]}\) denoting the radical restricted to primes
in that interval,

\[
 \boxed{
  \prod_{p\in T_n}p
  =
  \operatorname{rad}_{(n/2,n]}
  \gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right).
 }
\]

This exact identity was already isolated in Section 37 of
`Q32_CODEX_RESUME_2026-07-23.md`; an independent recurrence computation
rechecked it for every \(11\le n\le1000\). It is conceptually cleaner than
the interpolation carriers, but it is saturated:

\[
 \log C_n=(\log2)n+O(\log n),
\]

and every nonboundary candidate prime divides \(C_n\) independently of the
Apéry-zero condition. In fact its interval part is exactly the squarefree
primorial over all nonboundary candidate primes before intersecting with
\(b_n\); the binomial factor has no local selectivity. A subexponential
bound for the displayed gcd would
solve the top-half branch, but no known gcd theorem for fixed-dimensional
\(S\)-unit or constant-coefficient recurrence orbits applies to the pair
consisting of an Apéry \(P\)-recursive term and a factorial ratio. The
carrier is therefore an exact restatement, not the missing second
Casoratian observable.

#### Path D: the holonomic-gcd formulation is exact, but no theorem covers it

The carrier isolates a concise alternative target:

\[
 \log\operatorname{rad}_{(n/2,n]}
 \gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right)=o(n).
\tag{HD.1}
\]

The stronger assertion

\[
 \log\gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right)=o(n)
\tag{HD.2}
\]

would certainly suffice for the top-half branch, but (HD.1) is all that
branch asks from this carrier.  The extra small-prime valuations in (HD.2) are
automatically harmless. Kummer's theorem gives

\[
 v_p\binom n{\lfloor n/2\rfloor}
 \le \lfloor\log_p n\rfloor+1,
\]

and therefore

\[
 \sum_{p\le\sqrt n}
 v_p\binom n{\lfloor n/2\rfloor}\log p
 \le\pi(\sqrt n)\log(n+1)=O(\sqrt n).
\tag{HD.2a}
\]

For \(p>\sqrt n\), the binomial valuation is zero or one.  Thus the
difference between (HD.1) and (HD.2) is additional middle-prime support,
not a hidden linear contribution from powers of small primes.

In fact that additional support has an exact quotient-parity description.
Write \(n=qp+r\), \(0\le r<p\).  Legendre's formula gives

\[
 v_p\binom n{\lfloor n/2\rfloor}
 =
 \begin{cases}
  1,&q\ \text{odd and }r\le p-2,\\
  0,&q\ \text{even, or }r=p-1.
 \end{cases}
\tag{HD.2b}
\]

For \(q=2a\), both halves of \(n\) contain \(a\) copies of \(p\).  For
\(q=2a+1\), both contain only \(a\), except that the upper half contains
one more at the boundary \(r=p-1\).  Since that boundary is never an
Apéry target, the central carrier contains every odd-\(q\) middle target
simultaneously.  It may also contain primes coming from
\(p\mid b_q\), so this is a support inclusion, not an unconditional
equality for growing \(q\).  The even quotient arcs still require their
own carriers.

There is a useful fixed-denominator generalization.  For \(d\ge2\), put

\[
 B_{n,d}=\binom n{\lfloor n/d\rfloor}.
\]

Write \(q=dt+u\), \(0\le u<d\), and

\[
 \left\lfloor\frac nd\right\rfloor
 =tp+\left\lfloor\frac{up+r}{d}\right\rfloor.
\]

Legendre's formula then gives one carry precisely when

\[
 \left\lfloor\frac{up+r}{d}\right\rfloor>r.
\tag{HD.3}
\]

If \(q\equiv-1\pmod d\), so \(u=d-1\), condition (HD.3) holds for every
\(r\le p-2\) and fails only at \(r=p-1\).  Consequently every actual
target in every quotient arc

\[
 q\equiv-1\pmod d
\]

divides \(B_{n,d}\) exactly once.  The central-binomial parity law is the
case \(d=2\).

This yields a particularly clean sufficient Path D family:

\[
 \boxed{
 \text{for every prime }d,\qquad
 \log\gcd(b_n,B_{n,d})=o_d(n).
 }
\tag{HD.4}
\]

Indeed, a forbidden spike first localizes to one fixed quotient \(q\).
Choose any prime divisor \(d\mid q+1\); every target in that arc then
divides the gcd in (HD.4).  No uniformity in \(d\) is needed.  The carrier
height remains linear for fixed \(d\),

\[
 \log B_{n,d}=n\,h(1/d)+O_d(\log n),
\]

so (HD.4) is still a genuinely new gcd theorem, but it organizes all
quotient arcs into fixed factorial-ratio families rather than one carrier
per \(q\).

Equivalently, the fixed-\(q\) carriers in (10.32) give the
interval-restricted estimates

\[
 \log\operatorname{rad}_{(n/(q+1),\,n/q]}
 \gcd(b_n,C_{n,q})=o_q(n)
\tag{HD.5}
\]

for every fixed \(q\).  These turn the
all-middle problem into a family of interval-restricted gcd estimates
with no claim about irrelevant prime factors outside the indicated arc.

All sequences in (HD.1) and (HD.4) are holonomic: \(b_n\) is the Apéry
polynomial-coefficient recurrence, while for fixed \(d\) the \(d\)
residue-class subsequences of \(B_{n,d}\) are hypergeometric. This observation
does not put (HD.1) under an existing gcd theorem.  The moving-target
Subspace-Theorem results checked here, including Grieve--Wang, treat
polynomials evaluated on \(S\)-unit points and their applications to
algebraic **constant-coefficient** linear recurrences.  An Apéry
\(P\)-recursive orbit has a transition matrix depending on \(n\), and is
not a fixed-dimensional \(S\)-unit orbit.

Several tempting formal substitutes also fail:

1. coprimality of the two Ore annihilators is a characteristic-zero
   statement about solution spaces, not a bound for
   \(\gcd(b_n,C_n)\) at one moving index.  Even the sequences
   \(u_n=2^n\) and \(v_n=6^n\) have distinct coprime first-order
   annihilators and multiplicatively independent dominant bases, while
   \(\gcd(u_n,v_n)=2^n\);
2. multiplicative independence of the dominant Archimedean bases
   \(17+12\sqrt2\) and \(2\) gives no control of exact common prime
   divisors of their coefficient sequences;
3. the binomial sequence was chosen precisely because every top-half
   candidate prime divides it.  Its holonomic recurrence supplies no
   additional selective congruence at those primes;
4. a one-index common zero does not give the repeated shifted zeros needed
   for a recurrence resultant or Casoratian to accumulate valuation.

A targeted literature search on \(P\)-recursive, holonomic, D-finite, and
\(G\)-function coefficient gcds found definitions, representation and
asymptotic theorems, but no result implying (HD.1), (HD.4), or (HD.5). Path D is
therefore a legitimate new theorem interface, not an application waiting
only for citation.  To advance it one needs either a new arithmetic gcd
theorem for a \(G\)-function coefficient and a factorial ratio, or an
Apéry-specific low-height certificate which is selective inside the
universal binomial support.

The full-gcd strengthening (HD.2) is at least numerically plausible.  An
exact recurrence computation through \(n=10{,}000\) gives the following
maxima of

\[
 \frac1n\log\gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right)
\]

on successive blocks:

| \(n\)-block | maximum | attained at |
|---:|---:|---:|
| \(501\)--\(1{,}000\) | \(0.0326690\) | \(676\) |
| \(1{,}001\)--\(2{,}000\) | \(0.0249751\) | \(1{,}041\) |
| \(2{,}001\)--\(5{,}000\) | \(0.0133423\) | \(2{,}005\) |
| \(5{,}001\)--\(10{,}000\) | \(0.00547604\) | \(5{,}241\) |

At the last maximum the exact gcd is \(2{,}912{,}053{,}021{,}625\).
This is finite evidence only.  In particular it cannot substitute for a
theorem controlling exceptional \(n\), the precise pointwise issue of
Problem 3.2.

### 10.8 The corrected sufficient condition: cross-prime collision energy

There is a clean theorem-shaped replacement for (10.2).

On a dyadic shell, define

\[
 K_b(n)=
 \#\{\sqrt n<p\le n:n\bmod p\in Z_p\}
\]

as in (10.16), and its cross-prime falling second moment

\[
 \mathcal C(N)=
 \sum_{N<n\le2N}(K_b(n))_2,
\qquad (x)_2=x(x-1).
\tag{10.23}
\]

If

\[
 \mathcal C(N)=o\!\left(\frac{N^2}{\log^2N}\right),
\tag{10.24}
\]

then for every \(n\in(N,2N]\),

\[
 (K_b(n))_2\le\mathcal C(N),
\]

and therefore

\[
 \max_{N<n\le2N}K_b(n)=o(N/\log N).
\tag{10.25}
\]

Since every middle prime has logarithm \(O(\log N)\), (10.25) gives

\[
 \sum_{\substack{\sqrt n<p\le n\\n\bmod p\in Z_p}}\log p=o(N).
\tag{10.26}
\]

Together with the already controlled small-prime and companion channels,
(10.26) proves Problem 3.2: the Wronskian bound gives
\(v_p(G_n)\le6\) for \(p>\sqrt n\), so passing from the support radical to
the prime-power contribution costs only this absolute factor.

The probabilistic data predict the much stronger natural bound

\[
 \mathcal C(N)\ll N^{1+o(1)}.
\tag{10.27}
\]

Indeed, \(K_b(n)\) has measured mean and variance close to \(\log2\).
Bound (10.27) would give

\[
 K_b(n)\ll N^{1/2+o(1)}
\]

uniformly, already far more than needed.

The best immediate deterministic estimate from the proved vertical bound is
much larger. Dropping the moving endpoint restrictions,

\[
\begin{aligned}
 \sum_{N<n\le2N}K_b(n)
 &\le
 \sum_{\sqrt N<p\le2N}|Z_p|\left(\frac Np+1\right)\\
 &\ll \frac{N^{5/3}}{\log N}.
\end{aligned}
\]

Together with \(\max K_b(n)\ll N/\log N\), this gives only

\[
 \mathcal C(N)
 \ll\frac{N^{8/3}}{\log^2N},
\]

missing (10.24) by the factor \(N^{2/3}\). A black-box additive large sieve
is weaker: if
\(\mu_N=\sum_{\sqrt N<p\le2N}|Z_p|/p\), its natural output is

\[
 \mathcal C(N)\ll N\mu_N^2+N^2\mu_N
 \ll\frac{N^{8/3}}{\log N}.
\]

The \(N^2\) term reflects Farey spacing \(N^{-2}\) while the sampled
\(n\)-interval has length only \(N\). Even a hypothetical
\(|Z_p|=O(1)\) gives only \(\mu_N=O(1)\) and the black-box bound
\(\mathcal C(N)=O(N^2)\), still short of the required logarithmic
little-oh. The needed gain is therefore genuinely Apéry-specific, not a
formal consequence of ordinary additive large-sieve spacing.

Expanding (10.23) shows the exact arithmetic content. It counts solutions
with distinct primes \(p\ne q\) to

\[
 n=ap+r=bq+s,\qquad
 r\in Z_p,\quad s\in Z_q,\quad N<n\le2N,
\tag{10.28}
\]

under the middle-prime restrictions. For the top-half arc, \(a=b=1\), so
the core becomes the fixed-sum collision

\[
 p+r=q+s=n.
\tag{10.29}
\]

For later computations, the exact dyadic window must be retained. If
\(K_1(n)\) is the top-half count and \(p<\ell\), put

\[
 J_N(p,\ell)=
 \left[
  \max\{\ell-p,N-p+1,0\},\
  \min\{p-1,2N-p\}
 \right]\cap\mathbf Z.
\]

Then

\[
\begin{aligned}
 \sum_{N<n\le2N}\binom{K_1(n)}2
 =
 \sum_{\substack{N/2<p<\ell\le2N\\p,\ell\ {\rm prime}}}
 \sum_{r\in J_N(p,\ell)}
 \mathbf1_{Z_p}(r)
 \mathbf1_{Z_\ell}\bigl(r-(\ell-p)\bigr).
\end{aligned}
\tag{DE.1}
\]

Thus an independence benchmark for one prime pair is

\[
 |J_N(p,\ell)|\,\frac{|Z_p|}{p}\frac{|Z_\ell|}{\ell},
\tag{DE.2}
\]

not \(|Z_p||Z_\ell|/\ell\) with the overlap window silently replaced by
its maximal length. This correction changes finite pair statistics near
the shell boundary, but not the theorem-shaped target (10.24).

Equivalently, one needs additive-energy control for the family of translated
sets

\[
 p+Z_p\subset[p,2p).
\]

A row bound alone permits all these translated sets to pass through one
common \(n_0\), making \(\mathcal C(N)\) as large as the failure scale in
(10.24). Thus (10.24), unlike \(|Z_p|\le D\), explicitly excludes the
adversarial alignment.

For a dyadic prime block \(X<p\le2X\), put

\[
 K_X(m)=\#\{X<p\le2X:m\bmod p\in Z_p\}.
\]

On the full CRT representative range \(0\le m<4X^2\), the corresponding
higher-moment formulation is

\[
 \sum_{m<O(X^2)}(K_X(m))_k
 \ll X^{2+o(1)}\lambda_X^k,
\qquad
 \lambda_X=\sum_{X<p\le2X}\frac{|Z_p|}{p}.
\tag{10.30}
\]

With the proven \(|Z_p|\ll p^{2/3}\), any fixed \(k>6\) at the
independence scale in (10.30) gives a power-saving pointwise bound. If a
uniform constant row bound were also known, a much lower moment would
suffice. But in either case the horizontal moment, not the row bound, is
the decisive new input.

#### A weaker sufficient theorem after quotient localization

The all-middle estimate (10.24) is convenient but stronger than necessary.
For a fixed positive integer \(q\), define the quotient-arc count

\[
 K_q(n)=
 \#\left\{
  p\text{ prime}:\frac n{q+1}<p\le\frac nq,\quad n-qp\in Z_p
 \right\}
\]

There is an exact factorial carrier on every fixed quotient arc, extending
the central-binomial identity above. Define

\[
 C_{n,q}=\binom n{\lfloor n/(q+1)\rfloor}.
\]

If \(n=qp+r<(q+1)p\), \(p>\sqrt n\), and \(r\le p-2\), then
\[
 r<\left\lfloor\frac n{q+1}\right\rfloor<p.
\]
Kummer's theorem therefore gives exactly one base-\(p\) carry and

\[
 v_p(C_{n,q})=1.
\tag{10.31}
\]

At the only excluded boundary \(r=p-1\), there is no carry, but also
\(b_r\equiv1\pmod p\). On the other hand, Gessel--Lucas gives

\[
 b_n\equiv b_qb_r\pmod p.
\]

For fixed \(q\), all sufficiently large primes in the arc are larger than
the fixed integer \(b_q\). Hence, for all sufficiently large \(n\),

\[
 \prod_{\substack{n/(q+1)<p\le n/q\\b_{n-qp}\equiv0\pmod p}}p
 =
 \operatorname{rad}_{(n/(q+1),\,n/q]}
 \gcd(b_n,C_{n,q}).
\tag{10.32}
\]

This was independently checked in 197,070 admissible prime/index cases for
\(1\le q\le7\) and \(n\le2000\), after explicitly excluding
\(p\mid b_q\). The proof, not the computation, establishes (10.32).

The carrier height is

\[
 \log C_{n,q}
 =
 n\,h\!\left(\frac1{q+1}\right)+O(\log n),
\]

where \(h(x)=-x\log x-(1-x)\log(1-x)\). It has positive linear rate for
every fixed \(q\), but rate \(O((\log q)/q)\) as \(q\to\infty\). This is
another exact explanation for the quotient localization: growing-\(q\)
arcs are cheap, while every fixed-\(q\) gcd remains a genuinely new
horizontal theorem.

Now define the shell pair energy

\[
 \mathcal E_q(N)=
 \sum_{N<n\le2N}(K_q(n))_2.
\tag{10.33}
\]

Every collision counted here has the exact form

\[
 n=qp+r=q\ell+s,\qquad
 r\in Z_p,\quad s\in Z_\ell,
\]

or equivalently

\[
 r-s=q(\ell-p).
\tag{10.34}
\]

The following family of estimates suffices for the full lower channel:

\[
 \boxed{
  \text{for every fixed }q,\qquad
  \mathcal E_q(N)
  =o_q\!\left(\frac{N^2}{\log^2N}\right).
 }
\tag{10.35}
\]

Here the quantifiers are important: (10.35) is required on every
sufficiently large dyadic shell, separately for each fixed \(q\).

To prove sufficiency, suppose instead that along an infinite sequence

\[
 \sum_{\substack{\sqrt n<p\le n\\n\bmod p\in Z_p}}\log p
 \ge\varepsilon n.
\]

Choose \(\eta=\eta(\varepsilon)>0\) so small that Chebyshev's bound for
\(\vartheta(\eta n)\) makes all primes \(p\le\eta n\) contribute at most
\(\varepsilon n/2\), even without imposing the zero condition. The
remaining hits have \(q=\lfloor n/p\rfloor<1/\eta\), so they lie in only
finitely many fixed quotient arcs. They include
\(\gg_\varepsilon n/\log n\) distinct primes. By pigeonhole and then an
infinite subsequence, one fixed \(q\) satisfies

\[
 K_q(n)\gg_\varepsilon n/\log n.
\]

The single summand at that \(n\) then forces

\[
 \mathcal E_q(N)\gg_\varepsilon N^2/\log^2N
\]

on its dyadic shell, contradicting (10.35). In particular, the much rougher
bound

\[
 \mathcal E_q(N)\ll_q N^{2-\delta_q}
\]

for any \(\delta_q>0\) would be enough.

This argument does **not** assert
\(\mathcal C(N)=\sum_q\mathcal E_q(N)\); that identity is false because
\(\mathcal C(N)\) also contains cross-\(q\) pairs. Instead it gives a
different \(L^\infty\) route: first discard \(p\le\eta n\), then bound the
finitely many surviving quotient arcs separately. No uniformity in \(q\)
is required, since \(q\) is frozen only after \(\varepsilon\) and \(\eta\)
are fixed.

The exact same-\(q\) equation (10.34) does improve the completely
unrestricted CRT count: for fixed distinct \(p,\ell\), choosing
\(r\in Z_p\) determines \(s\) and \(n\), so there are at most
\(\min(|Z_p|,|Z_\ell|)\) collisions. With only
\(|Z_p|\ll p^{2/3}\), however, this yields merely

\[
 \mathcal E_q(N)\ll_q\frac{N^{8/3}}{\log^2N},
\]

again missing (10.35) by \(N^{2/3}\). Reflection changes only constants.
The saving must decide whether the one marked residue determined modulo
\(\ell\) is actually in \(Z_\ell\); present same-prime continuants do not.

The binary scan gives direct evidence for exactly (10.33), not just for
one-point sparsity. On \(1,000,000<n\le2,000,000\), with \(p\ge7\), the
first four fixed arcs give:

| \(q\) | \(N^{-1}\sum_nK_q(n)\) | \(N^{-1}\mathcal E_q(N)\) | square of the mean | \(\max K_q(n)\) |
|---:|---:|---:|---:|---:|
| 1 | 0.050072 | 0.002396 | 0.002507 | 3 |
| 2 | 0.030336 | 0.000946 | 0.000920 | 3 |
| 3 | 0.022023 | 0.000550 | 0.000485 | 3 |
| 4 | 0.017400 | 0.000328 | 0.000303 | 3 |

The natural prediction for fixed \(q\) is

\[
 \mathbf E K_q(n)
 \sim
 \log\!\left(
  \frac{\log(n/q)}{\log(n/(q+1))}
 \right)
 \sim\frac{\log(1+1/q)}{\log n},
\]

and hence \(\mathcal E_q(N)=N^{1+o(1)}/\log^2N\). Again this is evidence,
not the needed deterministic bound.

#### Moment-order bookkeeping

There are two different higher-moment thresholds. Formula (10.30) uses the
entire CRT representative range \(m<O(X^2)\). With
\(|Z_p|\ll p^{2/3}\), its interval length \(X^2\) explains the condition
\(k>6\).

For the actual dyadic \(n\)-shell, one may localize simultaneously in
\(N<n\le2N\) and \(X<p\le2X\). Put

\[
 K_{N,X}(n)=
 \#\{X<p\le2X:\sqrt n<p\le n,\ n\bmod p\in Z_p\}
\]

and

\[
 A_X=\sum_{X<p\le2X}\frac{|Z_p|}{p}
 \ll\frac{X^{2/3}}{\log X}.
\]

An independence-scale bound

\[
 \sum_{N<n\le2N}(K_{N,X}(n))_k
 \ll N^{1+o(1)}A_X^k
\tag{10.36}
\]

would give

\[
 \max K_X(n)
 \ll N^{1/k+o(1)}\frac{X^{2/3}}{\log X}.
\]

After restoring the logarithmic prime weight, the worst block is \(X\asymp
N\), with exponent \(2/3+1/k\). Thus \(k=4\) is the first integer moment
which closes using only the currently proved vertical exponent:

\[
 \frac23+\frac14=\frac{11}{12}<1.
\]

This fourth-moment route is stronger than the fixed-\(q\) criteria.  From
weakest to strongest, the clean sufficient hypotheses isolated above are:

1. nearly primitive, adjacent pure-cross logarithmic-gap depletion
   (PC.21);
2. nearly primitive fixed-\(q\) pair energy (KO.16);
3. full fixed-\(q\) pair energy (10.35);
4. with only \(|Z_p|\ll p^{2/3}\), the relative fourth moment (10.36).

### 10.9 Precise live route and corrected verdict

The two most focused genuinely live horizontal formulations exposed by the
computation are:

1. for every fixed quotient \(q\), prove (PC.21): after the pointwise
   Kummer-order pruning, the number of adjacent pure-cross target primes
   with prime gap \(O(\log n)\) is \(o(n/\log n)\);
2. exploit the extra common divisor in (SK.4) and prove (SK.5):
   same-kernel pure-cross pairs with subpolynomial gap are \(o(n)\).

The more standard but stronger alternatives are the residual pair energy
(KO.16), the full fixed-\(q\) pair energy (10.35), or the global
fixed-sum estimates (10.24)/(10.30).  By (KO.6), all of them only have to
treat the nearly primitive
sector

\[
 \gcd(p-1,n-q)\le\tau(n-q)\log n\,F(n),
 \qquad q=\lfloor n/p\rfloor,
\]

for any chosen subpolynomial \(F(n)\to\infty\). Thus fixed-order modular
forms and bounded-order Kummer motives can be removed unconditionally before
any collision estimate is attempted.

Path D, the holonomic-gcd family (HD.4), is an alternative formulation
rather than a collision hypothesis.  It is exact but presently requires a
new arithmetic theorem outside the scope of known constant-coefficient
recurrence gcd results.

This requires an Apéry-specific horizontal mechanism, for example:

- a reciprocity law turning \(p\mid b_{n-p}\) into divisibility of a
  bounded-height integer shared across different \(p\);
- a fixed-\((q,h)\) pure-cross certificate of sublinear height, extending
  (PC.9) beyond the degenerate continuant;
- a bounded-complexity algebraic parametrization of the off-center root
  position as \(p\) varies;
- a dispersion theorem for the moving character-index Hasse interpolation;
- the global sublinear-height certificate already isolated at the end of
  Section 8.

No such mechanism is presently supplied by the recurrence Casoratian,
finite-level defect dimensions, K3 generic ordinarity, or the central
eta-product.

The revised assessment is therefore:

1. the small zero fibers and target counts are high-quality evidence for
   the conjecture;
2. the claimed evenness has two central exceptions and is completely
   explained by elementary reflection;
3. the exact scan through two million raises the zero-fiber maximum from
   8 to 12 and the full-middle column maximum to 7;
4. the observed Poisson law predicts slow unbounded maxima, not absolute
   constants;
5. even a proved absolute bound \(|Z_p|\le D\) would strengthen averaged
   results but would **not** solve the pointwise problem without an
   additional horizontal collision estimate;
6. Routes A and B remain saturated for exactly the same reason as before:
   neither creates cross-prime low-height separation.
7. The divisor-sensitive order sieve is a genuine unconditional gain: it
   reduces the remaining Mellin problem to characters of order
   \(p^{1-o(1)}\), but leaves the same horizontal alignment obstruction in
   that sector.
8. Degenerate short-gap collisions are controlled pointwise by
   \(M_{q,h}^{\rm deg}(n)\ll_qh\), but all 2,764 observed top-half pairs
   are in the complementary pure-cross locus, where the continuant is a
   unit in both characteristics.
9. The exact remaining theorem can be weakened from global pair energy to
   either (PC.21), an adjacent logarithmic-gap statement inside the nearly
   primitive residual, or the structured same-kernel statement (SK.5).
   The prime-center reflection star (KO.17) shows that all currently known
   local constraints remain compatible with worst-case alignment there.
10. The divided-residue coherence (A.5) is exact but tautological for an
    arbitrary common quotient \(A\).  It becomes useful only if a new
    many-prime jet elimination turns it into a nonzero sublinear-height
    certificate.
11. The Lucas collapse gives the exact holonomic-gcd interfaces
    (HD.1), (HD.4), and (HD.5), but
    no checked theorem for \(P\)-recursive coefficients proves them.

This is a negative answer to the proposed shortcut, but a positive
sharpening of the frontier: the missing theorem is now the explicit
pure-cross estimate (PC.21) or (SK.5), or a new selective
holonomic-gcd mechanism such as (HD.4), not a vertical zero-fiber bound.

### 10.10 Sources checked for the geometry and congruences

- S. Ahlgren and K. Ono, *A Gaussian hypergeometric series evaluation and
  Apéry number congruences*, J. reine angew. Math. 518 (2000), 187--212,
  <https://doi.org/10.1515/crll.2000.004>.
- F. Bogomolov and Y. Zarhin, *Ordinary reduction of K3 surfaces*,
  <https://arxiv.org/abs/0902.1548>.
- A. Malik and A. Straub, *Divisibility properties of sporadic Apéry-like
  numbers*, <https://arxiv.org/abs/1508.00297>.
- X. Caruso, F. Fürnsinn, D. Vargas-Montoya, and W. Zudilin,
  *Galois groups of Apéry-like series modulo primes*, Theorems 1--2,
  <https://arxiv.org/abs/2510.23298>,
  <https://doi.org/10.1017/S0004972725100932>.
- N. Grieve and J. T.-Y. Wang, *Greatest common divisors with moving
  targets and consequences for linear recurrence sequences*, Theorem 1.2,
  <https://arxiv.org/abs/1902.09109>.

These sources support, respectively, the central modular congruence, the
fixed-K3 ordinary-reduction theorem, the Lucas/divisibility framework, the
square/quadratic-times-square factorization of the reduced generating
series, and the scope of moving-target gcd theorems. None states a uniform
bound for the coefficient-index zero fibers or a pointwise horizontal
collision theorem of the form (10.24) or (10.35).

### 10.11 Corrections required in earlier working notes

This audit supersedes several pointwise claims elsewhere in the research
directory.

1. Any table claiming that a bound on \(|Z_p|\) **alone**, even
   \(|Z_p|=O(1)\), implies the all-\(n\) version of Problem 3.2 is invalid.
   Such a bound supports first-moment and density statements, but the
   adversarial construction (10.11) disproves the pointwise implication.
2. In particular, the hierarchy in `Q3.2_paper_draft.md` which lists
   \(|Z_p|=O(p^{1-\delta})\), \(O(\log p)\), or \(O(1)\) as standalone
   pointwise sufficient conditions must be amended by a horizontal
   hypothesis such as (10.24) or (10.30).
3. The maximum-zero claim 16 through 80,000 in
   `Q3.2_density_theorem.tex`/`.md` must be replaced by 8. The number 16
   belongs to \(\max_aN_p(a)\), not \(N_p(0)\).
4. The statement in `Q32_MAIN_THEOREM_writeup.md` that the every-\(n\)
   residual is the density-zero problem for the **fixed** weight-four
   eta-product is not a valid reduction for the off-center roots. The
   later `Q32_CODEX_RESUME_2026-07-23.md` correctly repairs this: the fixed
   form controls only the central residue, whose pointwise weight is already
   \(O(\log n)\) by (10.9).
5. The abstract of `proof.tex` says that the Poisson model predicts
   \(\max_{p\le N}|Z_p|\sim2\log\log N\). This is not the extreme-value law
   for a fixed-mean Poisson sample. Since \(|Z_p|/2\) is modeled by
   \(\operatorname{Poisson}(1/2)\) over \(\pi(N)\) rows, the first-order
   prediction is
   \[
    \max_{p\le N}|Z_p|
    \sim\frac{2\log N}{\log\log N},
   \]
   up to the usual replacement of \(N\) by \(\pi(N)\) inside the
   logarithms.
6. Corollary `cor:divides` in `proof.tex` claims
   \(v_p(G_n)\le3\lfloor\log_p n\rfloor\) from
   \(G_n\mid d_n(a_nb_{n-1}-a_{n-1}b_n)\). The displayed divisibility is
   not justified: in the attempted Bézout combination the coefficient
   \(a_{n-1}\) is generally rational. The denominator-safe combination uses
   \(d_n^2W_{n-1}\) and proves
   \[
    v_p(G_n)\le6\lfloor\log_p n\rfloor-3v_p(n)
   \]
   for \(p\ge5\). Thus the safe middle-prime depth is 6. All density and
   pointwise-support reductions remain valid after multiplying their
   absolute constants by at most two.

The unconditional density theorem and the proved local
\(|Z_p|\ll p^{2/3}\) bound are not affected by these corrections. What is
removed is only the unsupported passage from vertical sparsity to a
pointwise horizontal bound.

## 11. Recovered endpoint tower: an unconditional fixed law through order seven

The source audit in Section 5 remains correct about the cited Liu and Sun
papers: those papers do not supply the proposed free Bernoulli defect
dimensions.  Subsequent work in a separate UIS research tree did, however,
produce a different and mathematically substantive result.  The complete
dependency chain has now been synchronized into this repository:

```text
Q32_TARGET_CUBIC_BLOCK_LAW_2026-07-28.md
Q32_ENDPOINT_BERNOULLI_RANK_ONE_2026-07-29.md
Q32_HARMONIC_LEMMA_PROOF_2026-07-29.md
Q32_QUARTIC_TARGET_BERNOULLI_ELIMINATION_2026-07-29.md
Q32_WEIGHT_FIVE_ENDPOINT_AND_SEXTIC_TARGET_2026-07-29.md
Q32_H6_STAR_PROOF_2026-07-29.md
Q32_WEIGHT_SEVEN_ENDPOINT_RANK_ONE_2026-07-29.md
Q32_ORDER_SEVEN_TARGET_FIXED_LAW_2026-07-29.md
```

The associated independent audit scripts are in `../scripts/` with the
prefix `q32_`.  They have been rerun through prime \(1000\); the cubic,
quartic/quintic, sextic, seventh-endpoint, H6, and final nonlinear target
checks all report zero failures.  This computation checks the displayed
identities but is not used in place of the proofs.

### 11.1 The exact chain of common quotient digits

Let \(p\ge7\), \(0\le r<p\), \(p\mid b_r\), and put

\[
 n=p+r,\qquad s=p-1-r,\qquad x_p=\frac{b_n}{p}.
\]

The shifted Apéry recurrence has exact rational fundamental solutions
\({\cal U}_r(X),{\cal V}_r(X)\).  Truncating their Taylor series at a fixed
order is legitimate because every denominator occurring for \(r<p\) is a
\(p\)-unit.  The direct and reflected decompositions are

\[
\begin{aligned}
 b_{p+r}
 &=b_p{\cal U}_r(p)-p^3b_{p-1}{\cal J}_r(p),\\
 b_{2p-1-s}
 &=b_{2p-1}{\cal U}_s(-2p)+8p^3b_{2p}{\cal J}_s(-2p),
\end{aligned}
\tag{11.1}
\]

where \({\cal J}_r(X)={\cal V}_r(X)/(1+X)^3\).  The endpoint defects in
\(b_p,b_{p-1},b_{2p},b_{2p-1}\) are then removed by fixed direct/reflected
combinations.

The resulting progression is:

| local precision for \(x_p\) | input |
|---|---|
| \(p,p^2,p^3\) | Gessel--Lucas and the universal quadratic/cubic shifted jets |
| \(p^4,p^5\) | the rank-one endpoint coordinate \(\Delta_p=b_{p-1}-1\), eliminated with weights \(8,7\) |
| \(p^6\) | the weight-five coordinate \(h_p=(b_p-5+7\Delta_p)/p^5\), eliminated with weights \(336,-5\) |
| \(p^7\) | the weight-seven endpoint scalar, reduced by H6 to a quadratic in \(\Delta_p\) and eliminated nonlinearly |

Thus, for a fixed set \(S\) of top-half target primes and
\[
 R_*=\prod_{p\in S\setminus\{7,331,769\}}p,\qquad
 A_*=\frac{b_n}{R_*},
\]
the local laws determine \(A_*\) modulo \(p^j\) for every \(p\mid R_*\)
and every \(1\le j\le7\), hence determine one compatible class modulo
\(R_*^7\).

This is not the speculative free defect algebra of Section 5.  At each
grade the endpoint calculation begins with concrete finite multiple
harmonic sums, and the direct/reflected rows remove the surviving
prime-dependent scalar.

### 11.2 The H6 identity is now proved

For \(p\ge11\), write

\[
 S_j=\sum_{k=1}^{p-1}k^{-j},\qquad
 D=\sum_{k=1}^{p-1}\frac{H_k^{(2)}}{k^3}.
\]

The missing finite-harmonic relation is

\[
 \boxed{S_2^2-5S_4-2pD\equiv0\pmod {p^3}.}
\tag{H6}
\]

The proof in `Q32_H6_STAR_PROOF_2026-07-29.md` has three independent
pieces.

First, Zhao's length-three formula and Proposition 3.13 give

\[
 H(1,4,1)\equiv\frac13B_{p-3}^2\pmod p.
\tag{11.2}
\]

A strict stuffle identity then yields

\[
 \sum_{k=1}^{p-1}\frac{H_k(1,1)}{k^4}
 \equiv-\frac16B_{p-3}^2\pmod p.
\tag{11.3}
\]

Second, an exact partial fraction on the additive triangle
\(u,v\ge1,\ u+v\le p-1\), evaluated both directly and after
\(v\mapsto p-v\), gives

\[
 D\equiv S_5-2H(1,4)\pmod {p^2}.
\tag{11.4}
\]

Third, the exact star-binomial identity

\[
 \sum_{k=1}^{p-1}\frac{(-1)^{k-1}\binom{p-1}{k}}{k^4}
 =H^\star(\{1\}^4;p-1)
\]

compares two expansions modulo \(p^3\).  Sun's second-order Kummer
congruence then supplies the required coordinate

\[
 \frac{S_4}{p}
 \equiv
 -8\frac{B_{p-5}}{p-5}
 +4\frac{B_{2p-6}}{2p-6}\pmod {p^2}.
\tag{11.5}
\]

Substitution proves H6.

The two non-elementary inputs have been checked against the primary
sources, not merely against later summaries.  Zhao's arXiv v6 source
prints the length-three display used in the proof of Proposition 3.13
and the specialization (11.2).  Sun's equation (1.1), with \(b=p-5\),
is exactly the affine second-order Kummer law used in (11.5).  The
remaining steps are explicit stuffle, partial-fraction, binomial, and
Faulhaber calculations.  See [Zhao](https://arxiv.org/abs/math/0301252)
and [Sun](https://doi.org/10.1016/S0166-218X(00)00184-0).

### 11.3 The apparent seventh coordinate is not free

Put

\[
 X=p^2S_2,\qquad Y=p^4S_4,\qquad Z=p^5D.
\]

The seventh endpoint block calculation gives

\[
\begin{aligned}
 \Delta_p&\equiv X+5Y-X^2+4Z,\\
 E_p&:=b_{2p}-73+824\Delta_p-\frac{752}{5}
       (b_p-5+7\Delta_p)\\
 &\equiv\frac{24}{5}(935X^2-830Y-332Z)
 \pmod {p^7}.
\end{aligned}
\tag{11.6}
\]

Let \(W_p=E_p/p^6\bmod p\).  H6 is
\[
 X^2-5Y-2Z\equiv0\pmod {p^7},
\]
so (11.6) collapses to

\[
 \boxed{
 p^6W_p\equiv\frac{24\cdot769}{5}\Delta_p^2\pmod {p^7}.
 }
\tag{11.7}
\]

Thus \(W_p\) is not an independent defect digit.  If \(D_7(p,r)\) and
\(Z_7(p,s)\) denote the explicit direct and reflected degree-seven
recurrence-jet residues, their two rows are

\[
\begin{aligned}
 D_7(p,r)&\equiv x_p\left(1-\frac15H_p\right),\\
 Z_7(p,s)&\equiv x_p\left(
 1-\frac{336}{25}H_p+\frac{103}{5\cdot769}p^6W_p
 \right)\pmod {p^7},
\end{aligned}
\tag{11.8}
\]

where \(H_p=b_p-5+7\Delta_p\).  Multiplying the first row by
\(\Delta_p^2\) and using (11.7) removes the nonlinear term:

\[
 \boxed{
 (1680+2472\Delta_p^2)D_7(p,r)-25Z_7(p,s)
 \equiv1655\,\frac{b_n}{p}\pmod {p^7}.
 }
\tag{11.9}
\]

The fixed exceptions are \(7\) (H6 range), \(331\) (right-hand unit), and
\(769\) (endpoint normalization).  Their total logarithmic cost is
\(O(1)\).  There is no hidden exponential coefficient in
\(\Delta_p^2\bmod p^7\): writing
\(\delta_p=(\Delta_p/p^3)\bmod p\) gives the representative
\[
 \Delta_p^2\equiv p^6\delta_p^2\pmod {p^7},
\]
of polynomial height for this fixed grade.

This is a genuine correction to the retracted linear-rank conclusion in
`Q32_ORDER_SEVEN_TARGET_BARRIER_2026-07-29.md`.

### 11.4 Exact global value and exact global limitation

Since \(\log b_n=cn+o(n)\), \(c=\log(17+12\sqrt2)\), the order-seven
modulus first crosses the top-half information threshold:

\[
 A_*<R_*^7
 \quad\Longleftrightarrow\quad
 \log R_*>\frac c8n+o(n)
 =0.440686\ldots n+o(n).
\tag{11.10}
\]

In that range, the CRT class supplied by (11.9) has least representative
equal to the actual integer \(A_*\).  This is exact reconstruction, but
not a contradiction.  It sharpens the local tower while leaving the
inverse-limit obstruction of Section 5.6 intact: all compatible residues
are residues of the diagonal integer \(A_*\), and their canonical
representatives must eventually stabilize to it.

Consequently two logically distinct inputs are still missing:

1. a local theorem extending the fixed elimination to every prescribed
   grade, or a proof identifying the first genuine obstruction; and
2. a global theorem producing a nonzero lift of sublinear height, or
   otherwise proving that the reconstructed diagonal value cannot
   stabilize in the forbidden large-radical regime.

The first input alone does not prove \(\log R_n=o(n)\).  Conversely, one
global small-lift theorem at arbitrary fixed grades would turn the local
tower into a pointwise argument.  The order-seven result therefore reopens
the endpoint-filtration problem, but it does not reopen Route B's
free-dimension count and does not repair Route A's tautological
Casoratian.

### 11.5 Source-safe geometric status

The finite Mellin identity in Section 10.6 is unconditional:
\[
 b_r=-\sum_{t\in\mathbf F_p^\times}{\cal A}_p(t)t^{-r}.
\]
For an integral Laurent-polynomial constant-term model \(\Phi\), there is
also an exact raw toric Kummer alternating trace
\[
 S_{p,r}=\sum_{x\in(\mathbf F_p^\times)^3}\omega_p^r(\Phi(x))
 \equiv-b_r\pmod{\mathfrak p_{p,r}}.
\tag{11.11}
\]

The stronger rank-two Frobenius package suggested by a K3/Mellin analogy
is not yet source-backed.  Peters proves over \(\mathbf C\) that the
rank-three Apéry variation is a symmetric square; the two finite
conifolds have rank-two exponents \(0,\tfrac12\), hence rank-three
semisimple reflection monodromy and middle-extension drop one.  If an
arithmetic rank-three companion with these tame local types is constructed,
Katz's Euler formula would give two-dimensional, pure weight-three Mellin
cohomology for a clean nontrivial character.

What is still missing is an arithmetic inversion self-duality with its
sign, the determinant/epsilon factor, and a trace comparison including
the conifold stalks.  Katz's dimension-two theorem assumes arithmetic
self-duality; his \(\operatorname{Sym}^2(\mathrm{Leg})\) theorem is a
specific example with that arithmetic input, not a general transfer
principle.  Therefore none of

```text
det(Frob)=p^3,
G_geom=G_arith=SL_2,
primitive Mellin trace == b_r,
b_r==0 iff the rank-two fibre is nonordinary
```

is used unconditionally here.  Even if all were later proved, first trace
divisibility would not imply torsion: reduction in the same characteristic
also kills a determinant of valuation three.  A horizontal argument would
still need a divided trace or another independent Frobenius-jet scalar.

## 12. Long-prefix saturation and the adjacent-minor family

There is one further exact reformulation worth recording because it
initially looks like a many-equation escape from Route A.

Put

\[
 t_k(n)=\binom nk^2\binom{n+k}k^2,\qquad
 S_K(n)=\sum_{k=0}^Kt_k(n),
\]

\[
 H=\left\lfloor\frac n3\right\rfloor,\qquad
 L=\left\lfloor\frac n2\right\rfloor,\qquad
 U_n=\prod_{n/2<p\le n}p.
\]

For \(p\in(n/2,n]\), write \(n=p+r\).  Since \(k<p\) for
\(k\le L\), exact valuation of the square root of the summand gives

\[
 \boxed{
 v_p(t_k(n))
 =2{\bf1}_{k>r}+2{\bf1}_{k>p-1-r}.
 }
\tag{12.1}
\]

The folded index
\[
 m_p=\min(r,p-1-r)
\]
always satisfies \(m_p\le H\).  Consequently

\[
 U_n^2\mid t_k(n)\quad(H<k\le L)
\tag{12.2}
\]

and all long prefixes have one common residue:

\[
 \boxed{
 S_K(n)\equiv S_H(n)\equiv b_{n-p}\pmod p
 \quad(H\le K\le L).
 }
\tag{12.3}
\]

Thus every target prime divides every prefix in the long block.  This is
not growing codimension: the block differs only by universally
\(U_n^2\)-divisible increments.

### 12.1 Exact Smith saturation

For a candidate prime put
\[
 r_p=n-p,\qquad s_p=2p-n-1,\qquad
 M_p=\max(r_p,s_p),
\]
and define
\[
 E_n^{(0)}
 =\prod_{\substack{n/2<p\le n\\M_p\le H}}p,\qquad
 E_n^{(1)}
 =\prod_{\substack{n/2<p\le n\\M_p\le H+1}}p.
\]
Equation (12.1) gives the exact top-prime part

\[
 \left(\gcd_{H<k\le L}t_k(n)\right)_{\{n/2<p\le n\}}
 =U_n^2\left(E_n^{(0)}\right)^2.
\tag{12.4}
\]

This corrects the earlier coarse condition \(p\le2H+1\), which is
false.  For example, at \(n=30,p=17\) one has
\[
 H=10,\quad (r_p,s_p)=(13,3),\quad
 v_{17}(t_{11})=2,
\]
although \(17\le2H+1\).  The second Kummer root contributes to the
term gcd exactly when \(M_p\le H\), not on the coarse interval.

After dividing the tail differences by their maximal universal
top-prime content while retaining one copy of each candidate prime, the
tail ideal in the semilocal ring at \(p\in(n/2,n]\) is exactly
\((U_n)\).  Adding the prefix coordinate gives

\[
 (S_H,\hbox{ saturated tail differences})
 =(S_H,U_n).
\tag{12.5}
\]

Its one Smith invariant is therefore

\[
 \prod_{\substack{n/2<p\le n\\p\mid b_{n-p}}}p.
\tag{12.6}
\]

There is no second selective elementary divisor.  Finite differences,
higher exterior powers, and Bézout combinations of the tail rows cannot
change this: after saturation the entire target condition is one additive
prefix coordinate.

The obstruction can also be realized exactly.  The prefix sequence obeys

\[
 (k+1)^4S_{k+1}
 -\bigl((n-k)^2(n+k+1)^2+(k+1)^4\bigr)S_k
 +(n-k)^2(n+k+1)^2S_{k-1}=0.
\tag{12.7}
\]

Its coefficient sum is zero, so \(S_k\mapsto S_k+C\) preserves the
recurrence and every tail difference.  CRT can choose the additive
constant so that an arbitrary subset of the candidate primes divides all
long prefixes, while retaining the actual \(t_k(n)\), all valuations, and
positivity.  Hence the tail package alone admits worst-case horizontal
alignment.

### 12.2 A genuine family of positive selective minors

Write the summand ratio in lowest terms:

\[
 \frac{t_{k+1}}{t_k}=\frac{a_k}{b_k},\qquad
 a_k=\frac{(n-k)^2(n+k+1)^2}{d_k},\quad
 b_k=\frac{(k+1)^4}{d_k},
\]

\[
 d_k=\gcd\bigl((n-k)^2(n+k+1)^2,(k+1)^4\bigr).
\]

For \(H<k<L\), define the normalized adjacent Hankel minor

\[
\begin{aligned}
 C_{n,k}
 &:=-\frac{S_{k-1}S_{k+1}-S_k^2}
          {\gcd(t_k,t_{k+1})}\\
 &=b_kS_k-a_kS_{k-1}.
\end{aligned}
\tag{12.8}
\]

The second expression proves integrality.  Strict decrease of the
positive ratios \(t_{k+1}/t_k\) proves \(C_{n,k}>0\).  At a target prime,
the prefix has one selective factor \(p\), while both adjacent summands
have at least \(p^2\); therefore

\[
 \boxed{R_n^{\rm top}\mid C_{n,k}\quad(H<k<L).}
\tag{12.9}
\]

For the first admissible \(k\),

\[
 \log C_{n,k}=n\log16+O(\log n).
\tag{12.10}
\]

Thus a single minor is a real characteristic-zero carrier, but it has
linear height.  Larger Hankel/Dodgson minors do not increase the selective
rank: writing \(S_{H+j}=X+U_n^2Y_j\), every normalized determinant is
linear in \(X=S_H\) after the universal factor is removed.

The family gcd

\[
 G_n^{\rm adj}=\gcd_{H<k<L}C_{n,k}
\tag{12.11}
\]

is a sharper exact interface.  Its complete top-prime Smith ideal can in
fact be computed.  Put
\[
 Q_n=U_n^2\left(E_n^{(1)}\right)^2,\qquad
 {\cal O}_n=\mathbf Z[a^{-1}:\gcd(a,U_n)=1].
\]
For \(n\ge60\),
\[
 \boxed{
 (C_{n,H+1},\ldots,C_{n,L-1}){\cal O}_n
 =(S_H,Q_n){\cal O}_n.
 }
\tag{12.11a}
\]
Equivalently,
\[
 \boxed{
 (G_n^{\rm adj})_{\{n/2<p\le n\}}
 =\gcd(S_H,Q_n).
 }
\tag{12.11b}
\]

Here is the essential determinant calculation.  Write
\[
 \delta_k=a_k-b_k,\qquad
 T_{k-1}=S_{k-1}-S_H,
\]
so
\[
 C_{n,k}=-\delta_kS_H+\gamma_k,\qquad
 \gamma_k=b_kt_k-\delta_kT_{k-1}.
\]
Then
\[
 \det
 \begin{pmatrix}
  -\delta_k&\gamma_k\\
  -\delta_\ell&\gamma_\ell
 \end{pmatrix}
 =\delta_\ell b_kt_k-\delta_kb_\ell t_\ell.
\tag{12.11c}
\]
For adjacent rows this factors as \(t_{k+1}\) times a ratio-difference
coefficient.  After writing \(N=n(n+1)\), \(x=k+1\), the numerator of that
coefficient is a nonzero polynomial of degree six in \(x\).  There are at
least seven adjacent positions when \(n\ge60\), so at every candidate
prime one determinant attains the exact universal valuation recorded in
\(Q_n\).
The first coefficient column is locally primitive by the analogous
degree-three calculation.  The two Smith invariants are therefore
\[
 1,\quad Q_n,
\]
which proves (12.11a).

In particular, for \(p\in(n/2,n]\),
\[
 v_p(G_n^{\rm adj})
 =\min\left(v_p(S_H),
 2+2{\bf1}_{M_p\le H+1}\right).
\tag{12.11d}
\]

The shift from \(H\) in (12.4) to \(H+1\) here is real: the
coefficient determinants acquire their second Kummer factor when the
outer folded root is at the first carrier row.  The set
\(\{p:M_p\le H+1\}\) contains at most two candidate moduli for each
\(n\).  This correction changes the higher-power Smith ledger but not
the radical or the linear-height obstruction.
The radical of its candidate-prime part is exactly the target radical,
but its higher powers are merely capped powers of the same prefix
coordinate.  Hence
\[
 \log R_n^{\rm top}
 \le\log (G_n^{\rm adj})_{\{n/2<p\le n\}}
 \le4\log R_n^{\rm top}.
\tag{12.11e}
\]
Sublinear height for this Smith generator is quantitatively equivalent to
the original target-radical theorem; the many-minor gcd has not created a
new transverse condition.

The ordinary integer gcd is nevertheless numerically small.  An
independent exact computation through \(n=1000\) gives the dyadic maxima

| interval | \(\max \log G_n^{\rm adj}/n\) |
|---:|---:|
| \((100,200]\) | \(0.0987381\) |
| \((200,400]\) | \(0.0507389\) |
| \((400,600]\) | \(0.0328444\) |
| \((600,800]\) | \(0.0315821\) |
| \((800,1000]\) | \(0.0168826\) |

This is strong evidence for the conjecture, not a new proof mechanism.  It
is on the same scale as the
earlier two-truncation gcd in
`Q32_CODEX_RESUME_2026-07-23.md`, Sections 42 and 46.  The additive-shift
countermodel is exact at every valuation allowed by (12.11d): CRT can
prescribe independently, for every candidate prime, any exponent between
zero and its Smith cap while retaining the actual tail differences,
recurrence, determinant data, and positivity.  Therefore no Smith,
finite-difference, or tail-only argument can prove that (12.11) has
sublinear height.

The one possible reopening is to exploit the distinguished boundary
\(S_0=1\).  Backward propagation from the long block to that boundary
must cross the \(p\)-singular folded index \(m_p\), which depends on the
candidate prime.  A useful theorem would have to bypass or jointly control
these moving singular steps; otherwise their product restores the
top-half primorial and linear height.

The exact reproducer is
`../scripts/q32_adjacent_minor_gcd.py`; it also checks the normalized
Hankel identity, positivity, target divisibility, both square-primorial
truncation congruences, and the gcd obtained by intersecting the adjacent
and two-truncation carrier families.

### 12.3 A first quotient bridge, and why it is still universal

Let
\[
 J=\left\lfloor\frac{n-1}{3}\right\rfloor,\qquad
 L_n=S_J(n),\qquad
 H_n^{\rm upper}=\sum_{k=\lceil n/2\rceil}^nt_k(n).
\]
The exact two-truncation calculation gives

\[
 U_n^2\mid H_n^{\rm upper}-4L_n.
\tag{12.12}
\]

Together with (12.2), this implies, with the floor adjustment at
\(3\mid n\) handled by its single boundary term,

\[
 b_n\equiv5L_n\pmod {U_n^2}.
\tag{12.13}
\]

At a target prime, division of (12.13) identifies the first quotient digit
\[
 \frac{b_n}{p}\equiv5\frac{L_n}{p}\pmod p.
\]
Likewise (12.8) identifies \(C_{n,k}/p\) with a fixed coefficient times
the same divided prefix digit.  This is a genuine bridge to the endpoint
quotient tower of Section 11.

However, the immediate cancellation
\[
 5C_{n,k}+(a_k-b_k)b_n
\]
lies in the universal \(U_n^2\)-tail ideal: it gains the second power for
every candidate prime, not selectively for targets.  To turn the bridge
into a proof one needs a higher corrected carrier which gains \(p^j\)
only on the target locus while keeping height \(O(n)\) with a constant
not growing linearly in \(j\).  No such global correction is presently
proved.

## 13. Precision eight: the target Casoratian continues

The first grade not visible in Section 11 is effective weight seven.
It contains new complete-block and depth finite-harmonic coordinates, so
the order-seven square law is not an induction step.  Nevertheless, the
direct/reflected target rows themselves can be extended unconditionally
one more digit.

For \(p\ge11\), put
\[
\begin{aligned}
 \Delta_p&=b_{p-1}-1,\\
 H_p&=b_p-5+7\Delta_p,\\
 {\cal E}_p&=b_{2p}-73+824\Delta_p-\frac{752}{5}H_p,\\
 {\cal F}_p&=b_{2p-1}-5-8\Delta_p-\frac{336}{5}H_p.
\end{aligned}
\tag{13.1}
\]
The proved raw endpoint laws and \(H6\) identity imply
\[
 {\cal E}_p,{\cal F}_p\in p^6\mathbf Z_{(p)},\qquad
 769{\cal F}_p+103{\cal E}_p\in p^7\mathbf Z_{(p)}.
\tag{13.1a}
\]

Let \(p\mid b_r\), set \(s=p-1-r\), and write
\(x=b_{p+r}/p\).  Extending the shifted fundamental solution
\({\cal U}\) through degree eight and the companion
\({\cal J}\) through degree five gives rows \(D_8,Z_8\) satisfying
\[
\begin{aligned}
 D_8&\equiv
 x\left(1-\frac15H_p\right)
 -\frac15p^2H_p{\cal J}_r(0),\\
 Z_8&\equiv
 x\left(1-\frac{336}{25}H_p-\frac15{\cal F}_p\right)
 +\frac{166144}{25}p^2H_p{\cal J}_s(0)
 \pmod {p^8}.
\end{aligned}
\tag{13.2}
\]

The companion values coincide at a reflected zero:
\[
 {\cal J}_s(0)\equiv{\cal J}_r(0)\pmod p.
\tag{13.3}
\]
Indeed, for the standard solution \(U_j(0)=b_j\) and the companion
\(V_0=0,V_1=1\),
\[
 U_kV_{k+1}-U_{k+1}V_k=\frac1{(k+1)^3}.
\]
At \(U_r=U_s=0\), this identity, the Apéry recurrence at \(r\), and
reflection identify \(V_r\) with \(V_s\).  Since
\({\cal J}_j(0)=V_j\), (13.3) follows.

The fixed combination
\[
\boxed{
\left(166144+33296H_p+{\cal F}_p\right)D_8
+5Z_8
\equiv166149\,x\pmod {p^8}.
}
\tag{13.4}
\]
therefore cancels the new target-dependent companion coordinate.
Here
\[
 166149=3^2\cdot18461,
\]
so the right coefficient is a fixed unit outside \(p=18461\) in the
stated range.  Exact symbolic reduction and all \(163\) target rows at
primes through \(1000\) pass independently, including both target rows
at \(p=769\).

This formulation is valid at \(p=769\).  Its integer-coefficient form is
\[
\begin{aligned}
 &[5b_{2p-1}+830695-40\Delta_p+166144H_p]D_8+25Z_8\\
 &\hspace{38mm}\equiv830745\,x\pmod {p^8}.
\end{aligned}
\tag{13.4a}
\]
For \(p\ne769\), the split coordinates
\[
 w_p=\frac{{\cal E}_p}{p^6},\qquad
 v_p=\frac{{\cal F}_p+(103/769){\cal E}_p}{p^7}
\tag{13.4b}
\]
are integral and recover the earlier \(w_p,v_p\) version of (13.4).
Thus \(769\) is only an artificial exception to that normalization;
the sole fixed inversion exception of the target law is \(18461\).

The proof and reproducer are

```text
Q32_ORDER_EIGHT_TARGET_CASORATIAN_2026-07-29.md
../scripts/q32_order_eight_target_audit.py
```

### 13.1 A conjectural all-\(m\) endpoint law

There is also a strong new computational pattern.  After retaining the
full next digit of \(w_p\), every endpoint residual through
\(p\le1000,\ m\le20\) is a fixed multiple of the single anchor \(v_p\).
Here \(w_p,v_p\) have the split normalization (13.4b), so this separate
all-\(m\) conjecture still excludes \(p=769\).
Define
\[
\begin{aligned}
 P_0(m)&=3845m^4-29268m^3+36974m^2-9112,\\
 Q_0(m)&=45371m^4-58102m^2+536,\\
 N_8&=305911296=2^9\,3^3\,22129,\\
 C_m&=-\frac{m^3(P_0(m)b_m+Q_0(m)b_{m-1})}{N_8},\\
 D_m&=\frac{m^3(Q_0(m)b_m+P_0(-m)b_{m-1})}{N_8}.
\end{aligned}
\tag{13.5}
\]
With the proved lower carriers \(E_m,F_m,P_m,Q_m,R_m,S_m\), the
observed law is
\[
\begin{aligned}
 b_{mp}-b_m
 &\equiv E_m\Delta_p+P_mH_p+p^6R_mw_p+p^7C_mv_p,\\
 b_{mp-1}-b_{m-1}
 &\equiv F_m\Delta_p+Q_mH_p+p^6S_mw_p+p^7D_mv_p
 \pmod {p^8}.
\end{aligned}
\tag{13.6}
\]
The normalization has
\[
 C_1=D_1=C_2=0,\qquad D_2=1.
\]

This was rationally reconstructed through \(p\le400\) and then checked
on disjoint larger primes through \(1000\), with \(6520\) exact
residual-divisibility checks and no failure.  It is not yet a theorem.
One part of the missing block reduction can now be proved.  For
\[
 \xi=\frac{H(6)}p,\quad \eta=\frac{H(2,4)}p,\quad
 A=H(2,2,3),\quad B=H(2,5)\pmod p,
\]
the exact relations
\[
 3\eta=2\xi,\qquad 3A=14\xi,\qquad 2B=-7\xi
\tag{13.6a}
\]
show that the effective weight-seven finite-MHS quotient is
one-dimensional.  The first identity uses the necessary lifted
reversal correction
\[
 H(2,4)-H(4,2)
 \equiv p\{2H(4,3)+4H(5,2)\}\pmod {p^2}.
\]
Thus ordinary reversal is no longer a gap.  What remains is the full
change from the block basis to \((\Delta_p,H_p,w_p,v_p)\), including
lifted lower-coordinate terms, followed by the termwise
identity/Gosper certificate.  Projecting only the primitive MHS vector
does not retain those lower lifts.

The precise conjecture and audit are

```text
Q32_WEIGHT_EIGHT_ENDPOINT_RANK_ONE_2026-07-29.md
Q32_WEIGHT_SEVEN_MHS_RANK_ONE_2026-07-29.md
../scripts/q32_weight_eight_endpoint_rank_audit.py
../scripts/q32_weight_seven_mhs_rank_one_audit.py
```

### 13.2 Why this is still local reconstruction

Equation (13.4) determines \(x=b_n/p\pmod {p^8}\), hence
\(b_n\pmod {p^9}\), for each target prime outside a fixed set.  Across
all such primes with radical \(R\), the unique CRT carrier is
\[
 b_n-(b_n\bmod R^9)
 =R^9\left\lfloor\frac{b_n}{R^9}\right\rfloor.
\tag{13.7}
\]
It is positive only while \(R^9\le b_n\), where it gives a fixed-order
linear bound; once \(R^9>b_n\), it is zero because the local data have
reconstructed the actual diagonal integer.

The same saturation already occurs at order seven, with \(R^8\).
Pairing target primes merely factors the same remainder into pairwise
CRT coordinates and increases height faster than selective exponent.
Thus the new local digit is genuine progress on the endpoint algebra,
but it does not supply the independent nonzero lift required for
\(\log R_n=o(n)\).

## 14. The fixed boundary does not rescue the long-prefix route

The additive-shift countermodel in Section 12 does not by itself use the
distinguished value \(S_0=1\).  The actual boundary can be propagated
exactly, but the propagation stops at the same folded Kummer edge.

For the reduced ratio \(t_{k+1}/t_k=a_k/b_k\), the prefix recurrence is
\[
 b_kS_{k+1}-(a_k+b_k)S_k+a_kS_{k-1}=0.
\tag{14.1}
\]
With \(C_k=b_kS_k-a_kS_{k-1}\), repeated Euclidean substitution gives
\[
\left(\prod_{i=1}^Ka_i\right)S_0
=\left(\prod_{i=1}^Kb_i\right)S_K
-\sum_{j=1}^K
 \left(\prod_{i=1}^{j-1}b_i\right)
 \left(\prod_{i=j+1}^Ka_i\right)C_j.
\tag{14.2}
\]

For a target \(p=n-r\), put
\[
 m=\min(r,p-1-r).
\]
The exact Kummer support law gives
\[
 S_k\equiv C_k\equiv0\pmod p\quad(k\ge m),
\qquad
 S_{m-1},C_{m-1}\not\equiv0\pmod p.
\tag{14.3}
\]
The transition coefficient at the fold satisfies
\[
 v_p(a_m)=
 \begin{cases}
 2,&r\ne p-1-r,\\
 4,&r=p-1-r.
 \end{cases}
\tag{14.4}
\]
All terms in (14.2) which reach below the fold contain this missing
factor.  Thus \(S_0=1\) is remembered only after paying \(p^2\), or
\(p^4\) centrally.  In DVR language the folded transfer has Smith form
\(\operatorname{diag}(1,p^2)\), respectively
\(\operatorname{diag}(1,p^4)\); unit transfers on either side cannot
remove that invariant factor.

Reflection exchanges the two Kummer roots and fixes \(m\).  The high
suffix is a scalar copy of the same folded prefix modulo \(p\), and a
different cutoff detects the target only if it remains on the far side
of \(m\).  For two different target characteristics the two folded DVR
modules form a direct CRT product; the pure-cross continuant is a unit
in both and does not couple their singular factors.

Therefore the distinguished boundary closes the formal
additive-constant objection but yields a sharper route-specific no-go:
every identity generated by the one-dimensional summand transfer must
pay the moving folded determinant.  Clearing those determinants
simultaneously restores a linear-height Kummer product.  Reopening this
route requires an identity outside that transfer module which couples
folded edges in different characteristics.

## 15. Correction to the short endpoint-transfer audit

A proposed follow-up attempted to combine the two precision-eight target
laws by transporting the recurrence state from \(2p-1\) to \(2q-1\),
where \(p<q=p+h<2p\).  Its main local calculation was not correct.

Put
\[
 S_k=
 \begin{pmatrix}
  P(k)/(k+1)^3&-k^3/(k+1)^3\\
  1&0
 \end{pmatrix},
 \qquad
 T_{p,q}=S_{2q-2}\cdots S_{2p-1}.
\tag{15.1}
\]
Then
\[
 y_{2q-1}=T_{p,q}y_{2p-1},
 \qquad
 \det T_{p,q}=\left(\frac{2p-1}{2q-1}\right)^3.
\tag{15.2}
\]
In the short-gap range \(0<h<p/2\), the index interval in (15.1)
contains the two \(p\)-singular steps \(2p-1,2p\), but it contains no
\(q\)-singular step.  The singular pair has integral numerator
\[
 M_{2p}M_{2p-1}\equiv
 \begin{pmatrix}-25&5\\-5&1\end{pmatrix}\pmod p.
\tag{15.3}
\]
It has rank one and a unit entry, while its determinant has valuation
six and the transfer denominator has valuation three.  Hence
\[
 \operatorname{SmithExp}_p(T_{p,q})=(-3,3).
\tag{15.4}
\]

At \(q\), every \(k\) and \(k+1\) occurring in (15.1) is a \(q\)-unit:
for \(h\le p-2\),
\[
 q<2p-1\le k<k+1\le2q-1<2q.
\]
Consequently every individual \(S_k\) belongs to
\(\operatorname{GL}_2(\mathbb Z_q)\), and so do their product and its
inverse:
\[
 \boxed{\operatorname{SmithExp}_q(T_{p,q})=(0,0).}
\tag{15.5}
\]
In particular, the symmetric claim
\[
 \operatorname{SmithExp}_q(T_{p,q}^{-1})=(-3,3)
\]
is false: reversing the same transfer does not introduce the absent
steps \(2q-1,2q\).

There are two full-range endpoint effects which are irrelevant to a
logarithmic prime gap but useful for auditing the indices.  Exact
fraction arithmetic for every checked pair \(5\le p<100\) gives
\[
\begin{array}{c|c|c}
\text{condition}&\operatorname{SmithExp}_p(T_{p,q})
 &\operatorname{SmithExp}_q(T_{p,q})\\ \hline
2h=p+1&(-3,0)&(0,0),\\
h=p-1&(-3,3)&(0,3),\\
\text{otherwise}&(-3,3)&(0,0).
\end{array}
\tag{15.6}
\]
The short row (15.4)--(15.5) is proved above.  Table (15.6) is recorded
as an exact finite audit, not as a proved all-\(p\) classification of the
additional \(p\)-singular block.
The last \(q\)-endpoint case cannot be an actual top-half target pair
for \(q>5\): it gives \(q=2p-1\), and \(q\le n<2p\) then forces
\(q=n\), whereas \(b_q\equiv b_1=5\pmod q\).
The exact reproducer checks all \(212\) prime pairs with
\(5\le p<100\), including \(104\) short-gap pairs:

```text
../scripts/q32_endpoint_transfer_smith_audit.py
```

This correction invalidates the claimed two-sided Smith/Fitting proof of
a no-go theorem.  It does **not** by itself produce a positive carrier.
On the corrected lattice, the \(q\)-transport is unimodular and the only
moving denominator clearer comes from \(p\).  Thus any elimination that
uses that clearer obtains its \(p\)-factor from a coefficient, not from
the target condition.  A rigorous conclusion about the primitive scalar
elimination ideal still requires the actual presentation of all local
target coordinates; the schematic two-row argument does not supply it.

### 15.1 Corrected nested-reflection check

The related nested-Casoratian proposal also needs its endpoints oriented
correctly.  For the three-target example \(n=321\),
\[
 321=193+128=179+142,
\]
and reflection gives
\[
 128\longleftrightarrow64\pmod {193},
 \qquad
 142\longleftrightarrow36\pmod {179}.
\]
The actual nested intervals are therefore
\[
 [64,128]\subset[36,142],
\tag{15.7}
\]
not the intervals used in the proposed numerical check.

For the four ordered endpoints \(36,64,128,142\), the Plücker identity is
\[
\begin{aligned}
 &K(36,64)K(128,142)
 -K(36,128)K(64,142)\\
 &\hspace{31mm}
 +K(36,142)K(64,128)=0.
\end{aligned}
\tag{15.8}
\]
The first product has total boundary gap \(28+14=42=3h\).
The last product contains the two target reflection determinants and is
divisible by \(179\cdot193\) after denominator clearing.  But the middle
product still contains gaps \(92\) and \(78\), and each of its two factors
is a unit in both relevant characteristics.  Exact continuant evaluation
confirms all four unit statements.

There is a more basic selectivity obstruction.  The prime factor in each
reflection continuant is the universal factor
\(2x+d+1\) imposed by reflection of the recurrence; it is present whether
or not the distinguished Apéry solution vanishes at the endpoints.
After this universal factor is removed, the reflection determinant has
no remaining target-prime divisibility in the example.  Hence nested
Plücker condensation of the ordinary Casoratians cannot be the missing
pure-cross carrier.  A viable nested construction would have to use a
divided target quotient or another observable beyond the universal
reflection transfer.

## 16. Path D as a coefficientwise inverse-denominator problem

The fixed factorial-ratio carriers admit a useful characteristic-zero
reformulation which is stronger than merely observing that both factors are
holonomic.  Fix \(d\ge2\), put
\[
 B_{n,d}=\binom n{\lfloor n/d\rfloor},
 \qquad
 q_{n,d}=\frac{B_{n,d}}{\gcd(b_n,B_{n,d})}.
\tag{16.1}
\]
Thus \(q_{n,d}\) is the reduced denominator of the rational number
\(b_n/B_{n,d}\).  The sufficient gcd theorem (HD.4) is equivalently a
pointwise **near-maximal denominator** theorem:
\[
 \log q_{n,d}
 =\log B_{n,d}-o_d(n)
 =n\,h(1/d)+o_d(n).
\tag{16.2}
\]
This reverses the direction of the usual denominator theorems for
\(G\)-functions, which give a common upper denominator for a prefix.

For \(0\le a<d\), define
\[
\begin{aligned}
 A_{d,a}(x)&=\sum_{m\ge0}b_{dm+a}x^m,\\
 R_{d,a}(x)&=\sum_{m\ge0}\binom{dm+a}{m}^{-1}x^m,\\
 H_{d,a}(x)&=A_{d,a}\mathbin{\star}R_{d,a}
 =\sum_{m\ge0}\frac{b_{dm+a}}{\binom{dm+a}{m}}x^m ,
\end{aligned}
\tag{16.3}
\]
where \(\star\) denotes Hadamard product.  The coefficient ratio in the
second line is
\[
 \frac{[x^{m+1}]R_{d,a}}{[x^m]R_{d,a}}
 =
 \frac{(m+1)\prod_{i=1}^{d-1}((d-1)m+a+i)}
      {\prod_{i=1}^{d}(dm+a+i)}.
\tag{16.4}
\]
Hence \(R_{d,a}\) is hypergeometric.  The exact row-lcm identity
\[
 \operatorname {lcm}_{0\le k\le n}\binom nk
 =\frac{\operatorname {lcm}(1,\ldots,n+1)}{n+1}
\tag{16.5}
\]
shows that the common denominator of its first \(M+1\) coefficients is
\(\exp(O_d(M))\).  Therefore \(R_{d,a}\) is a \(G\)-function directly
from the definition.  Sections and Hadamard products preserve
\(G\)-functions, so \(H_{d,a}\) is a \(G\)-function whose \(m\)-th reduced
coefficient denominator is exactly
\[
 q_{dm+a,d}.
\tag{16.6}
\]
For \(d=2\), the reciprocal factors are explicitly
\[
 R_{2,0}(x)={}_2F_1(1,1;1/2;x/4),
 \qquad
 R_{2,1}(x)={}_2F_1(1,2;3/2;x/4).
\tag{16.7}
\]
Thus even the clean central case asks for a lower denominator theorem for
the individual coefficients of an explicit Hadamard \(G\)-function, not
an upper denominator theorem.

### 16.1 Every residue section has asymptotically maximal prefix lcm

There is an unconditional theorem at the prefix level, and it holds
separately in every residue section.  Put
\[
 Q_{d,a}(M)=
 \operatorname {lcm}_{0\le m\le M}q_{dm+a,d}.
\tag{16.8}
\]
Then, for every fixed \(d\ge2\) and \(0\le a<d\),
\[
 \boxed{\log Q_{d,a}(M)=dM+o_d(M).}
\tag{16.9}
\]

For the upper bound, (16.5) gives
\[
 q_{dm+a,d}\mid B_{dm+a,d}
 \mid\operatorname {lcm}(1,\ldots,dm+a+1),
\]
and hence
\[
 Q_{d,a}(M)\mid\operatorname {lcm}(1,\ldots,dM+a+1).
\tag{16.10}
\]
The prime number theorem gives
\(\log Q_{d,a}(M)\le dM+o_d(M)\).

For the matching lower bound, let \(p\) be a sufficiently large prime and
choose the unique \(c\in\{0,\ldots,d-1\}\) with
\[
 p+c\equiv a\pmod d.
\]
Set \(n=p+c=dm+a\).  Once \(p\) exceeds a constant depending only on
\((d,a)\), one has \(c<m<p\).  The numerator of
\(\binom{p+c}{m}\) contains \(p\) once and neither denominator factorial
contains \(p\), so
\[
 v_p(B_{n,d})=1.
\tag{16.11}
\]
Gessel--Lucas gives
\[
 b_{p+c}\equiv b_1b_c=5b_c\not\equiv0\pmod p
\tag{16.12}
\]
after enlarging the same fixed threshold.  Consequently
\(p\mid q_{dm+a,d}\).  Every sufficiently large prime
\[
 p\le dM+a-(d-1)
\]
is obtained in this way with \(m\le M\), and therefore divides
\(Q_{d,a}(M)\).  The lower bound in (16.9) follows from
\(\vartheta(dM-O_d(1))=dM+o_d(M)\).

This strengthens the parity-section observation: generic \(G\)-function
prefix denominator bounds are already saturated inside each fixed
section, not only after mixing residue classes.  It is nevertheless not
the desired pointwise theorem.  The lower bound is supplied by the sparse
indices \(n=p+c\); it is compatible with arbitrarily small denominators at
other individual indices.

### 16.2 A target is an isolated denominator hole

The inverse-denominator formulation has one further exact local feature.
Take \(d=2\), write
\[
 C_n=\binom n{\lfloor n/2\rfloor},
 \qquad q_n=q_{n,2},
\]
and let \(p\in(n/2,n]\), \(p\ge7\), be a target.  Write \(n=p+r\).
The endpoint congruences \(b_0=1\) and
\(b_{p-2}\equiv b_1=5\pmod p\) show that
\[
 1\le r\le p-3.
\tag{16.13}
\]
The central carry law gives
\[
 v_p(C_{n-1})=v_p(C_n)=v_p(C_{n+1})=1.
\tag{16.14}
\]
On the other hand, Lucas and nonconsecutivity of the Apéry zero set give
\[
\begin{aligned}
 b_n&\equiv5b_r\equiv0\pmod p,\\
 b_{n-1}&\equiv5b_{r-1}\not\equiv0\pmod p,\\
 b_{n+1}&\equiv5b_{r+1}\not\equiv0\pmod p.
\end{aligned}
\tag{16.15}
\]
It follows that the reduced coefficient denominators have the exact
valuation pattern
\[
 \boxed{
 (v_p(q_{n-1}),v_p(q_n),v_p(q_{n+1}))=(1,0,1).
 }
\tag{16.16}
\]
Thus a top-half target is not merely a common divisor in (HD.1): it is an
isolated missing prime in the denominator sequence of the central
Hadamard \(G\)-function.

This is genuine additional organization, but not yet horizontal
separation.  Reduced denominators depend nonlinearly on the coefficient
recurrence, and the order-two recurrence over \(\mathbb F_p\) permits an
isolated zero with two nonzero neighbors.  Different target primes again
live in a direct CRT product.  To turn (16.16) into P3.2 one would need a
new theorem excluding linearly many simultaneous large-prime holes at one
coefficient index; neither prefix maximality (16.9) nor a standard
\(G\)-function denominator bound supplies that pointwise assertion.

### 16.3 The three-coefficient hole carrier

The pattern (16.16) can be packaged as one integer which cancels the
denominator support common to all three coefficients.  Define
\[
 \mathfrak d_n=
 \frac{\gcd(q_{n-1},q_{n+1})}
      {\gcd(q_{n-1},q_n,q_{n+1})}.
\tag{16.17}
\]
Equivalently, for every prime \(\ell\),
\[
 v_\ell(\mathfrak d_n)=
 \max\{0,\min(v_\ell(q_{n-1}),v_\ell(q_{n+1}))
                   -v_\ell(q_n)\}.
\tag{16.18}
\]
Thus \(\mathfrak d_n\) is the positive local denominator curvature at
the middle coefficient.

For \(n>10\), it has the exact top-half support
\[
 \boxed{
 \operatorname {rad}_{(n/2,n]}\mathfrak d_n
 =\prod_{p\in T_n}p.
 }
\tag{16.19}
\]
Indeed, write \(n=p+r\).  The three central binomial coefficients all
contain \(p\) exactly once precisely in the interior range
\(1\le r\le p-3\).  On this range Lucas gives
\[
 v_p(q_{n+j})=
 \begin{cases}
 0,&p\mid b_{r+j},\\
 1,&p\nmid b_{r+j},
 \end{cases}
 \qquad j=-1,0,1.
\tag{16.20}
\]
If \(p\mid b_r\), nonconsecutivity gives the pattern \(1,0,1\), hence
\(p\mid\mathfrak d_n\).  Conversely, (16.18) can be positive only for
that pattern, and therefore forces \(p\mid b_r\).  The three excluded
endpoint residues are nontargets, so no boundary exception remains.

This is a more selective characteristic-zero carrier than \(C_n\):
candidate primes which occur in all three reduced denominators cancel
from (16.17).  Its exact values are strikingly small.  Through
\(n=10{,}000\), its largest bit length is \(45\), at \(n=6792\);
on \(5001\le n<10{,}000\),
\[
 \max\frac{\log\mathfrak d_n}{n}
 =0.0045749844\ldots .
\tag{16.21}
\]
This is finite evidence only.  Formula (16.19) means that a proof of
\[
 \log\mathfrak d_n=o(n)
\tag{16.22}
\]
would settle the top-half branch, but no denominator-curvature theorem
for rational \(G\)-function coefficients currently gives (16.22).

There is an exact comparison which prevents overinterpreting the small
data.  Put \(g_n=\gcd(b_n,C_n)\) and
\[
 U_n=
 \frac{\gcd(C_{n-1},C_{n+1})}
      {\gcd(C_{n-1},C_n,C_{n+1})}.
\tag{16.23}
\]
Prime by prime, writing
\[
 c_j=v_\ell(C_{n+j}),\qquad
 e_j=v_\ell(\gcd(b_{n+j},C_{n+j})),
\]
formula (16.18) gives
\[
\begin{aligned}
 v_\ell(\mathfrak d_n)
 &=
 \max\{0,\min(c_{-1}-e_{-1},c_1-e_1)
                  -(c_0-e_0)\}\\
 &\le e_0+\max\{0,\min(c_{-1},c_1)-c_0\}.
\end{aligned}
\]
Consequently
\[
 \mathfrak d_n\mid g_nU_n.
\tag{16.24}
\]
The adjacent central coefficients make \(U_n\) explicit.  If \(n=2m\),
then \(C_{n-1}=C_n/2\), so \(U_n=1\).  If \(n=2m+1\), write
\[
 C_{2m}=(m+1)K,\quad
 C_{2m+1}=(2m+1)K,\quad
 C_{2m+2}=2(2m+1)K;
\]
then \(U_n=\gcd(m+1,2)\).  Hence
\[
 \boxed{\mathfrak d_n\mid2\gcd(b_n,C_n),}
\qquad
 \mathfrak d_n\mid\gcd(b_n,C_n)\quad(n\ {\rm even}).
\tag{16.25}
\]
The hole carrier therefore removes irrelevant common support and gives
the exact local curvature interpretation, but its observed small height
is still evidence about the original Apéry gcd rather than an independent
proved compression.

### 16.4 The normalized recurrence permits arbitrary simultaneous holes

The recurrence itself cannot prove a height bound for \(\mathfrak d_n\).
This can be made exact rather than inferred from the isolated-zero
pattern.  Put
\[
 y_n=\frac{b_n}{C_n}.
\]
The adjacent central-binomial ratios give a primitive
parity-dependent row
\[
 \alpha_ny_{n+1}-\beta_ny_n+\epsilon_ny_{n-1}=0,
\tag{16.26}
\]
where
\[
\begin{array}{c|ccc}
 &\alpha_n&\beta_n&\epsilon_n\\ \hline
 n\ {\rm even}
 &(n+1)^4&\dfrac{n+2}{2}P(n)&
      \dfrac{n^3(n+2)}4\\[3mm]
 n\ {\rm odd}
 &\dfrac{2(n+1)^3}{\gamma_n}&
      \dfrac{P(n)}{\gamma_n}&
      \dfrac{n^2(n+1)}{2\gamma_n}
\end{array},
\qquad
 \gamma_n=\gcd(n+1,5).
\tag{16.27}
\]
Indeed
\(\gcd(n+1,P(n))=\gcd(n+1,5)\) on the odd rows, and the even polar
coefficients are coprime.  Thus
\[
 \gcd(\alpha_n,\beta_n,\epsilon_n)=1
\]
in both parities.  Omitting \(\gamma_n\) would incorrectly call the odd
row primitive when \(n\equiv4\pmod5\).

Take an interior top-half candidate \(p\), so that
\[
 n=p+r,\qquad1\le r\le p-3,\qquad p\ge7.
\]
Both \(\alpha_n\) and \(\epsilon_n\) are \(p\)-units.  Therefore, if
\[
 v_p(y_n)\ge0,\qquad v_p(y_{n-1})=-1,
\]
the recurrence forces \(v_p(y_{n+1})=-1\).  For the actual Apéry
solution, writing
\[
 Y_\pm=p\,y_{n\pm1}\pmod p
\]
gives only the automatic polar cancellation
\[
\alpha_nY_++\epsilon_nY_-\equiv0\pmod p.
\tag{16.28}
\]
No singular coefficient or second zero condition occurs.

More strongly, let \(R\) be any squarefree product of interior candidate
primes at the same \(n\), and choose integers \(u,t\) coprime to \(R\).
Prescribe
\[
 y_{n-1}=\frac{\alpha_nu}{R},\qquad
 y_n=\alpha_nt,\qquad
 y_{n+1}=\beta_nt-\frac{\epsilon_nu}{R}.
\tag{16.29}
\]
This triple satisfies (16.26) identically.  For every \(p\mid R\), the
two polar reduced numerators are units, so
\[
 v_p(\operatorname {den}y_{n-1},
     \operatorname {den}y_n,
     \operatorname {den}y_{n+1})=(1,0,1).
\tag{16.30}
\]
The two adjacent values determine a unique global rational solution in
both directions.  Its reduced-denominator curvature at \(n\) is divisible
by the arbitrarily prescribed product \(R\).

Thus even simultaneous holes at all candidate primes are compatible with
the exact primitive normalized row.  Since \(\alpha_n\) is a \(p\)-unit,
the coefficient map
\[
 (X_-,X_0,X_+)\longmapsto
 \epsilon_nX_- -\beta_nX_0+\alpha_nX_+
\]
is surjective over \(\mathbb Z_p\); its cokernel and zeroth Fitting ideal
carry no target torsion.  A coefficient discriminant or any
state-independent bounded-window Fitting construction therefore cannot
detect the regular hole.

The solution in (16.29), however, depends
on \(n\) and \(R\).  The construction rules out a local estimate uniform
over rational initial states; it does **not** rule out a global theorem
whose constants or proof use the one fixed initial state
\((b_0,b_1)=(1,5)\).

This quantifier distinction can itself be made sharp.  The model can be
propagated through any fixed window of the exact Apéry rows.  It can also
be joined to any protected initial prefix, with the exact initial values,
by altering two intervening steering rows and then restoring the Apéry
operator before the target window.  Hence generic D-finiteness, the two
initial constants viewed in isolation, and exact local recurrence geometry
still permit arbitrary holes.  What the model cannot preserve is both the
exact Apéry row at **every** intervening index and the exact initial vector:
those data uniquely determine the actual solution.  A valid positive
theorem may use precisely that full global connection.

There is a separate exact warning against invoking only the abstract
\(G\)-function class.  Define
\[
 h_k=
 \begin{cases}
  1,&k\ \text{even},\\
  2^{-k},&k\ \text{odd}.
 \end{cases}
\]
Its generating series is rational:
\[
 \sum_{k\ge0}h_kx^k
 =\frac1{1-x^2}+\frac{x/2}{1-x^2/4}.
\tag{16.31}
\]
If \(Q_k\) is the reduced denominator of \(h_k\), then at \(k=2m\)
\[
 (Q_{k-1},Q_k,Q_{k+1})
 =(2^{2m-1},1,2^{2m+1}),
\]
so
\[
 \frac{\gcd(Q_{k-1},Q_{k+1})}
      {\gcd(Q_{k-1},Q_k,Q_{k+1})}
 =2^{2m-1}.
\tag{16.32}
\]
Even a rational \(G\)-function can therefore have exponential local
denominator curvature at infinitely many indices.  This example is
reducible and does not share the Apéry operator; its role is only to show
that a class-level denominator theorem needs substantial extra hypotheses.

Any proof of (16.22) must consequently use the distinguished initial state
or an arithmetic property equivalent in strength.  The normalized
recurrence and local denominator propagation alone do not suffice, while a
global theorem specific to the explicit Hadamard \(G\)-function remains a
logically valid opening.

### 16.5 Complete primitive denominator clearing is still exponentially tall

The recurrence does give a genuine state-dependent integer divisible by
the whole curvature carrier, but its evaluated height is already linear.
Write the three actual reduced fractions as
\[
 y_{n+j}=\frac{u_j}{q_{n+j}},
 \qquad \gcd(u_j,q_{n+j})=1,\qquad j=-1,0,1,
\]
and put
\[
 g=\gcd(q_{n-1},q_n,q_{n+1}),\qquad d=\mathfrak d_n.
\]
There are positive integers \(x_-,x_0,x_+\) such that
\[
 q_{n-1}=gdx_-,\qquad
 q_n=gx_0,\qquad
 q_{n+1}=gdx_+,
\qquad
 \gcd(x_-,x_+)=\gcd(d,x_0)=1.
\tag{16.33}
\]
Clearing (16.26) by exactly these reduced denominators gives
\[
 x_0K_n^{\rm hole}
 =\beta_nu_0\,d\,x_-x_+,
 \qquad
 K_n^{\rm hole}
 =\alpha_nu_{+1}x_-+\epsilon_nu_{-1}x_+.
\tag{16.34}
\]
Hence
\[
 \boxed{\mathfrak d_n\mid K_n^{\rm hole}.}
\]
This has removed the common triple denominator, the common hole factor
from both polar denominators, their coprime residual cofactors, and the
primitive row content.  There is no hidden factorial gauge in (16.34).

Nevertheless all terms in \(K_n^{\rm hole}\) are positive, and
\[
 K_n^{\rm hole}\ge\alpha_nu_{+1}
 \ge\alpha_ny_{n+1}.
\]
Since
\[
 y_n=\frac{b_n}{C_n}
 =\Theta\!\left(
   \left(\frac{17+12\sqrt2}{2}\right)^n n^{-1}\right),
\]
one obtains the exact route-specific height obstruction
\[
 \log K_n^{\rm hole}
 \ge
 n\log\!\left(\frac{17+12\sqrt2}{2}\right)-O(\log n).
\tag{16.35}
\]
The raw factorial transfer has \(n\log n\) height; primitive
central-binomial removal lowers it to the linear scale (16.35), not to
\(o(n)\).  This does not prove that every nonlinear global carrier is
large.  It proves that the canonical fully reduced one-row carrier cannot
close Path D.

The recurrence permits the two polar terms of valuation \(-1\) to
cancel and produce a \(p\)-integral middle term; after primitive
numerator-denominator clearing this is again the congruence
\(p\mid b_r\), not a polynomial singular-factor condition.

The exact reproducer checks the hypergeometric coefficient ratio, the
residue-section prime anchors, and (16.19):

```text
../scripts/q32_inverse_denominator_audit.py --limit 2500
```

Its default audit passes 219,903 candidate incidences and finds 248
exact holes.  It also checks 2,500 primitive recurrence rows, 2,500
identities (16.34), and 1,000 prescribed simultaneous-hole valuations.

## 17. Higher tensor contraction: an exact bounded-class no-go

There is a precise algebraic conservation law behind the failure of
higher Casoratian and Plücker forms on a pure-cross pair.  It is useful
because it closes all bounded-degree constructions whose smallness is
forced only by the universal dominant rank-one limit.

Use the top-half notation (PC.1): \(p<\ell=p+h\),
\[
 p\mid b_{s+h},\qquad \ell\mid b_s,
\]
and put
\[
 (A,B)=(a_s,b_s),\qquad(C,D)=(a_{s+h},b_{s+h}),
 \qquad\Delta=AD-BC.
\tag{17.1}
\]
On the pure-cross locus,
\[
\begin{array}{c|c|c}
 &\bmod\ell&\bmod p\\ \hline
 B&0&\ne0\\
 D&\ne0&0\\
 \Delta&AD\ne0&-BC\ne0.
\end{array}
\tag{17.2}
\]
After the standard gap-continuant denominator clearing, an integral
multiple \(\delta_{s,h}\) of \(\Delta\) has
\[
 \log|\delta_{s,h}|=O(h\log n),
\tag{17.3}
\]
but (17.2) says that it is a unit in both target characteristics.

The first primitive cross-divisible contraction occurs in degree three:
\[
 \boxed{
 F_3=\delta_{s,h}\{pB-\ell D\},\qquad p\ell\mid F_3.
 }
\tag{17.4}
\]
Writing \(B=\ell\beta\), \(D=p\gamma\) proves the divisibility.  The
coefficient pair \((p,-\ell)\) is primitive, so the factor \(p\ell\)
has not simply been inserted into every coefficient.  Nevertheless
\(F_3\) is exponentially large in the base index.  The elementary
monotonic estimate \(b_{k+1}\ge5b_k\) gives
\[
 |pB-\ell D|
 =\ell D-pB
 \ge(\ell5^h-p)B,
\tag{17.5}
\]
so one dominant Apéry slot remains after the single determinant
contraction.  For completeness, the estimate follows inductively from
\(b_k\ge b_{k-1}\) and
\[
 P(k)-k^3-5(k+1)^3
 =28k^3+36k^2+12k\ge0.
\]

The algebraic reason persists in every degree.  Work in
\[
 R=\mathbb Z[A,B,C,D].
\]
Formal divisibility under the two pure-cross specializations is the
ideal
\[
 {\cal C}_{p,\ell}
 =(\ell,B)\cap(p,D)
 =(\ell,B)(p,D).
\tag{17.6}
\]
Universal dominant-diagonal cancellation is divisibility by
\[
 \Delta=AD-BC,
\]
the defining equation of the rank-one \(2\times2\) matrix locus.
Modulo \((\ell,B)\), \(\Delta\) becomes \(AD\); modulo \((p,D)\), it
becomes \(-BC\).  Each is a non-zero-divisor in the corresponding
polynomial domain.  Therefore, for every \(m\ge0\),
\[
 \boxed{
 ({\cal C}_{p,\ell}:\Delta^m)={\cal C}_{p,\ell}.
 }
\tag{17.7}
\]

Suppose a homogeneous tensor form has been reduced to the two endpoint
states, all recurrence-only Smith content has been saturated, and the
largest completely contracted factor is \(\Delta^m\):
\[
 F=\Delta^mG.
\tag{17.8}
\]
Equation (17.7) says that pure-cross divisibility must still be carried
by \(G\).  If \(\deg G>0\), at least one unpaired dominant slot remains.
If \(\deg G=0\), then
\[
 {\cal C}_{p,\ell}\cap\mathbb Z=p\ell\,\mathbb Z,
\tag{17.9}
\]
so the fully contracted coefficient already has content \(p\ell\).
For example,
\[
 \Delta BD
\tag{17.10}
\]
is a primitive cross-divisible quartic but retains two dominant slots,
whereas the completely contracted quartic \(\Delta^2\) is a unit modulo
both primes; its cross-divisible line is only
\[
 p\ell\,\mathbb Z\Delta^2.
\tag{17.11}
\]

Nearby states and additional Plücker brackets do not enlarge the
completely contracted algebra.  Regular recurrence transport reduces
every nearby state to the same two endpoint states, and every bracket
is a recurrence scalar times \(\Delta\).  Plücker relations only relate
those scalars.  Singular-step factors must first be divided out when
they occur for every initial state; after this primitive saturation they
do not become target-selective.

Hence bounded-degree higher tensors obey the exact conservation law:

> Primitive pure-cross selectivity leaves a dominant tensor slot;
> complete universal dominant cancellation leaves coefficient content
> \(p\ell\).

This no-go has a deliberate scope boundary.  It does not cover an
identity whose leading cancellation uses the distinguished Apéry initial
state in a way not valid on the full rank-one diagonal, a genuinely new
opposite-quotient equation, or a nonlocal arithmetic relation.  It does
close the proposal that degree three, degree four, or further universal
symmetric-power/Plücker contractions alone will repair Route A.

### 17.1 What an initial-state-specific bounded tensor would have to prove

The preceding scope boundary can be made arithmetic rather than left as a
formal exception.  Put
\[
 \lambda=17+12\sqrt2.
\]
For each fixed shift \(j\), the two distinguished solutions have the
Poincaré asymptotics
\[
\begin{aligned}
 b_{s+j}
   &=\kappa\,\lambda^{s+j}s^{-3/2}
       \{1+O_j(s^{-1})\},\\
 a_{s+j}
   &=\zeta(3)b_{s+j}
       +O_j(\lambda^{-s}s^{-3/2}),
\end{aligned}
\tag{17.12}
\]
with analogous full inverse-power expansions.  The same conclusion holds
uniformly for \(j=O(\log s)\), after allowing subexponential coefficient
height.

Consider a fixed-degree homogeneous tensor expression in finitely many
shifted endpoint states, with rational recurrence coefficients, after all
universal determinant factors have been removed.  Before specializing the
relative endpoint scales, its restriction to collinear state vectors is a
polynomial
\[
 {\cal P}(X;t_1,\ldots,t_e),
\tag{17.13}
\]
where \(X\) is the common projective direction and the \(t_i\) are scalar
ratios between endpoint magnitudes.  Substitute the dominant part
\[
 (a_{s+j},b_{s+j})
 \rightsquigarrow
 (\zeta(3),1)\,\kappa\lambda^{s+j}s^{-3/2}
\tag{17.14}
\]
and divide by the common exponential and power of \(s\).  Its first
nonzero asymptotic symbol has the form
\[
 P_h(\zeta(3)),
 \qquad
 P_h(X)\in\mathbb Q(\sqrt2)[X],
 \qquad
 \deg P_h\le d,
\tag{17.15}
\]
where \(d\) is the tensor degree and \(h\) denotes the finitely many
endpoint gaps.  The quadratic coefficient field occurs because every
\(\lambda^h\) lies in \(\mathbb Q(\sqrt2)\); explicitly,
\[
 P_h(X)={\cal P}(X;\lambda^{h_1},\ldots,\lambda^{h_e}).
\]

If \({\cal P}(X;t_1,\ldots,t_e)\) vanishes identically before the
\(t_i\) are specialized, the cancellation holds on the entire rank-one
matrix locus.  Algebraically it is therefore in the determinantal ideal
already saturated in (17.7).  If only \(P_h(X)\) vanishes identically,
the cancellation may instead use the characteristic root
\(\lambda\).  A bounded template then removes only finitely many terms of
the inverse-power asymptotic expansion and leaves the same exponential
slot.  Eliminating that slot to every order would require either an exact
recurrence identity, whose universal transfer factors must again be
saturated, or a growing-order construction.

If \(P_h\) is nonzero but the distinguished leading term vanishes, then
\[
 P_h(\zeta(3))=0.
\tag{17.16}
\]
Such a bounded tensor would consequently prove that \(\zeta(3)\) is
algebraic over \(\mathbb Q(\sqrt2)\), of degree at most \(d\).  Already the
linear case would prove
\(\zeta(3)\in\mathbb Q(\sqrt2)\), which is not excluded by Apéry's
rational irrationality theorem.  For \(d\ge2\), no algebraicity result of
the required kind is known either.

There is a second possibility when the tensor degree, the gaps, or the
coefficient template grows with \(s\): the polynomials \(P_s\) may vary and
\(P_s(\zeta(3))\) may be very small without being zero.  To turn such a
construction into a subexponential-height contraction one needs a lower
bound of the shape
\[
 |P(\zeta(3))|
 \ge \exp\{-o(s)\}
\tag{17.17}
\]
for the particular polynomials of degree and height produced by the
template.  The known rational irrationality measures concern degree-one
polynomials over \(\mathbb Q\); they do not give (17.17) over
\(\mathbb Q(\sqrt2)\), much less at growing degree.

Thus “use the distinguished initial state” is not a routine loose end in
the tensor route.  At bounded degree it confronts an unknown algebraicity
statement for \(\zeta(3)\); at growing degree it confronts an unknown
polynomial-approximation measure.  An escape which avoids both must obtain
smallness from a new arithmetic mechanism rather than from cancellation of
the dominant Apéry line alone.

## 18. Audit of the \(S(p)\)/Katz handoff (2026-07-30)

The six-page handoff `P32_proof_v7.tex` proposes to finish Part B from a
uniform bound on
\[
 {\cal S}_a(p)=
 \{1\le m<p:\ p\mid P_{mp},\ p\nmid b_{mp}\}.
\tag{18.1}
\]
Its matrix statistic is
\[
 \Phi_m=M_{m-1}\cdots M_1,\qquad
 a_m\equiv6(\Phi_m)_{11}\pmod p.
\tag{18.2}
\]
Computation through \(p<10^5\) gives a maximum of \(11\) for this set and a
clean paired-Poisson histogram.  None of this closes the pointwise problem.
There are three independent fatal gaps.

### 18.1 The proposed multiplier set omits the lower Apéry channel

For \(\sqrt n<p\le n\), write
\[
 n=qp+r,\qquad0\le r<p.
\]
The audited block identity is
\[
 p^3a_{qp+r}\equiv a_qb_r\pmod p,
\tag{18.3}
\]
and hence, apart from the already separated harmless small primes,
\[
 p\mid G_n
 \quad\Longleftrightarrow\quad
 a_qb_r\equiv0\pmod p.
\tag{18.4}
\]
Thus there are two channels:
\[
 p\mid a_q
 \qquad\hbox{or}\qquad
 p\mid b_r.
\tag{18.5}
\]
The first is the companion channel.  The height/short-interval split at
\(q=n^{1/3}\) already proves that its total logarithmic weight is
\(O(n^{2/3})\).  Indeed, for small \(q\) use
\(\log|\operatorname {num}(a_q)|=O(q)\); for large \(q\), every relevant
prime is at most \(n^{2/3}\), so Chebyshev gives \(O(n^{2/3})\).

The statistic (18.1)--(18.2) concerns precisely this already controlled
companion channel.  It says nothing about the second channel in (18.5),
which is the actual obstruction.

There is a minimal exact counterexample to the v7 Part-B classification:
\[
 n=16,\qquad p=11,\qquad q=1,\qquad r=5.
\]
Here
\[
 b_5=819005\equiv0\pmod {11},
\]
and exact recurrence arithmetic gives
\[
 P_{16}\equiv Q_{16}\equiv0\pmod {11},
 \qquad v_{11}(G_{16})=1.
\tag{18.6}
\]
On the other hand,
\[
 P_{11}\equiv6\pmod {11},\qquad
 b_{11}\equiv5\pmod {11},
\tag{18.7}
\]
so \(1=\lfloor16/11\rfloor\notin{\cal S}_a(11)\).
Consequently a Part-B prime occurs outside the proposed bad-multiplier
streaks.  This refutes the Part-B proposition before one asks how large
\({\cal S}_a(p)\) is.

For \(p>n/2\), (18.4) has \(q=1\) and \(a_1=6\), so the exact remaining
support is instead
\[
 p\mid b_{n-p}.
\tag{18.8}
\]
Lucas gives, away from the small exceptional characteristic,
\[
 \{p\in(n/2,n]:p\mid G_n\}
 =
 \{p\in(n/2,n]:p\mid b_n\}.
\tag{18.9}
\]
This is the top-half large-prime radical problem isolated earlier.

### 18.2 Even an absolute row bound does not imply a pointwise column bound

Suppose, contrary to the preceding channel mismatch, that one had found
sets \(S(p)\) which really parametrized all bad residues and had proved
\(|S(p)|\le C\) for every prime.  The v7 proof still replaces a
deterministic count at a fixed \(n\) by the expectation
\(\sum_p(C-1)/p\).  This is invalid.

For a clean abstract countermodel, choose a rapidly increasing sequence
\(N_j\) such that the intervals \((N_j/2,N_j]\) are disjoint.  For each
prime in the \(j\)-th interval put
\[
 S(p)=\{1,p-2,p-1\}.
\tag{18.10}
\]
This has constant size, contains the claimed forced endpoint, and respects
the reflection pair \(1\leftrightarrow p-2\).  Yet at \(n=N_j\),
\(\lfloor n/p\rfloor=1\) for every prime in the top half, so
\[
 \sum_{\substack{n/2<p\le n\\
                  \lfloor n/p\rfloor\in S(p)}}\log p
 =\vartheta(n)-\vartheta(n/2)
 \sim\frac n2.
\tag{18.11}
\]
The actual Apéry sets have extra arithmetic constraints, but using them
requires a cross-prime theorem.  A one-prime cardinality estimate, even an
absolute one, cannot supply that theorem.

The reported Poisson\((1/2)\) fit points in the same direction.  It predicts
bounded mean and excellent tails for an individual prime, but the maximum
over the primes up to \(x\) grows slowly, on the scale
\[
 \frac{\log\pi(x)}{\log\log\pi(x)},
\tag{18.12}
\]
and is therefore unbounded.  A maximum of \(11\) below \(10^5\) is
compatible with this model; it is not evidence for a universal constant.

### 18.3 The claimed Katz theorem has the wrong object and the wrong output

The Beukers--Peters monodromy describes the K3 Picard--Fuchs local system as
its geometric base parameter varies.  In contrast, \(\Phi_m\) in (18.2) is
a variable-length prefix product in the coefficient index of one
recurrence solution.  It is not a Frobenius conjugacy class of that K3
sheaf at the base point \(m\).

No bounded-conductor lisse sheaf representing
\[
 m\longmapsto\Phi_m
\tag{18.13}
\]
has been constructed.  In fact \((\Phi_m)_{11}=a_m/6\), and direct
interpolation over \(\mathbb F_p\) has degree of order \(p\).  This is
exactly the growing-complexity behavior excluded by the claimed
bounded-monodromy application.

There is also a quantitative error.  Even if one independently proved,
for all nontrivial additive characters,
\[
 \left|\sum_{m\in\mathbb F_p}
   \psi\big((\Phi_m)_{11}\big)\right|
 \le C\sqrt p,
\tag{18.14}
\]
Fourier inversion would give only
\[
 \#\{m:(\Phi_m)_{11}=0\}=1+O(\sqrt p),
\tag{18.15}
\]
not \(O(1)\).  A genuinely bounded-degree algebraic map would have \(O(1)\)
zeros for the elementary reason that a nonzero bounded-degree function has
only boundedly many roots; the measured degree \(\asymp p\) is precisely
why that argument is unavailable.

Accordingly, the phrase “Katz equidistribution for products of
algebraically varying matrices” is not a citable proof.  It would require
an exact theorem, an actual construction matching (18.13), a conductor
bound uniform in \(p\), and then a separate cross-prime argument.  The
usual Deligne/Katz square-root theorem supplies none of the last two
steps merely from the Picard--Fuchs monodromy group.

### 18.4 The self-similarity congruence is proved elsewhere, but not by v7

Small-prime exact computation supports
\[
 p^3\frac{a_{mp}}{b_{mp}}
 \equiv\frac{a_m}{b_m}\pmod p
\tag{18.16}
\]
when the displayed ratios are \(p\)-integral after the indicated scaling.
The canonical `Q3.2_density_theorem.md` in fact proves it.  Its explicit
companion formula, together with the two Kummer valuation tables, proves
\[
 p^3a_{mp}\equiv a_m\pmod p.
\tag{18.16a}
\]
Gessel--Lucas gives \(b_{mp}\equiv b_m\pmod p\), so division by the latter
unit gives (18.16).  The same argument followed by recurrence propagation
in \(r\) proves the complete block law (18.3).

The v7 telescoping proof is nevertheless invalid as written.  Terms for
which an intermediate \(b_j\) vanishes have more negative \(p\)-adic
valuation, not less; saying they are “dominated” by the generic valuation
reverses the valuation comparison.  The asserted integrality of the
complementary block sum is not proved either.

There is a second overstatement: congruence (18.16) is not formally
equivalent to exact equality of valuations when its right-hand side
vanishes modulo \(p\).  What follows rigorously is the mod-\(p\) set
identity needed to recognize the companion channel.  This correction still
cannot rescue v7, because Sections 18.1 and 18.2 refute its claimed Part-B
consequence independently.

### 18.5 Correct surviving target

Put
\[
 T_n=\{p\in(n/2,n]:p\mid b_n\},\qquad K_n=|T_n|.
\tag{18.17}
\]
A genuine sufficient theorem is the shell-uniform pair-energy estimate
\[
 H_2(N):=
 \sum_{N<n\le2N}(K_n)_2
 =o\!\left(\frac{N^2}{\log^2N}\right).
\tag{18.18}
\]
Indeed,
\[
 \max_{N<n\le2N}K_n
 \le1+\sqrt{H_2(N)}
 =o(N/\log N),
\tag{18.19}
\]
so every top-half logarithmic radical is \(o(n)\).  The dyadic-shell
quantifier in (18.18) is essential; a cumulative estimate along a sparse
sequence of endpoints need not control every shell.

For the full lower channel, the sharp formulation is quotientwise.  For a
fixed \(q\ge1\), put
\[
 K_q(n)=
 \#\left\{p:
   \frac n{q+1}<p\le\frac nq,\quad
   p\mid b_{n-qp}\right\},
\tag{18.19a}
\]
and
\[
 H_{2,q}(N)=\sum_{N<n\le2N}(K_q(n))_2.
\tag{18.19b}
\]
It is sufficient to prove, separately for every fixed \(q\),
\[
 H_{2,q}(N)=o_q\!\left(\frac{N^2}{\log^2N}\right).
\tag{18.19c}
\]
No uniformity in growing \(q\) is needed.  Indeed, first discard
\(p\le\eta n\), whose total logarithmic mass is \(O(\eta n)\).  The
remaining primes have \(q<1/\eta\), only finitely many quotient classes.
For each of them, nonnegativity and (18.19c) give
\[
 \max_{N<n\le2N}K_q(n)=o_q(N/\log N).
\]
Sum over the finite set of \(q\), then let \(\eta\downarrow0\).  Together
with the \(O(n^{2/3})\) companion channel and the small-prime estimate, this
proves the full \(\log G_n=o(n)\).

The fixed-gap Selberg carriers prove sublinear height for each fixed gap,
but charging all \(h\le A\log n\) separately has total height \(O(n)\).
Radicalizing, harmonic formal exponents, or dyadic bookkeeping cannot
remove this cost: every target prime incident to a short edge must still
occur with an integer exponent somewhere.  A positive continuation needs
an Apéry-specific inverse-gap saving or a common horizontal rank condition,
not another local equidistribution claim.

### 18.6 Irrationality measures cannot remove the denominator-clearing rate

There is also no escape through the real quality of the Apéry
approximants.  Put
\[
 \lambda=17+12\sqrt2,\qquad L=\log\lambda.
\]
The standard asymptotics give
\[
 \log Q_n=(3+L)n+o(n),\qquad
 \left|\zeta(3)-\frac{P_n}{Q_n}\right|
 =\exp\{-2Ln+o(n)\}.
\tag{18.20}
\]
If \(\mu_0\) is any proved finite upper bound for the irrationality
exponent of \(\zeta(3)\), apply the corresponding lower bound to the
reduced denominator \(q_n=Q_n/G_n\).  One obtains exactly
\[
 \log G_n
 \le
 \left(3+L-\frac{2L}{\mu_0}\right)n+o(n).
\tag{18.21}
\]
With the Rhin--Viola value
\(\mu_0=5.51389062\ldots\), the coefficient in (18.21) is
\[
 5.2467260\ldots.
\tag{18.22}
\]
Even the conjecturally optimal value \(\mu_0=2\) leaves
\[
 \log G_n\le3n+o(n).
\tag{18.23}
\]
The irreducible loss is exactly the \(d_n^3\) denominator-clearing factor.
Thus no improvement of the irrationality exponent alone can yield
subexponential content.

The adjacent Casoratian gives the complementary pair budget
\[
 \log G_n+\log G_{n-1}\le6n+o(n),
\tag{18.24}
\]
but permits alternating exponential contents and gives no pointwise
estimate.  This limitation is formal: multiply the numerator and
denominator of the Fibonacci convergents by \(M^n\).  The rational
approximants remain exponentially accurate and satisfy a second-order
Poincaré recurrence; their adjacent determinants remain explicit and
supported on the fixed primes dividing \(M\), while their content is
\(M^n\).  Approximation quality plus a smooth nonzero Casoratian therefore
cannot control content without a separate primitive-lattice theorem.

The unconditional frontier is therefore the fixed-\(q\) family (18.19c);
(18.18) is its necessary \(q=1\) core.  Equivalently, one needs a
cross-prime large-radical theorem for the selected residues.  The
\(S(p)\)/Katz handoff does not advance that frontier.

### 18.7 All-degree saturation of the direct moving-zero carrier

There is an exact algebraic reason why repeated attempts to manufacture a
small integer directly from the individual target congruences return the
unknown radical itself.

Fix \(n>10\), and put
\[
 {\cal P}_n=\{p\text{ prime}:n/2<p\le n\},\qquad
 r_p=n-p,\qquad
 R_n=\prod_{\substack{p\in{\cal P}_n\\p\mid b_{r_p}}}p.
\tag{18.25}
\]
In
\[
 {\cal R}_n=\mathbb Z[X_0,\ldots,X_{\lfloor(n-1)/2\rfloor}]
\]
consider the moving-zero ideals
\[
 {\mathfrak j}_p=(p,X_{r_p}),\qquad
 {\mathfrak I}_n=\bigcap_{p\in{\cal P}_n}{\mathfrak j}_p.
\tag{18.26}
\]
The ideals \({\mathfrak j}_p\) are pairwise comaximal, since distinct
candidate primes generate the unit ideal in \(\mathbb Z\).  Consequently
\[
 {\mathfrak I}_n=\prod_{p\in{\cal P}_n}{\mathfrak j}_p.
\tag{18.27}
\]
Under evaluation \(X_j\mapsto b_j\), one has
\[
 \operatorname {ev}_b({\mathfrak j}_p)
   =(p,b_{r_p})
   =
   \begin{cases}
      p\mathbb Z,&p\mid b_{r_p},\\
      \mathbb Z,&p\nmid b_{r_p}.
   \end{cases}
\tag{18.28}
\]
Evaluation commutes with products of ideals, and hence
\[
 \boxed{\operatorname {ev}_b({\mathfrak I}_n)=R_n\mathbb Z.}
\tag{18.29}
\]
There is no degree restriction in (18.29).  Thus arbitrary polynomial
combinations, determinants, resultants, Smith or Fitting ideals, and
growing tensor or exterior powers cannot produce an independent small
multiple of \(R_n\) when their divisibility proof uses only the separate
formal implications
\[
 p\mid b_{r_p}.
\]
After complete saturation, the evaluation ideal is exactly the radical one
is trying to bound.

The same obstruction is visible without polynomial language.  If
\[
 \Pi_n=\prod_{p\in{\cal P}_n}p,
\]
then
\[
 \boxed{
 R_n=
 \gcd\!\left(
 \Pi_n,\
 \left\{\frac{\Pi_n}{p}\,b_{n-p}:p\in{\cal P}_n\right\}
 \right).}
\tag{18.30}
\]
Indeed, for a candidate \(\ell\), every entry except the one indexed by
\(\ell\) already contains \(\ell\), while the exceptional entry contains
\(\ell\) exactly when \(\ell\mid b_{n-\ell}\).  Formula (18.30) is an exact
universal carrier, but every input has logarithmic height
\(\Theta(n)\), and its saturated gcd is precisely \(R_n\).  It therefore
does not prove \(o(n)\).

One can also quotient a finite recurrence window by all universal Apéry
recurrence relations, reducing it to the two initial-state coordinates.
The evaluation image remains (18.29): recurrence identities merely change
the presentation of the same local target ideal.  This extends the earlier
bounded-template no-go to arbitrary degree and arbitrary recurrence
window.  Its scope must nevertheless be stated accurately.  It does not
exclude a second target-selective congruence special to the distinguished
Apéry initial vector, nor an independently proved Archimedean
near-cancellation.

### 18.8 Signed Newton carriers and exact central-alias cancellation

The constant-term model gives a particularly tempting way to seek such a
cancellation.  Write
\[
 d_k=\Delta^k b_0=\operatorname {CT}(\Lambda-1)^k,\qquad
 b_r=\sum_{k=0}^r\binom rk d_k.
\tag{18.31}
\]
The constant monomial of \(\Lambda-1\) is \(4\).  Set
\[
 m=\left\lfloor\frac{n-1}{2}\right\rfloor,\qquad
 Z_4(c)=\sum_{k=0}^m c_k4^k.
\]
Suppose a rational linear carrier, after clearing a denominator prime to
every \(p\in{\cal P}_n\), has coefficient vector \(c=(c_0,\ldots,c_m)\)
and its modular proof is the formal implication
\[
 \sum_{k=0}^{r_p}\binom{r_p}{k}d_k=0
 \quad\Longrightarrow\quad
 \sum_{k=0}^m c_kd_k=0
 \pmod p.
\tag{18.32}
\]
The kernel of the first nonzero linear form in (18.32) has codimension
one.  Therefore the second form must be a scalar multiple of it, giving
\[
 c_k\equiv c_0\binom{r_p}{k}
      \equiv c_0\binom nk\pmod p
\quad(0\le k\le m).
\tag{18.33}
\]
Chinese remaindering over all candidates yields
\[
 c_k=c_0\binom nk+\Pi_nt_k,\qquad t_k\in\mathbb Z.
\tag{18.34}
\]
Moreover,
\[
 Z_4(c)\equiv
 c_0\sum_{k=0}^{r_p}\binom{r_p}{k}4^k
 =c_0\,5^{r_p}\pmod p.
\tag{18.35}
\]
All asymptotically relevant candidates are at least \(7\).  Hence exact
central-alias cancellation \(Z_4(c)=0\) forces \(p\mid c_0\) for every
candidate \(p\), and then (18.33) forces
\[
 \Pi_n\mid c_k\qquad(0\le k\le m).
\tag{18.36}
\]
Every nonzero integer supplied by such a \(p\)-safe linear carrier
therefore has
\[
 \log|C|\ge\log\Pi_n
 =\vartheta(n)-\vartheta(n/2)
 =\left(\frac12+o(1)\right)n.
\tag{18.37}
\]
Signed coefficients do not remove the linear height.

There is an all-degree version of the same transversality.  Replace
\({\mathfrak j}_p\) in (18.26) by
\[
 {\mathfrak J}_p=
 \left(p,\sum_{k=0}^{r_p}\binom{r_p}{k}X_k\right),
\qquad
 {\mathfrak I}^{N}_n=\bigcap_p{\mathfrak J}_p,
\tag{18.38}
\]
and let \({\mathfrak M}_4\) be the kernel of \(X_k\mapsto4^k\).
Since the Newton form evaluates at this alias to
\[
 \sum_{k=0}^{r_p}\binom{r_p}{k}4^k=5^{r_p}\not\equiv0\pmod p,
\]
one has
\[
 {\mathfrak I}^{N}_n+{\mathfrak M}_4=(1),\qquad
 {\mathfrak I}^{N}_n\cap{\mathfrak M}_4
   ={\mathfrak I}^{N}_n{\mathfrak M}_4.
\tag{18.39}
\]
After evaluation \(X_k\mapsto d_k\),
\[
 \operatorname {ev}_d(
 {\mathfrak I}^{N}_n\cap{\mathfrak M}_4)
 =R_n g_m\mathbb Z,\qquad
 g_m=\gcd_{0\le k\le m}(d_k-4^k).
\tag{18.40}
\]
For \(m\ge3\),
\[
 d_2-4^2=48,\qquad d_3-4^3=1176,
\]
so \(g_m\mid24\).  Exact central-alias cancellation changes the unknown
radical by at most this fixed factor; it gives no independent height
contraction.  A positive toric route must use a quantitative
near-cancellation among the actual Apéry aliases, or a genuinely second
congruence, rather than an exact formal alias identity.

### 18.9 The factorial-gcd reformulation is stronger, not standard

A natural reformulation is
\[
 \gcd(b_n,n!).
\tag{18.41}
\]
For squarefree support, this adds no help.  Every prime \(p\le\sqrt n\)
has total logarithmic weight \(O(\sqrt n)\), while for \(p>\sqrt n\)
membership in the radical of (18.41) is exactly the same moving
large-prime support already isolated above.  The central-binomial
restricted radical is therefore the cleaner interface.

The full gcd in (18.41) is substantially harder because it retains
\[
 \min\{v_p(b_n),v_p(n!)\}
\]
at small primes.  The \(p\)-Lucas law controls only support modulo \(p\).
Lifting from \(p\) to \(p^2\) introduces a new derivative or jet
coordinate, and higher powers introduce further data.  A sufficient local
estimate such as
\[
 v_p(b_n)=O(\log_p n)
\tag{18.42}
\]
uniformly in \(p,n\) would make the small-prime valuation contribution
\(O(\sqrt n)\), but no result audited here proves (18.42).  Known
factorial-ratio valuation theorems generally certify divisibility rather
than upper-bound these singular depths.

Nor is the desired content estimate a hidden step in the standard proofs
of Apéry's theorem.  Van der Poorten's exposition clears denominators with
\(d_n^3\), estimates the resulting integral numerator and denominator, and
applies the elementary irrationality criterion to those unreduced
integers.  Beukers replaces the recurrence construction by a triple
integral with the same lcm-denominator arithmetic.  Zudilin's
hypergeometric approach controls denominators of linear forms.  None of
these arguments needs the fraction to be in lowest terms, and the sources
audited here do not state a subexponential bound for
\[
 \gcd(d_n^3a_n,d_n^3b_n).
\tag{18.43}
\]
Thus there is no standard integral-, congruence-, or
hypergeometric proof of (18.43) available to import.  More cautiously:
the literature search performed for this audit located no published
standalone theorem proving it.  The pointwise cross-prime estimate
(18.19c), not a conventional denominator lemma, remains the missing
input.

### 18.10 Why the observed bound \(K_n\le3\) is not a credible lemma

The scan through \(n\le2\cdot10^6\) finds no top-half column with four
targets.  This is useful experimental information, but the proved local
package has no finite-multiplicity content.

Indeed, fix \(K\).  For arbitrarily large \(n\), choose \(K\) distinct
primes \(p_i\in(n/2,n]\), put \(r_i=n-p_i\), and prescribe the abstract
reflected doublets
\[
 Z_{p_i}=\{r_i,p_i-1-r_i\}.
\tag{18.44}
\]
The primes can be chosen greedily to avoid:

1. the three exceptional relations
   \(2r_i=p_i-2,p_i-1,p_i\), which make the reflected points consecutive
   or equal; and
2. for every earlier \(p_j\), the two relations
   \[
   p_i=2n+1-2p_j,\qquad p_j=2n+1-2p_i,
   \tag{18.45}
   \]
   which make a cross-pair degenerate.

Only finitely many candidates are excluded at each step, whereas the top
half contains arbitrarily many primes.  The resulting data satisfy:

- exact reflection and two zeros per active row;
- nonconsecutivity and every currently proved vertical zero bound;
- the target condition \(r_i\in Z_{p_i}\) for all \(i\);
- pure-cross status for every distinct pair; and
- only the universal reflected same-row continuant factor.

This is not a construction of Apéry zero sets.  It is a logical model of
all presently used local predicates, and it shows that those predicates
alone are compatible with arbitrary fixed \(K\).  Any theorem
\(K_n=O(1)\), if true for the actual Apéry sequence, must use new
cross-characteristic information about the distinguished initial state.
The natural candidates are divided quotient jets, a bounded global
eliminant, or a pointwise old-prime theorem; none is currently available.

The theorem actually needed is much weaker:
\[
 K_n=o(n/\log n).
\tag{18.46}
\]
The pure-cross amplification in PC.20--PC.21 is calibrated to (18.46).
A forbidden spike forces many adjacent pairs with gaps
\(h\ll\log n\).  The fixed-\(h\) Selberg estimate has the needed saving,
but summing over this logarithmic range loses it.  An averaged
Apéry-specific depletion over \(h\), rather than a universal bound of
three, is therefore the sharper positive target.

### 18.11 Exact mixed-characteristic dispersion barrier

The fixed-\(q\) criterion (18.19c) can be written without any probabilistic
language.  Let
\[
 z_p(x)=1_{\{b_x\equiv0\pmod p\}},\qquad 0\le x<p,
\]
and extend \(z_p\) periodically only for the purpose of Fourier expansion.
For an ordered collision \(p<\ell\), put
\[
 h=\ell-p,\qquad s=n-q\ell.
\]
Then the two lower digits are \(s\) and \(s+qh\), and the exact target
conditions are
\[
 z_\ell(s)z_p(s+qh)=1.
\tag{18.47}
\]
Consequently
\[
 H_{2,q}(N)
 =2\sum_{p<\ell}
   \sum_{\substack{s:\ N<q\ell+s\le2N\\
                   0\le s<\ell\\
                   0\le s+q(\ell-p)<p}}
 z_\ell(s)z_p\!\left(s+q(\ell-p)\right).
\tag{18.48}
\]
All primes in (18.48) lie in a \(q\)-dependent compact multiple of \(N\).
This is the correct bilinear incidence form; the two factors live in
different residue characteristics.

Put
\[
 \widehat z_p(a)=
 \sum_{x\bmod p}z_p(x)e_p(-ax).
\tag{18.49}
\]
Parseval gives
\[
 \sum_{a\bmod p}|\widehat z_p(a)|^2=p|Z_p|.
\tag{18.50}
\]
Substituting (18.49) in (18.48) produces a double prime-modulus Fourier
kernel.  A standard \(L^2\) large sieve can use only the row masses
(18.50).  With the proved \(|Z_p|\ll p^{2/3}\), its resulting scale is
\[
 O_q\!\left(\frac{N^{8/3}}{\log^2N}\right),
\tag{18.51}
\]
far above (18.19c).  More decisively, even the hypothetical local bound
\(|Z_p|\le C\) yields only
\[
 O_q\!\left(\frac{N^2}{\log^2N}\right),
\tag{18.52}
\]
the critical scale rather than little-\(o\).  The reflected-star model of
Section 18.10 saturates this use of row norms, so (18.52) is not repaired
by sharpening constants or by ordinary Parseval bookkeeping.

What would suffice is a genuinely Apéry-specific estimate
\[
 \sum_{p<\ell}\sum_s
 z_\ell(s)z_p(s+q(\ell-p))
 =o_q\!\left(\frac{N^2}{\log^2N}\right)
\tag{18.53}
\]
with precisely the arc and shell restrictions in (18.48), or a
power-saving smoothed version.  This is a mixed-characteristic local-limit
dispersion theorem, not a classical Barban--Davenport--Halberstam
statement: the tested sequence changes with the modulus.

Nor is \(z_p\) a known bounded-conductor trace function.  Detecting exact
vanishing uses
\[
 z_p(x)=\frac1p\sum_{u\bmod p}e_p(ub_x),
\tag{18.54}
\]
and the interpolation degree of \(x\mapsto b_x\) is of order \(p\).
The corresponding Artin--Schreier pullback therefore has growing
complexity.  Mellin descriptions of the underlying K3 trace do not fix
this: in the surviving nearly primitive sector the Kummer character order
also grows like \(p\).  Existing fixed-field bounded-conductor bilinear
trace estimates and large sieves for a fixed compatible family do not
apply to (18.53).

Thus (18.53) is a precise new theorem interface, not a theorem presently
available from Katz, BDH, or the standard prime-moduli large sieve.  Any
analytic claim closing P3.2 must exhibit cancellation among the actual
coefficients \(\widehat z_p(a)\) as the characteristic varies; cardinality
bounds for \(Z_p\) cannot substitute for it.

### 18.12 Fixed recurrence polynomials do not explain the targets

The validated binary scan through \(n\le2\cdot10^6\) contains \(106039\)
top-half incidences.  An exact parse tested divisibility of the natural
recurrence factors
\[
 n,\quad n+1,\quad 2n+1,\quad 17n^2+17n+5
\tag{18.55}
\]
at every target prime.  The first, second, and fourth factors cover no
target at all.  The factor \(2n+1\) covers only
\[
 (n,p,r)=(16,11,5),\qquad(4705,3137,1568),
\tag{18.56}
\]
the two exceptional self-reflected fibers \(r=(p-1)/2\).  In particular,
none of the \(159\) incidences in the \(53\) triple columns divides any
product of the factors in (18.55).

This is only a finite exclusion, not a theorem about arbitrary
polynomials.  It nevertheless rules out the most immediate conjecture
that the observed polynomial-size radicals are forced by a singular factor
of the recurrence.  The generic target support is invisible to all of
(18.55).

### 18.13 The first \(p\)-adic lift is still rank one

The most natural source of a second target equation is the first lift of
the two normalized gcd coordinates.  It does not work.  For a prime
\(p\ge5\) and \(0\le r<p\), define
\[
 D_r=
 2\sum_{k=0}^r
 \binom rk^2\binom{r+k}k^2
 \bigl(H_{r+k}-H_{r-k}\bigr).
\tag{18.57}
\]
Every term in (18.57) is \(p\)-integral.  If \(r+k\ge p\), the squared
binomial supplies \(p^2\), while the harmonic difference loses at most one
factor \(p\).

Gessel's shifted congruence, specialized to outer digit one, is
\[
 \boxed{b_{p+r}\equiv5(b_r+pD_r)\pmod {p^2}.}
\tag{18.58}
\]
The companion coordinate satisfies the parallel formula
\[
 \boxed{p^3a_{p+r}\equiv6(b_r+pD_r)\pmod {p^2}.}
\tag{18.59}
\]
Formula (18.59) can be proved without a new supercongruence.  The
differentiated Apéry recurrence shows that
\(X_r=b_r+pD_r\) satisfies the recurrence shifted by \(p\), modulo
\(p^2\).  The exact Wronskian gives
\[
 p^3a_p\equiv6\pmod {p^2},
\]
and the recurrence at \(p\) gives
\[
 p^3a_{p+1}\equiv6(5+12p)\pmod {p^2}.
\]
These two initial values propagate (18.59) throughout \(0\le r<p\),
because the relevant leading coefficients are \(p\)-units.

On a target write
\[
 b_r=pq_{p,r},\qquad
 \xi_{p,r}=q_{p,r}+D_r\pmod p.
\tag{18.60}
\]
Dividing (18.58)--(18.59) by \(p\) gives
\[
 \boxed{
 \left(\frac{b_{p+r}}p,\frac{p^3a_{p+r}}p\right)
 \equiv \xi_{p,r}(5,6)\pmod p.}
\tag{18.61}
\]
Thus the first lift has only the universal relation
\[
 6b_{p+r}-5p^3a_{p+r}\equiv0\pmod {p^2};
\tag{18.62}
\]
it supplies no second target-selective equation.  Formally, the map
\[
 \mathbb F_p[X,Y]\longrightarrow\mathbb F_p[\xi],
 \qquad X\mapsto5\xi,\quad Y\mapsto6\xi
\]
has kernel \((6X-5Y)\).  This is an arbitrary-degree rank-one statement,
not merely a failure of one proposed linear combination.

The remaining scalar is genuinely unconstrained by targetness.  Take
\[
 p=73,\qquad r=2,\qquad b_2=73,\qquad D_2=210.
\]
Then
\[
 \xi_{73,2}=1+210\equiv65\not\equiv0\pmod {73}.
\tag{18.63}
\]
Consequently \(73\) is a target at \(n=75\), but
\[
 v_{73}(G_{75})=1.
\tag{18.64}
\]
The first lift detects only the deeper event:
\[
 v_p(G_{p+r})\ge2
 \quad\Longleftrightarrow\quad
 \xi_{p,r}=0
 \quad\Longleftrightarrow\quad
 p^2\mid b_{p+r}.
\tag{18.65}
\]
A positive jet route must therefore discover a new relation between the
actual values \(b_r/p\) and \(D_r\), or a genuinely independent higher
coordinate.  Neither the support condition, the recurrence, the
Wronskian, nor the first Gessel lift provides it.

### 18.14 Pointwise short-gap dispersion: the exact analytic target

Section 18.11 records the shell pair-energy formulation.  The
pure-cross amplification gives a weaker, but more sharply pointwise,
interface.  Put
\[
 f_p(x)=1_{\{b_{x\bmod p}\equiv0\pmod p\}},\qquad z_p=|Z_p|,
\]
and for fixed \(A>0\), \(H=A\log N\).  A sufficient theorem is
\[
 \sup_{N<n\le2N}
 \left|
 \sum_{\substack{2\le h\le H\\2\mid h}}
 \sum_{\substack{n/2<p\le n-h\\p,\ p+h\ {\rm prime}}}
 \left(
 f_p(n)f_{p+h}(n)-\frac{z_pz_{p+h}}{p(p+h)}
 \right)
 \right|
 =o_A(N/\log N).
\tag{18.66}
\]
The mean term in (18.66) is already harmless.  The Selberg prime-pair
bound gives \(O_A(N/\log N)\) candidate pairs, and
\(z_pz_{p+h}/(p(p+h))=O(N^{-2/3})\), so their total is
\[
 O_A(N^{1/3}/\log N).
\tag{18.67}
\]
After subtracting the \(O_A(\log^2N)\) degenerate pairs, (18.66) implies
that the number of pure-cross pairs of gap at most \(A\log N\) is
\(o_A(N/\log N)\), uniformly in \(n\).

If \(K_n\ge\varepsilon N/\log N\), order the target primes.  Their total
span is at most \(N\), so at most \(N/(A\log N)\) adjacent gaps exceed
\(A\log N\).  Taking \(A>4/\varepsilon\) leaves
\(\gg_\varepsilon N/\log N\) short adjacent pure-cross pairs, contradicting
(18.66).  Thus (18.66) proves the top-half theorem.

Fourier inversion makes the missing cancellation explicit:
\[
 f_p(n)=\frac{z_p}{p}
 +\frac1p\sum_{a\ne0}\widehat f_p(a)e_p(an).
\tag{18.68}
\]
For \(\ell=p+h\), the centered product is the sum of
\[
 \frac{z_\ell}{p\ell}
 \sum_{a\ne0}\widehat f_p(a)e_p(an),
\qquad
 \frac{z_p}{p\ell}
 \sum_{b\ne0}\widehat f_\ell(b)e_\ell(bn),
\tag{18.69}
\]
and
\[
 \frac1{p\ell}
 \sum_{\substack{a\ne0\\b\ne0}}
 \widehat f_p(a)\widehat f_\ell(b)
 e_p(an)e_\ell(bn).
\tag{18.70}
\]
All three terms must cancel after summing nearby prime pairs, uniformly at
one fixed \(n\).

No collection of one-row Fourier estimates can prove this.  The reflected
two-point adversary has
\[
 |\widehat f_p(a)|\le2
\]
at every nonzero frequency, yet one can select a positive proportion of
the top-half primes through one column, discard the finite-degree
degeneracy graph, and retain \(\gg_A N/\log N\) adjacent pure-cross gaps
below \(A\log N\).  It satisfies stronger individual Fourier bounds than
Deligne would provide while violating (18.66).

The missing theorem is therefore genuinely two-characteristic and
actual-state-specific: cancellation in (18.69)--(18.70) as both \(p\) and
\(p+h\) vary.  Standard trace-function theorems fix one finite field and a
bounded-conductor sheaf; standard large sieves average the external
integer \(n\).  Neither matches the supremum and moving-modulus structure
of (18.66).

## 19. Synchronization with the two uisai2 research sessions

The two live P3.2 sessions on uisai2 were audited read-only on
2026-07-30.  Both effective sessions, `zinan:0` (`dm`) and `zinan:3`
(`family`), were running `gpt-5.6-sol max`.  The separate `dm-` Claude
window had only repeated HTTP 529/500 failures and no usable new result.
The uisai2 notes live in the separate `zinan-memory` repository, so the
proof-relevant conclusions are recorded here rather than treating that
repository as the canonical Ramanujan Challenge checkout.

### 19.1 CFVZ, marked traces, and growing exterior systems

The `dm` session gives a useful strengthening of the Route B no-go.  For
the CFVZ factorization
\[
 A_p(t)=g_p(t)B_p(t)^2\quad\text{in }\mathbb F_p[t],
\]
the coefficient map from \(B_p\) to \(A_p\) is triangular with diagonal
2.  After adjoining all formal Hensel grades, reflected quotient digits,
differentiated recurrence rows, and the moving-prefix forcing row, the
local target algebra still has the form
\[
 \mathcal A_{p,r}^{(Q,J)}
 \simeq(\mathbb Z/p^Q\mathbb Z)[X],
 \qquad
 {\rm SNF}=\operatorname {diag}(1,\ldots,1,0).
\tag{19.1}
\]
Here \(X=b_r/p\) is one cyclic nuisance coordinate.  The apparent new
digits are its successive digits, not independent target equations.  The
moving-prefix forcing coefficient is a unit; the actual target
\((n,p,r)=(20,17,3)\) has nonzero forcing residue \(2\bmod17\).  Thus the
entire adjacent CFVZ/Hensel tower remains a graph over one free coordinate.
It cannot supply the growing defect codimension requested in Route B.

A marked cyclotomic trace does produce, for each individual candidate
prime, a rational integer \(G_{p,d,e}\) of height \(O(\log n)\) such that
\[
 p\mid G_{p,d,e}\quad\Longleftrightarrow\quad p\mid b_{n-p}.
\tag{19.2}
\]
This is genuinely prime-selective locally.  It does not batch: retaining
the marked prime above \(p\) makes the scalar depend on \(p\), and rational
CRT recombination has full candidate-primorial cost.  After normalizing the
first residual trace, its centered CRT batch is exactly the class
\[
 b_n\bmod \prod_{n/2<p<n}p,
\tag{19.3}
\]
so its gcd with the candidate primorial is the original target radical.
This is an especially clean instance of the Route B quantifier
obstruction: a cheap private coordinate for every prime does not give one
cheap common coordinate.

The rank-12 exterior construction reaches the same boundary.  An
uncrossed top prime sees the two-coordinate state
\[
 (b_r,\eta_r),\qquad
 \eta_r=[z^r]\,B(z)\log(1+z),
\tag{19.4}
\]
and targetness forces only \(b_r=0\).  A window
\(h=o(n/\log n)\) either captures only the \(o(n)\)-weight boundary band
or requires the false implication \(b_r=0\Rightarrow\eta_r=0\).  Reaching
below every candidate needs a linear window, after which the full exterior
Fitting ideal returns the original prime-free gcd.  Hence neither growing
rank nor formal precision by itself escapes the Casoratian/defect
saturation.

### 19.2 A genuine universal second-layer congruence

The `family` session proved the most substantial new positive identity.
Set
\[
 \Theta_n=
 \sum_{0\le k\le\lfloor(n-1)/3\rfloor}
 \binom nk^2\binom{n+k}{k}^2
\tag{19.5}
\]
and
\[
 U_n=\prod_{\substack{n/2<p\le n\\p\ge7}}p.
\]
For every such prime,
\[
 b_n\equiv5\Theta_n\pmod {p^2},
\qquad\text{hence}\qquad
 U_n^2\mid b_n-5\Theta_n.
\tag{19.6}
\]
The congruence was independently rechecked locally through \(n=300\).
It gives the exact top-half support identity
\[
 p\mid b_n\quad\Longleftrightarrow\quad p\mid\Theta_n
 \qquad(n/2<p\le n),
\tag{19.7}
\]
and therefore a new holonomic-gcd interface
\[
 \prod_{p\in T_n}p
 =\operatorname {rad}_{(n/2,n]}
   \gcd(b_n,\Theta_n).
\tag{19.8}
\]

This does not yet give a subexponential carrier.  The congruence (19.6)
holds for every candidate prime, not only for a target, and
\[
 \log\left|\frac{b_n-5\Theta_n}{U_n^2}\right|
 =\bigl(\log(17+12\sqrt2)-1+o(1)\bigr)n.
\tag{19.9}
\]
The divided value is not target-divisible: at
\((n,p,r)=(20,17,3)\),
\[
 \frac{b_n-5\Theta_n}{p^2}\equiv3\pmod p.
\tag{19.10}
\]
Fixed-order neighboring determinants retain the same universal
candidate-square content and lose the target factor after saturation.
The session is testing a growing collection of third divided digits; the
natural cutoffs already have many target counterexamples, so no theorem
has emerged from that experiment.

The same session also supplied a corrected proof of
\(|Z_p|=O(p^{2/3})\) using nonzero gap continuants.  This agrees with
Section 10 and remains a one-characteristic row bound; it does not control
the moving column at a fixed \(n\).

### 19.3 Net effect on the live frontier

The uisai2 work is useful, but it does not close P3.2:

1. (19.1)--(19.4) close broad CFVZ, marked-trace, Hensel-jet, and growing
   exterior variants of Route B unless an identity outside their proved
   graph modules is found.
2. (19.6) is a real new second-layer congruence and (19.8) is a clean new
   exact carrier, but its known height remains linear.
3. Neither session proves the fixed-\(q\) pair energy (18.19c) or the
   two-characteristic dispersion estimate (18.66).  Those remain the two
   honest positive interfaces.

## 20. Audit of the proposed \(L^1+L^2\) decomposition

A subsequent handoff proposed splitting the top-half problem into
\[
 {\rm (L1)}\qquad
 \sum_{p\le x}z_p\ll\pi(x)
\tag{20.1}
\]
and
\[
 {\rm (L2)}\qquad
 K_n=o(n/\log n)\quad\text{pointwise},
\tag{20.2}
\]
where \(z_p=|Z_p|\) and \(K_n=|T_n|\).  This is useful bookkeeping, but
not a new reduction of the hard pointwise assertion.

Indeed,
\[
 \sum_{N<n\le2N}K_n
 \le\sum_{p\le2N}z_p.
\tag{20.3}
\]
Thus (20.1) gives a shell first moment \(O(N/\log N)\), average
\(O(1/\log N)\), and
\[
 \#\{N<n\le2N:K_n\ge\varepsilon N/\log N\}
 =O(1/\varepsilon).
\tag{20.4}
\]
The missing step that this finite exceptional set is empty is exactly
(20.2), hence exactly the top-half pointwise theorem.

There is a sharp abstract obstruction.  Choose rapidly increasing \(N_j\)
whose intervals \((N_j/2,N_j]\) are disjoint.  For every prime \(p\) in
that interval, put \(N_j-p\) and its reflected mate into \(Z_p\), with no
other optional zeros.  Then \(z_p\le2\), so (20.1) holds even pointwise,
but
\[
 K_{N_j}=\pi(N_j)-\pi(N_j/2)
 \sim\frac{N_j}{2\log N_j}.
\tag{20.5}
\]
The spikes have disjoint prime support.  Therefore a statement about
shared primes or spacing between two spikes cannot exclude one sparse
spike in every dyadic shell.

The proposed Hasse--Witt Mellin formula also requires calibration.  If
\[
 A_p(t)=\sum_{0\le j<p}b_jt^j,
\]
then for \(1\le r\le p-2\),
\[
 b_r\equiv-\sum_{t\in\mathbb F_p^\times}A_p(t)t^{-r}\pmod p
\tag{20.6}
\]
is the elementary power-sum coefficient-extraction identity.  Identifying
\(A_p\) with a Hasse--Witt polynomial supplies geometric meaning, but not
an estimate for how many coefficient-extraction sums vanish modulo their
varying characteristic.  To prove (20.1) from (20.6) one still needs a
precise integral or \(p\)-adic trace object, an average-over-\(p\) moment
theorem, and a valid passage from that moment to divisibility by \(p\).
Ordinary \(\ell\)-adic equidistribution of fiber traces does not supply
these steps.

Finally, a \(p^2\) congruence for a full transfer matrix is not, merely by
being matrix-valued, the two-characteristic relation in (18.66).  The
determinant is prescribed, the columns obey the same order-two recurrence,
and the all-grade CFVZ calculation (19.1) leaves one cyclic target
coordinate.  Higher precision within one characteristic must first
produce an explicit second target-selective invariant; even then it must
be batched across primes at sublinear height.  No such invariant is
currently known.  Similarly, geometric \(\mathrm{SL}_2\) monodromy of the
Picard--Fuchs family does not imply Haar equidistribution for the
variable-length coefficient-index walk
\(\Phi_r=M_r\cdots M_1\).

The honest value of (20.1) is therefore as a potentially strong average
theorem.  It would reduce possible counterexamples to \(O_\varepsilon(1)\)
per shell, but it does not replace either the fixed-\(q\) pair-energy
target or the pointwise two-characteristic dispersion target.

The computation through all \(17982\) primes \(p\le200000\) gives
\[
 \frac1{\pi(200000)}\sum_{p\le200000}z_p=1.00801,
 \qquad\max z_p=12.
\tag{20.7}
\]
The dyadic means show no visible upward drift.  This is unusually clean
evidence for (20.1), but the proved continuant bound gives only
\[
 \sum_{p\le x}z_p\ll\frac{x^{5/3}}{\log x}.
\tag{20.8}
\]
Indeed, if \(\rho_h(p)\) denotes the number of roots modulo \(p\) of the
gap continuant \(N_h\), the packing proof gives
\[
 z_p\le1+\frac pH+\sum_{2\le h\le H}\rho_h(p).
\tag{20.9}
\]
The known degree estimate sums to \(O(H^2)\), optimized at
\(H=p^{1/3}\).  A uniform growing-family estimate with \(O(1)\) average
roots per \(h\le p^{1/2}\) improves (20.8) only to the \(x^{3/2+o(1)}/
\log x\) scale.  Reaching (20.1) requires a global constraint on an
entire zero fibre, not separate Chebotarev estimates for each gap.

## 21. The modular central zero: exact structure and exact logical scope

Let \(n=p+r\) with \(0\le r<p\).  The fixed recurrence factor satisfies
\[
 p\mid 2n+1
 \quad\Longleftrightarrow\quad
 p\mid 2r+1
 \quad\Longleftrightarrow\quad
 r=\frac{p-1}{2}.
\tag{21.1}
\]
In the last case \(n=(3p-1)/2\).  Beukers' congruence, subsequently placed
in the hypergeometric modular framework of Ahlgren--Ono, is
\[
 b_{(p-1)/2}\equiv a_p(f)\pmod p,
 \qquad
 f(z)=\eta(2z)^4\eta(4z)^4.
\tag{21.2}
\]
Consequently the targets caught by the factor \(2n+1\) are exactly the
central nonordinary targets.  The two scan rows
\[
 (n,p,r)=(16,11,5),\qquad(4705,3137,1568)
\tag{21.3}
\]
match the two primes \(11,3137\) for which the modular coefficient was
found divisible by \(p\) in the stated finite range.  This is an exact
explanation of the coincidence, not a proof that these are the only such
primes.

The same two inputs give the exact parity identity
\[
 z_p\equiv{\bf1}_{p\mid a_p(f)}\pmod2.
\tag{21.4}
\]
This is unconditional for the primes in the congruence range; the finite
scan merely checks it.

The consequence drawn in the handoff must nevertheless be weakened.
Writing
\[
 z_p=|Z_p|
 =2\,\#\{\text{noncentral reflected zero pairs}\}
  +{\bf1}_{p\mid a_p(f)},
\tag{21.5}
\]
a bound \(z_p\le C\) is compatible with an arbitrary set of nonordinary
primes: it simply allows the last summand to be zero or one.  Similarly,
\[
 \sum_{p\le x}{\bf1}_{p\mid a_p(f)}\le\pi(x)
\tag{21.6}
\]
is already enough for its contribution to (20.1).  Neither a
classification nor a density theorem for nonordinary primes is needed for
these upper bounds.

Thus (21.2) rules out any argument asserting eventual central
nonvanishing, any claim that the factor \(2n+1\) covers all targets, and
any exact classification of empty zero fibres which ignores the modular
event.  It does not rule out a constant bound for the noncentral pairs,
and it does not make an average route logically unique.  In particular,
the statement that Hypothesis \(z_p=O(1)\) ``contains the nonordinary-prime
problem as a special case'' confuses an upper bound with an exact
classification.

References for (21.2) are F. Beukers, *Another congruence for the Apéry
numbers*, J. Number Theory 25 (1987), 201--210, and S. Ahlgren--K. Ono,
*A Gaussian hypergeometric series evaluation and Apéry number
congruences*, J. reine angew. Math. 518 (2000), 187--212.

## 22. Latest uisai2 file-channel results

Only the `dm` and `family` session files were synchronized for this
update.

### 22.1 Two Laurent models remain selectively rank one

The `dm` session compares two primitive Laurent periods
\(\Lambda_1,\Lambda_2\), both with
\(\operatorname {CT}\Lambda_i^m=b_m\).  Their natural fixed-\(n\)
prefixes \(S_n,D_n\) satisfy, for every top candidate,
\[
 S_n\equiv D_n\equiv b_{n-p}\pmod p,
 \qquad b_n\equiv5S_n\equiv5D_n\pmod p,
\tag{22.1}
\]
and
\[
 U_n\mid D_n-S_n,\qquad
 {\rm SNF}(S_n,D_n)={\rm diag}(1,U_n)
\tag{22.2}
\]
in the indicated cross-model basis.  The primitive quotient has
\[
 \log\frac{D_n-S_n}{U_n}
   =(\log16-\tfrac12+o(1))n.
\tag{22.3}
\]
The explicit birational mutation commutator adds a nonzero row but its
last primitive invariant is a unit at the target
\((20,17,3)\).  Hence the second model produces full-candidate content,
not a second target coordinate.

The same session corrected the selected Kummer-trace lift.  After the
necessary normalization, its second digit is related to the first Witt
coordinates by an invertible coordinate change (determinant \(-10\) for
\(p>5\)).  The target rows at \(p=17,r=3\) and the central
\(p=11,r=5\) have nonzero divided defects.  Thus the first divided marked
trace remains a free coordinate rather than a target zero.

### 22.2 A universal third layer on one quarter of the top band

The `family` session gives a genuine new congruence.  Put
\[
 T_n=\sum_{0\le k\le\lfloor(n-1)/2\rfloor}
       \binom nk^2\binom{n+k}k^2.
\tag{22.4}
\]
For every prime \(p\ge7\) in the top half,
\[
 b_n\equiv5T_n\pmod {p^2},
\tag{22.5}
\]
and, in the lower subband,
\[
 b_n\equiv5T_n\pmod {p^3}
 \qquad
 \left(\frac n2<p\le\frac{3n+1}{4}\right).
\tag{22.6}
\]
If \(U_n\) is the full top candidate product and \(V_n\) the subband
product, then
\[
 U_n^2V_n\mid b_n-5T_n.
\tag{22.7}
\]
The proof uses the complete low block and the Apéry WZ certificate; the
cutoff beyond both Kummer thresholds contributes \(p^4\).  This harvests
one extra quarter-unit of universal prime-number-theorem content.

It remains nonselective.  The residual height bound is
\[
 \log\left|\frac{b_n-5T_n}{U_n^2V_n}\right|
 \le\bigl(\log(17+12\sqrt2)-\tfrac54+o(1)\bigr)n,
\tag{22.8}
\]
and the actual central target \((16,11,5)\) has a nonzero fourth divided
digit.  Thus (22.6) is positive filtration progress but not a proof of
the target radical bound.

## 23. Full transfer frames and the remaining local scalar

For
\[
 M_k=
 \begin{pmatrix}
 c_k/(k+1)^3&-k^3/(k+1)^3\\
 1&0
 \end{pmatrix},
 \qquad
 c_k=(2k+1)(17k^2+17k+5),
\tag{23.1}
\]
one has \(\det M_k=k^3/(k+1)^3\).  The integral block
\[
 L_k=(k+1)^3M_k
\tag{23.2}
\]
at the singular crossing is
\[
 L_{p-1}\equiv
 \begin{pmatrix}-5&1\\0&0\end{pmatrix}\pmod p.
\tag{23.3}
\]
It therefore has rank one.  This is the matrix form of the scalar Lucas
collapse.

More precisely,
\[
 \Theta_r=\Phi_{p+r}\Phi_p^{-1}
 \quad\Longrightarrow\quad
 \det\Theta_r=\frac{p^3}{(p+r)^3}.
\tag{23.4}
\]
At a target, the first column is primitive because consecutive Apéry
zeros are excluded.  The local Smith valuations are therefore exactly
\((0,3)\).  It follows that the frame is cyclic modulo \(p^s\) for
\(s\le3\).

On \(p\mid b_r\), write
\[
 \xi_r=\frac{b_r}{p}+D_r\pmod p.
\]
The verified scalar lift gives
\[
 \frac{b_{p+r}}p\equiv5\xi_r,\qquad
 \frac{p^3a_{p+r}}p\equiv6\xi_r\pmod p.
\tag{23.5}
\]
Every subsequent step before the next singular row is invertible modulo
\(p\), so neighboring first-layer coordinates are affine functions of
the same \(\xi_r\).  The Wronskian is the determinant of the same transfer
product, not another independent row.

There is also a precise limitation on this no-go.  Right-multiplying a
local fundamental frame by \(I+ptE_{21}\) preserves its recurrence,
determinant, and reduction modulo \(p\), while varying a divided
coordinate.  This proves that the listed **local** identities cannot fix
\(\xi_r\).  The deformation changes the distinguished initial frame, so
it does not prove that the actual Apéry connection data can never supply
a new equation.  In fact the second Smith direction first becomes visible
modulo \(p^4\), where it is represented by the companion column
\[
 -p^3\begin{pmatrix}a_r/6\\a_{r-1}/6\end{pmatrix}\pmod {p^4}.
\tag{23.6}
\]
Thus the determinant argument closes precisions at most three, not every
finite lift.  The new fourth-grade direction is one-characteristic
companion data; it still needs both a target-selective implication and a
cross-characteristic batching mechanism to address (18.66).

## 24. Reduced companion numerators and the double-Casoratian test

Write
\[
 a_n=\frac{A_n}{C_n},\qquad
 \gcd(A_n,C_n)=1,\quad C_n>0.
\tag{24.1}
\]
The identity
\[
 G_n=\frac{d_n^3}{C_n}\gcd(A_n,b_n)
\tag{24.2}
\]
is immediate.  Exact recurrence computation found
\[
 \gcd(A_n,b_n)=1\qquad(0\le n\le10000).
\tag{24.3}
\]
If (24.3) held for all \(n\), P3.2 would become the lower-denominator
statement
\[
 \log(d_n^3/C_n)=o(n).
\tag{24.4}
\]
Neither assertion is presently proved.

There is an exact denominator identity which shows that this reformulation
does not simplify the hard support.  Put \(D_n=d_n^3\) and
\(P_n=D_na_n\).  Since \(C_n\mid D_n\),
\[
 \frac{D_n}{C_n}=\gcd(D_n,P_n),
\qquad
 v_p(C_n)=3\lfloor\log_pn\rfloor
  -\min\{3\lfloor\log_pn\rfloor,v_p(P_n)\}.
\tag{24.5}
\]
For \(p>\sqrt n\), write \(n=qp+r\).  The block law gives
\[
 p\mid D_n/C_n
 \quad\Longleftrightarrow\quad
 p\mid a_qb_r.
\tag{24.6}
\]
In particular, for \(n/2<p\le n\) and \(p\ge7\),
\[
 p\mid D_n/C_n
 \quad\Longleftrightarrow\quad
 p\mid b_{n-p}
 \quad\Longleftrightarrow\quad
 p\mid b_n.
\tag{24.7}
\]
Hence (24.4) already contains the top-half radical theorem.  Proving
coprimality in (24.3) would identify \(G_n\) with this denominator defect,
but would not prove that the defect is subexponential.

The adjacent Casoratians give a tempting necessary condition.  Clearing
the reduced denominators in
\[
 a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}
\]
shows that \(g_n=\gcd(A_n,b_n)\) divides both
\[
 \Delta_n=\frac{6C_{n-1}C_n}{n^3},
 \qquad
 \Delta_{n+1}=\frac{6C_nC_{n+1}}{(n+1)^3}.
\tag{24.8}
\]
Hence
\[
 g_n\mid\gcd(b_n,\Delta_n,\Delta_{n+1}).
\tag{24.9}
\]
The right side is far from one.  An exact scan for \(1\le n<2000\)
found 1952 nontrivial values; at \(n=20\) it retains the prime \(17\).
Thus two neighboring determinants do not explain (24.3).  A proof would
have to use how the distinguished companion numerator selects a primitive
line inside this much larger denominator-supported carrier, rather than
another rank-two determinant alone.

## 25. Exact companion-denominator curvature

The reduced-denominator experiment has a proof-grade strengthening.  It is
important because it both validates the isolated-hole picture and prevents
the cubic-divisibility observation from being mistaken for a solution.

Let \(p\ge7\), \(n=p+r\), \(0<r<p\), and suppose \(p\mid b_n\).  Lucas,
reflection, and the endpoint values imply
\[
 2\le r\le p-3,\qquad
 p\nmid b_{n-1}b_{n+1}.
\tag{25.1}
\]
At the singular row,
\[
 p^3a_p\equiv6,\qquad p^3a_{p+1}\equiv30\pmod p.
\tag{25.2}
\]
The first congruence follows from the Casoratian at \(p\), using
\(b_{p-1}\equiv1\pmod {p^2}\); the second follows from the recurrence
row \(p\).  Since \(p\ge7\), both displayed residues are units.  Every
transfer row from \(p+1\) through \(n-1\) lies in
\({\rm GL}_2(\mathbb Z_p)\), so
\[
 \min\{v_p(a_{n-1}),v_p(a_n)\}=-3.
\tag{25.3}
\]

Put
\[
 t=v_p(b_n)\ge1,\qquad
 \alpha=v_p(a_n),\qquad\beta=v_p(a_{n-1}).
\]
If \(\beta>-3\), (25.3) gives \(\alpha=-3\), and the two terms in
\[
 a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}
\tag{25.4}
\]
would have valuations \(-3\) and \(>-3\), contradicting the unit
right-hand side.  Hence \(\beta=-3\).  Comparing the two term valuations
\(\alpha\) and \(t-3\) in (25.4) gives
\[
 \boxed{
 \alpha=
 \begin{cases}
 t-3,&t=1,2,\\
 \ge0,&t=3,\\
 0,&t\ge4.
 \end{cases}}
\tag{25.5}
\]
One further invertible row gives \(v_p(a_{n+1})=-3\).  Since a reduced
denominator records the negative part of the valuation, this proves
\[
 \boxed{
 (v_p(C_{n-1}),v_p(C_n),v_p(C_{n+1}))
 =(3,\,3-\min(t,3),\,3).}
\tag{25.6}
\]

Define the companion-denominator curvature
\[
 \kappa_n=
 \frac{\gcd(C_{n-1},C_{n+1})}
      {\gcd(C_{n-1},C_n,C_{n+1})}.
\tag{25.7}
\]
Then
\[
 \boxed{v_p(\kappa_n)=\min\{v_p(b_n),3\}}
 \qquad(n/2<p\le n,\ p\ge7).
\tag{25.8}
\]
If \(p\nmid b_n\), the block congruence
\(p^3a_{p+r}\equiv6b_r\pmod p\) makes the middle denominator exponent
three, while the neighboring exponents are at most three, so the
curvature is zero.  Thus (25.8) also has the exact converse at the level
of support.

There are two immediate corollaries.

First, if \(p\mid A_n\) as well as \(p\mid b_n\), then
\(\alpha\ge1\).  Formula (25.5) forces
\[
 \boxed{v_p(b_n)=3.}
\tag{25.9}
\]
This is an equality, not only a lower bound.  Conversely a cubic target
need not divide \(A_n\), so (25.9) does not characterize the intersection.
An exact scan through \(p\le50000\) found 5170 shifted targets:
5165 have valuation one and the five
\[
 (p,r)=(41,30),(97,25),(151,14),(1453,1180),(6781,3974)
\tag{25.10}
\]
have valuation two; none has valuation at least three.  This is evidence,
not an exclusion theorem.

Second, if \(H_n=\gcd(A_n,b_n)\) and \(p^e\mid H_n\), clearing (25.4)
gives the stronger neighboring bounds
\[
 v_p(C_{n-1})\ge e+3v_p(n),\qquad
 v_p(C_{n+1})\ge e+3v_p(n+1).
\tag{25.11}
\]
There is no conflict with (24.8): the integer
\[
 A_nC_{n-1}b_{n-1}-A_{n-1}C_nb_n
 =\frac{6C_nC_{n-1}}{n^3}
\tag{25.12}
\]
is divisible by \(H_n\).  In particular, the division by \(n^3\) in
(24.8) is valid; an intermediate audit which objected to it had omitted
the integrality and was withdrawn.

Delaygue's Theorem 1 does not sharpen (25.5).  It proves only
\[
 v_p(b_m)\ge\alpha_p(b,m),
\tag{25.13}
\]
where \(\alpha_p\) counts base-\(p\) digits lying in \(Z_p\).  For
\(m=p+r\), \(p\ge7\), the digits are \((r,1)\) and \(b_1=5\), so a target
gives only \(v_p(b_{p+r})\ge1\), exactly the Lucas lower bound.  The
Landau--Legendre equality in that paper applies to the individual
factorial ratios before summation, not to the Apéry multisum.

The conceptual conclusion is exact: \(\kappa_n\) is a characteristic-zero
carrier for the complete top-half target set, but (25.8) makes a height
bound for it equivalent to the missing denominator-defect estimate.

The `uisai2` `dm` session proved that this curvature is an exact carrier at
every index, not only in the top half.  Put
\[
 T_n=\frac{C_{n-1}}{\gcd(C_{n-1},C_n)},\qquad
 U_n=\frac{C_{n+1}}{\gcd(C_n,C_{n+1})}.
\tag{25.14}
\]
Primewise, \(\kappa_n=\gcd(T_n,U_n)\).  Clearing the Casoratian and reducing
modulo the coprime denominator drop gives
\[
 T_n\mid n^3b_n,\qquad U_n\mid(n+1)^3b_n.
\tag{25.15}
\]
Since consecutive integers are coprime,
\[
 \boxed{\kappa_n\mid b_n}
\tag{25.16}
\]
for every \(n\).

There is also a reverse comparison outside a sublinear exceptional factor.
Let
\[
 R_n=\operatorname {rad}\gcd(b_n,d_n).
\tag{25.17}
\]
Choose \(y=\lceil n^{2/3}\rceil\), \(A=\lfloor n/y\rfloor\), and put
\[
 {\rm Err}_n=
 \operatorname {rad}(d_y n(n+1))
 \prod_{1\le q\le A}
   \operatorname {rad}\!\left(b_q\,d_q^3a_q\right).
\tag{25.18}
\]
The factors in (25.18) are integers and
\[
 \log{\rm Err}_n=O(n^{2/3}).
\tag{25.19}
\]
If a prime \(p\mid R_n\) does not divide \({\rm Err}_n\), then
\(p>y\), \(n=qp+r\) lies away from both block boundaries, and both upper
coordinates \(b_q,a_q\) are \(p\)-units.  Lucas and the companion block
law give
\[
 p\mid b_r,\qquad
 p^3a_n\equiv0,\qquad
 p^3a_{n-1},p^3a_{n+1}\not\equiv0\pmod p.
\tag{25.20}
\]
Thus the denominator exponents at \(n-1,n,n+1\) have a strict valley and
\(p\mid\kappa_n\).  Together with (25.16),
\[
 \operatorname {rad}(\kappa_n)\mid R_n
 \mid\operatorname {rad}(\kappa_n)\,{\rm Err}_n.
\tag{25.21}
\]
The denominator theorem also bounds the small-prime multiplicities, while
every exponent above \(\sqrt n\) is at most three.  Consequently
\[
 \boxed{
 \log R_n=o(n)\quad\Longleftrightarrow\quad
 \log\kappa_n=o(n).}
\tag{25.22}
\]
Combined with the already proved small-prime and companion-channel
reductions, (25.22) is equivalent to P3.2.  This is the cleanest
target-free statement presently available:
\[
 \log\gcd\!\left(
 \frac{C_{n-1}}{\gcd(C_{n-1},C_n)},
 \frac{C_{n+1}}{\gcd(C_n,C_{n+1})}
 \right)=o(n).
\tag{25.23}
\]
It remains a pointwise strict-valley anti-concentration theorem; ordinary
telescoping of denominator rises and drops does not control the overlap in
(25.23).

## 26. The single common coefficient carrier

The moving-remainder formulation can be replaced, up to a proved
sublinear error, by a large-prime-factor question for the one integer
\(b_n\).  Put
\[
 B(n)=\{p>\sqrt n:p\mid b_{n\bmod p}\},\qquad
 C(n)=\{p>\sqrt n:p\mid b_n\}.
\tag{26.1}
\]
If \(n=qp+r\), Lucas gives
\[
 b_n\equiv b_qb_r\pmod p,
\tag{26.2}
\]
and hence \(B(n)\subset C(n)\).  A prime in \(C(n)\setminus B(n)\)
divides \(b_q\).

Fix \(1\le Q\le\sqrt n\).  For the false positives with \(q\le Q\),
distinct primes are all factors of \(\prod_{q\le Q}b_q\), so
\[
 \sum_{\substack{p\in C(n)\setminus B(n)\\q\le Q}}\log p
 \le\sum_{q\le Q}\log b_q=O(Q^2).
\tag{26.3}
\]
For \(q>Q\), one has \(p<n/Q\), and therefore
\[
 \sum_{\substack{p\in C(n)\setminus B(n)\\q>Q}}\log p
 \le\vartheta(n/Q)=O(n/Q).
\tag{26.4}
\]
Choosing \(Q=n^{1/3}\) proves
\[
 \boxed{
 \sum_{p\in C(n)}\log p
 =\sum_{p\in B(n)}\log p+O(n^{2/3}).}
\tag{26.5}
\]

Thus the unresolved channel is equivalent to
\[
 \boxed{
 \log\operatorname {rad}_{\sqrt n<p\le n,\ p\mid b_n}b_n=o(n).}
\tag{26.6}
\]
This is more than the one-way observation that channel-\(B\) primes divide
\(b_n\): (26.5) proves that the extra large prime factors of \(b_n\) have
sublinear total weight.

The reformulation does not supply its own bound.  The exponential height
\[
 \log b_n=n\log(17+12\sqrt2)+O(\log n)
\tag{26.7}
\]
allows a linear large-prime radical.  Theorems of Luca--Shparlinski give
density-one lower bounds on prime factors; Delaygue gives lower
\(p\)-adic valuations; standard \(G\)-function denominator results say
nothing about numerator prime factors.  No checked theorem gives the
all-\(n\) upper bound (26.6).

For comparison, the ambient class cannot suffice: the central binomial
coefficient is hypergeometric, holonomic, globally bounded, and a rational
diagonal, yet every prime in \((n,2n]\) divides \(\binom{2n}{n}\).  Any
successful theorem must use a feature specific to the Apéry initial state,
not merely membership in one of these classes.

## 27. Mellin normalization, exact-zero complexity, and density one

Let
\[
 \Phi(x,y,z)=(1+x)(1+y)(1+z)(1+1/(xyz)).
\tag{27.1}
\]
Then \(b_m=\operatorname {CT}\Phi^m\).  For the toric pencil
\(1-u\Phi=0\), the first Hasse--Witt scalar in the standard residue
trivialization is
\[
 F_{<p}(u)=\sum_{0\le m<p}b_mu^m,
\tag{27.2}
\]
because
\[
 \operatorname {CT}(1-u\Phi)^{p-1}
 \equiv\sum_{m<p}b_mu^m\pmod p.
\tag{27.3}
\]
If the family is written as \(\Phi=t\), the period parameter is
\(u=t^{-1}\); a change of Hodge-line trivialization may add a monomial
unit.  Consequently the literal identity \(F_{<p}(t)=\alpha_p(t)\) in
the \(\Phi=t\) coordinate is not correct.

Finite Fourier inversion gives, for \(1\le r\le p-2\),
\[
 b_r=-\sum_{t\in\mathbb F_p^\times}
       F_{<p}(t)t^{-r}\pmod p.
\tag{27.4}
\]
At the boundary the sum is
\(-b_0-b_{p-1}=-2\), so no zero mode was lost.

The obstacle to the proposed trace-formula proof is exact.  Expanding
\({\bf1}_{S=0}\) additively gives
\[
 \frac1p\sum_{\lambda\in\mathbb F_p}
  e_p\!\left(\lambda\sum_t
    F_{<p}(t)t^{-r}\right)
 =
 \frac1p\sum_\lambda\prod_t
  e_p\!\left(\lambda F_{<p}(t)t^{-r}\right).
\tag{27.5}
\]
The product involves all \(p-1\) fibers.  Equivalently,
\({\bf1}_{S=0}=1-S^{p-1}\) uses a moment whose order grows with \(p\).
Thus ordinary fixed-order Katz--Deligne correlations do not prove the
bounded first moment for \(z_p\).  This is a no-go for that expansion,
not for every possible use of the geometry.

There is nevertheless an unconditional density-one theorem from the
proved vertical estimate.  Let \(W_B(n)\) be the logarithmic weight of
the channel-\(B\) primes.  On \(N<n\le2N\),
\[
\begin{aligned}
 \sum_{N<n\le2N}W_B(n)
 &\le
 \sum_{\sqrt N<p\le2N}
   \left(\frac Np+1\right)|Z_p|\log p\\
 &=O(N^{5/3}),
\end{aligned}
\tag{27.6}
\]
using \(|Z_p|\ll p^{2/3}\).  Hence for every
\(2/3<\alpha<1\),
\[
 \#\{N<n\le2N:W_B(n)>N^\alpha\}
 =O(N^{5/3-\alpha})=o(N).
\tag{27.7}
\]
The small-prime and companion channels are already \(o(n)\) pointwise.
It follows that
\[
 \boxed{\log G_n=o(n)}
\tag{27.8}
\]
outside a set of natural density zero.  The all-\(n\) problem is precisely
the removal of the sparse exceptional columns left by (27.7).

## 28. Exact obstruction for plain binomial windows

For \(p>n/2\), Kummer's carry criterion simplifies to
\[
 p\mid\binom nk
 \quad\Longleftrightarrow\quad
 p>\max(k,n-k).
\tag{28.1}
\]
Indeed \(n=p+r\), \(0\le r<p\), and the unique base-\(p\) carry occurs
exactly when \(r<k<p\).

Let
\[
 {\cal B}(n,H)=
 \gcd_{\lvert k-n/2\rvert\le H}\binom nk.
\tag{28.2}
\]
Every prime \(p>n/2+H\) divides every coefficient in the window, whence
\[
 \log{\cal B}(n,H)
 \ge\vartheta(n)-\vartheta(n/2+H)
 =n/2-H+o(n)
\tag{28.3}
\]
whenever the prime-number-theorem error is uniform in the chosen range.
Therefore a subexponential window gcd requires
\(H=n/2-o(n)\), but then it captures only
\(p>n-o(n)\).  At the other endpoint, the single central coefficient
captures all of \((n/2,n]\) and has height \(n\log2+O(\log n)\).

This proves a rigid coverage-height tradeoff for the whole plain-window
family.  It does not exclude signed factorial ratios or an
Apéry-dependent multiscale carrier, but any such construction must leave
this one-parameter family.
