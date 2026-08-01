# Ramanujan Challenge — Problem 3.2

**Optimality of Apéry's irrationality-measure bound for ζ(3)**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

---

## Honest statement of what is proved

The challenge asks for a proof that

```
G_n = gcd(d_n a_n, d_n b_n) = e^{o(n)},      d_n = lcm(1,…,n)^3.
```

**This package does not contain a proof of that statement for all n.** It contains:

* an unconditional theorem giving `G_n = e^{o(n)}` for a set of n of **natural density 1**,
  with a polylogarithmic exceptional set;
* several new unconditional congruences and structural theorems (listed below), one of which
  makes the top-window step of the reduction exact rather than empirical;
* a precise identification of the single remaining obstruction, with the reasons every
  standard technique fails;
* large-scale computations confirming the conjecture in the range tested.

Nothing in this package is asserted beyond what is proved. Every numerical claim has a
runnable script in `scripts/`.

## New unconditional results in this submission

**1. Vector Lucas property (new; the numerator half appears to be unrecorded).**
For a prime p ≥ 7 and n = qp + r with 1 ≤ q < p, 0 ≤ r < p,

```
( p^3 a_{qp+r} ,  b_{qp+r} )  ≡  b_r · ( a_q , b_q )   (mod p).
```

The second coordinate is the classical Apéry–Lucas congruence `b_{qp+r} ≡ b_q b_r`
(Gessel 1982; Malik–Straub 2016). The first coordinate is the numerator counterpart; the
factor p^3 appears because the leading coefficient (n+1)^3 of the recurrence vanishes at
n = p−1. Both cases are proved in the paper. For q = 1 there is a direct five-step argument; for
general q the proof splits into a boundary computation `p^3 a_{qp} ≡ a_q (mod p)` (Kummer:
v_p(C(qp,k)) = 0 exactly when p | k, and the correction terms with p ∤ m have denominator
valuation ≤ 2, so only m = tp survives) and a propagation step (modulo p the recurrence in
r is exactly the b-recurrence, with u_{qp+1} ≡ 5u_{qp} since P(0) = 5), giving
u_{qp+r} ≡ u_{qp}·b_r for both coordinates.
The ζ(2) Apéry pair satisfies the same law with p^2 (3198 checks); the r-propagation is formal
once p^κ A_n is p-integral, but the boundary identity is separate input for each family.
A literature check (Gessel 1982; Beukers 1985, 1987; Coster; Malik–Straub 2016; Delaygue 2018;
Rowland–Yassawi–Krattenthaler 2021; Straub 2024) found no numerator counterpart, so the claim
made here is "apparently new in this normalisation", not absolute novelty.

*Consequence.* If p ∤ b_n then p^3 a_n is a p-unit, so v_p(D_n) = 3 and p does not divide
G_n. (Strictly, v_p(G_n) = e_p(n) + min(v_p(A_n), v_p(b_n)); in the top window with
v_p(b_n) ≤ 2 the second term vanishes because v_p(a_n) < 0.) Hence **only primes dividing b_n contribute to G_n**, and the top-window step of the
reduction is unconditional. It also unifies the trichotomy of lower-digit / leading-digit /
companion-block primes: modulo p the pair at index n is b_r times the pair at index q (a
Cartier/Frobenius block law — the action is b_r·I, so rank two, not a rank-one collapse), and
the right-hand side vanishes in exactly three ways.

**2. The marked scalar is a moment of point counts.** For every prime p ≥ 5 and every M,

```
C_M(p−1) ≡ − Σ_{x,y,z ∈ F_p^*} Λ(x,y,z)^M ≡ − Σ_t t^r N_p(t)   (mod p),
```

with Λ the Apéry Laurent polynomial and N_p its fibre point count. So the target condition is
the vanishing of a moment of the point-count function of the Apéry family.

**3. Parity of the zero count.** `|Z_p|` is **odd if and only if p is non-ordinary** for the
weight-4 level-8 newform 8.4.a.a — a corollary of the palindromy `b_{p−1−r} ≡ b_r (mod p)`
(Malik–Straub, Lemma 6.2) and the Ahlgren–Ono supercongruence. Verified for all 2260 primes
5 ≤ p ≤ 20000: odd exactly at p = 11 and p = 3137, exactly the non-ordinary primes.

**4. Unconditional discharge below n/log n.** Since D_n | d_n, one has
e_p(n) ≤ 3⌊log_p n⌋, so `Σ_{p ≤ n/log n} e_p(n) log p = O(n/log n) = o(n)` by Chebyshev
alone. Combined with the O(log²n) high-digit contribution, the conjecture is *equivalent* to

```
Σ_{p ∈ (n/log n, n], (n mod p) ∈ Z_p} log p = o(n).
```

## The remaining obstruction, stated precisely

For p in the top window the bad condition is a pure congruence `n mod p ∈ Z_p`, so the problem
is: one residue class per prime, modulus ≍ range. This is the **occupancy-one** barrier, not
the level-of-distribution barrier that BFI/Fouvry–Iwaniec/Zhang cross. Equivalently, in
∏_p F_p the intersection ⋂ B_p = ∏ Z_p is automatically nonempty, and the question is whether
one of its CRT representatives lies in the short interval — a *short-interval CRT discrepancy*
question, on which covering-system technology is silent by design.

The paper documents why each standard route fails: rank-one is forced by the geometry
(Beukers–Vlasenko), the mod-p² digit is not target-selective, the large sieve meets the Q²
barrier, abc points the wrong way, and only the first moment is supported by the range.

## A concrete proof path: two-prime Weil correlation (new, this session)

**5. Vertical Weil bound.** The complete exponential sum `C_p(h) = Σ_{a=0}^{p-1} e(h b_a/p)`
satisfies `|C_p(1)| ≤ 3.27 √p` for all 166 primes p ≤ 1000, and average ratio 1.26.
This is a trace-function / Deligne-type phenomenon. Verified by `scripts/q32_vertical_weil_audit.py`.

**6. Two-prime shifted correlation (the decisive gate).** For distinct primes p, q and
shift d = |p−q|, the correlation `Σ_m e(b_m/p − b_{m+d}/q)` satisfies the Weil-scale bound
`|Corr| ≤ 2.09 √M` across **all 127 tested prime pairs** (p, q ≤ 709, 102 pairs in [200,600]).
Verified by `scripts/q32_two_prime_correlation_audit.py`.

**Proof chain (if the two-prime bound is proved):** two-prime Weil ⟹ 4th-moment bound
`Σ |S_h|^4 ≤ C N³/log²N` ⟹ pointwise `max |S_h| = o(P_n)` ⟹ Fejér ⟹ the conjecture.

**The precise missing theorem:** construct a bounded-conductor geometrically irreducible
ℓ-adic sheaf G_{p,q,d} on A¹ whose Frobenius trace is `ψ_p(b_m) ψ_q(−b_{m+d})`, then apply
Deligne's Riemann Hypothesis. The Apéry Picard-Fuchs/Gauss-Manin sheaf provides the natural
candidate; the gap is showing that the *product* sheaf (two different prime reductions of the
same global motive, with shift) has no geometrically trivial component.

## Why the conjecture is nevertheless true

* `K(n) = #{p ∈ (n/2,n] : p | b_n} ≤ 3` for **all n ≤ 1,000,000**.
* `log rad_{p≤n}(b_n)/n` ≤ 0.109, 0.025, 0.0041, 0.00059 on [10²,10³], [10³,10⁴], [10⁴,10⁵], [10⁵,10⁶];
  the data fit `O(log n log log n)` against a required o(n).
* ω(b_n) ∈ {3,…,7} for n ≤ 26 — b_n barely has prime factors, and they are huge.
* The conjecture **fails** for Cooper's sporadic sequences s7, s10, s18 and for C(2n,n),
  whose zero sets contain macroscopic intervals; so no proof can rest on D-finiteness,
  modularity, Lucas congruences or polytope reflexivity alone. The separating property is the
  size and geometry of Z_p.

## Reproducing

```bash
python3 scripts/q32_numerator_congruence_audit.py    # the new vector Lucas congruence
python3 scripts/q32_marked_scalar_character_sum.py   # moment identity, palindromy, parity
python3 scripts/q32_shift_operator_audit.py          # 27-term criterion, shift operator
python3 scripts/q32_actual_Gn_audit.py               # G_n = d_n/D_n, targets are its large primes
python3 scripts/q32_seam_ray_split_audit.py          # seam ray split
python3 scripts/q32_family_compare.py                # the family dichotomy
python3 scripts/q32_top_window_target_counts.py      # K(n) over a range
python3 scripts/q32_vertical_weil_audit.py           # vertical Weil bound
python3 scripts/q32_two_prime_correlation_audit.py   # two-prime Weil correlation
gcc -O3 -pthread -o big_scan scripts/q32_big_scan.c -lm && ./big_scan 200000 200000 8
```

All scripts use only the Python standard library (the family comparison additionally reports
Cooper's sequences; `q32_big_scan.c` needs a C compiler and pthreads).
