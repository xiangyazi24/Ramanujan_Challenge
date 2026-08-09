# Codex Task 004: Two-flip identity + Line-section statistics

## BACKGROUND
The CED proof sketch has two fatal flaws identified by oracle audit.
The honest reduction is: (HM)_3 ⟸ (MC) + (AT), both named.

## TASK A: Verify the two-flip reciprocity identity

For CRT lift m = m(r,s,t) of (r mod p, s mod q, t mod ℓ) with n = pqℓ:

e(km/n) = e_p( k(r−s)·(qℓ)⁻¹ ) · e_ℓ( k(t−s)·(pq)⁻¹ ) · e(ks/n)

Verify this symbolically for small examples (p,q,ℓ = 5,7,11 and various k,r,s,t).
Write a verification script.

## TASK B: Line-section statistics of the incidence cloud

The additive incidence cloud is:
  C_p = {(r, h) ∈ F_p² : r ∈ Z_p, r+h ∈ Z_p, 2 ≤ h ≤ p}
  
(equivalently: N_h(r) ≡ 0 mod p AND r ∈ Z_p)

For the CED argument, M_p(k,k') equals the number of points of C_p
on the line L_{k,k'}: k'h ≡ (k-k')r (mod p), summed over small c.

Compute for primes p up to 2000 with Z(p) ≥ 2:
1. The size of C_p
2. The distribution of |C_p ∩ L| over a pencil of lines through the origin
3. The second moment Σ_L |C_p ∩ L|² (the "line-section energy")
4. Compare with the random prediction for a cloud of that size

## TASK C: Formulate (AT) precisely and test numerically

(AT): max_{m < X²} K_X(m) ≪ X^{2/3+o(1)} λ_X

Currently proved: max K ≤ X λ (from (HM)_2 + Chebyshev)
Needed: X^{1/3} improvement

For dyadic blocks X = 2^j, j = 4,...,12:
1. Compute K_X(m) = #{p ∈ (X,2X] : m mod p ∈ Z_p} for all m < X²
2. Report max K_X(m) and compare with X^{2/3} λ_X and X λ_X
3. Is the actual max K much smaller than the trivial bound?

Output all results to problems/3.2/ORACLE_COMM/codex_result_004.md
