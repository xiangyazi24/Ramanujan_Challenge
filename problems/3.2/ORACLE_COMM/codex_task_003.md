# Codex Task 003: Multiplicative Correlation M_p(k,k') Numerics (W-D)

## Definition
M_p(k,k') = #{(r,r') ∈ Z_p × Z_p : kr ≡ k'r' (mod p)}

For Z_p = zero set of Apéry numbers mod p.

## Task
Write a script that computes, for all primes p up to 5000:

1. For K₀ = floor(√p), compute:
   Σ_{k≠k', 1≤k,k'≤K₀} M_p(k,k')
   
   Compare with random prediction: K₀²·Z(p)²/p + K₀·Z(p)

2. For K₀ = floor(p^{1/3}), same computation.

3. Report the ratio (actual / random prediction) for each p with Z(p) ≥ 2.

4. Also compute the "multiplicative energy":
   E_mult(p) = Σ_{k=1}^{p-1} M_p(k,1)²
   Compare with p·Z(p)² (random prediction).

## Output
Write results + analysis to problems/3.2/ORACLE_COMM/codex_result_003.md
Write the script to problems/3.2/ORACLE_COMM/mp_numerics.py

This is time-sensitive — the campaign hinges on whether M_p is random-like.
