# Codex Task 002: CED Verification + Lemma 1 Proof

## STRATEGIC BREAKTHROUGH (from Fable R3)

The X-factor in (HM)_3 CRT error dissolves structurally at k=3:
- Fix output prime p, set v = qℓ (E₂-number in (X², 4X²])
- After reciprocity shuffle, Cauchy-Schwarz in v produces complete sums mod pp'
- Weil bound on the completed sum harvests the missing factor X = √(pp')

## TASK A: Verify the trilinear identity

For distinct primes p, q, ℓ ∈ (X, 2X], n = pqℓ:

#{m < X² : m mod p ∈ Z_p, m mod q ∈ Z_q, m mod ℓ ∈ Z_ℓ}
= (X²/n)·Z(p)Z(q)Z(ℓ) + E_{p,q,ℓ}

where E is the Fourier error involving F_p(a) = Σ_{r∈Z_p} e_p(ar).

Write a Python/Sage script to:
1. Compute E_{p,q,ℓ} exactly for all triples of primes in (X, 2X] with X = 50
2. Verify the Fourier expansion formula
3. Verify the exact orthogonality: Σ_{v mod pp'} F_p(kv⁻¹) F̄_{p'}(k'v⁻¹) = Z(p)Z(p') for p ≠ p'

## TASK B: Prove Lemma 1 (Fourier non-concentration)

Lemma: Σ_{1≤|k|≤K₀} |F_p(k)|² ≤ C(K₀·Z(p) + p²/K₀)

Proof route:
1. Majorize the window indicator by Fejér kernel of width p/K₀
2. The resulting expression: 2K₀[Z(p) + Σ_{0<h≤p/K₀} A_p(h)(1-hK₀/p)]
3. Insert A_p(h) ≤ 3(h-1) (from gap polynomial degree + nonvanishing)
4. Sum: Σ 3h over h ≤ p/K₀ gives p²/K₀²

Write this up as a formal proof in problems/3.2/ORACLE_COMM/lemma1_proof.tex

## TASK C: Numerical measurement

For primes p up to 10000:
1. Compute Σ_{k≠k', |k|,|k'|≤K₀} M_p(k,k') where M_p(k,k') = #{(r,r')∈Z_p²: kr ≡ k'r' (mod p)}
2. Compare with random prediction K₀²Z²/p + K₀Z
3. Test at K₀ = √p and K₀ = p^{1/3}

Output: problems/3.2/ORACLE_COMM/codex_result_002.md + scripts
