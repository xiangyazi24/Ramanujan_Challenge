# Codex Task 005: Write new paper sections

Write LaTeX for the following new results to be added to proof.tex.
Output to: problems/3.2/new_sections.tex

## Section: The diagonal Cartier representation (Theorem A)

State and prove:
For every prime p ≥ 5 and 0 ≤ r ≤ p-1:
  b_r ≡ [x^r y^r z^r w^r] F(x,y,z,w)^{p-1} (mod p)
where F = (1-x-y)(1-z-w) - xyzw.

Proof: Straub 2014 diagonal + Cartier operator (5 lines).
Cite: Straub, "Multivariate Apéry numbers and supercongruences", 
Algebra Number Theory 8(8):1985-2008, 2014.

## Section: Fourier non-concentration of zero sets (Lemma 1)

State and prove:
For every prime p ≥ 7 and 1 ≤ K₀ ≤ p:
  Σ_{1≤|k|≤K₀} |F_p(k)|² ≤ C(K₀ Z(p) + p²/K₀)

Proof: Fejér kernel majorization + A_p(h) ≤ 3(h-1) from gap polynomial.

## Section: The moment collapse

State and prove the one-line reduction:
  Σ (K_X(m))_k ≤ (max K)^{k-2} · Σ (K_X(m))_2 ≤ 5X²λ² · (max K)^{k-2}

Corollary: max K_X(m) ≪ λ_X · X^{o(1)} implies (HM)_k for all k,
hence G_n = e^{o(n)} for all n.

State (AT″) as a hypothesis and show it implies the full conjecture.

## Section: The honest reduction of (HM)_3

State as a theorem:
(HM)_3 ⟸ (MC) + (AT), where:
- (MC): multiplicative correlation bound on Z_p
- (AT): max K_X(m) ≪ X^{2/3+o(1)} λ_X

Explain the CED architecture (sketch, citing Linnik dispersion method).
Identify the anchored star as the sharp countermodel.

## Section: Computational evidence for (AT″)

Table: X, max K, X^{2/3}λ, Xλ, ratio for X = 16,...,4096.
Poisson extreme-value prediction.

Write professional mathematical prose suitable for a research paper.
