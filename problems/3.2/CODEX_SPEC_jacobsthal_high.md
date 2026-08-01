# SPEC: Jacobsthal parametrization — FAST RECON (gpt-5.6 high)

Repo: ~/repos/Ramanujan_Challenge, workdir problems/3.2/. You may run python3 freely.
Context files (read first): research/working_notes/FABLE_NOTES_energy_bootstrap.md sections 28-30 (tower theorem, structure theorem, quarter-point law); scripts research/scripts/quarterpoint_check.py and sp_spotcheck.py (working code for the series).

ESTABLISHED (verified): F(x) = sum b_k x^k (Apery zeta(3)); tau = sqrt(F), sigma = sqrt(F/q), q = 1-34x+x^2, as power series mod p. Quarter-point law: p ≡ 5 mod 24 => tau_{(p-1)/4} ≡ 0 mod p; p ≡ 23 mod 24 => sigma_{(p-3)/4} ≡ 0; the vanishing classes are exactly the primes represented by 2x^2+3y^2 (disc -24 non-principal genus). Half-integer recurrences: 4(j+2)^2 tau_{j+2} = 2(68j^2+170j+107) tau_{j+1} - (2j+1)^2 tau_j (sigma analogous — derive it).

MISSION (fast, empirical, report within your session):
1. For every prime p < 3000 in the vanishing classes (p ≡ 5, 23 mod 24): compute the FULL zero set Y_p = {j in truncation window: tau_j ≡ 0 (resp sigma_j ≡ 0) mod p}, and the representation p = 2x^2 + 3y^2 (x, y > 0; enumerate).
2. Hunt the parametrization: test whether zero POSITIONS beyond the forced quarter point are functions of (x, y) — candidates: j = (p-1)/4 ± linear combos of x, y; positions where j/p ≈ {x/(x+y)}-type ratios; #zeros vs class-group data.
3. Jacobsthal-value hunt: compute the NEAR-quarter values tau_{(p-1)/4 ± 1}, tau_{(p-1)/4 ± 2} mod p and test classical-style laws: ≡ c·x, c·y, c·xy, c·x/y mod p (fit the constant c across primes; also try 4x^2 mod p = (p-3y^2)-related forms). The classical template: p = a^2+b^2 => C((p-1)/2,(p-1)/4) ≡ 2a mod p (Gauss); analogues for disc -24 exist in Hudson-Williams-style literature (binomial coefficients and quartic/sextic residues).
4. Also test the NON-vanishing classes (p ≡ 1, 19, principal form x^2+6y^2): is tau_{(p-1)/4} ≡ c·(something in x,y) mod p there (nonzero but structured)?
5. Report: exact empirical laws found (with fitted constants and verification counts), zero-set size statistics per class, and the sharpest conjecture the data supports. Write findings to problems/3.2/CODEX_JACOBSTHAL_RECON.md and commit with message prefix "codex-high:".
Speed over depth. No proofs needed — laws + evidence.
