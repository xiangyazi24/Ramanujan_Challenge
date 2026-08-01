# DOCTRINE — Final Strike (2026-08-01 evening, automode)

## Goal (one sentence)

Prove the first unconditional exponent improvement |Z_p| ≪ p^{2/3−δ} (any δ>0)
— or the prime-averaged Poisson-mean theorem — by cracking the family
compatibility barrier at its sharpest point: the gap-polynomial family
{N_m mod p}_{m ≤ H} whose root counts control zeros via the restart identity.

## Standing facts (banked today, ledger §64–87)

- Restart identity: after every zero x_{s}=∞ the orbit restarts at 0 exactly;
  every zero is a root s of the gap polynomial N_m (m = gap from previous
  restart). |Z_p| ≤ 1 + Σ_{m≤H} #roots(N_m mod p) + p/H. Worst-case root
  bound deg N_m = O(m) gives the known 2/3. THE question: is the true
  Σ_{m≤H} #roots ~ H (Poisson) — and can any part of that be proved?
- Unconditional annealed Poisson (O(1/p)); LEMMA-CODEGREE (audited);
  SLIDING-WEIL-L; short-gap lemma (typical primes); Gal(M_h)=S_{m_h}
  certified h≤11 (cron); N̂_h irreducible over Q h≤14 (cron);
  joint Mellin group = SL₂×Sp₄.

## Avenues

(a) **Two-parameter surface point count** (= cron's top route, our side).
    The set {(s,m): N_m(s) ≡ 0 mod p, m ≤ H} is the F_p-point set of a
    two-parameter recurrence family (3-term recurrence in m, polynomial
    in s). Attack: (a1) empirical root statistics FIRST (is the truth ~H?);
    (a2) pair-correlation via Res(N_m, N_{m'}) and gcd structure (fixed
    pairs are Weil-able; sum over pairs may beat H² without full family
    theorem); (a3) algebraic model of the m-direction (continuant
    functional equations; Ostafe–Shparlinski orbit bounds applicability).
    Terminal: proof of Σ ≪ H^{2−δ} (⟹ |Z_p| ≪ p^{2/3−δ'}), or written
    proof-of-failure of each sub-vector.

(b) **Prime-averaged theorem, bankable now.** For fixed gap m, Chebotarev
    with the certified S_{m_h} groups gives the exact density of primes
    with a given root count; compute E_{p≤X} Σ_{m≤M₀} #roots for fixed M₀
    (exact main term), combine with the 2/3 tail ⟹ an unconditional
    average statement E_p[|Z_p| restricted to small gaps] = O(1)·M₀-terms.
    Terminal: theorem + proof written into the ledger, or the precise
    obstruction (uniformity in m of the Chebotarev error) named with its
    exact dependence.

(c) **Q6462 two-scale + four-mechanism harvest** (in-tab now): merge its
    verdicts into (a); if its two-scale decomposition gives an
    unconditional exponent via fixed-h Weil alone, execute immediately.

(d) **Fallback:** consolidate/write-up today's banked theorems (codegree
    repairs list, annealed chain) while probes run — never idle.

## Fallbacks if all avenues stall

Coordinate with cron's (R)/(C) five-piece package; the coloring-lemma BDH
variance route (their second-ranked) as joint target.

## Terminal conditions for the run

Success: an unconditional exponent < 2/3 proof-sketch machine-checked, or
the average theorem banked. Proof-of-failure: each sub-vector of (a) has a
written why-it-fails verdict AND (b)'s obstruction is precisely named.
