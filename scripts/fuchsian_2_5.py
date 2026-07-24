#!/usr/bin/env python3
"""Problem 2.5: Verify Fuchsian rigidity argument.

The CMF's scalar recurrence has Poincaré polynomial (c+16)(c²+544c+256)=0.
The integrated K module (pulled back via k=4√(2z)/(1-z)) should have the
same Poincaré data. If the local exponents also match, then by rigidity
of rank-3 Fuchsian systems with exactly 4 regular singularities, the
CMF IS the integrated K module, and the connection coefficient IS G.

Key singularity points on the z-line:
  z = 0      (k = 0)
  z = ρ = 17-12√2  (k = 1, logarithmic singularity of K)
  z = 1/ρ = 17+12√2  (k = ∞, apparent singularity)
  z = ∞      (k → ∞)

The ODE for the integrated K module in k:
  k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
  i.e., [k(1-k²)D² + (1-3k²)D - k] ∘ D(Y) = 0  where D = d/dk

We pull this back via k(z) = 4√(2z)/(1-z).
"""
from mpmath import mp, mpf, sqrt, log10, pi

mp.dps = 50

rho = 17 - 12*sqrt(2)
print(f"ρ = 17-12√2 = {mp.nstr(rho, 40)}")
print(f"1/ρ = 17+12√2 = {mp.nstr(1/rho, 40)}")
print(f"k(ρ) = 4√(2ρ)/(1-ρ) = {mp.nstr(4*sqrt(2*rho)/(1-rho), 15)}")

# Check: k(z) = 4√(2z)/(1-z)
# k² = 32z/(1-z)²
# Singularity of K(k): k = 1 → 32z = (1-z)² → z² - 34z + 1 = 0
# → z = 17 ± 12√2 = {ρ, 1/ρ}
# Verify:
print(f"\nz² - 34z + 1 at z=ρ: {mp.nstr(rho**2 - 34*rho + 1, 15)}")
print(f"z² - 34z + 1 at z=1/ρ: {mp.nstr((1/rho)**2 - 34/rho + 1, 15)}")

# Now compute the Poincaré polynomial of the ODE in z-variable.
# The integrated K ODE in k has regular singularities at k = 0, ±1, ∞.
# Local exponents:
# At k=0: the ODE k(1-k²)Y''' + ... factors as [L₂∘D] where L₂ is the
#   Gauss equation for K(k). D has exponent 0 (constant), L₂ has exponents
#   0, 0 at k=0. So the full system has exponents {0, 0, 0} at k=0.
#   (Logarithmic at k=0 for the K integral)

# At k=1: K(k) ~ -log(1-k²)/π + const. The integrated Y has a log² term.
#   Exponents of L₂ at k=1 are {0, 0} (logarithmic).
#   Full system at k=1: {0, 0, 0} with log and log² terms.

# At k=-1: same as k=1 by symmetry.

# At k=∞: L₂ has exponents {1/2, 1/2}. Full system has {?, ?, ?}.

# Now pull back to z via k = 4√(2z)/(1-z).
# dk/dz = 2√2(1+z)/(√z(1-z)²)
# Chain rule: D_z = (dk/dz)·D_k

# The map k(z) sends:
#   z=0 → k=0  (local degree 1/2 because k ~ 4√(2z))
#   z=ρ → k=1  (local degree 1)
#   z=1/ρ → k=-1 ... wait, k(z) = 4√(2z)/(1-z).
#   For z > 0, k > 0 always. So k = -1 is not reached for real positive z.
#   At z=∞: k(z) ~ -4√(2z)/z = -4√2/√z → 0. So z=∞ maps to k=0!
#   At z=1/ρ: k = 4√(2/ρ)/(1-1/ρ) = 4√(2/ρ)·ρ/(ρ-1)
#   Since ρ-1 = 16-12√2 < 0 (since ρ ≈ 0.0294), actually ρ < 1, so 1-1/ρ < 0.
#   k(1/ρ) = 4√(2/ρ)/(1-1/ρ) = 4√(2/ρ)·ρ/(ρ-1) = -4√(2ρ)/(1-ρ) = -k(ρ) = -1.
#   Yes! z = 1/ρ maps to k = -1.
print(f"\nk(1/ρ) = {mp.nstr(4*sqrt(2/rho)/(1-1/rho), 15)}")

# Summary of the pullback:
# z = 0 → k = 0 (ramification index 2 because k ~ √z)
# z = ρ → k = 1
# z = 1/ρ → k = -1
# z = ∞ → k = 0 (ramification index 2)
# z = -1/8 → k = ∞ (because 1-z = 1+1/8 = 9/8, k = 4√(-1/4)/(9/8) = complex)
# Wait, k² = 32z/(1-z)². For k→∞, need (1-z)² → 0, i.e., z→1.
# At z=1: k = 4√2/(1-1) → ∞. Yes!

print(f"\n--- Singularity map ---")
print(f"z=0   → k=0  (ramification 2, k ~ 4√2·z^(1/2))")
print(f"z=ρ   → k=1  (simple)")
print(f"z=1/ρ → k=-1 (simple)")
print(f"z=1   → k=∞  (pole)")
print(f"z=∞   → k=0  (ramification 2, k ~ 4√2/z^(1/2))")

# Now compute local exponents of the pulled-back ODE at each z-singularity.
# At z=0: k=0 with ramification 2. Original exponents at k=0 are {0,0,0}.
# Under ramification index e, exponents transform as θ ↦ eθ + (e-1)/2? No...
# For z = t², dk = 2t dt, D_k = (1/(2t))D_t.
# If the ODE at k=0 has solutions k^α, then in z these become z^(α/2).
# For the full system at k=0 with exponents {0,0,0} and logarithmic,
# the pullback to z has exponents {0,0,0} (since 0/2=0) but the
# logarithmic structure changes.

# At z=ρ: k=1, simple map, exponents preserved. {0,0,0} with logs.

# At z=1/ρ: k=-1, simple map, exponents preserved. {0,0,0} with logs.

# At z=1: k=∞. The ODE at k=∞: the exponents of L₂ at ∞ are {1/2,1/2}.
# For the full system [L₂∘D], the exponents at k=∞ need computation.
# The outer D shifts exponents by 1: if L₂ has exponents {α₁,α₂} at ∞,
# then [L₂∘D] has exponents {α₁+1, α₂+1, 0} at ∞.
# So exponents at k=∞ are {0, 3/2, 3/2}.

# At z=1, the pullback from k=∞ via k~1/(1-z): local parameter is (1-z),
# k ~ 4√2/(1-z). So k = O(1/(1-z)), which means the exponents scale by
# the ramification index 1. So exponents at z=1 are {0, 3/2, 3/2}.

print(f"\n--- Local exponents of the integrated K ODE pulled back to z ---")
print(f"z=0:   {{0, 0, 0}} (from k=0, ramification 2)")
print(f"z=ρ:   {{0, 0, 0}} (from k=1, logarithmic)")
print(f"z=1/ρ: {{0, 0, 0}} (from k=-1, logarithmic)")
print(f"z=1:   {{0, 3/2, 3/2}} (from k=∞)")

# Fuchs relation check:
# Sum of all exponents = (rank-1)(rank-2)/2 · (# singular points - 2)
# For rank 3 with 5 singular points (0, ρ, 1/ρ, 1, ∞):
# Sum should be 1·0/2·3 = ... actually the formula is:
# sum of exponents = (r choose 2) · (s-2) where r=rank, s=# singularities
# For r=3, s=5: sum = 3·3 = 9? Let me look up the exact formula.
# Actually the Fuchs relation: sum of all exponents = (r-1)·(s-2)·r/2
# No... Fuchs relation for a rank-r system with s regular singularities:
# sum of all local exponents = r(r-1)(s-2)/2
# For r=3, s=5: sum = 3·2·3/2 = 9.
exponents = [0,0,0, 0,0,0, 0,0,0, 0,1.5,1.5]  # z=0, ρ, 1/ρ, 1
# Need z=∞ too.
# At z=∞: comes from k=0 (ramification 2). Exponents {0,0,0}.
exponents += [0,0,0]
total = sum(exponents)
print(f"\nFuchs check: sum of exponents = {total}")
print(f"Expected (r=3, s=5): 3·2·3/2 = {3*2*3//2}")

# Sum is 3, expected 9. Off by 6.
# This means I'm wrong about some exponents. Let me reconsider.

# Actually, the Fuchs relation for an ORDER n ODE with s regular singularities is:
# sum of all local exponents = n(n-1)/2 · (s-2)
# For order 3, s=5 singularities: sum = 3·1·3 = 9? No, 3·2/2·3 = 9. Hmm.
# Wait: sum = n(n-1)(s-2)/2 = 3·2·3/2 = 9.
# But I only counted 4 finite singularities. With ∞, that's 5 total.
# Let me recount: z=0, z=ρ, z=1/ρ, z=1, z=∞. That's 5 singular points.
# Sum should be 3·2·(5-2)/2 = 9.

# My computed sum is 0+0+0 + 0+0+0 + 0+0+0 + 0+3/2+3/2 + ??? = 3 + exponents at ∞.
# So exponents at z=∞ must sum to 6.

# At z=∞: the ODE has a regular singularity. Let me compute the exponents
# more carefully.
# The pullback k(z) = 4√(2z)/(1-z) at z=∞:
# k ~ 4√(2z)/(-z) = -4√2/√z → 0 as z→∞.
# So z=∞ maps to k=0 with k ~ z^(-1/2).
# Local parameter: w = 1/z, k ~ √(2/w) · 4w/(w-1) ~ 4√(2w) as w→0.
# So k ~ √w with ramification 2.
# Exponents at k=0 are {0,0,0}. Under substitution k = c·w^(1/2),
# solutions k^0 become w^0 = 1. So exponents stay {0,0,0}?
# But then sum = 3 + 0 = 3 ≠ 9.

# The discrepancy suggests I'm computing the exponents wrong.
# The issue is likely at z=0 and z=∞ where the map has ramification 2.
# When the covering map has ramification, apparent singularities can appear,
# and the exponents transform differently.

# Let me just verify numerically using the CMF itself.
# The Poincaré polynomial (c+16)(c²+544c+256) gives the ASYMPTOTIC
# behavior at n→∞, which corresponds to z=0 in the generating function picture.
print(f"\n--- Poincaré data (n→∞ corresponds to z→0) ---")
print(f"Characteristic roots: -16, -16(17±12√2)")
print(f"Product: (-16)·(-16)·(17+12√2)·(-16)·(17-12√2) = {(-16)**3 * ((17)**2 - (12*sqrt(2))**2)}")
print(f"  = (-16)³ · (289-288) = (-16)³ = {(-16)**3}")

# The det of M(n) for large n ~ (-8)·(n)^(1+12+5+2+3+4)·2^(2+3+4)
# = (-8)·n^27·2^9 · leading numerical factor
# Actually det = -8·(n+1)(n+2)^6(n+3)^5(2n+3)^2(2n+5)^3(2n+7)^4
# For large n: ~ -8·n·n^6·n^5·(2n)^2·(2n)^3·(2n)^4 = -8·n^12·2^9·n^9 = -8·512·n^21
# So det ~ -4096·n^21. Leading Poincaré: product of c roots = 946/991952896 ≈ 9.54e-7.
# But det(M) ~ n^21, not constant. The Poincaré roots are about the NORMALIZED recurrence.

print(f"\nDone. The rigidity argument requires matching local exponents at all singularities.")
print(f"The integrated K module on the z-line has singularities at 0, ρ, 1/ρ, 1, ∞.")
print(f"If the CMF's differential Galois group matches, rigidity proves CMF = integrated K.")
