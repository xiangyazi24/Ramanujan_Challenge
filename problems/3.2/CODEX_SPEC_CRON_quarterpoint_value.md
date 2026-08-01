# SPEC: computational discovery — closed form for tau_{(p-1)/4} mod p (P3.2 quarter-point value)

## Context
tau: 4(j+2)^2 tau_{j+2} = 2(68j^2+170j+107) tau_{j+1} - (2j+1)^2 tau_j, tau_0=1, tau_1=5/2 (rationals; mod p use inverse of 2).
Proven: for chi(p)=(-6|p)=+1, p ≡ 1 mod 4: tau_{(p-1)/4} ≡ 0 mod p iff p ≡ 5 mod 24. OPEN: closed form of v_p := tau_{(p-1)/4} mod p for p ≡ 1 (mod 24).
Known failures (do not re-test alone): v^2 vs 4x^2-2p, 2x^2, x^2, 6y^2, 24y^2, 4a^2-2p; and v/C((p-1)/2,(p-1)/4) squared vs same list. Reference sequences implementation: problems/3.2/CRON_sqrt_recurrences.py.

## Task
1. Compute v_p for ALL p ≡ 1 (mod 24), p < 20000 (Python ok; a few hundred primes). For each p also compute: representations p = x^2+6y^2 (normalize x>0,y>0), p = a^2+b^2 (a odd, a ≡ 1 mod 4), p = c^2+2d^2 if applicable, binomials G1 = C((p-1)/2,(p-1)/4) mod p, G2 = C((p-1)/4, (p-1)/8) if defined, Jacobsthal sums phi(k) = sum_t legendre(t(t^2+k),p) for k in {1,2,3,6,-1,-2,-3,-6}, quartic residue symbols of 2,3,6, and the sqrt-of-b midpoint tau_{(p-1)/2} (≡ (-2|p), sanity).
2. Systematic fit: search for identities of the form v ≡ u * G1^e1 * (algebraic monomial in x,y,a,b,phi-values) with u a bounded rational constant and e1 in {-1,0,1}, monomial degree <= 2. Also test v^2 and v^4 against degree-<=2 polynomials in (x,y,a,b,p) with small coefficients (LLL/exhaustive small-coeff search over a few primes, then verify on all). Also test whether v/G1 or v*G1 is +-x, +-2x, +-a, +-2a, +-(x+-...) TIMES a quartic character of a small number (i.e., allow a 4th-root-of-unity-valued factor depending on p mod 16 / quartic class of 2 or 3).
3. Any identity that fits >= 10 primes must be verified on ALL computed primes; report exact hit/miss counts. If nothing closes, report the three best near-misses with their exact failure sets AND the empirical distribution of v/(2a*quartic-sign-guess) to guide theory.
4. Write problems/3.2/CRON_quarterpoint_value_report.md with the verified identity (or negative report), the code as CRON_quarterpoint_value.py, git add + commit.

## Discipline
- New files only (CRON_ prefix); no modification of existing files; machine-verify all claims.
