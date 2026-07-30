#!/usr/bin/env python3
"""End-to-end + stepwise verification of H2 == 2p*Xi_p (mod p^3), p>=7.
Steps: (A) H2 == p[6B_{2p-4}-8B_{3p-5}+3B_{4p-6}] mod p^3 (Faulhaber);
(B) strong Kummer Delta^2(B_m/m) == 0 mod p^2 along m=(p-3)+k(p-1);
(C) end-to-end identity."""
PMAX = 600
primes = [p for p in range(7, PMAX) if all(p % q for q in range(2, int(p**.5)+1))]
gA = gB = gC = 0
for p in primes:
    p3, p4 = p**3, p**4
    SB = lambda k: sum(pow(s, k, p4) for s in range(1, p)) % p4
    H2 = sum(pow(s, -2, p3) for s in range(1, p)) % p3
    assert H2 == (6*SB(2*p-4) - 8*SB(3*p-5) + 3*SB(4*p-6)) % p3; gA += 1
    m = [p-3+k*(p-1) for k in range(4)]
    beta = [SB(mk)*pow(mk, -1, p3) % p3 for mk in m]      # p * B_m/m mod p^3
    assert (beta[2]-2*beta[1]+beta[0]) % p3 == 0
    assert (beta[3]-2*beta[2]+beta[1]) % p3 == 0; gB += 1
    rhs = (2*SB(2*p-4)*pow(2*p-4, -1, p3) - 4*SB(p-3)*pow(p-3, -1, p3)) % p3
    assert H2 == rhs; gC += 1
print(f"primes {len(primes)}: step A {gA}, step B {gB}, end-to-end {gC}, all pass")
