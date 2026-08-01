#!/usr/bin/env python3
"""CRON fresh-eyes: rank-2 recurrences for tau=sqrt(F), sigma=sqrt(F/q) and the
quarter-point law.

  4(j+2)^2 tau_{j+2}   = 2(68j^2+170j+107) tau_{j+1}   - (2j+1)^2 tau_j,   tau_0=1, tau_1=5/2
  4(j+2)^2 sigma_{j+2} = 2(68j^2+238j+209) sigma_{j+1} - (2j+3)^2 sigma_j, sigma_0=1, sigma_1=39/2

Quarter-point law (verified p < 1000, zero failures):
  p = 5 (24)  => tau_{(p-1)/4} = 0 (mod p);   p = 1 (24)  => nonzero.
  p = 23 (24) => sigma_{(p-3)/4} = 0 (mod p); p = 19 (24) => nonzero.
Genus-theoretic reading (disc -24): vanishing class = p represented by 2x^2+3y^2.
"""
from sympy import primerange

def tau_mod(p, top):
    T = [1 % p, 5 * pow(2, -1, p) % p]
    for j in range(0, top-1):
        inv = pow(4*(j+2)*(j+2) % p, -1, p)
        T.append((2*(68*j*j+170*j+107)*T[j+1] - (2*j+1)**2*T[j]) * inv % p)
    return T

def sigma_mod(p, top):
    S = [1 % p, 39 * pow(2, -1, p) % p]
    for j in range(0, top-1):
        inv = pow(4*(j+2)*(j+2) % p, -1, p)
        S.append((2*(68*j*j+238*j+209)*S[j+1] - (2*j+3)**2*S[j]) * inv % p)
    return S

if __name__ == "__main__":
    stats = {1: [0,0], 5: [0,0], 19: [0,0], 23: [0,0]}
    for p in primerange(13, 1000):
        c = p % 24
        if c in (1, 5):
            k = (p-1)//4
            z = tau_mod(p, k+1)[k] == 0
        elif c in (19, 23):
            k = (p-3)//4
            z = sigma_mod(p, k+1)[k] == 0
        else:
            continue
        stats[c][0] += z; stats[c][1] += 1
    for c in (5, 23, 1, 19):
        z, t = stats[c]
        print(f"p={c:2d} (mod 24): quarter-point zero {z}/{t}"
              + ("  (LAW: always)" if c in (5, 23) else "  (LAW: never)"))
