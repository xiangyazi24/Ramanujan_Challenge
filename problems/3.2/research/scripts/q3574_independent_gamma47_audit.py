#!/usr/bin/env python3
"""Independent exact audit of the Q3573 Gamma_47 carrier.

This file is intentionally standalone: it imports no Ramanujan_Challenge
computation code and does not reuse the Fable implementation.  It reconstructs
Gamma_47 directly from the definitions recorded in Q3573/Q3445:

  * Apéry b_n from the three-term recurrence;
  * sampled slope-s forward differences d^(s)_K;
  * U_(s,K)=binom(-X_s,K) and full legal-family gcd G_s;
  * T_(s,K)=sum_{r<K} U_(s,r)d^(s)_r;
  * canonical primitive three-row C*_(s,K)=a V_K+b V_(K+1);
  * Gamma_47 = gcd of all C* over s=4,7 and legal K,K+1,K+2.

The q=6 legal windows are reconstructed from the exact Q3445 fixture formula.
All arithmetic determining b, selectors, carriers, gcds, and target products is
integer-exact.  No third-party packages are used.
"""

from math import comb, gcd, isqrt, prod

H = 0
DELTAS = (0, 1, 2)
FIXTURES = (
    (19, 8),
    (97, 25),
    (139, 61),
    (181, 19),
    (293, 47),
)
COLLISION_M = 2932
COLLISION_PAIRS = ((439, 298), (443, 274))


def apery_numbers(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    b = [1, 5]
    for n in range(1, limit):
        num = (34*n**3 + 51*n**2 + 27*n + 5) * b[n] - n**3 * b[n-1]
        den = (n + 1)**3
        assert num % den == 0
        b.append(num // den)
    return b[:limit+1]


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0:1] = b"\x00"
    if limit >= 1:
        sieve[1:2] = b"\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00" * (((limit - p*p)//p) + 1)
    return [p for p in range(2, limit+1) if sieve[p]]


def sigma(m: int, h: int, s: int) -> int:
    hits = [x for x in range(h+1, h+s+1) if (m + 1 + x) % s == 0]
    assert len(hits) == 1
    return hits[0]


def legal_windows(m: int, h: int, delta2: int):
    s2 = sigma(m, h, 2)
    F2 = 3*s2 + 2*h - 3 + 2*delta2
    s4 = sigma(m, h, 4)
    s7 = sigma(m, h, 7)
    X4 = (m + 1 + s4)//4
    X7 = (m + 1 + s7)//7
    Phi = (m - F2)//6
    A4 = max(1, 2*Phi - X4 + 1)
    B4 = min(X7 - 1, 1 + (m - s4)//4)
    A7 = max(1, Phi - X7 + 1)
    B7 = min(X7 - 1, 1 + (m - s7)//7)
    return F2, s4, X4, s7, X7, A4, B4, A7, B7


def neg_binom(X: int, k: int) -> int:
    assert X >= 1 and k >= 0
    return (-1 if k & 1 else 1) * comb(X + k - 1, k)


def first_forward_differences(values: list[int]) -> list[int]:
    row = values[:]
    out = []
    while row:
        out.append(row[0])
        row = [row[i+1] - row[i] for i in range(len(row)-1)]
    return out


def slope_carriers(m: int, s: int, sigma_s: int, X: int, A: int, B: int,
                   apery: list[int]):
    assert A <= B - 2, (m, s, A, B)
    # C*_K for K<=B-2 uses d_K,d_(K+1), hence samples through y_(B-1).
    sampled = [apery[sigma_s + s*t] for t in range(B)]
    d = first_forward_differences(sampled)
    U = [neg_binom(X, k) for k in range(B+1)]

    G = 0
    for k in range(A, B+1):
        G = gcd(G, abs(U[k]))
    assert G > 0
    e = [u // G for u in U]
    assert all(U[k] % G == 0 for k in range(A, B+1))
    egcd = 0
    for k in range(A, B+1):
        egcd = gcd(egcd, abs(e[k]))
    assert egcd == 1

    # T_K = sum_{r=0}^{K-1} U_r d_r.  We need T through B.
    T = [0] * (B + 1)
    for k in range(B):
        T[k+1] = T[k] + U[k] * d[k]

    carriers = []
    for k in range(A, B-1):
        Ak = X + 2*k + 1
        if Ak & 1:
            a = (Ak + 1)//2
            bb = -(Ak - 1)//2
            gs = 1
        else:
            a = -1
            bb = 1
            gs = 2
        assert a*Ak + bb*(Ak+2) == gs

        Vk = (X+k)*T[k] + (k+1)*T[k+1]
        Vkp1 = (X+k+1)*T[k+1] + (k+2)*T[k+2]
        C = a*Vk + bb*Vkp1

        # Independent normal-form consistency check from Q3573 (11)-(12).
        Z = ((a*(k+1) + bb*(Ak+2)) * e[k] * d[k]
             + bb*(k+2) * e[k+1] * d[k+1])
        assert C == gs*T[k] + G*Z
        carriers.append((k, C))

    return G, carriers


def gamma47(m: int, delta2: int, apery: list[int]):
    F2,s4,X4,s7,X7,A4,B4,A7,B7 = legal_windows(m, H, delta2)
    G4, C4 = slope_carriers(m, 4, s4, X4, A4, B4, apery)
    G7, C7 = slope_carriers(m, 7, s7, X7, A7, B7, apery)
    gamma = 0
    for _, value in C4 + C7:
        gamma = gcd(gamma, abs(value))
    return {
        "m": m, "delta2": delta2, "F2": F2,
        "s4": s4, "X4": X4, "A4": A4, "B4": B4, "G4": G4,
        "s7": s7, "X7": X7, "A7": A7, "B7": B7, "G7": G7,
        "count4": len(C4), "count7": len(C7), "gamma": gamma,
    }


def q6_targets(m: int, apery: list[int], primes: list[int]):
    targets = []
    for p in primes:
        q, j = divmod(m, p)
        if q == 6 and apery[j] % p == 0:
            targets.append((p, j))
    return targets


def vp(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def main():
    fixture_ms = [6*p+j for p,j in FIXTURES]
    max_m = max(fixture_ms + [COLLISION_M])
    apery = apery_numbers(max_m)
    primes = primes_upto(max_m)

    print("Q3574 INDEPENDENT GAMMA47 AUDIT")
    print("implementation=standalone exact recurrence/Newton/Bézout; no project imports")
    print(f"H={H} deltas={DELTAS}")

    for p,j in FIXTURES:
        m = 6*p+j
        assert m//p == 6 and m%p == j
        targets = q6_targets(m, apery, primes)
        target_product = prod(q for q,_ in targets)
        print(f"FIXTURE p={p} j={j} m={m} b_j_mod_p={apery[j] % p} vp_bj={vp(apery[j],p)}")
        print(f"  q6_targets={targets} target_product={target_product}")
        for delta2 in DELTAS:
            rec = gamma47(m, delta2, apery)
            gamma = rec["gamma"]
            assert gamma % target_product == 0, (m, delta2, gamma, target_product)
            eps = gamma // target_product
            print(
                "  delta2={delta2} F2={F2} "
                "I4=[{A4},{B4}] I7=[{A7},{B7}] "
                "count=({count4},{count7}) gamma={gamma} epsilon={eps}".format(
                    eps=eps, **rec
                )
            )

    print(f"COLLISION m={COLLISION_M}")
    for p,j in COLLISION_PAIRS:
        assert COLLISION_M == 6*p+j
        print(
            f"  p={p} j={j} b_j_mod_p={apery[j] % p} "
            f"vp_bj={vp(apery[j],p)}"
        )
    targets = q6_targets(COLLISION_M, apery, primes)
    target_product = prod(q for q,_ in targets)
    print(f"  q6_targets={targets} target_product={target_product}")
    for delta2 in DELTAS:
        rec = gamma47(COLLISION_M, delta2, apery)
        gamma = rec["gamma"]
        assert gamma % target_product == 0
        eps = gamma // target_product
        print(
            "  delta2={delta2} F2={F2} "
            "I4=[{A4},{B4}] I7=[{A7},{B7}] "
            "count=({count4},{count7}) gamma={gamma} epsilon={eps}".format(
                eps=eps, **rec
            )
        )

    claimed = 5*439*443
    print(f"CLAIMED_COLLISION_VALUE 5*439*443={claimed}")


if __name__ == "__main__":
    main()
