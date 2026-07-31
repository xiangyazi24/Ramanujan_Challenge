#!/usr/bin/env python3
"""
Q2714 same-kernel quotient-digit rigidity experiment (v2, corrected).

For each n, find target primes p in (n/2, n] with p | b_n,
compute c_p = (b_n / p) mod p exactly,
then test for cross-prime structure grouped by kernel g = gcd(p-1, n-1).

Tests:
  1. CRT consistency: is there a small integer X such that X = c_p mod p for all targets?
  2. Affine: c_p = A*u + B mod p where p = 1 + g*u?
  3. Cross-prime pairs: for same-kernel primes p,q, does c_p*c_q or c_p+c_q have structure?
"""

import sys
from math import gcd, isqrt, log
from collections import defaultdict


def sieve_primes(limit):
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = 0
    return [i for i in range(2, limit + 1) if is_prime[i]]


def apery_b_exact(n):
    """Compute b_n exactly."""
    if n == 0:
        return 1
    if n == 1:
        return 5
    b_prev, b_curr = 1, 5
    for r in range(1, n):
        num = (34*r**3 + 51*r**2 + 27*r + 5) * b_curr - r**3 * b_prev
        b_next = num // (r + 1)**3
        b_prev, b_curr = b_curr, b_next
    return b_curr


def find_targets(n, bn, primes):
    """Find target primes and compute quotient digits exactly."""
    targets = []
    for p in primes:
        if p <= n // 2 or p > n:
            continue
        if bn % p != 0:
            continue
        c_p = (bn // p) % p
        vp = 0
        temp = bn
        while temp % p == 0:
            temp //= p
            vp += 1
        r_p = n - p
        g = gcd(p - 1, n - 1)
        u = (p - 1) // g
        v = (n - 1) // g
        targets.append({
            'p': p, 'c_p': c_p, 'r_p': r_p,
            'g': g, 'u': u, 'v': v, 'vp': vp
        })
    return targets


def crt_two(a1, m1, a2, m2):
    g = gcd(m1, m2)
    if (a1 - a2) % g != 0:
        return None
    lcm = m1 * m2 // g
    _, x, _ = extended_gcd(m1 // g, m2 // g)
    diff = (a2 - a1) // g
    t = (diff * x) % (m2 // g)
    return (a1 + m1 * t) % lcm


def multi_crt(pairs):
    if not pairs:
        return None, 0
    X, M = pairs[0]
    for a, m in pairs[1:]:
        X_new = crt_two(X, M, a, m)
        if X_new is None:
            return None, 0
        M = M * m // gcd(M, m)
        X = X_new
    return X, M


def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x1, y1 = extended_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def test_crt_lift(targets_by_kernel, label=""):
    """Test if c_p values in each kernel class have a common small CRT lift."""
    results = []
    for g, tgts in sorted(targets_by_kernel.items()):
        if len(tgts) < 2:
            continue
        # Only use targets with vp = 1 (c_p nonzero)
        good = [t for t in tgts if t['vp'] == 1]
        if len(good) < 2:
            continue
        X, M = multi_crt([(t['c_p'], t['p']) for t in good])
        if X is not None:
            if X > M // 2:
                X -= M
            ratio = abs(X) / M if M > 0 else 1
            log_X = log(abs(X) + 1) if X != 0 else 0
            log_M = log(M) if M > 0 else 0
            results.append((g, len(good), X, log_X, log_M, ratio))
    return results


def test_affine_residue(targets_by_kernel):
    """For each kernel class, test if c_p = A*u + B mod p (lifted to integer).
    Use first two to determine (A,B), check rest."""
    results = []
    for g, tgts in sorted(targets_by_kernel.items()):
        good = [t for t in tgts if t['vp'] == 1]
        if len(good) < 3:
            continue
        # Try all pairs as basis, check remaining
        best_match = 0
        for i in range(min(len(good), 5)):
            for j in range(i+1, min(len(good), 5)):
                t1, t2 = good[i], good[j]
                p1, p2 = t1['p'], t2['p']
                # c = A*u + B mod p
                # c1 = A*u1 + B mod p1
                # c2 = A*u2 + B mod p2
                # Need to lift to integers: try small A, check B
                matched = 0
                for A_try in range(-50, 51):
                    B1 = (t1['c_p'] - A_try * t1['u']) % p1
                    B2 = (t2['c_p'] - A_try * t2['u']) % p2
                    # Check consistency via CRT
                    B_crt = crt_two(B1, p1, B2, p2)
                    if B_crt is None:
                        continue
                    # Check remaining targets
                    cnt = 2
                    for k in range(len(good)):
                        if k == i or k == j:
                            continue
                        tk = good[k]
                        predicted = (A_try * tk['u'] + B_crt) % tk['p']
                        if predicted == tk['c_p']:
                            cnt += 1
                    if cnt > matched:
                        matched = cnt
                best_match = max(best_match, matched)
        results.append((g, len(good), best_match))
    return results


def main():
    max_n = 2000
    if len(sys.argv) > 1:
        max_n = int(sys.argv[1])

    primes = sieve_primes(max_n + 1)

    print(f"Same-kernel quotient-digit experiment, n up to {max_n}")
    print(f"{'n':>6} {'#tgt':>5} {'#v1':>4} {'#kern':>5}")
    print("-" * 35)

    all_by_n = {}
    hit_ns = []

    for n in range(10, max_n + 1):
        bn = apery_b_exact(n)
        targets = find_targets(n, bn, primes)
        if len(targets) >= 2:
            by_kernel = defaultdict(list)
            for t in targets:
                by_kernel[t['g']].append(t)
            n_v1 = sum(1 for t in targets if t['vp'] == 1)
            print(f"{n:>6} {len(targets):>5} {n_v1:>4} {len(by_kernel):>5}", end="")

            # CRT test
            crt_res = test_crt_lift(by_kernel)
            for g, cnt, X, lx, lm, ratio in crt_res:
                if ratio < 0.1:
                    print(f"  [g={g}: {cnt}tgt, X={X}, lnX/lnM={lx/lm:.3f}]", end="")
            print()

            all_by_n[n] = targets
            hit_ns.append(n)

    # Summary statistics
    print(f"\n{'='*60}")
    print(f"Summary: {len(hit_ns)} values of n with >= 2 target primes")

    # Global CRT test: across ALL target primes for a given n
    print(f"\nGlobal CRT lift analysis (all targets for each n):")
    print(f"{'n':>6} {'#tgt':>5} {'ln|X|':>8} {'lnM':>8} {'ratio':>8}")
    print("-" * 45)
    for n in hit_ns[:30]:
        targets = all_by_n[n]
        good = [t for t in targets if t['vp'] == 1]
        if len(good) < 2:
            continue
        X, M = multi_crt([(t['c_p'], t['p']) for t in good])
        if X is not None:
            if X > M // 2:
                X -= M
            lx = log(abs(X) + 1)
            lm = log(M)
            ratio = lx / lm
            print(f"{n:>6} {len(good):>5} {lx:>8.2f} {lm:>8.2f} {ratio:>8.4f}")

    # Cross-n kernel structure: for a fixed kernel g,
    # does the c_p pattern repeat across different n?
    print(f"\nKernel-1 (g=1) targets across n values:")
    g1_data = []
    for n in hit_ns:
        for t in all_by_n[n]:
            if t['g'] == 1 and t['vp'] == 1:
                g1_data.append((n, t['p'], t['c_p'], t['r_p'], t['u']))
    print(f"Total g=1, v_p=1 targets: {len(g1_data)}")

    # For the same p appearing as target of different n values:
    by_prime = defaultdict(list)
    for n, p, c, r, u in g1_data:
        by_prime[p].append((n, c, r))
    multi_n = {p: vs for p, vs in by_prime.items() if len(vs) >= 2}
    if multi_n:
        print(f"\nPrimes appearing as targets for multiple n ({len(multi_n)} primes):")
        for p, vs in sorted(multi_n.items())[:20]:
            print(f"  p={p}: ", end="")
            for n, c, r in vs:
                print(f"n={n}(r={r},c={c}) ", end="")
            print()


if __name__ == "__main__":
    main()
