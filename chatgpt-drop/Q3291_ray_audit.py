#!/usr/bin/env python3
from __future__ import annotations

# Triggered by the temporary Q3291 GitHub Actions runner.
import hashlib
import sys
from collections import defaultdict
from math import gcd

sys.set_int_max_str_digits(0)

NMAX = 5000


def P(n: int) -> int:
    return 34*n**3 + 51*n**2 + 27*n + 5


def main() -> None:
    E = [0] * (NMAX + 2)
    b = [0] * (NMAX + 2)
    b[0], b[1] = 1, 5

    Y_prev, Y_cur = 0, 6
    F = 1
    E[0] = 1
    E[1] = F // gcd(F, Y_cur)

    for n in range(1, NMAX + 1):
        if n >= 2:
            F *= n**3
            E[n] = F // gcd(F, Y_cur)

        den = (n + 1)**3
        nb = P(n)*b[n] - n**3*b[n-1]
        assert nb % den == 0
        b[n+1] = nb // den

        Y_next = P(n)*Y_cur - n**6*Y_prev
        Y_prev, Y_cur = Y_cur, Y_next

    F *= (NMAX + 1)**3
    E[NMAX + 1] = F // gcd(F, Y_cur)

    groups: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    digest = hashlib.sha256()

    for n in range(1, NMAX + 1):
        gm = gcd(E[n-1], E[n])
        gp = gcd(E[n], E[n+1])
        Tm = E[n-1] // gm
        Tp = E[n+1] // gp
        mum = gcd(6*Tm, n**3)
        mup = gcd(6*Tp, (n+1)**3)
        Xm = 6*Tm // mum
        Xp = 6*Tp // mup
        W = gcd(Xm, Xp)
        ray = (Xm // W, Xp // W)
        assert gcd(*ray) == 1
        assert Xm == ray[0]*W and Xp == ray[1]*W
        assert (6*b[n]) % W == 0
        groups[ray].append((n, W))
        digest.update(f"{n}:{ray[0]}:{ray[1]}:{W}\n".encode())

    repeated = [(ray, occ) for ray, occ in groups.items() if len(occ) >= 2]
    repeated.sort(key=lambda item: item[1][0][0])
    repeated_w = [
        (ray, [(n, W) for n, W in occ if W > 1])
        for ray, occ in repeated
        if any(W > 1 for _, W in occ)
    ]

    max_all = max((len(occ), ray, occ) for ray, occ in groups.items())
    max_w = max(
        ((sum(W > 1 for _, W in occ), ray, occ) for ray, occ in groups.items()),
        default=(0, (0, 0), []),
    )

    print(f"NMAX={NMAX}")
    print(f"record_sha256={digest.hexdigest()}")
    print(f"distinct_rays={len(groups)}")
    print(f"repeated_ray_count={len(repeated)}")
    print(f"repeated_rays_with_W_gt_1_count={len(repeated_w)}")
    print(f"max_total_repetition={max_all[0]} ray={max_all[1]}")
    print("max_total_occurrences=" + ",".join(f"{n}:{W}" for n, W in max_all[2]))
    print(f"max_W_gt_1_repetition={max_w[0]} ray={max_w[1]}")
    print("max_W_gt_1_occurrences=" + ",".join(f"{n}:{W}" for n, W in max_w[2] if W > 1))

    print("\n=== EVERY REPEATED RAY ===")
    for ray, occ in repeated:
        print(f"ray={ray[0]}/{ray[1]} count={len(occ)} W_gt_1={sum(W > 1 for _, W in occ)}")
        print("occurrences=" + ",".join(f"{n}:{W}" for n, W in occ))

    print("\n=== EVERY W>1 OCCURRENCE ON A REPEATED RAY ===")
    for ray, occ in repeated_w:
        print(f"ray={ray[0]}/{ray[1]} count={len(occ)}")
        print("occurrences=" + ",".join(f"{n}:{W}" for n, W in occ))


if __name__ == "__main__":
    main()
