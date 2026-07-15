#!/usr/bin/env python3
"""Exact modular screen for polynomial Ore intertwiners R.

Tests R = sum_{i=0}^d r_i(n) S^i, deg r_i <= D, in
    L27 * R = Q * LZ
by right-reducing L27*R modulo LZ.  All arithmetic is exact in F_p;
full column rank at one prime proves that no nonzero polynomial R exists
over Q for the stated (d,D) ansatz.

Ore rule: S f(n) = f(n+1) S.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Sequence, Tuple


def A(n: int) -> int:
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)


def B(n: int) -> int:
    return 128*(2*n+7)**3*(2*n+9)**3*(
        104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3
        +92943995*n*n+102256019*n+46709052
    )


def C(n: int) -> int:
    return 16*(n+3)**4*(2*n+9)**3*(
        3784*n**5+57792*n**4+351019*n**3+1059230*n*n+1587211*n+944620
    )


def Dcoef(n: int) -> int:
    return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)


def zp3(n: int) -> int:
    return 2*(946*n*n-731*n+153)*(2*n+1)*(n+1)**3


def zp2(n: int) -> int:
    return -2*(
        104060*n**6+127710*n**5+12788*n**4-34525*n**3
        -8482*n*n+3298*n+1071
    )


def zp1(n: int) -> int:
    return 2*(3784*n**5-1032*n**4-1925*n**3+853*n*n+328*n-184)*n


def zp0(n: int) -> int:
    return -(946*n*n+1161*n+368)*n*(n-1)**3


def l27_at(n: int, p: int) -> List[int]:
    return [(-Dcoef(n)) % p, C(n+1) % p, (-B(n+2)) % p, A(n+2) % p]


def lz_at(n: int, p: int) -> List[int]:
    # LZ = sum_j zpj(n+2) S^j
    return [zp0(n+2) % p, zp1(n+2) % p, zp2(n+2) % p, zp3(n+2) % p]


def basis_remainder(i: int, k: int, n0: int, p: int) -> Optional[Tuple[int,int,int]]:
    """Evaluate right remainder of L27 * (n^k S^i) mod LZ at n=n0."""
    l27 = l27_at(n0, p)
    coeff = [0] * (i + 4)
    # (a_j(n) S^j)(n^k S^i) = a_j(n) (n+j)^k S^(j+i)
    for j in range(4):
        coeff[i+j] = l27[j] * pow((n0+j) % p, k, p) % p

    for m in range(i+3, 2, -1):
        shift = m - 3
        z = lz_at(n0 + shift, p)
        lead = z[3]
        if lead == 0:
            return None
        q = coeff[m] * pow(lead, p-2, p) % p
        if q:
            for j in range(4):
                coeff[shift+j] = (coeff[shift+j] - q*z[j]) % p
    return coeff[0] % p, coeff[1] % p, coeff[2] % p


def rank_mod(rows: Sequence[Sequence[int]], p: int) -> int:
    if not rows:
        return 0
    a = [list(map(lambda x: x % p, row)) for row in rows]
    nr, nc = len(a), len(a[0])
    rank = 0
    for col in range(nc):
        pivot = None
        for r in range(rank, nr):
            if a[r][col] % p:
                pivot = r
                break
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], p-2, p)
        a[rank] = [(x*inv) % p for x in a[rank]]
        for r in range(nr):
            if r != rank and a[r][col]:
                c = a[r][col]
                a[r] = [(a[r][j] - c*a[rank][j]) % p for j in range(nc)]
        rank += 1
        if rank == nc:
            return rank
    return rank


def matrix_for(d: int, D: int, p: int, extra_blocks: int = 10) -> Tuple[List[List[int]], List[int]]:
    cols = [(i,k) for i in range(d+1) for k in range(D+1)]
    need_blocks = (len(cols) + 2)//3 + extra_blocks
    rows: List[List[int]] = []
    ns: List[int] = []
    n0 = 0
    while len(ns) < need_blocks and n0 < 20*need_blocks + 100:
        rems = []
        valid = True
        for i,k in cols:
            r = basis_remainder(i,k,n0,p)
            if r is None:
                valid = False
                break
            rems.append(r)
        if valid:
            for comp in range(3):
                rows.append([r[comp] for r in rems])
            ns.append(n0)
        n0 += 1
    if len(ns) < need_blocks:
        raise RuntimeError(f"not enough safe evaluation points mod {p}")
    return rows, ns


def screen(max_d: int = 5, max_D: int = 30) -> None:
    primes = [1000003, 1000033, 1000037]
    print("Q5025 polynomial intertwiner modular rank screen")
    print("Convention: L27*R = Q*LZ, right remainder by LZ")
    print("primes:", primes)
    print()
    survivors = []
    for d in range(max_d+1):
        print(f"=== order d={d} ===")
        for DD in range(max_D+1):
            u = (d+1)*(DD+1)
            ranks = []
            full_at = None
            for p in primes:
                rows, ns = matrix_for(d, DD, p)
                rr = rank_mod(rows, p)
                ranks.append(rr)
                if rr == u:
                    full_at = p
                    break
            if full_at is not None:
                print(f"d={d:2d} D={DD:2d} unknowns={u:3d}: FULL rank mod {full_at}")
            else:
                print(f"d={d:2d} D={DD:2d} unknowns={u:3d}: DEFICIENT ranks={ranks}")
                survivors.append((d,DD,u,ranks))
        print()
    print("SURVIVORS", survivors)


if __name__ == "__main__":
    screen(max_d=5, max_D=30)
