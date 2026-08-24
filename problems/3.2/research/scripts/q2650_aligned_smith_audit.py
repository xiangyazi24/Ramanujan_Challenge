#!/usr/bin/env python3
"""Dependency-free exact p-adic Smith audit for Q2650.

Computes the first singular aligned Apéry frame p=31, H=93, m=94,
j=32 (d=j-1=31), under the primitive Newton normalization used by
q32_pade_total_positivity.py / Q2638.

The certificate is valuation-theoretic, not a raw-rank test:
  * T_(m,j) = tau_(H,d) is computed exactly by Bareiss;
  * v_p Delta_(H,d) is the minimum v_p among *all* maximal minors of
    R_(H,d);
  * likewise for Delta_(H,d+1);
  * v_p K_(m,j) follows from the primitive Wronskian formula.

No third-party packages are used.
"""

from math import comb

P = 31
H = 93
M = H + 1
J = 32
D = J - 1


def apery_values(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    b = [1, 5]
    for n in range(1, limit):
        num = (34*n**3 + 51*n**2 + 27*n + 5) * b[n] - n**3 * b[n-1]
        den = (n + 1)**3
        assert num % den == 0
        b.append(num // den)
    return b[:limit+1]


def forward_differences(values: list[int]) -> list[int]:
    row = values[:]
    out = []
    while row:
        out.append(row[0])
        row = [row[i+1] - row[i] for i in range(len(row)-1)]
    return out


def M_entry(k: int, ell: int, c: list[int]) -> int:
    if ell < 0 or ell > k:
        return 0
    return comb(k, ell) * sum(
        comb(ell, t) * c[k-ell+t]
        for t in range(ell+1)
    )


def bareiss_det(A: list[list[int]]) -> int:
    n = len(A)
    assert all(len(row) == n for row in A)
    if n == 0:
        return 1
    if n == 1:
        return A[0][0]
    B = [row[:] for row in A]
    sign = 1
    prev = 1
    for k in range(n-1):
        if B[k][k] == 0:
            swap = next((i for i in range(k+1, n) if B[i][k] != 0), None)
            if swap is None:
                return 0
            B[k], B[swap] = B[swap], B[k]
            sign = -sign
        pivot = B[k][k]
        for i in range(k+1, n):
            bik = B[i][k]
            for j in range(k+1, n):
                num = B[i][j] * pivot - bik * B[k][j]
                assert num % prev == 0
                B[i][j] = num // prev
            B[i][k] = 0
        prev = pivot
    return sign * B[-1][-1]


def vp(n: int, p: int = P):
    if n == 0:
        return None
    n = abs(n)
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def unit_residue(n: int, e: int, p: int = P) -> int:
    assert n != 0
    return (n // (p**e)) % p


def frame(H: int, d: int, c: list[int]) -> list[list[int]]:
    # R_(H,d) = M[H-d+1,...,H ; 0,...,d]
    return [
        [M_entry(k, ell, c) for ell in range(d+1)]
        for k in range(H-d+1, H+1)
    ]


def maximal_minor_certificate(A: list[list[int]], p: int = P):
    rows = len(A)
    cols = len(A[0])
    assert cols == rows + 1
    vals = []
    residues = []
    determinants = []
    for deleted in range(cols):
        minor = [row[:deleted] + row[deleted+1:] for row in A]
        det = bareiss_det(minor)
        determinants.append(det)
        e = vp(det, p)
        vals.append(e)
        residues.append(None if e is None else unit_residue(det, e, p))
    finite = [e for e in vals if e is not None]
    assert finite
    minimum = min(finite)
    witnesses = [i for i, e in enumerate(vals) if e == minimum]
    return minimum, witnesses, vals, residues, determinants


def tau(H: int, d: int, c: list[int]) -> int:
    return bareiss_det([
        [M_entry(k, ell, c) for ell in range(d+1)]
        for k in range(H-d, H+1)
    ])


def fmt(vals):
    return '[' + ', '.join('inf' if x is None else str(x) for x in vals) + ']'


def main() -> None:
    b = apery_values(H)
    c = forward_differences(b)

    R31 = frame(H, D, c)
    R32 = frame(H, D+1, c)

    vD31, w31, vals31, res31, dets31 = maximal_minor_certificate(R31)
    vD32, w32, vals32, res32, dets32 = maximal_minor_certificate(R32)

    T = tau(H, D, c)
    vT = vp(T)
    assert vT is not None

    # T_(94,32) is exactly the d=32 maximal minor deleting column 32.
    assert T == dets32[D+1]
    assert vals32[D+1] == vT

    vbin = vp(comb(M, J))
    assert vbin is not None
    vU = vT - vD31
    vV = vT - vD32
    vK = vbin + vU + vV

    print('Q2650 exact aligned Smith certificate')
    print(f'p={P} H={H} m={M} j={J} d={D}')
    print(f'v31(binomial({M},{J}))={vbin}')
    print(f'v31(T_({M},{J}))={vT}')
    print(f'T normalized unit residue mod31={unit_residue(T, vT)}')
    print(f'v31(Delta_({H},{D}))={vD31}')
    print(f'Delta_({H},{D}) witness deleted columns={w31}')
    print(f'Delta_({H},{D}) maximal-minor valuations={fmt(vals31)}')
    print('Delta_({},{}) witness unit residues={}'.format(
        H, D, {i: res31[i] for i in w31}))
    print(f'v31(Delta_({H},{D+1}))={vD32}')
    print(f'Delta_({H},{D+1}) witness deleted columns={w32}')
    print(f'Delta_({H},{D+1}) maximal-minor valuations={fmt(vals32)}')
    print('Delta_({},{}) witness unit residues={}'.format(
        H, D+1, {i: res32[i] for i in w32}))
    print(f'check T is delete-column-{D+1} minor in Delta_({H},{D+1}) frame: yes')
    print(f'v31(U_({M},{J}))=vT-vDelta31={vU}')
    print(f'v31(V_({M},{J}))=vT-vDelta32={vV}')
    print(f'v31(K_({M},{J}))={vK}')
    print(f'K_({M},{J}) is a 31-unit: {vK == 0}')

    # Strong sanity checks: determinantal-divisor valuation cannot exceed any
    # maximal-minor valuation, and primitive boundary valuations are nonnegative.
    assert all(e is None or e >= vD31 for e in vals31)
    assert all(e is None or e >= vD32 for e in vals32)
    assert vU >= 0 and vV >= 0


if __name__ == '__main__':
    main()
