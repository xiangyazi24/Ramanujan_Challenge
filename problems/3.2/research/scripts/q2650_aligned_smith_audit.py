#!/usr/bin/env python3
"""Q2650: exact p-adic Smith/minimal-maximal-minor audit.

Dependency-free.  No raw-rank inference is used for the certificate.
For an integer matrix A and p=31, `smith_exponents_mod` performs Smith
reduction over the DVR Z_p, represented modulo p^PREC.  At each stage it
extracts the least p-power from the whole residual block, pivots on a unit,
and splits off one invariant factor.  The accumulated pivot exponents are
the p-adic Smith exponents.  Their sum is the valuation of the last
determinantal divisor for a full-row-rank rectangular matrix, and the
valuation of the determinant for a square matrix.

For R_(H,d), we additionally compute the Smith valuation of every maximal
minor obtained by deleting one column.  Thus Delta_(H,d) is certified by
an explicit vector of maximal-minor valuations and its minimum.
"""

from math import comb

P = 31
PREC = 20
H = 93
M = 94


def apery_values(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    b = [1, 5]
    for n in range(1, limit):
        num = (34*n**3 + 51*n**2 + 27*n + 5)*b[n] - n**3*b[n-1]
        den = (n+1)**3
        assert num % den == 0
        b.append(num // den)
    return b[:limit+1]


def forward_differences(values: list[int]) -> list[int]:
    row = values[:]
    out = []
    while row:
        out.append(row[0])
        row = [row[i+1]-row[i] for i in range(len(row)-1)]
    return out


def M_entry(k: int, ell: int, c: list[int]) -> int:
    if ell < 0 or ell > k:
        return 0
    return comb(k, ell)*sum(
        comb(ell, t)*c[k-ell+t]
        for t in range(ell+1)
    )


def vp_int(n: int, p: int = P):
    if n == 0:
        return None
    n = abs(n)
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def vp_residue(x: int, p: int, precision: int) -> int:
    """Valuation of a residue mod p^precision, capped at precision."""
    modulus = p**precision
    x %= modulus
    if x == 0:
        return precision
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def smith_exponents_mod(A, p: int = P, precision: int = PREC):
    """p-adic Smith exponents, exact provided precision is not exhausted.

    Row/column operations use only p-adic units.  If the residual block has
    minimum valuation e, factoring p^e from that entire block shifts every
    remaining invariant factor by e; `offset` records this shift.
    """
    if not A:
        return []
    rows = len(A)
    cols = len(A[0])
    assert all(len(row) == cols for row in A)
    modulus = p**precision
    B = [[x % modulus for x in row] for row in A]
    offset = 0
    exponents = []
    current_precision = precision

    while B and B[0]:
        r = len(B)
        c = len(B[0])
        # Full row rank is expected for the rectangular frames; square minors
        # may be singular over Q only in principle.  Precision exhaustion is
        # treated as a hard failure, never as rank information.
        best = current_precision
        bi = bj = None
        for i in range(r):
            for j in range(c):
                e = vp_residue(B[i][j], p, current_precision)
                if e < best:
                    best, bi, bj = e, i, j
                    if best == 0:
                        break
            if best == 0:
                break
        if bi is None:
            raise RuntimeError(
                f'precision exhausted: residual {r}x{c} block is 0 mod p^{current_precision}'
            )

        if best:
            pe = p**best
            new_precision = current_precision - best
            if new_precision <= 0:
                raise RuntimeError('precision exhausted while extracting Smith factor')
            new_modulus = p**new_precision
            for i in range(r):
                for j in range(c):
                    # Representatives are divisible by p^best because all
                    # entries have valuation at least best.
                    assert B[i][j] % pe == 0
                    B[i][j] = (B[i][j] // pe) % new_modulus
            offset += best
            current_precision = new_precision
            modulus = new_modulus

        # Move a unit to (0,0).
        if bi != 0:
            B[0], B[bi] = B[bi], B[0]
        if bj != 0:
            for row in B:
                row[0], row[bj] = row[bj], row[0]
        pivot = B[0][0] % modulus
        assert pivot % p != 0
        inv = pow(pivot, -1, modulus)

        # Clear first column by p-adic-unimodular row additions.
        for i in range(1, r):
            factor = (B[i][0] * inv) % modulus
            if factor:
                B[i] = [
                    (B[i][j] - factor*B[0][j]) % modulus
                    for j in range(c)
                ]
        # Clear first row by column additions.  Column 0 is now zero below
        # the pivot, so this does not alter the residual block rows.
        for j in range(1, c):
            factor = (B[0][j] * inv) % modulus
            if factor:
                for i in range(r):
                    B[i][j] = (B[i][j] - factor*B[i][0]) % modulus

        exponents.append(offset)
        B = [row[1:] for row in B[1:]]

    return exponents


def frame(H: int, d: int, c):
    return [
        [M_entry(k, ell, c) for ell in range(d+1)]
        for k in range(H-d+1, H+1)
    ]


def solid_T(H: int, d: int, c):
    return [
        [M_entry(k, ell, c) for ell in range(d+1)]
        for k in range(H-d, H+1)
    ]


def maximal_minor_valuations(A):
    r = len(A)
    c = len(A[0])
    assert c == r+1
    vals = []
    exponent_lists = []
    for deleted in range(c):
        minor = [row[:deleted] + row[deleted+1:] for row in A]
        exps = smith_exponents_mod(minor)
        assert len(exps) == r
        vals.append(sum(exps))
        exponent_lists.append(exps)
    minimum = min(vals)
    witnesses = [i for i, v in enumerate(vals) if v == minimum]
    return minimum, witnesses, vals, exponent_lists


def edge_certificate(j: int, c):
    d = j-1
    Aleft = frame(H, d, c)
    Aright = frame(H, d+1, c)
    Tmat = solid_T(H, d, c)

    expsT = smith_exponents_mod(Tmat)
    assert len(expsT) == j
    vT = sum(expsT)

    vDl, wl, valsl, expl = maximal_minor_valuations(Aleft)
    vDr, wr, valsr, expr = maximal_minor_valuations(Aright)

    # T_(m,j) is literally the right-frame maximal minor deleting column j.
    assert valsr[j] == vT

    vbin = vp_int(comb(M, j))
    assert vbin is not None
    vU = vT-vDl
    vV = vT-vDr
    vK = vbin+vU+vV
    assert vU >= 0 and vV >= 0

    print(f'EDGE j={j} d={d}')
    print(f'  v31(binomial({M},{j}))={vbin}')
    print(f'  T Smith exponents={expsT}')
    print(f'  v31(T_({M},{j}))={vT}')
    print(f'  v31(Delta_({H},{d}))={vDl}')
    print(f'  Delta_({H},{d}) witness deleted columns={wl}')
    print(f'  Delta_({H},{d}) maximal-minor valuations={valsl}')
    for w in wl:
        print(f'    witness left delete {w} Smith exponents={expl[w]}')
    print(f'  v31(Delta_({H},{d+1}))={vDr}')
    print(f'  Delta_({H},{d+1}) witness deleted columns={wr}')
    print(f'  Delta_({H},{d+1}) maximal-minor valuations={valsr}')
    for w in wr:
        print(f'    witness right delete {w} Smith exponents={expr[w]}')
    print(f'  T is right-frame delete-column-{j} minor: yes')
    print(f'  v31(U_({M},{j}))={vU}')
    print(f'  v31(V_({M},{j}))={vV}')
    print(f'  v31(K_({M},{j}))={vK}')
    print(f'  K_({M},{j}) is a 31-unit: {vK == 0}')
    return vK


def self_test():
    # Known Smith examples, including a case where raw mod-p rank sees only
    # one defect but its invariant exponent is 2.
    assert smith_exponents_mod([[1, 0], [0, P]]) == [0, 1]
    assert smith_exponents_mod([[1, 0], [0, P**2]]) == [0, 2]
    assert smith_exponents_mod([[P, 0], [0, P**2]]) == [1, 2]
    assert smith_exponents_mod([[P, P], [0, P**2]]) == [1, 2]


def main():
    self_test()
    b = apery_values(H)
    c = forward_differences(b)
    print('Q2650 exact p-adic Smith/minimal-maximal-minor certificate')
    print(f'p={P} H={H} m={M} precision={PREC}')
    # Sanity data for the first singular base frame.
    zeros = [r for r in range(P) if b[r] % P == 0]
    print(f'base Apéry zeros mod31={zeros}')
    print(f'v31(b_8)={vp_int(b[8])}, v31(b_22)={vp_int(b[22])}')
    print(f'(b_8/31) mod31={(b[8]//P)%P}, (b_22/31) mod31={(b[22]//P)%P}')

    vK32 = edge_certificate(32, c)
    if vK32 != 0:
        print('j=32 is not a unit; computing requested aligned fallback j=63')
        edge_certificate(63, c)


if __name__ == '__main__':
    main()
