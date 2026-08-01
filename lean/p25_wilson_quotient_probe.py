#!/usr/bin/env python3
"""Search for a polynomial 3-to-2 quotient gauge from P2.5 to Wilson--Pade."""

P = 2**61 - 1


def inv(x):
    return pow(x % P, P - 2, P)


def positive_matrix(n):
    return [
        [(2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         480*n**4+4980*n**3+19210*n**2+32690*n+20730],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240)],
    ]


def normalized_matrix(n):
    d = 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    di = inv(d)
    return [[x % P * di % P for x in row] for row in positive_matrix(n)]


def wilson_matrix(m):
    den0 = 8*(m+1)**2*(4*m+3)**2
    den1 = den0*(4*m+5)**2
    return [
        [((4*m+1)**2*(40*m*m+68*m+29) % P)*inv(den0)%P,
         ((1536*m**4+4224*m**3+4288*m**2+1904*m+313)%P)*inv(den0)%P],
        [((4*m+1)**2*(1536*m**4+6144*m**3+9152*m**2+6016*m+1473)%P)*inv(den1)%P,
         ((59392*m**6+300032*m**5+620544*m**4+672000*m**3+401640*m**2+125612*m+16077)%P)*inv(den1)%P],
    ]


def rank_mod(rows, ncols):
    rows = [[x % P for x in row] for row in rows]
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = inv(rows[rank][col])
        rows[rank] = [x*scale % P for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                scale = rows[r][col]
                rows[r] = [(rows[r][c] - scale*rows[rank][c]) % P
                           for c in range(ncols)]
        rank += 1
        if rank == len(rows):
            break
    return rank


def polynomial_gauge_nullity(degree, shift=2, holdout=5):
    # Unknown order: s[i,j,r] for a 3x2 matrix, r=0..degree.
    width = degree + 1
    ncols = 6*width
    rows = []
    samples = degree + 2 + holdout
    for n in range(samples):
        A = normalized_matrix(n)
        W = wilson_matrix(n+shift)
        # C = W^T, and A(n) S(n+1) = S(n) C(n).
        npows = [pow(n, r, P) for r in range(width)]
        spows = [pow(n+1, r, P) for r in range(width)]
        for i in range(3):
            for j in range(2):
                row = [0]*ncols
                for k in range(3):
                    for r in range(width):
                        idx = (k*2+j)*width+r
                        row[idx] = (row[idx] + A[i][k]*spows[r]) % P
                for ell in range(2):
                    # C[ell,j] = W[j,ell].
                    c = W[j][ell]
                    for r in range(width):
                        idx = (i*2+ell)*width+r
                        row[idx] = (row[idx] - c*npows[r]) % P
                rows.append(row)
    return ncols-rank_mod(rows, ncols)


for shift in range(0, 6):
    print(f"shift={shift}", flush=True)
    for degree in range(0, 26):
        nullity = polynomial_gauge_nullity(degree, shift)
        if nullity:
            print(f"  degree={degree}: nullity={nullity}", flush=True)
            break
        if degree in (0, 5, 10, 15, 20, 25):
            print(f"  degree={degree}: 0", flush=True)
