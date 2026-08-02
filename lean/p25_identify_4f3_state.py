#!/usr/bin/env python3
"""Identify the CMF product with the terminating 4F3 theta state."""
from fractions import Fraction as Q


def raw(n):
    n = Q(n)
    return [
        [(136*n**5+1672*n**4+8134*n**3+19808*n**2+24469*n+12456)/((n+4)**2*(2*n+7)**2*(2*n+9)),
         (n+1)*(2*n+3)*(96*n**4+972*n**3+3596*n**2+5851*n+3600)/(2*(n+4)**2*(2*n+7)**2*(2*n+9)),
         (n+1)*(2*n+3)*(128*n**5+1336*n**4+5404*n**3+11194*n**2+12973*n+7200)/(4*(n+4)**2*(2*n+7)**2*(2*n+9))],
        [(384*n**6+6096*n**5+40680*n**4+146720*n**3+302710*n**2+339941*n+162969)/(2*(n+4)**4*(2*n+7)**2*(2*n+9)),
         (n+1)*(2*n+3)*(272*n**5+3720*n**4+20484*n**3+57010*n**2+80258*n+45681)/(4*(n+4)**4*(2*n+7)**2*(2*n+9)),
         (n+1)*(2*n+3)*(384*n**6+5552*n**5+33880*n**4+113572*n**3+224502*n**2+252676*n+128259)/(8*(n+4)**4*(2*n+7)**2*(2*n+9))],
        [3*(96*n**4+1036*n**3+4136*n**2+7311*n+4871)/(2*(n+4)**4*(2*n+7)**2*(2*n+9)),
         3*(n+1)*(2*n+3)*(64*n**3+540*n**2+1460*n+1271)/(4*(n+4)**4*(2*n+7)**2*(2*n+9)),
         (n+1)*(2*n+3)*(32*n**5+816*n**4+6080*n**3+19656*n**2+29798*n+18063)/(8*(n+4)**4*(2*n+7)**2*(2*n+9))],
    ]


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mv(a, v):
    return [sum(a[i][j]*v[j] for j in range(3)) for i in range(3)]


def inv(a):
    # Gauss-Jordan, tiny exact matrices.
    x = [row[:] + [Q(i == j) for j in range(3)] for i, row in enumerate(a)]
    for j in range(3):
        p = next(i for i in range(j, 3) if x[i][j])
        x[j], x[p] = x[p], x[j]
        q = x[j][j]
        x[j] = [z/q for z in x[j]]
        for i in range(3):
            if i != j:
                q = x[i][j]
                x[i] = [x[i][k]-q*x[j][k] for k in range(6)]
    return [row[3:] for row in x]


def tr(a): return [list(x) for x in zip(*a)]


def poch(a, k):
    out = Q(1)
    for j in range(k): out *= a+j
    return out


def state(N):
    # F_N = 4F3(-N,-N-1/2,N+4,N+4;7/2,7/2,4;z).
    out = [Q(0), Q(0), Q(0)]
    for k in range(N+1):
        t = poch(Q(-N), k)*poch(Q(-N)-Q(1,2), k)*poch(Q(N+4), k)**2
        t /= poch(Q(7,2), k)**2*poch(Q(4), k)*poch(Q(1), k)
        for r in range(3): out[r] += Q(k**r)*t
    return out


I = [[Q(i == j) for j in range(3)] for i in range(3)]
D = I
for N in range(6):
    s = state(N)
    candidates = {
        "D e0": mv(D, [Q(1), Q(0), Q(0)]),
        "D^-1 e0": mv(inv(D), [Q(1), Q(0), Q(0)]),
        "D^T e0": mv(tr(D), [Q(1), Q(0), Q(0)]),
        "D^-T e0": mv(tr(inv(D)), [Q(1), Q(0), Q(0)]),
    }
    print("N", N, "state", s)
    for name, value in candidates.items():
        if value == s:
            print(" MATCH", name)
        else:
            # Try projective equality as conventions may rescale.
            scale = next((value[i]/s[i] for i in range(3) if s[i]), None)
            if scale is not None and all(value[i] == scale*s[i] for i in range(3)):
                print(" PROJECTIVE", name, scale)
    D = mm(D, raw(N))
