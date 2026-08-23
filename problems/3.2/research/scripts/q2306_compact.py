#!/usr/bin/env python3
"""Compact dependency-free reproducer for Q2306."""
from math import comb

p, r, s = 181, 19, 47


def tr(a):
    a = [x % p for x in a] or [0]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return tr([(a[i] if i < len(a) else 0) +
               (b[i] if i < len(b) else 0)
               for i in range(max(len(a), len(b)))])


def sub(a, b):
    return tr([(a[i] if i < len(a) else 0) -
               (b[i] if i < len(b) else 0)
               for i in range(max(len(a), len(b)))])


def sc(a, c):
    return tr([c * x for x in a])


def mul(a, b):
    z = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            z[i + j] += x * y
    return tr(z)


def dm(a, b):
    a, b = tr(a), tr(b)
    q = [0] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        d = len(a) - len(b)
        c = a[-1] * pow(b[-1], -1, p) % p
        q[d] = c
        for j, y in enumerate(b):
            a[d + j] = (a[d + j] - c * y) % p
        a = tr(a)
    return tr(q), a


def monic(a):
    a = tr(a)
    return sc(a, pow(a[-1], -1, p))


def gcdp(a, b):
    while tr(b) != [0]:
        a, b = b, dm(a, b)[1]
    return monic(a)


def ev(a, x):
    z = 0
    for c in reversed(a):
        z = (z * x + c) % p
    return z


def resultant(f, g):
    f, g = tr(f), tr(g)
    m, n = len(f) - 1, len(g) - 1
    if n == 0:
        return pow(g[0], m, p)
    if m < n:
        return ((-1) ** (m * n) * resultant(g, f)) % p
    _, h = dm(f, g)
    if h == [0]:
        return 0
    k = len(h) - 1
    return ((-1) ** (m * n) * pow(g[-1], m - k, p) *
            resultant(g, h)) % p


def det(a):
    a = [[x % p for x in row] for row in a]
    z = 1
    for j in range(len(a)):
        i = next((i for i in range(j, len(a)) if a[i][j]), None)
        if i is None:
            return 0
        if i != j:
            a[i], a[j] = a[j], a[i]
            z = -z
        z = z * a[j][j] % p
        u = pow(a[j][j], -1, p)
        for i in range(j + 1, len(a)):
            c = a[i][j] * u % p
            for k in range(j, len(a)):
                a[i][k] = (a[i][k] - c * a[j][k]) % p
    return z % p


def sylvester(f, g):
    m, n = len(f) - 1, len(g) - 1
    fh, gh, N = f[::-1], g[::-1], m + n
    rows = []
    for d in range(n):
        row = [0] * N
        row[d:d + m + 1] = fh
        rows.append(row)
    for d in range(m):
        row = [0] * N
        row[d:d + n + 1] = gh
        rows.append(row)
    return det(rows)


def lam(n):
    return n * (n + 1) % p


def D(n, k):
    return comb(n, k) * comb(n + k, k) % p


q = [[1]]
for k in range(s):
    q.append(sc(mul(q[-1], [-lam(k), 1]), pow((k + 1) ** 2, -1, p)))


def A(n):
    z = [0]
    for k in range(n + 1):
        z = add(z, sc(q[k], D(n, k)))
    return z


def K(a, b):
    return sum(D(a, k) * D(b, k) for k in range(min(a, b) + 1)) % p


Ar, As = A(r), A(s)
assert ev(Ar, lam(r)) == ev(As, lam(s)) == 0
assert ev(Ar, lam(s)) == ev(As, lam(r)) == K(r, s) == 134
assert gcdp(Ar, As) == [1]
assert resultant(Ar, As) == sylvester(Ar, As) == 34

# Monic Racah recurrence and associated block.
def aa(n):
    return (n * n + n + 1) * pow(2, -1, p) % p


def bb(n):
    return pow(n, 6, p) * pow(4 * (4 * n * n - 1), -1, p) % p


P = [monic(A(n)) for n in range(s + 1)]
for n in range(s):
    rhs = mul([aa(n), 1], P[n])
    if n:
        rhs = sub(rhs, sc(P[n - 1], bb(n)))
    assert rhs == P[n + 1]

h = s - r
Sblk = [[1], [aa(r + 1), 1]]
while len(Sblk) <= h - 1:
    j = len(Sblk) - 1
    n = r + j + 1
    Sblk.append(sub(mul([aa(n), 1], Sblk[j]), sc(Sblk[j - 1], bb(n))))
Sblk = Sblk[h - 1]
assert dm(sub(P[s], mul(P[r + 1], Sblk)), P[r])[1] == [0]

adj = resultant(P[r], P[r + 1])
block = resultant(P[r], Sblk)
monic_res = resultant(P[r], P[s])
scale = pow(Ar[-1], s, p) * pow(As[-1], r, p) % p
assert (Ar[-1], As[-1], scale, adj, block, monic_res) == (56, 38, 119, 98, 19, 52)
assert adj * block % p == monic_res
assert scale * monic_res % p == 34

print({
    "lambda_r": lam(r), "lambda_s": lam(s),
    "b_r": ev(Ar, lam(r)), "b_s": ev(As, lam(s)),
    "K_rs": K(r, s), "gcd": gcdp(Ar, As),
    "Res(A_r,A_s)": resultant(Ar, As),
    "lc_r": Ar[-1], "lc_s": As[-1], "scale": scale,
    "adjacent": adj, "block": block, "monic": monic_res,
    "associated_block": Sblk,
})
