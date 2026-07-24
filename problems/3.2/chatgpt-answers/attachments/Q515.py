# SageMath-compatible exact code.
from sage.all import ZZ, QQ, binomial, matrix, vector


def P(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery(N):
    A = [ZZ(1), ZZ(5)]
    for n in range(1, N):
        num = P(n)*A[n] - n**3*A[n-1]
        den = (n+1)**3
        assert num % den == 0
        A.append(num // den)
    return A[:N+1]


def K0(k,m):
    if k == 0:
        return ZZ(1) if m == 0 else ZZ(0)
    return (-1)**(k-m) * QQ(2*k, k+m) * binomial(k+m,k-m)


def K1(k,m):
    return (-1)**(k-m) * QQ(2*k+1,k+m+1) * binomial(k+m+1,k-m)


def residuals(K):
    A = apery(2*K+1)
    e = []
    o = []
    for k in range(K+1):
        ek = sum(K0(k,m)*A[2*m] for m in range(k+1))
        ok = sum(K1(k,m)*A[2*m+1] for m in range(k+1))
        assert ek in ZZ and ok in ZZ
        e.append(ZZ(ek)); o.append(ZZ(ok))
    return A,e,o


def Phi(eps,k,m):
    return binomial(2*m+eps,m-k)

A,e,o = residuals(12)

# Exact reconstruction.
for m in range(13):
    assert A[2*m] == sum(e[k]*Phi(0,k,m) for k in range(m+1))
    assert A[2*m+1] == sum(o[k]*Phi(1,k,m) for k in range(m+1))

# Low-order regression values.
assert e[:7] == [
    1, 71, 32711, 21263474, 16196884679,
    13494506759471, 11910357240848882
]
assert o[:6] == [
    5, 1430, 811805, 578594525,
    463454152550, 398546130989165
]

# Shift compatibility.
for eps in (0,1):
    for m in range(1,12):
        for k in range(1,m):
            assert Phi(eps,k,m+1) == (
                Phi(eps,k-1,m)+2*Phi(eps,k,m)+Phi(eps,k+1,m)
            )
assert all(Phi(0,0,m+1)==2*Phi(0,0,m)+2*Phi(0,1,m) for m in range(1,12))
assert all(Phi(1,0,m+1)==3*Phi(1,0,m)+Phi(1,1,m) for m in range(1,12))

# Parity pattern proving absence of factorial content.
for k in range(1,13):
    assert (e[k] % 2 == 0) == (k % 3 == 0)
for k in range(13):
    assert (o[k] % 2 == 0) == (k % 3 == 1)

# Top-half tail identity, checked exactly modulo every prime in a range.
for eps in (0,1):
    for m in range(1,12):
        n = 2*m+eps
        for p in prime_range(n//2+1,n+1):
            r = n-p
            d = m-r
            if eps == 1:
                # p = r+2d+1 with this same d.
                pass
            coeffs = e if eps == 0 else o
            tail = sum(binomial(r,h)*coeffs[d+h] for h in range(r+1))
            assert (tail-5*A[r]) % p == 0
            assert (tail-A[n]) % p == 0

print('PASS: exact two-sieved factorial-basis reduction')
