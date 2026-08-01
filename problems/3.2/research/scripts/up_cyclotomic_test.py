"""Test cron's value-separation: tau_{(p-1)/4} = 2A * U_p (p = A^2+B^2, 2A = C(e,e/2) mod p).
Sharp hypothesis test: is U_p a root of unity of order dividing 24 (resp. gcd(24,p-1)) in F_p*?
Also record multiplicative order of U_p exactly.
"""
from sympy import isprime

def branch_quarter(p):
    # b_n mod p via recurrence, then tau = sqrt(F) or sigma = sqrt(F/q) truncation coefficients
    e = (p-1)//2
    b = [1 % p, 5 % p]
    for n in range(1, e+2):
        nxt = ((34*n**3+51*n**2+27*n+5)*b[n] - n**3*b[n-1]) * pow((n+1)**3, p-2, p) % p
        b.append(nxt)
    chi = pow(-6 % p, e, p)
    chi = 1 if chi == 1 else -1
    # series to invert: F (chi=+1) or F/q (chi=-1)
    m = e//2 + 2
    if chi == 1:
        c = b[:m+1]
    else:
        # F/q = F * (1-34t+t^2)^{-1}: q^{-1} coefficients satisfy g_k = 34 g_{k-1} - g_{k-2}
        g = [1, 34 % p]
        for k in range(2, m+1):
            g.append((34*g[k-1] - g[k-2]) % p)
        c = [sum(b[i]*g[k-i] for i in range(k+1)) % p for k in range(m+1)]
    # sqrt via s_0=1, s_k = (c_k - sum_{i=1}^{k-1} s_i s_{k-i}) * inv(2) ... standard
    s = [1]
    inv2 = pow(2, p-2, p)
    for k in range(1, m+1):
        acc = sum(s[i]*s[k-i] for i in range(1, k)) % p
        s.append((c[k] - acc) * inv2 % p)
    idx = (p-1)//4 if chi == 1 else (p-3)//4
    return chi, s[idx] if idx <= m else None

def two_A(p):
    e = (p-1)//2
    # C(e, e//2) mod p
    num = 1
    for i in range(e//2):
        num = num * (e-i) % p * pow(i+1, p-2, p) % p
    return num

def mult_order(x, p):
    from sympy import factorint
    o = p-1
    for q, k in factorint(p-1).items():
        for _ in range(k):
            if pow(x, o//q, p) == 1: o //= q
            else: break
    return o

rows = []
for p in range(29, 4000):
    if not isprime(p) or p % 4 != 1: continue
    if p % 24 == 5: continue  # vanishing class
    chi, tq = branch_quarter(p)
    if tq is None or tq == 0: continue
    A2 = two_A(p)
    U = tq * pow(A2, p-2, p) % p
    o = mult_order(U, p)
    rows.append((p, p % 24, chi, U, o, o <= 24, (p-1) % o == 0))

n24 = sum(1 for r in rows if r[5])
print(f"tested {len(rows)} primes; U_p order<=24 in {n24} cases")
from collections import Counter
print("order histogram:", Counter(r[4] for r in rows).most_common(12))
print("by class mod 24:", {c: Counter(r[4] for r in rows if r[1]==c).most_common(5) for c in sorted(set(r[1] for r in rows))})
