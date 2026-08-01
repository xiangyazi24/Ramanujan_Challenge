#!/usr/bin/env python3
"""Q6360 verification: v_p = tau_{(p-1)/4} == [Q^p] F(Q) mod p, F = Q*D(Q^4)*prod[...]^{1/2}.

In-sample: the 12 supplied primes. OUT-OF-SAMPLE (decisive): 769, 937, 1009, 1033, 1153.
[Q^p]F = F_J where D(q)H(q) = q^{1/4} sum F_n q^n, J=(p-1)/4.
H(q) = q^{1/4} * P(q)^{1/2}, P = prod (1-q^n)(1-q^2n)(1-q^3n)(1-q^6n).
D(q) = (1/2)[E2(q) - 2 E2(q^2) - 3 E2(q^3) + 6 E2(q^6)].
tau_n: sqrt of Apery generating series mod p (tau_0=1), v_p = tau_{(p-1)/4} mod p.
"""
from fractions import Fraction
import sys

PRIMES_IN = [73, 97, 193, 241, 313, 337, 409, 433, 457, 577, 601, 673]
V_SUPPLIED = {73:19, 97:1, 193:187, 241:128, 313:200, 337:175, 409:71, 433:432,
              457:338, 577:386, 601:30, 673:432}
PRIMES_OUT = [769, 937, 1009, 1033, 1153]
NMAX = max(PRIMES_IN + PRIMES_OUT) // 4 + 2   # need F_n for n <= (p-1)/4

def eta_quotient_P(N):
    """P = prod_{n>=1} (1-q^n)(1-q^2n)(1-q^3n)(1-q^6n), integer series to q^N."""
    P = [0]*(N+1); P[0] = 1
    for scale in (1, 2, 3, 6):
        for n in range(1, N//scale + 1):
            e = scale*n
            # multiply by (1 - q^e)
            for i in range(N, e-1, -1):
                P[i] -= P[i-e]
    return P

def sqrt_series(P, N):
    """S with S^2 = P, S[0]=1, rational coefficients."""
    S = [Fraction(0)]*(N+1); S[0] = Fraction(1)
    for n in range(1, N+1):
        acc = Fraction(P[n])
        for k in range(1, n):
            acc -= S[k]*S[n-k]
        S[n] = acc / 2
    return S

def E2(N, scale):
    """E2(q^scale) = 1 - 24 sum sigma1(n) q^{scale*n}, to q^N."""
    from sympy import divisor_sigma
    E = [Fraction(0)]*(N+1); E[0] = Fraction(1)
    for n in range(1, N//scale + 1):
        E[scale*n] = Fraction(-24*int(divisor_sigma(n, 1)))
    return E

def mul(A, B, N):
    C = [Fraction(0)]*(N+1)
    for i in range(N+1):
        if A[i] == 0: continue
        for j in range(N+1-i):
            if B[j]: C[i+j] += A[i]*B[j]
    return C

def tau_quarter_mod(p):
    """tau_{(p-1)/4} mod p via sqrt of Apery series mod p."""
    J = (p-1)//4
    b = [1 % p, 5 % p]
    for n in range(1, J+1):
        num = ((34*n**3+51*n**2+27*n+5)*b[n] - n**3*b[n-1]) % p
        b.append(num * pow(pow(n+1, 3, p), p-2, p) % p)
    inv2 = pow(2, p-2, p)
    t = [1] + [0]*J
    for n in range(1, J+1):
        acc = b[n]
        for k in range(1, n):
            acc -= t[k]*t[n-k]
        t[n] = acc % p * inv2 % p
    return t[J]

def main():
    N = NMAX
    print(f"building series to q^{N} ...", flush=True)
    P = eta_quotient_P(N)
    S = sqrt_series(P, N)             # H = q^{1/4} S(q)
    D = [Fraction(0)]*(N+1)
    for sc, coef in ((1, Fraction(1,2)), (2, Fraction(-1)), (3, Fraction(-3,2)), (6, Fraction(3))):
        E = E2(N, sc)
        for i in range(N+1): D[i] += coef*E[i]
    F = mul(D, S, N)                  # D*H = q^{1/4} sum F_n q^n
    print("series ready; testing", flush=True)
    ok_in = ok_out = 0
    for p in PRIMES_IN + PRIMES_OUT:
        J = (p-1)//4
        f = F[J]
        num, den = f.numerator, f.denominator
        if den % p == 0:
            print(f"p={p}: DENOMINATOR HIT (den divisible by p) — cannot reduce"); continue
        lhs = num % p * pow(den % p, p-2, p) % p
        rhs = tau_quarter_mod(p)
        tag = 'in ' if p in PRIMES_IN else 'OUT'
        match = lhs == rhs
        extra = ''
        if p in V_SUPPLIED:
            extra = f" supplied={V_SUPPLIED[p]}{'✓' if rhs==V_SUPPLIED[p] else '✗ SUPPLIED MISMATCH'}"
        print(f"[{tag}] p={p}: F_J mod p = {lhs}  tau_J mod p = {rhs}  {'MATCH' if match else 'MISMATCH'}{extra}")
        if match:
            if p in PRIMES_IN: ok_in += 1
            else: ok_out += 1
    print(f"\nin-sample {ok_in}/{len(PRIMES_IN)}  OUT-OF-SAMPLE {ok_out}/{len(PRIMES_OUT)}")

main()
