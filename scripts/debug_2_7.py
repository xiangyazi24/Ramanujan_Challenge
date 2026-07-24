#!/usr/bin/env python3
"""Debug Problem 2.7: check recurrence indexing and initial conditions."""
from mpmath import mp, mpf, zeta, log10, fabs

mp.dps = 100

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    n = mpf(n)
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Initial conditions from PDF
p0 = mpf(-612218384750)
p1 = mpf(-9525021973931919) / mpf(18100)
p2 = mpf(-29561828382772029) / mpf(65380)

q0 = mpf(-215040420000)
q1 = mpf(-1672822650043404) / mpf(905)
q2 = mpf(-964185327658080) / mpf(6071)

target = zeta(2) + zeta(3)

# Test DIFFERENT recurrence interpretations

# Version 1: As written in PDF
# u_{n+1} = (B_n/A_n)u_n - (C_{n-1}/A_{n-1})u_{n-1} + (D_{n-2}/A_{n-2})u_{n-2}, n>=2
def gen_v1(u0, u1, u2, N):
    u = [u0, u1, u2]
    for n in range(2, N):
        u_next = B(n)/A(n)*u[n] - C(n-1)/A(n-1)*u[n-1] + D(n-2)/A(n-2)*u[n-2]
        u.append(u_next)
    return u

# Version 2: Standard companion - all coefficients at same index
# A_n u_{n+1} - B_n u_n + C_n u_{n-1} - D_n u_{n-2} = 0
def gen_v2(u0, u1, u2, N):
    u = [u0, u1, u2]
    for n in range(2, N):
        u_next = (B(n)*u[n] - C(n)*u[n-1] + D(n)*u[n-2]) / A(n)
        u.append(u_next)
    return u

# Version 3: Different sign convention
# A_n u_{n+1} = B_n u_n + C_{n-1} u_{n-1} - D_{n-2} u_{n-2}
def gen_v3(u0, u1, u2, N):
    u = [u0, u1, u2]
    for n in range(2, N):
        u_next = (B(n)*u[n] + C(n-1)*u[n-1] - D(n-2)*u[n-2]) / A(n)
        u.append(u_next)
    return u

# Version 4: Maybe the PDF formula means something else:
# 0 = -A_n u_{n+1} + B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}
# i.e., same as v1 but let me double-check the recurrence form in the PDF

# The PDF says: "Let (p_n)_{n>=0} and (q_n)_{n>=0} be the two solutions of the recurrence"
# u_{n+1} = B_n/A_n u_n - C_{n-1}/A_{n-1} u_{n-1} + D_{n-2}/A_{n-2} u_{n-2}, n >= 2

# This is: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# Multiplying by A_n:
# A_n u_{n+1} = B_n u_n - (A_n C_{n-1}/A_{n-1}) u_{n-1} + (A_n D_{n-2}/A_{n-2}) u_{n-2}

# But this is NOT a standard 4-term recurrence with polynomial coefficients!
# Unless A_n/A_{n-1} and A_n/A_{n-2} are rational functions of n.

# Wait: A_n = 1024(2n+5)^4(2n+7)^3(2n+9)^3(946n^2+6407n+10860)
# A_{n-1} = 1024(2n+3)^4(2n+5)^3(2n+7)^3(946(n-1)^2+6407(n-1)+10860)
# So A_n/A_{n-1} = [(2n+5)^4/(2n+3)^4] * [(2n+9)^3/(2n+5)^3] * [Q(n)/Q(n-1)]
# = (2n+5)(2n+9)^3/[(2n+3)^4] * Q(n)/Q(n-1)

# This is rational. So the equivalent 4-term recurrence with polynomial coefficients:
# A_n A_{n-1} A_{n-2} u_{n+1} - B_n A_{n-1} A_{n-2} u_n + C_{n-1} A_n A_{n-2} u_{n-1} - D_{n-2} A_n A_{n-1} u_{n-2} = 0
# But this has huge degree (~36).

# Actually wait, I think the PDF is just using a shorthand where the recurrence is:
# A_n u_{n+1} = B_n u_n - C_n u_{n-1} + D_n u_{n-2}  (standard 4-term)
# But they WRITE it as u_{n+1} = (B/A)u_n - (C/A)u_{n-1} + (D/A)u_{n-2}
# and shift n appropriately so the second term has C_{n-1}/A_{n-1} etc.

# Hmm, looking at the PDF again: it literally says
# u_{n+1} = B_n/A_n u_n - C_{n-1}/A_{n-1} u_{n-1} + D_{n-2}/A_{n-2} u_{n-2}

# This is a specific choice of where to evaluate the coefficients.
# It's NOT the same as A_n u_{n+1} = B_n u_n - C_n u_{n-1} + D_n u_{n-2}

# Let me try the STANDARD 4-term form: homogeneous with all coefficients at same n
N = 60

for label, gen in [("V1: PDF literal", gen_v1), ("V2: standard companion", gen_v2), ("V3: + C - D", gen_v3)]:
    p = gen(p0, p1, p2, N)
    q = gen(q0, q1, q2, N)
    ratio = p[-1]/q[-1]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
    else:
        digits = mpf('inf')
    print(f"{label}: p/q = {mp.nstr(ratio, 20)}, digits = {mp.nstr(digits, 5)}")

# Also try: maybe the indices in the PDF notation mean the recurrence should be
# evaluated at n -> n+1 (shifted)?
# I.e., for n>=2: u_{n+1} = (B_n/A_n)u_n - (C_{n-1}/A_{n-1})u_{n-1} + (D_{n-2}/A_{n-2})u_{n-2}
# But maybe the STANDARD 4-term recurrence is:
# A(n) u(n+3) - B(n) u(n+2) + C(n) u(n+1) - D(n) u(n) = 0, n>=0
# Then with a shift m = n+2: A(m-2)u(m+1) - B(m-2)u(m) + C(m-2)u(m-1) - D(m-2)u(m-2) = 0
# => u(m+1) = (B(m-2)/A(m-2))u(m) - (C(m-2)/A(m-2))u(m-1) + (D(m-2)/A(m-2))u(m-2)
# This would mean: for m>=2, use coefficients at n=m-2.

print("\n--- Shifted versions ---")
def gen_shift(u0, u1, u2, N, shift):
    u = [u0, u1, u2]
    for n in range(2, N):
        ns = n + shift
        u_next = B(ns)/A(ns)*u[n] - C(ns)/A(ns)*u[n-1] + D(ns)/A(ns)*u[n-2]
        u.append(u_next)
    return u

for shift in [-2, -1, 0, 1, 2]:
    try:
        p_s = gen_shift(p0, p1, p2, N, shift)
        q_s = gen_shift(q0, q1, q2, N, shift)
        ratio = p_s[-1]/q_s[-1]
        err = fabs(ratio - target)
        if err > 0:
            digits = -log10(err)
        else:
            digits = mpf('inf')
        print(f"  shift={shift:+d}: p/q = {mp.nstr(ratio, 20)}, digits = {mp.nstr(digits, 5)}")
    except Exception as e:
        print(f"  shift={shift:+d}: ERROR {e}")
