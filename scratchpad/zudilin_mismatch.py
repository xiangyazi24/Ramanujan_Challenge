#!/usr/bin/env python3
"""
Analyze the mismatch between Cauchy-twisted Zudilin basis and P2.7 sequences.
Goal: find the correction operator T that maps one to the other.
"""
from fractions import Fraction as Q
from math import comb

N = 20

# Zudilin Theorem 2 recurrence (for the inner sum a_n)
def z3(n):
    n = Q(n)
    return 2*(946*n*n-731*n+153)*(2*n+1)*(n+1)**3

def z2(n):
    n = Q(n)
    return -2*(104060*n**6+127710*n**5+12788*n**4
               -34525*n**3-8482*n**2+3298*n+1071)

def z1(n):
    n = Q(n)
    return 2*n*(3784*n**5-1032*n**4-1925*n**3
                +853*n**2+328*n-184)

def z0(n):
    n = Q(n)
    return -(946*n*n+1161*n+368)*n*(n-1)**3

def zudilin(init):
    u = list(map(Q, init))
    for n in range(2, N):
        u.append(-(z2(n)*u[n] + z1(n)*u[n-1]
                   + z0(n)*u[n-2]) / z3(n))
    return u

a   = zudilin([1, 7, 163])
b2  = zudilin([0, Q(23,2), Q(2145,8)])
b3  = zudilin([0, Q(17,2), Q(3135,16)])

# c_j(1) coefficients
c = [Q(1)]
for n in range(N+5):
    cm1 = c[n-1] if n >= 1 else Q(0)
    cm2 = c[n-2] if n >= 2 else Q(0)
    c.append(((440*n+660)*c[n] + (160-16*n)*cm1
              + (2*n-43)*cm2) / (8*(n+1)))

def twist(u):
    return [sum(c[n-m]*u[m] for m in range(n+1)) / Q(64)**n
            for n in range(N)]

ta  = twist(a)
tb2 = twist(b2)
tb3 = twist(b3)

# P2.7 recurrence
def A(n):
    n=Q(n)
    return Q(1024)*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3 \
        *(946*n*n+6407*n+10860)

def B(n):
    n=Q(n)
    return Q(128)*(2*n+7)**3*(2*n+9)**3 * \
      (104060*n**6+1745370*n**5+12145238*n**4
       +44886481*n**3+92943995*n**2+102256019*n+46709052)

def C(n):
    n=Q(n)
    return Q(16)*(n+3)**4*(2*n+9)**3 * \
      (3784*n**5+57792*n**4+351019*n**3
       +1059230*n**2+1587211*n+944620)

def D(n):
    n=Q(n)
    return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

q = [Q(-215040420000),
     Q(-167282265043404,905),
     Q(-964185327658080,6071)]
for n in range(2, N):
    q.append(B(n)/A(n)*q[n] - C(n-1)/A(n-1)*q[n-1]
             + D(n-2)/A(n-2)*q[n-2])

# Fit twisted basis to match q[0:3]
# q_n = alpha * ta[n] + beta * tb2[n] + gamma * tb3[n]
# At n=0,1,2: 3x3 linear system
from sympy import Matrix, Rational as SR

M_mat = Matrix([[SR(ta[i].numerator, ta[i].denominator),
                 SR(tb2[i].numerator, tb2[i].denominator),
                 SR(tb3[i].numerator, tb3[i].denominator)]
                for i in range(3)])
y_vec = Matrix([SR(q[i].numerator, q[i].denominator) for i in range(3)])
coef = M_mat.LUsolve(y_vec)

alpha = Q(int(coef[0].p), int(coef[0].q))
beta  = Q(int(coef[1].p), int(coef[1].q))
gamma = Q(int(coef[2].p), int(coef[2].q))

qstar = [alpha*ta[n] + beta*tb2[n] + gamma*tb3[n] for n in range(N)]

print("=== Mismatch analysis ===")
print(f"Basis coefficients: alpha={float(alpha):.6e}, beta={float(beta):.6e}, gamma={float(gamma):.6e}")
print()

delta = [qstar[n] - q[n] for n in range(N)]
print("delta[n] = qstar[n] - q[n]:")
for n in range(min(N, 12)):
    print(f"  n={n}: delta = {float(delta[n]):.15e}  (rel to q: {float(delta[n]/q[n]) if q[n] != 0 else 'inf':.15e})")

print()
print("Ratio delta[n]/delta[n-1]:")
for n in range(4, min(N, 12)):
    if delta[n-1] != 0:
        print(f"  n={n}: {float(delta[n]/delta[n-1]):.15f}")

print()
print("Growth: |delta[n]|^(1/n):")
for n in range(3, min(N, 12)):
    if delta[n] != 0:
        import math
        val = abs(float(delta[n]))
        if val > 0:
            print(f"  n={n}: {val**(1/n):.15f}")

# Check if delta satisfies the P2.7 recurrence
print()
print("Does delta satisfy P2.7 recurrence?")
for n in range(2, min(N-1, 10)):
    res = delta[n+1] - B(n)/A(n)*delta[n] + C(n-1)/A(n-1)*delta[n-1] - D(n-2)/A(n-2)*delta[n-2]
    print(f"  n={n}: residual = {float(res):.6e}")

# Check if delta satisfies Zudilin recurrence (twisted)
print()
print("Does delta satisfy Zudilin recurrence?")
for n in range(2, min(N-1, 8)):
    # qstar satisfies the twisted Zudilin recurrence
    # q satisfies P2.7
    # delta = qstar - q satisfies neither individually
    pass

# Try to find the recurrence that delta satisfies
# Since qstar satisfies the twisted-Zudilin rec and q satisfies P2.7,
# delta satisfies the LCLM of both recurrences. The LCLM has order 6.
# But delta is a SPECIFIC solution, so it might satisfy a lower-order recurrence.
print()
print("Attempting to guess a recurrence for delta using 4-term ansatz...")
# delta[0] = delta[1] = delta[2] = 0, delta[3] != 0
# So effectively delta is a sequence starting at n=3
# shifted: d_n = delta[n+3] for n = 0, 1, 2, ...
d = delta[3:]  # length N-3

# Try order 3, varying degree d
for deg in range(1, 15):
    # c3(n) d[n+3] + c2(n) d[n+2] + c1(n) d[n+1] + c0(n) d[n] = 0
    # where ci(n) = sum_j aij * n^j for j=0..deg
    # Total unknowns: 4*(deg+1)
    # Equations: N-3-3 = N-6
    num_unk = 4*(deg+1)
    num_eq = len(d) - 3
    if num_eq < num_unk:
        continue

    from sympy import symbols, Eq, solve, Rational
    # Build system over Q
    rows = []
    for n in range(num_eq):
        row = []
        for i in range(4):  # c0, c1, c2, c3
            for j in range(deg+1):
                # coefficient of a_{i,j} is n^j * d[n+i]
                val = Q(n)**j * d[n+i]
                row.append(SR(val.numerator, val.denominator))
        rows.append(row)

    M_sys = Matrix(rows[:num_unk])  # square system
    ker = M_sys.nullspace()
    if ker:
        print(f"  degree {deg}: kernel dim {len(ker)}")
        if len(ker) == 1:
            # Verify with remaining equations
            v = ker[0]
            ok = True
            for n in range(num_unk, num_eq):
                row = []
                for i in range(4):
                    for j in range(deg+1):
                        val = Q(n)**j * d[n+i]
                        row.append(SR(val.numerator, val.denominator))
                test = sum(row[k]*v[k] for k in range(len(v)))
                if test != 0:
                    ok = False
                    break
            if ok:
                print(f"    *** VERIFIED: delta satisfies order-3 degree-{deg} recurrence! ***")
                break
            else:
                print(f"    (kernel found but verification failed)")
    else:
        print(f"  degree {deg}: no kernel")
