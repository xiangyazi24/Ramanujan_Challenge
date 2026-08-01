#!/usr/bin/env python3
"""R13's decisive experiment:
A) verify phi_2, phi_3, phi_4 expansions + Green recursion;
B) curve nullity: does the Apery evaluation orbit {P_n} lie on a curve of degree <= 12?"""
import sympy as sp
import numpy as np
x = sp.symbols('x')

# Part A
def Pp(z): return 34*z**3 + 51*z**2 + 27*z + 5
K = [sp.Integer(0), sp.Integer(1)]
for d in range(1, 6):
    K.append(sp.expand(Pp(x+d)*K[d] - (x+d)**6*K[d-1]))
def Den(d):
    r = sp.Integer(1)
    for j in range(1, d+1): r *= (x+j)**3
    return sp.expand(r)
phi = {d: sp.expand(K[d] - (x+d)**3*K[d-1] - Den(d-1)) for d in range(1, 6)}
ok_rec = all(sp.expand(phi[d+1] - ((x+d)**3*phi[d] + 4*(2*x+2*d+1)**3*K[d])) == 0 for d in range(1, 5))
print("Green recursion phi_{d+1} = (x+d)^3 phi_d + 4(2x+2d+1)^3 K_d:", "VERIFIED" if ok_rec else "FAIL")
print("phi_2 =", sp.factor(phi[2]))
z = sp.symbols('z')
Q5 = 140*z**5 + 64*z**4 - 61*z**3 - 25*z**2 + 13*z + 5
print("phi_3 == 4(2z-1)Q5(z), z=x+2:", sp.expand(phi[3] - 4*(2*(x+2)-1)*Q5.subs(z, x+2)) == 0)
R8 = 2380*z**8 + 9412*z**7 + 9511*z**6 - 3920*z**5 - 8051*z**4 + 832*z**3 + 2528*z**2 - 4*z - 340
print("phi_4 == 8(2z+1)R8(z), z=x+2:", sp.expand(phi[4] - 8*(2*(x+2)+1)*R8.subs(z, x+2)) == 0)

# Part B: curve nullity mod p
def nullity_test(p, DMAX=12):
    lammu = lambda n: (pow(pow(n+1,3,p), p-2, p), 4*pow(2*n+1,3,p) % p)
    C = [[1,0],[0,1]]
    pts = {(1, 0)}
    for n in range(1, p-1):
        lam, mu = lammu(n)
        S = [[(1 + lam*mu) % p, lam], [mu % p, 1]]
        C = [[(S[0][0]*C[0][0]+S[0][1]*C[1][0]) % p, (S[0][0]*C[0][1]+S[0][1]*C[1][1]) % p],
             [(S[1][0]*C[0][0]+S[1][1]*C[1][0]) % p, (S[1][0]*C[0][1]+S[1][1]*C[1][1]) % p]]
        pts.add((C[0][0], C[0][1]))
    pts = list(pts)
    print(f"p={p}: distinct evaluation points = {len(pts)}")
    def rank_mod(M, p):
        M = [row[:] for row in M]; rows = len(M); cols = len(M[0]); r = 0
        for c in range(cols):
            piv = next((i for i in range(r, rows) if M[i][c] % p), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]
            inv = pow(M[r][c], p-2, p)
            M[r] = [v*inv % p for v in M[r]]
            for i in range(rows):
                if i != r and M[i][c]:
                    f = M[i][c]
                    M[i] = [(M[i][j] - f*M[r][j]) % p for j in range(cols)]
            r += 1
            if r == rows: break
        return r
    for D in (2, 3, 4, 6, 9, 12):
        mons = [(i, j) for t in range(D+1) for i in range(t+1) for j in [t-i]]
        # subsample points for speed if huge: rank needs at most len(mons)+1 independent pts
        sample = pts[: min(len(pts), 3*len(mons) + 20)]
        M = [[pow(a, i, p)*pow(b, j, p) % p for (i, j) in mons] for (a, b) in sample]
        rk = rank_mod(M, p)
        nul = len(mons) - rk
        print(f"  degree <= {D}: monomials={len(mons)} rank={rk} nullity={nul}")
for p in (1009, 5003):
    nullity_test(p)
