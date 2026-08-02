#!/usr/bin/env python3
"""Numerical signs of the exact Meijer-G moving terminal vector."""

import mpmath as mp

mp.mp.dps = 80
G = mp.catalan
L2 = mp.log(2)
F = [
    mp.sqrt(mp.pi)*(150*G-128*L2-mp.mpf(146)/3),
    mp.sqrt(mp.pi)*(mp.mpf(24745)/4*G-mp.mpf(14624)/3*L2-mp.mpf(823511)/360),
    mp.sqrt(mp.pi)*(mp.mpf(7225281)/32*G-mp.mpf(886784)/5*L2-mp.mpf(2818419551)/33600),
]


def T(n):
    n = mp.mpf(n)
    return mp.matrix([
        [4*(n+2)*(17*n**3+111*n**2+240*n+171)/((n+1)*(n+3)*(2*n+3)*(2*n+5)),
         (n+2)*(24*n**2+101*n+102)/((n+1)*(2*n+3)),
         (n+2)*(2*n+5)*(16*n**2+81*n+90)/(2*(n+1)*(2*n+3))],
        [(96*n**4+780*n**3+2384*n**2+3273*n+1723)/((n+1)*(n+2)*(n+3)*(2*n+3)*(2*n+5)),
         (68*n**3+398*n**2+778*n+523)/(2*(n+1)*(n+2)*(2*n+3)),
         (96*n**4+884*n**3+2970*n**2+4360*n+2403)/(4*(n+1)*(n+2)*(2*n+3))],
        [-5*(24*n**2+117*n+143)/((n+1)*(n+2)*(n+3)*(2*n+3)*(2*n+5)),
         -5*(16*n+41)/(2*(n+1)*(n+2)*(2*n+3)),
         (8*n**3-44*n**2-478*n-801)/(4*(n+1)*(n+2)*(2*n+3))],
    ])


U0 = mp.matrix([
    [1, mp.mpf(152)/5, mp.mpf(195477)/175],
    [0, mp.mpf(1723)/90, mp.mpf(1963751)/2800],
    [0, -mp.mpf(143)/18, -mp.mpf(165201)/560],
])
y = mp.matrix([F])*U0**-1
for n in range(20):
    print(n, [mp.nstr(y[0, j], 16) for j in range(3)],
          [mp.sign(y[0, j]) for j in range(3)])
    y = y*T(n)
