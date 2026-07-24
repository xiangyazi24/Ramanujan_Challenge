#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp

mp.mp.dps = 100
TOL = mp.mpf('1e-70')


# ------------------------------------------------------------------
# A-polynomial of 7_2, abelian factor removed
# ------------------------------------------------------------------

def A_coefficients(M):
    c4 = M**14 - M**12 + 3*M**4 + 4*M**2 - 2
    c3 = (
        -2*M**18 + 5*M**16 + M**14 - 4*M**12
        + 6*M**8 + 5*M**6 + 2*M**4 - 4*M**2 + 1
    )
    c2 = (
        M**22 - 4*M**20 + 2*M**18 + 5*M**16 + 6*M**14
        - 4*M**10 + M**8 + 5*M**6 - 2*M**4
    )
    c1 = -2*M**22 + 4*M**20 + 3*M**18 - M**10 + M**8
    c0 = M**22
    return c4, c3, c2, c1, c0


def Apoly(M, L):
    c4, c3, c2, c1, c0 = A_coefficients(M)
    return ((((L + c4)*L + c3)*L + c2)*L + c1)*L + c0


def A_L(M, L):
    c4, c3, c2, c1, _ = A_coefficients(M)
    return (((5*L + 4*c4)*L + 3*c3)*L + 2*c2)*L + c1


def A_M(M, L):
    # Numerical differentiation is sufficient for this verification.
    return mp.diff(lambda x: Apoly(x, L), M)


# ------------------------------------------------------------------
# Endpoint roots and deformation chart
# ------------------------------------------------------------------

def endpoint_roots():
    sa = mp.findroot(lambda s: Apoly(s*s, s), (mp.mpf('0.58'), mp.mpf('0.60')))
    sb = mp.findroot(lambda s: Apoly(s, s), (mp.mpf('0.40'), mp.mpf('0.42')))
    return sa, sb


def chart_shapes(M, L):
    X = M*M
    u = (L + X**3)/(X*(L + X))
    r = -(1 + mp.sqrt(1 + 4*u*u))/(2*u)
    tau = 1 - r*r
    t = tau
    v = u/X
    w = 1/(1 - u*X)

    Aaux = X*tau - r
    Baux = tau - r*X
    Caux = tau - r

    return {
        'X': X, 'r': r, 't': t, 'u': u, 'v': v, 'w': w,
        'A': Aaux, 'B': Baux, 'C': Caux,
    }


def zp(z):
    return 1/(1-z)


def zpp(z):
    return (z-1)/z


def check_shapes(M, L, S):
    X, r = S['X'], S['r']
    t, u, v, w = S['t'], S['u'], S['v'], S['w']
    Aaux, Baux, Caux = S['A'], S['B'], S['C']

    residuals = {
        'edge_s': u**2*zp(v)*zpp(v)*w*zp(w) - 1,
        'edge_0': t**2*zp(t)*v*zpp(w) - 1,
        'edge_1': zpp(t)**2*zp(u)*zp(v)*w - 1,
        'edge_2': zp(t)*zp(u)*zpp(u)**2*v*zpp(v)*zp(w)*zpp(w) - 1,
        'reduced_gluing': X*r**4*(1-r*r) - Aaux*Baux*Caux,
        'meridian': u/v - X,
        'longitude': -X**2*Aaux/Baux - L,
        'A_polynomial': Apoly(M, L),
    }
    return residuals


# ------------------------------------------------------------------
# Upper-half-plane limiting logarithms and flattenings
# ------------------------------------------------------------------

def upper_logs_real(z):
    z = mp.mpf(z)
    if z < 0:
        return mp.log(-z) + mp.pi*1j, mp.log(1-z)
    if z > 1:
        return mp.log(z), mp.log(z-1) - mp.pi*1j
    return mp.log(z), mp.log(1-z)


def lifted_triple(z):
    lz, l1 = upper_logs_real(z)
    lp = -l1
    lpp = l1 - lz + mp.pi*1j
    return lz, lp, lpp


def check_lifted_edges(S):
    t, u, v, w = S['t'], S['u'], S['v'], S['w']
    Lt, Ltp, Ltpp = lifted_triple(t)
    Lu, Lup, Lupp = lifted_triple(u)
    Lv, Lvp, Lvpp = lifted_triple(v)
    Lw, Lwp, Lwpp = lifted_triple(w)

    two_pi_i = 2*mp.pi*1j
    return {
        'E_s': 2*Lu + Lvp + Lvpp + Lw + Lwp - two_pi_i,
        'E_0': 2*Lt + Ltp + Lv + Lwpp - two_pi_i,
        'E_1': 2*Ltpp + Lup + Lvp + Lw - two_pi_i,
        'E_2': Ltp + Lup + 2*Lupp + Lv + Lvpp + Lwp + Lwpp - two_pi_i,
    }


# ------------------------------------------------------------------
# Extended Rogers values on the upper real boundary
# ------------------------------------------------------------------

def Li2_upper_real(x):
    x = mp.mpf(x)
    if x <= 1:
        return mp.polylog(2, x)
    return (
        mp.pi**2/6
        - mp.polylog(2, 1-x)
        - mp.log(x)*mp.log(x-1)
        + mp.pi*1j*mp.log(x)
    )


def Rhat_upper_real(x):
    lz, l1 = upper_logs_real(x)
    return Li2_upper_real(x) + mp.mpf('0.5')*lz*l1 - mp.pi**2/6


def regulator(S, eps=(1, 1, 1, 1)):
    vals = [S['t'], S['u'], S['v'], S['w']]
    return sum(e*Rhat_upper_real(z) for e, z in zip(eps, vals))


# ------------------------------------------------------------------
# Follow the A-polynomial branch to check the differential index
# ------------------------------------------------------------------

def track_L(M0, L0, M1, steps=80):
    M, L = mp.mpf(M0), mp.mpf(L0)
    for j in range(1, steps + 1):
        Mn = M0 + (M1-M0)*mp.mpf(j)/steps
        L = mp.findroot(lambda y: Apoly(Mn, y), L)
        M = Mn
    return L


def differential_index(M, L):
    # The tangent to A(M,L)=0.
    dL = -A_M(M, L)/A_L(M, L)

    h = mp.mpf('1e-30')
    Mp, Mm = M+h, M-h
    Lp = mp.findroot(lambda y: Apoly(Mp, y), L + dL*h)
    Lm = mp.findroot(lambda y: Apoly(Mm, y), L - dL*h)

    Sp, Sm = chart_shapes(Mp, Lp), chart_shapes(Mm, Lm)
    S0 = chart_shapes(M, L)

    Omega = 0
    for key in ('t', 'u', 'v', 'w'):
        lz0, l10 = upper_logs_real(S0[key])
        lzp, l1p = upper_logs_real(Sp[key])
        lzm, l1m = upper_logs_real(Sm[key])
        dlz = (lzp-lzm)/(2*h)
        dl1 = (l1p-l1m)/(2*h)
        Omega += lz0*dl1 - l10*dlz

    omega_M = mp.log(M)*(dL/L) - mp.log(L)/M
    return Omega/omega_M


def show_dict(title, D):
    print(title)
    for k, v in D.items():
        print(f'  {k:20s} {mp.nstr(v, 25)}')
    print()


def main():
    sa, sb = endpoint_roots()
    Ma, La = sa**2, sa
    Mb, Lb = sb, sb

    Sa = chart_shapes(Ma, La)
    Sb = chart_shapes(Mb, Lb)

    print('endpoint roots')
    print('s_alpha =', mp.nstr(sa, 60))
    print('s_beta  =', mp.nstr(sb, 60))
    print()

    for name, M, L, S in (
        ('alpha', Ma, La, Sa),
        ('beta', Mb, Lb, Sb),
    ):
        print(name, 'shapes')
        for k in ('t', 'u', 'v', 'w'):
            print(f'  {k} = {mp.nstr(S[k], 60)}')
        print()
        show_dict(name + ' multiplicative residuals', check_shapes(M, L, S))
        show_dict(name + ' lifted-log residuals', check_lifted_edges(S))

    Ra = regulator(Sa)
    Rb = regulator(Sb)
    DeltaR = Rb - Ra
    target = -4*mp.pi**2/85

    print('endpoint Rogers sums')
    print('R(alpha) =', mp.nstr(Ra, 70))
    print('R(beta)  =', mp.nstr(Rb, 70))
    print('difference=', mp.nstr(DeltaR, 70))
    print('target    =', mp.nstr(target, 70))
    print('difference-target =', mp.nstr(DeltaR-target, 30))
    print('(difference-target)/pi^2 =',
          mp.nstr((DeltaR-target)/mp.pi**2, 30))
    print()

    Mmid = (Ma + Mb)/2
    Lmid = track_L(Ma, La, Mmid)
    print('differential index Omega/omega_M =',
          mp.nstr(differential_index(Mmid, Lmid), 40))
    print('expected N_M = -2')


if __name__ == '__main__':
    main()
