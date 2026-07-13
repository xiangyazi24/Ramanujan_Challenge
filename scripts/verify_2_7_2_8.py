#!/usr/bin/env python3
"""Numerical verification for Problems 2.7 and 2.8 (large coefficients).
Run on uisai2 with high precision."""
from mpmath import mp, mpf, zeta, pi, sqrt

mp.dps = 200  # 200 digits for serious verification

def verify_2_7():
    """Problem 2.7: Efficient four-term recurrence for ζ(2)+ζ(3)."""
    print("=== Problem 2.7: 4-term recurrence for ζ(2)+ζ(3) ===")
    target = zeta(2) + zeta(3)

    def A(n):
        return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

    def B(n):
        return 128*(2*n+7)**3*(2*n+9)**3*(
            104060*n**6 + 1745370*n**5 + 12145238*n**4
            + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)

    def C(n):
        return 16*(n+3)**4*(2*n+9)**3*(
            3784*n**5 + 57792*n**4 + 351019*n**3
            + 1059230*n**2 + 1587211*n + 944620)

    def D(n):
        return (n+3)**4*(n+4)**6*(946*n**2 + 4515*n + 5399)

    # Recurrence: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
    # Initial conditions (VERY large integers)
    p = [mpf(-612218384750),
         mpf(-9525021973931919) / 18100,
         mpf(-29561828382772029) / 65380]
    q = [mpf(-215040420000),
         mpf(-167282265043404) / 905,
         mpf(-964185327658080) / 6071]

    N = 50
    for n in range(2, N):
        pn = (B(n)/A(n))*p[-1] - (C(n-1)/A(n-1))*p[-2] + (D(n-2)/A(n-2))*p[-3]
        qn = (B(n)/A(n))*q[-1] - (C(n-1)/A(n-1))*q[-2] + (D(n-2)/A(n-2))*q[-3]
        p.append(pn)
        q.append(qn)

    ratio = p[-1] / q[-1]
    print(f"Target ζ(2)+ζ(3) = {target}")
    print(f"p_n/q_n (N={N}) = {ratio}")
    print(f"Difference      = {ratio - target}")
    print()

def verify_2_8():
    """Problem 2.8: Very fast rational approximation of √10005/π."""
    print("=== Problem 2.8: √10005/π ===")
    target = sqrt(mpf(10005)) / pi
    print(f"Target √10005/π = {target}")

    R = mpf(151931373056001)

    def make_M(n):
        u = 2*n + 3
        w = u * (3*u - 2) * (3*u + 2)

        a1 = (144*R - 99)*u**5 - (288*R - 333)*u**4 + (144*R - 229)*u**3 - 114*u**2 + 40*u + 64
        a2 = (432*R - 243)*u**4 - (864*R - 909)*u**3 + (432*R - 868)*u**2 - 80*u + 272
        a3 = (432*R - 153)*u**3 - (864*R - 648)*u**2 + (432*R - 860)*u + 360

        b1 = 9*u**4 - (144*R - 63)*u**3 + 158*u**2 + 168*u + 64
        b2 = 36*u**3 + (216*R - 189)*u**2 - 316*u - 168
        b3 = 54*u**2 + (108*R - 189)*u - 158

        c1 = 18*u**5 + (54 + 45)*u**4 - (288*R**2 - 378*R + 251)*u**3 + (948*R - 1086)*u**2 + (1008*R - 1384)*u + (384*R - 576)
        c2 = (153*R - 72)*u**4 - (657*R - 702)*u**3 - (432*R**2 - 1292*R + 1069)*u**2 + (2064*R - 2508)*u + (1072*R - 1512)
        c3 = (180*R - 108)*u**3 - (891*R - 864)*u**2 - (216*R**2 - 1450*R + 1385)*u + (1116*R - 1422)
        c4 = (6*R - 4)*u**2 - (33*R - 32)*u - (4*R**2 - 58*R - 236337691420383)

        # Wait, the paper has c4 = (6R-4)u^2 - (33R-32)u - (4R^2 - 58R - 236337691420383)
        # Let me re-read the paper more carefully...
        # Actually the exact c4 formula is:
        # c4 = (6R-4)u^2 - (33R-32)u - (4R^2 - 58R - 236337691420383)
        # This doesn't look right. Let me use the raw definition from the paper.
        # The issue is the paper gives very specific polynomial formulas that I'd need to copy exactly.

        # For now, let me just check if the matrix product works by computing a few terms
        # and checking the ratio against sqrt(10005)/pi

        # Actually, the M(n) matrix is 4x4 with specific structure.
        # The exact transcription from the PDF is complex. Skip for now.
        pass

    print("NOTE: Problem 2.8 requires exact transcription of 4x4 matrix entries")
    print("from the challenge PDF. Deferred to manual entry or OCR on uisai2.")
    print()

if __name__ == "__main__":
    verify_2_7()
    verify_2_8()
