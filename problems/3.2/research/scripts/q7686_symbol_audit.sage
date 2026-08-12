from sage.all import *

# Q7686 exact audit: R-dependent top-symbol cancellation for the level-6
# principal-part basis Omega_m = Psi * t^{-m}.  No floating point is used.
ELL = 5
LM = ELL**3 + 1  # E(T_ell f) eigenvalue used only as a consistency label.


def P(n):
    n = ZZ(n)
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_and_kappa(N):
    b = [ZZ(1), ZZ(5)]
    for n in range(1, N):
        num = P(n)*b[n] - ZZ(n)**3*b[n-1]
        den = ZZ(n+1)**3
        assert num % den == 0
        b.append(num // den)

    Rt = PowerSeriesRing(QQ, 'x', default_prec=N+2)
    x = Rt.gen()
    F = Rt(sum(QQ(b[n])*x**n for n in range(N+1))).add_bigoh(N+2)
    Delta = (1 - 34*x + x**2).add_bigoh(N+2)
    G = (1/(F**2 * Delta.sqrt())).add_bigoh(N+2)
    g = [ZZ(G[n]) for n in range(N+1)]
    assert g[:6] == [1, 7, 192, 5520, 165168, 5068320]

    kap = [QQ(0), QQ(-36)]
    for r in range(2, N+1):
        kap.append((P(r-1)*kap[r-1] - ZZ(r-1)**3*kap[r-2] - 5*g[r]) / ZZ(r)**3)
    return b, g, kap


BSEQ, GSEQ, KAP = apery_and_kappa(520)


def build_A(prec):
    S = PowerSeriesRing(QQ, 'q', default_prec=prec)
    q = S.gen()
    A = S(1)
    # q/t = prod_{3 not| n} (1+q^n)^12.
    for n in range(1, prec):
        if n % 3 != 0:
            A = (A * (1 + q**n)**12).add_bigoh(prec)
    assert [ZZ(A[i]) for i in range(4)] == [1, 12, 78, 364]
    return S, q, A


def build_psi(prec, b):
    S, q, A = build_A(prec)
    t = (q / A).add_bigoh(prec)
    # E(q)=F(t(q)); composition is legal because v_q(t)=1.
    Fq = S(sum(QQ(b[n])*q**n for n in range(min(len(b), prec)))).add_bigoh(prec)
    E = Fq(t).add_bigoh(prec)
    H = (1 - q*A.derivative()/A).add_bigoh(prec)
    Psi = (E*H).add_bigoh(prec)
    assert [ZZ(E[i]) for i in range(5)] == [1, 5, 13, 23, 29]
    assert [ZZ(H[i]) for i in range(5)] == [1, -12, -12, -12, -12]
    # CT_q Omega_m = [q^m] Psi*A^m = b_m.
    for m in range(0, min(20, prec-1)):
        ph = (Psi*(A**m)).add_bigoh(prec)
        assert ZZ(ph[m]) == b[m]
    return S, q, A, Psi


# ------------------------------------------------------------
# Normalized finite difference
#   delta_l(A)=A S - S^l A, S=t^{-1},
#   A_{l,k}=S^{-kl} delta_l^k(T_l).
# Binomially,
#   A_{l,k}=sum_{h=0}^k (-1)^(k-h) C(k,h) t^{lh} T_l t^{-h}.
# Its V_l branch is
#   l^3 B_l(q)^k V_l,
#   B_l(q)=t(q)^l/t(q^l)-1=A(q^l)/A(q)^l-1
#          = -12 l q + O(q^2).
# Therefore on Omega_m the top surviving pole is l*m-k and has coefficient
#   d_{l,k}=l^3(-12l)^k.
# ------------------------------------------------------------


def subst_qell(f, ell, S, q, prec):
    out = S(0)
    for n in range((prec-1)//ell + 1):
        out += f[n] * q**(ell*n)
    return out.add_bigoh(prec)


def selected_matrix(R, k, ell=ELL):
    need = ell*(R-1) + k + 5
    S, q, A, Psi = build_psi(need, BSEQ)
    Aell = subst_qell(A, ell, S, q, need)
    Psiell = subst_qell(Psi, ell, S, q, need)
    B = (Aell/(A**ell) - 1).add_bigoh(need)
    assert B[0] == 0
    assert ZZ(B[1]) == -12*ell
    Bk = (B**k).add_bigoh(need)
    # Divide by q^k exactly; B has valuation one.
    Bnorm = S([Bk[n+k] for n in range(need-k)]).add_bigoh(need-k)
    d = ZZ(ell**3) * ZZ(-12*ell)**k
    assert ZZ(Bnorm[0]) == ZZ(-12*ell)**k

    # Rows are j_i=ell*(R+i)-k, columns m_j=R+j, with i,j=1..R.
    # The U_l part is far below these rows; the selected block is exactly the V_l part.
    M = matrix(ZZ, R, R)
    base = (Bnorm * Psiell).add_bigoh(need-k)
    Apow = (Aell**(R+1)).add_bigoh(need-k)
    for jj in range(R):
        core = (base * Apow).add_bigoh(need-k)
        for ii in range(jj+1):
            deg = ell*(jj-ii)
            M[ii,jj] = ZZ(ell**3 * core[deg])
        Apow = (Apow*Aell).add_bigoh(need-k)

    assert all(M[i,j] == 0 for i in range(R) for j in range(i))
    assert all(M[i,i] == d for i in range(R))
    # U_l pole bound for every binomial summand h: after t^{ell h},
    # pole <= floor((m+h)/ell)-ell*h <= floor(m/ell).
    # This is strictly below ell*m-k for the tested window.
    assert ell*(R+1)-k > (2*R+k)//ell
    return M, d


print('Q7686 SYMBOL AUDIT ell', ELL, 'lambda_E', LM)
for R in [20, 50, 100]:
    k = floor(sqrt(R))
    M, d = selected_matrix(R, k)
    ptest = next_prime(2*R)
    if ptest == ELL:
        ptest = next_prime(ptest)
    rankp = M.change_ring(GF(ptest)).rank()
    assert rankp == R
    assert gcd(d, ptest) == 1
    # Determinant is d^R by exact triangularity; do not expand the huge integer.
    print('MATRIX', 'R', R, 'k', k, 'shape', M.dimensions(),
          'diag', d, 'diag_factor', factor(abs(d)),
          'det_bits', (abs(d)**R).nbits(),
          'fresh_test_p', ptest, 'rank_mod_p', rankp,
          'top_row_first', ELL*(R+1)-k,
          'top_row_last', ELL*(2*R)-k)

# ------------------------------------------------------------
# Locked BFH targets.  A failure of the first evaluation E0=CT_q already
# proves that the t-shifted finite-difference family is not target-safe.
# We compute CT_q(A_{ell,k} Omega_r) modulo p using principal parts only.
# ------------------------------------------------------------

BIGPREC = ELL*492 + 40
SbigQ, qbigQ, AbigQ = build_A(BIGPREC)
PSIPREC = 540
SpsiQ, qpsiQ, ApsiQ, PsismallQ = build_psi(PSIPREC, BSEQ)


def coeffs_mod(seriesQ, upto, K):
    vals = []
    for i in range(upto):
        z = QQ(seriesQ[i])
        assert z.denominator() == 1
        vals.append(K(ZZ(z)))
    return vals


def target_E0_mod(p, r, k, ell=ELL):
    K = GF(p)
    bigprec = ell*r + 2
    Sb = PowerSeriesRing(K, 'q', default_prec=bigprec)
    qb = Sb.gen()
    A = Sb(coeffs_mod(AbigQ, bigprec, K)).add_bigoh(bigprec)

    smallprec = r + k + 3
    Ss = PowerSeriesRing(K, 'z', default_prec=smallprec)
    z = Ss.gen()
    As = Ss(coeffs_mod(ApsiQ, smallprec, K)).add_bigoh(smallprec)
    Psi = Ss(coeffs_mod(PsismallQ, smallprec, K)).add_bigoh(smallprec)

    invstep = (A**(-ell)).add_bigoh(bigprec)
    invpow = Sb(1)
    ans = K(0)

    for h in range(k+1):
        m = r + h
        phi = (Psi * (As**m)).add_bigoh(smallprec)

        def c_pole(j):
            if j < 0 or j > m:
                return K(0)
            return phi[m-j]

        def T_pole(j):
            # [q^{-j}] T_ell Omega_m.
            out = c_pole(ell*j)
            if j % ell == 0:
                out += K(ell**3) * c_pole(j//ell)
            return out

        # CT(t^{ell h} T_ell Omega_m): t^{ell h}=q^{ell h}A^{-ell h}.
        ct = K(0)
        maxs = ell*r
        for s in range(maxs+1):
            ct += invpow[s] * T_pole(s + ell*h)
        ans += K((-1)**(k-h) * binomial(k,h)) * ct
        invpow = (invpow * invstep).add_bigoh(bigprec)

    return ans


LOCKED = [
    # r=13 lies in (8,16], and 17>16 is fresh. sqrt(8) floor = 2.
    (17, 13, 8, [1, 2]),
    # r=492 lies in (300,600], and 2237>600 is fresh. sqrt(300) floor = 17.
    (2237, 492, 300, [1, 5, 17]),
]
for p, r, R, ks in LOCKED:
    assert BSEQ[r] % p == 0
    assert KAP[r].denominator() % p != 0
    assert KAP[r].numerator() % p == 0
    print('LOCKED_TARGET', 'p', p, 'r', r, 'R', R,
          'bmod', BSEQ[r] % p, 'kappanum_mod', KAP[r].numerator() % p)
    for k in ks:
        val = target_E0_mod(p, r, k)
        print('LOCKED_E0', 'p', p, 'r', r, 'k', k, 'E0_A_mod_p', ZZ(val))

print('Q7686_AUDIT_SUCCESS')
