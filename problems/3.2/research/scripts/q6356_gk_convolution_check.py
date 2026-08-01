# Third-party check of Q6356 GK convolution (cron verified p=13,17; here fresh primes)
# b_r = -(1/N) sum_t J(phi psi^{t+r}, psi^{N/2-t})^2 J(phi psi^{t-r}, psi^{N/2-t})^2 mod p
# mod-p model: psi = identity character, J(m,n) = sum_{x!=0,1} x^m (1-x)^n
def check(p, rs):
    N = p-1; half = N//2
    # precompute J table lazily
    Jcache = {}
    def J(m, n):
        m %= N; n %= N
        if (m,n) not in Jcache:
            Jcache[(m,n)] = sum(pow(x,m,p)*pow((1-x)%p,n,p) for x in range(2,p)) % p
        return Jcache[(m,n)]
    b=[1,5]
    for k in range(1,p):
        b.append(((34*k**3+51*k**2+27*k+5)*b[k]-k**3*b[k-1])*pow((k+1)**3,p-2,p)%p)
    bad=[]
    for r in rs:
        S = sum(J(half+t+r, half-t)**2 * J(half+t-r, half-t)**2 for t in range(N)) % p
        # -(1/N) = 1 mod p since N = -1
        if S != b[r] % p: bad.append((r, S, b[r]%p))
    print(f"p={p}: {'ALL MATCH' if not bad else 'MISMATCH'} on {len(rs)} values of r", bad[:4])
check(29, list(range(1,28)))
check(37, list(range(1,36)))
check(41, list(range(1,40)))
