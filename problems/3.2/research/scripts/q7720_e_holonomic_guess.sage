from sage.all import *

N = 100

H1 = [QQ(0)]*(2*N+3)
H2 = [QQ(0)]*(2*N+3)
for n in range(1, len(H1)):
    H1[n] = H1[n-1] + QQ(1)/n
    H2[n] = H2[n-1] + QQ(1)/(n*n)


def e_value(s):
    ans = QQ(0)
    for a in range(s+1):
        T = ZZ(binomial(s,a))**2 * ZZ(binomial(s+a,a))**2
        q = H1[s+a] + H1[s-a] - 2*H1[a]
        d = H1[a] - H1[s-a]
        ans += T * (-H2[a] - H2[s-a] + 2*q*d)
    return ans


def D_value(s):
    ans = QQ(0)
    for a in range(s+1):
        T = ZZ(binomial(s,a))**2 * ZZ(binomial(s+a,a))**2
        ans += T*(H1[s+a] - H1[s-a])
    return ans

seq = [e_value(s) for s in range(N+1)]
b = [sum(ZZ(binomial(s,a))**2 * ZZ(binomial(s+a,a))**2 for a in range(s+1)) for s in range(N+1)]
D = [D_value(s) for s in range(N+1)]

print('FIRST_E', seq[:10])
print('FIRST_CLEARED_FACT2', [seq[s]*factorial(s)^2 for s in range(10)])
print('FIRST_RATIO_DEC', [RR(seq[s]/b[s]) for s in range(1,10)])
print('MINUS_2_ZETA2', RR(-2*zeta(2)))

def P(x):
    return 34*x^3 + 51*x^2 + 27*x + 5
S = [None,None]
for n in range(2,N+1):
    S.append(n^3*seq[n] - P(n-1)*seq[n-1] + (n-1)^3*seq[n-2])
print('APERY_DEFECT_FIRST', S[2:12])


def guess_rec(sequence, fit_end=70, verify_end=100, max_order=10, max_deg=18):
    for r in range(1,max_order+1):
        for d in range(0,max_deg+1):
            cols=(r+1)*(d+1)
            rows=[]
            for n in range(0,fit_end-r+1):
                row=[]
                for j in range(r+1):
                    for k in range(d+1):
                        row.append(QQ(n)^k * sequence[n+j])
                rows.append(row)
            if len(rows) < cols-1:
                continue
            M=matrix(QQ,rows)
            ker=M.right_kernel()
            if ker.dimension()==0:
                continue
            for v in ker.basis():
                ok=True
                for n in range(fit_end-r+1,verify_end-r+1):
                    z=QQ(0); idx=0
                    for j in range(r+1):
                        for k in range(d+1):
                            z += v[idx]*QQ(n)^k*sequence[n+j]; idx+=1
                    if z != 0:
                        ok=False; break
                if ok:
                    den=lcm([x.denominator() for x in v])
                    w=vector(ZZ,[ZZ(x*den) for x in v])
                    g=gcd([abs(x) for x in w if x])
                    if g: w=vector(ZZ,[x//g for x in w])
                    nz=[x for x in w if x]
                    if nz and nz[-1] < 0: w=-w
                    R.<x>=PolynomialRing(ZZ)
                    pol=[]; idx=0
                    for j in range(r+1):
                        pol.append(sum(w[idx+k]*x^k for k in range(d+1)))
                        idx += d+1
                    return r,d,pol
    return None

rec = guess_rec(seq)
print('REC_E', rec)
if rec:
    print('REC_E_FACTORS')
    for j,Q in enumerate(rec[2]): print(j, factor(Q))

seqclr=[seq[s]*factorial(s)^2 for s in range(N+1)]
recclr=guess_rec(seqclr,max_order=10,max_deg=18)
print('REC_FACT2_E',recclr)
if recclr:
    print('REC_FACT2_E_FACTORS')
    for j,Q in enumerate(recclr[2]): print(j, factor(Q))

# Defect starts at n=2; shift it to t_k=S_{k+2} and guess separately.
Sshift=[S[k+2] for k in range(N-1)]
recS=guess_rec(Sshift,fit_end=65,verify_end=98,max_order=8,max_deg=14)
print('REC_APERY_DEFECT_SHIFTED',recS)
if recS:
    print('REC_APERY_DEFECT_SHIFTED_FACTORS')
    for j,Q in enumerate(recS[2]): print(j, factor(Q))
    # characteristic polynomial at infinity
    Rz.<z>=PolynomialRing(QQ)
    lead=[]
    d=recS[1]
    for Q in recS[2]: lead.append(Q[d])
    print('DEFECT_CHAR_INF',factor(sum(lead[j]*z^j for j in range(len(lead)))))

# Search low-degree expression for Apéry defect using b_n,b_{n-1},D_n,D_{n-1}.
def fit_source(maxdeg=8):
    for d in range(maxdeg+1):
        rows=[]; rhs=[]
        for n in range(2,60):
            row=[]
            for arr,shift in [(b,0),(b,-1),(D,0),(D,-1)]:
                for k in range(d+1): row.append(QQ(n)^k*arr[n+shift])
            rows.append(row); rhs.append(S[n])
        M=matrix(QQ,rows); y=vector(QQ,rhs)
        try: sol=M.solve_right(y)
        except ValueError: continue
        ok=True
        for n in range(60,N+1):
            z=QQ(0); idx=0
            for arr,shift in [(b,0),(b,-1),(D,0),(D,-1)]:
                for k in range(d+1): z += sol[idx]*QQ(n)^k*arr[n+shift]; idx+=1
            if z != S[n]: ok=False; break
        if ok: return d,sol
    return None
print('SOURCE_FIT_bD',fit_source())

W=[]
for s in range(N): W.append(b[s]*seq[s+1]-b[s+1]*seq[s])
print('FIRST_W',W[:10])
print('FIRST_W_RATIO_DEC',[RR(W[s]/(b[s]*b[s+1])) for s in range(1,10)])
print('Q7720_GUESS_DONE')
