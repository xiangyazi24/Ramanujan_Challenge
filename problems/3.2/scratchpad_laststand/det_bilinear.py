import numpy as np, math
# B(t) = sum_{r,u} W(u-r) e_p(t*(b_r c_u - b_u c_r)), W = indicator of strip (D/2, D]
def orbit(p):
    N=p-2; b=[1,5]; c=[0,6]
    for n in range(1,N+1):
        A=(34*n**3+51*n**2+27*n+5)%p; Bc=(n**3)%p; Dn=pow((n+1)**3%p,p-2,p)
        b.append(((A*b[n]-Bc*b[n-1])*Dn)%p); c.append(((A*c[n]-Bc*c[n-1])*Dn)%p)
    return np.array(b[:N+1]),np.array(c[:N+1]),N

for p in [1009, 2003, 4003]:
    b,c,N = orbit(p); D=int(math.isqrt(N)*2)
    rs=np.arange(1,N+1-D)
    maxB=0; tot=0; cnt=0
    for t in range(1,min(p,60)):   # sample of frequencies
        Bt=0
        for dd in range(D//2+1, D+1):
            u=rs+dd
            det=(t*(b[rs]*c[u]-b[u]*c[rs]))%p
            Bt+=np.exp(2j*np.pi*det/p).sum()
        a=abs(Bt); maxB=max(maxB,a); tot+=a; cnt+=1
    npairs=len(rs)*(D-D//2)
    print(f"p={p} N={N} D={D} #strip pairs={npairs}  max|B(t)|={maxB:.0f}  mean|B|={tot/cnt:.0f}  sqrt(pairs)={math.sqrt(npairs):.0f}  needed<=p^1={p}")
