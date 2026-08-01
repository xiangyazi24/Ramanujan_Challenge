import numpy as np, math
def orbit(p):
    N=p-2; b=[1,5]; c=[0,6]
    for n in range(1,N+1):
        A=(34*n**3+51*n**2+27*n+5)%p; Bc=(n**3)%p; Dn=pow((n+1)**3%p,p-2,p)
        b.append(((A*b[n]-Bc*b[n-1])*Dn)%p); c.append(((A*c[n]-Bc*c[n-1])*Dn)%p)
    return np.array(b[:N+1]),np.array(c[:N+1]),N
# S(a,b) = sum_r e_p(a b_r + b c_r): compute full 2D DFT via 2D histogram
for p in [1009, 2003]:
    b,c,N = orbit(p)
    H=np.zeros((p,p))
    for r in range(1,N+1): H[b[r],c[r]]+=1
    # injectivity of r -> v_r ?
    dup = int((H*(H-1)).sum()/2)
    S=np.fft.fft2(H); A=np.abs(S)
    A0=A.copy(); A0[0,0]=0
    print(f"p={p} N={N}: vector-collisions(v_r=v_r')={dup}  max|S(a,b)|={A0.max():.0f} = {A0.max()/math.sqrt(N):.2f} sqrt(N)  mean|S|={A0.mean():.1f}  mean|S|^2={ (A0**2).mean():.0f} (Parseval predicts N={N})")
