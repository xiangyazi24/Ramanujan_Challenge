# Experiment 3 (Q6762): does a low-bidegree A(r,d) vanish on ALL collision points (r,d)?
# Compare kernel dimension vs random point sets of the same size.
import random, math

def orbit_collisions(p, Dmax):
    N = p-2
    b=[1,5]; c=[0,6]
    for n in range(1,N+1):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; Dn=pow((n+1)**3%p,p-2,p)
        b.append(((A*b[n]-B*b[n-1])*Dn)%p); c.append(((A*c[n]-B*c[n-1])*Dn)%p)
    xi=[]
    for r in range(N+1):
        xi.append(('INF',) if c[r]%p==0 else (b[r]*pow(c[r],p-2,p))%p)
    pts=[]
    for d in range(1,Dmax+1):
        for r in range(1,N+1-d):
            if xi[r]==xi[r+d]: pts.append((r,d))
    return pts

def kernel_dim(pts, p, degr, degd):
    # monomials r^i d^j, i<=degr, j<=degd ; Gaussian elimination mod p
    rows=[]
    for (r,d) in pts:
        row=[]
        ri=1
        for i in range(degr+1):
            dj=1
            for j in range(degd+1):
                row.append((ri*dj)%p); dj=(dj*d)%p
            ri=(ri*r)%p
        rows.append(row)
    ncols=(degr+1)*(degd+1)
    # rank mod p
    rank=0; rows=[list(x) for x in rows]; col=0; nrows=len(rows)
    for c0 in range(ncols):
        piv=None
        for i in range(rank,nrows):
            if rows[i][c0]%p!=0: piv=i; break
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        inv=pow(rows[rank][c0],p-2,p)
        rows[rank]=[(x*inv)%p for x in rows[rank]]
        for i in range(nrows):
            if i!=rank and rows[i][c0]%p!=0:
                f=rows[i][c0]
                rows[i]=[(a-f*bb)%p for a,bb in zip(rows[i],rows[rank])]
        rank+=1
        if rank==min(nrows,ncols): break
    return ncols-rank

for p in [997, 1999]:
    Dmax=int(math.isqrt(p)*2)
    pts=orbit_collisions(p,Dmax)
    npts=len(pts)
    rnd=[(random.randrange(1,p-2), random.randrange(1,Dmax+1)) for _ in range(npts)]
    for (dr,dd) in [(6,6),(10,10),(14,8)]:
        ncols=(dr+1)*(dd+1)
        if ncols>npts+5:  # need enough points to constrain
            continue
        kc=kernel_dim(pts,p,dr,dd); kr=kernel_dim(rnd,p,dr,dd)
        print(f"p={p} D={Dmax} #pts={npts} bideg=({dr},{dd}) ncols={ncols}: kernel(collision)={kc} kernel(random)={kr}")
