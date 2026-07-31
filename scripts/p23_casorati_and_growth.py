import math
N=25
A={-1:1,0:1}; B={-1:0,0:1}
for m in range(1,N+1):
    A[m]=(2*m+1)*A[m-1]+m*m*A[m-2]; B[m]=(2*m+1)*B[m-1]+m*m*B[m-2]
D={0:1,1:0}
for m in range(2,N+1): D[m]=(m-1)*(D[m-1]+D[m-2])
F={m: math.factorial(m) for m in range(0,N+1)}

print("Casorati of Lambert: A_m B_{m-1} - A_{m-1} B_m  vs  (-1)^{m+1} (m!)^2")
for m in range(0,7):
    lhs=A[m]*B[m-1]-A[m-1]*B[m]; rhs=(-1)**(m+1)*math.factorial(m)**2
    print(f"  m={m}: {lhs}  {rhs}  {'OK' if lhs==rhs else 'MISMATCH'}")

print("Casorati of derangement system: D_m F_{m-1} - D_{m-1} F_m vs (-1)^m (m-1)!")
for m in range(1,8):
    lhs=D[m]*F[m-1]-D[m-1]*F[m]; rhs=(-1)**m*math.factorial(m-1)
    print(f"  m={m}: {lhs}  {rhs}  {'OK' if lhs==rhs else 'MISMATCH'}")

# growth: A_m / m! should behave like (1+sqrt2)^m
import mpmath as mp
mp.mp.dps=40
print("A_m/(m! (1+sqrt2)^m):")
for m in [5,10,15,20,25]:
    print("  m=%2d  %s" % (m, mp.nstr(mp.mpf(A[m])/(mp.factorial(m)*(1+mp.sqrt(2))**m),12)))
print("remainder R_m = A_m*pi/4 - B_m, and R_m/A_m:")
for m in [1,3,5,10,15,20]:
    R=mp.mpf(A[m])*mp.pi/4-B[m]
    print("  m=%2d  R=%s   R/A=%s   (3-2sqrt2)^m=%s" % (m, mp.nstr(R,10), mp.nstr(R/A[m],10), mp.nstr((3-2*mp.sqrt(2))**m,10)))
