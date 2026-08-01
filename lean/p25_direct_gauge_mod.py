#!/usr/bin/env python3
"""Fast finite-field polynomial gauge scan for the direct Catalan 3F2 orbit."""

P = 2147483647


def inv(x):
    return pow(x % P, P - 2, P)


def div(a, b):
    return a % P * inv(b) % P


def source(n):
    return [
        [div(2*n**3+25*n**2+24*n+5,(n+1)**2*(2*n+1)),
         div(2*(8*n**2+9*n+2),(n+1)**2),
         div(48*n**4+112*n**3+97*n**2+37*n+5,(n+1)**2*(2*n+1))],
        [div(8*(3*n+2),(n+1)*(2*n+1)), div(17*n+13,n+1),
         div(2*(24*n**3+48*n**2+33*n+8),(n+1)*(2*n+1))],
        [div(4*(4*n+3),(n+1)**2*(2*n+1)), div(2*(6*n+5),(n+1)**2),
         div(34*n**3+69*n**2+48*n+12,(n+1)**2*(2*n+1))],
    ]


def target(n):
    g=2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    e=[
      (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
      384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
      480*n**4+4980*n**3+19210*n**2+32690*n+20730,
      (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
      (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
      (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
      (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
      (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
      (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240),
    ]
    return [[div(e[3*i+j],g) for j in range(3)] for i in range(3)]


def transpose(a): return [list(x) for x in zip(*a)]


def equations(degree, variant='plain', offset=0, samples_extra=5, denominator='one'):
    rows=[]
    for n in range(1, degree+samples_extra+1):
        b=source(n+offset);t=target(n)
        if variant in ('bt','both'): b=transpose(b)
        if variant in ('tt','both'): t=transpose(t)
        xp=[pow(n,k,P) for k in range(degree+1)]
        yp=[pow(n+1,k,P) for k in range(degree+1)]
        def den(x):
            if denominator == 'det':
                return (x+1)*(x+2)*(2*x+3)**2*(2*x+5)**2 % P
            if denominator == 'detn':
                return x*(x+1)*(x+2)*(2*x+1)**2*(2*x+3)**2*(2*x+5)**2 % P
            return 1
        dn,dnp=den(n),den(n+1)
        for i in range(3):
          for j in range(3):
            row=[0]*(9*(degree+1))
            for k in range(degree+1):
              for a in range(3):
                row[k*9+a*3+j]=(row[k*9+a*3+j]+dn*b[i][a]*yp[k])%P
              for c in range(3):
                row[k*9+i*3+c]=(row[k*9+i*3+c]-dnp*xp[k]*t[c][j])%P
            rows.append(row)
    return rows


def rank(rows):
    if not rows:return 0
    m=len(rows);n=len(rows[0]);r=0
    for c in range(n):
        pivot=next((i for i in range(r,m) if rows[i][c]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        z=inv(rows[r][c])
        # Echelon form is enough for rank; work only on the active suffix.
        for j in range(c,n): rows[r][j]=rows[r][j]*z%P
        for i in range(r+1,m):
            if rows[i][c]:
                z=rows[i][c]
                for j in range(c,n): rows[i][j]=(rows[i][j]-z*rows[r][j])%P
        r+=1
        if r==m:return r
    return r


if __name__=='__main__':
  import sys
  variants=(sys.argv[1],) if len(sys.argv)>1 else ('plain','bt','tt','both')
  offsets=(int(sys.argv[2]),) if len(sys.argv)>2 else range(-3,4)
  maximum=int(sys.argv[3]) if len(sys.argv)>3 else 24
  denominator=sys.argv[4] if len(sys.argv)>4 else 'one'
  for variant in variants:
    for offset in offsets:
      print('variant',variant,'offset',offset,flush=True)
      for d in range(maximum+1):
        rows=equations(d,variant,offset,denominator=denominator)
        nullity=9*(d+1)-rank(rows)
        if d%4==0 or nullity:print(d,nullity,flush=True)
        if nullity:break
