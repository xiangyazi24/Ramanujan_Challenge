#!/usr/bin/env python3
"""Compare challenge brackets with Catalan partial sums at c*6^N."""

import mpmath as mp
mp.mp.dps=100

def mat(n):
 d=mp.mpf(2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
 e=[
 (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
 480*n**4+4980*n**3+19210*n**2+32690*n+20730,
 (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
 (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
 (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
 (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
 (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
 (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240)]
 return [[mp.mpf(e[3*i+j])/d for j in range(3)] for i in range(3)]

def tail(K):
 a=mp.mpf(K)+mp.mpf('.5')
 return (-1)**K/16*(mp.polygamma(1,a/2)-mp.polygamma(1,(a+1)/2))

import sys
base=int(sys.argv[1]) if len(sys.argv)>1 else 6
limit=int(sys.argv[2]) if len(sys.argv)>2 else 100
G=mp.catalan
for c in range(1,limit+1):
 p=list(map(mp.mpf,[30921,32972,8240]));q=list(map(mp.mpf,[33750,36000,9000]))
 ok=True;worst=mp.inf;rows=[]
 for N in range(30):
  K=2*c*base**N
  lo=G-tail(K);hi=G-tail(K+1)
  ratios=[p[j]/q[j] for j in range(3)]
  margins=(lo-ratios[2],ratios[0]-hi)
  worst=min(worst,*margins)
  if min(margins)<=0:ok=False;break
  A=mat(N)
  p=[sum(p[i]*A[i][j] for i in range(3)) for j in range(3)]
  q=[sum(q[i]*A[i][j] for i in range(3)) for j in range(3)]
 if ok:
  print('works',c,'worst',mp.nstr(worst,12))
  break
 else:
  print('fails',c,'at',N,'margins',*[mp.nstr(x,8) for x in margins])
