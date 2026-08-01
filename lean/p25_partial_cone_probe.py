#!/usr/bin/env python3
import mpmath as mp
mp.mp.dps=500
exec(open('p25_accelerated_partial_probe.py').read().split('import sys')[0])

base=7;c=53;G=mp.catalan
p=list(map(mp.mpf,[30921,32972,8240]));q=list(map(mp.mpf,[33750,36000,9000]))
for N in range(81):
 K=2*c*base**N
 lo=G-tail(K);hi=G-tail(K+1)
 low=[lo*q[j]-p[j] for j in range(3)]
 up=[hi*q[j]-p[j] for j in range(3)]
 mid=[(low[j]+up[j])/2 for j in range(3)]
 if N<15 or N in (20,30,40,50,60,70,80):
  def coords(e):
   return [mp.nstr(v/(-e[0]),12) for v in e]
  print(N,'low signs',[mp.sign(v) for v in low],coords(low),
        'up signs',[mp.sign(v) for v in up],coords(up),
        'unc/mid',mp.nstr((hi-lo)*max(q)/max(abs(x) for x in mid),8))
 if N<80:
  A=mat(N)
  p=[sum(p[i]*A[i][j] for i in range(3)) for j in range(3)]
  q=[sum(q[i]*A[i][j] for i in range(3)) for j in range(3)]
