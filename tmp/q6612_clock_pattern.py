#!/usr/bin/env python3
import sys,time
sys.setrecursionlimit(10000000)
def bit(x,i):return (x>>i)&1
def moves(n,c):
 out=[]
 def a(v):return bit(c,2*v)
 def b(v):return bit(c,2*v+1)
 for v in range(n):
  pred=(v-1)%n;succ=(v+1)%n; root=(v==0); succroot=(succ==0)
  xl=a(pred); xr=b(succ) if succroot else a(succ); upr=0 if succroot else b(succ)
  if root:
   if a(v)==xr and upr==0:out.append((v,'B',c^(1<<(2*v))))
   if b(v)!=xl:out.append((v,'T',c^(1<<(2*v+1))))
  else:
   if a(v)!=xl:out.append((v,'L',(c^(1<<(2*v)))|(1<<(2*v+1))))
   if a(v)==xr and b(v)==1 and upr==0:out.append((v,'U',c&~(1<<(2*v+1))))
 return out
def worst(n):
 V=4**n; status=bytearray(V); dist=[0]*V; choice=[None]*V
 def dfs(v):
  mm=moves(n,v)
  if len(mm)==1:return 0
  if status[v]==1:raise RuntimeError('cycle')
  if status[v]==2:return dist[v]
  status[v]=1; best=-1;bc=None
  for site,k,w in mm:
   z=1+dfs(w)
   if z>best:best=z;bc=(site,k,w)
  status[v]=2;dist[v]=best;choice[v]=bc;return best
 best=-1;arg=None
 for v in range(V):
  z=dfs(v)
  if z>best:best=z;arg=v
 path=[];v=arg
 while choice[v] is not None:
  site,k,w=choice[v];path.append((site,k));v=w
 def fields(c):return [(bit(c,2*i),bit(c,2*i+1)) for i in range(n)]
 return best,arg,fields(arg),path
for n in range(3,11):
 t=time.time();best,arg,f,p=worst(n);print('N',n,'BEST',best,'ARG',arg,'FIELDS',f,'PATH',p,'SEC',round(time.time()-t,3),flush=True)
