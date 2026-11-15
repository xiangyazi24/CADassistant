#!/usr/bin/env python3
from fractions import Fraction
from math import gcd
from functools import reduce
import json

N=180

def mul_factor(a,m,e):
    if e>=0:
        for _ in range(e):
            for n in range(N,m-1,-1): a[n]-=a[n-m]
    else:
        for _ in range(-e):
            for n in range(m,N+1): a[n]+=a[n-m]

def series(exps,shift):
    a=[0]*(N+1); a[0]=1
    for m in range(1,N+1):
        e=exps.get(m%25,0)
        if e: mul_factor(a,m,e)
    return [0]*shift+a[:N+1-shift]

def conv(a,b):
    c=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(min(len(b),N+1-i)):
                if b[j]: c[i+j]+=x*b[j]
    return c

def addE(d,g,e):
    d[g%25]=d.get(g%25,0)+e; d[(-g)%25]=d.get((-g)%25,0)+e

eu={}
for g in (1,7):addE(eu,g,1)
for g in (2,11):addE(eu,g,-1)
ev={}
for g in (1,7,3,4):addE(ev,g,1)
for g in (2,11):addE(ev,g,-2)
U=series(eu,1); V=series(ev,2)
print('series ready',flush=True)
Up=[[1]+[0]*N]; Vp=[[1]+[0]*N]
for i in range(7):Up.append(conv(Up[-1],U))
for j in range(4):Vp.append(conv(Vp[-1],V))
mons=[]; ij=[]
for j in range(5):
  for i in range(8):
    ij.append((i,j)); mons.append(Up[i] if j==0 else Vp[j] if i==0 else conv(Up[i],Vp[j]))
print('monomials ready',flush=True)
rows=[[m[n] for m in mons] for n in range(100)]
# choose 39 independent rows modulo p
p=1000000007; basis={}; chosen=[]
for idx,row in enumerate(rows):
  v=[x%p for x in row]
  while True:
    q=next((k for k,x in enumerate(v) if x),None)
    if q is None: break
    if q not in basis:
      inv=pow(v[q],p-2,p); v=[x*inv%p for x in v]; basis[q]=v; chosen.append(idx); break
    t=v[q]; w=basis[q]; v=[(x-t*y)%p for x,y in zip(v,w)]
  if len(chosen)==39: break
print('chosen',chosen[-1],chosen,flush=True)
A=[[Fraction(x) for x in rows[r]] for r in chosen]
# Rational RREF of 39x40.
rr=0; piv=[]
for col in range(40):
  k=next((k for k in range(rr,39) if A[k][col]),None)
  if k is None: continue
  A[rr],A[k]=A[k],A[rr]
  z=A[rr][col]; A[rr]=[x/z for x in A[rr]]
  for k in range(39):
    if k!=rr and A[k][col]:
      z=A[k][col]; A[k]=[x-z*y for x,y in zip(A[k],A[rr])]
  piv.append(col); rr+=1
  if rr==39: break
free=next(c for c in range(40) if c not in piv)
v=[Fraction(0) for _ in range(40)]; v[free]=1
for r in range(38,-1,-1):
  c=piv[r]; v[c]=-sum(A[r][j]*v[j] for j in range(c+1,40))
lcm=1
for x in v: lcm=lcm*x.denominator//gcd(lcm,x.denominator)
z=[int(x*lcm) for x in v]
g=reduce(gcd,(abs(x) for x in z if x),0); z=[x//g for x in z]
if next(x for x in reversed(z) if x)<0:z=[-x for x in z]
# exact verification to N
bad=[]
for n in range(N+1):
  s=sum(c*mons[k][n] for k,c in enumerate(z))
  if s: bad.append((n,s)); break
print('kernel',z,flush=True); print('bad',bad,flush=True)
open('q6154/mini-relation.json','w').write(json.dumps({'ij':ij,'coeffs':z,'chosen_last':chosen[-1],'chosen':chosen,'verified':N,'bad':bad},indent=2))
if bad: raise SystemExit(1)
print('PASS MINI RELATION',flush=True)
