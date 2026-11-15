#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
N=5000

def mul_factor(a,m,e):
    if e>=0:
        for _ in range(e):
            for n in range(N,m-1,-1): a[n]-=a[n-m]
    else:
        for _ in range(-e):
            for n in range(m,N+1): a[n]+=a[n-m]

def build(exps,shift):
    a=[0]*(N+1); a[0]=1
    for m in range(1,N+1):
        e=exps.get(m%25,0)
        if e: mul_factor(a,m,e)
    return [0]*shift+a[:N+1-shift]

def addE(d,g,e):
    d[g%25]=d.get(g%25,0)+e; d[(-g)%25]=d.get((-g)%25,0)+e

def pack(v,w):
    z=0
    for x in reversed(v): z=(z<<w)+x
    return z

def conv(a,b):
    la=len(a); lb=len(b)
    while la and not a[la-1]: la-=1
    while lb and not b[lb-1]: lb-=1
    if not la or not lb:return [0]*(N+1)
    aa=a[:la];bb=b[:lb]
    ap=[max(x,0) for x in aa];an=[max(-x,0) for x in aa]
    bp=[max(x,0) for x in bb];bn=[max(-x,0) for x in bb]
    bound=min(la,lb,N+1)*max(max(ap),max(an))*max(max(bp),max(bn))
    w=max(2,bound.bit_length()+1); mask=(1<<w)-1
    zz=[pack(ap,w)*pack(bp,w),pack(an,w)*pack(bn,w),pack(ap,w)*pack(bn,w),pack(an,w)*pack(bp,w)]
    c=[0]*(N+1)
    for k in range(N+1):
        c[k]=(zz[0]&mask)+(zz[1]&mask)-(zz[2]&mask)-(zz[3]&mask)
        zz=[z>>w for z in zz]
    return c

def add(*items):
    c=[0]*(N+1)
    for scale,a in items:
        for i,x in enumerate(a):c[i]+=scale*x
    return c

def digest(a):return hashlib.sha256(','.join(map(str,a)).encode()).hexdigest()

eu={}
for g in (1,7):addE(eu,g,1)
for g in (2,11):addE(eu,g,-1)
ev={}
for g in (1,7,3,4):addE(ev,g,1)
for g in (2,11):addE(ev,g,-2)
print('building U,V',flush=True)
U=build(eu,1);V=build(ev,2)
print('powers of U',flush=True)
U2=conv(U,U);U3=conv(U2,U);U4=conv(U3,U);U5=conv(U4,U);U6=conv(U5,U);U7=conv(U6,U)
A=add((1,U),(-1,U3),(-2,U4))
B=add((-1,U3),(-1,U4),(-2,U5))
print('Horner verification',flush=True)
t=add((1,V),(-1,U3))
t=add((1,conv(t,V)),(1,A))
t=add((1,conv(t,V)),(1,B))
t=add((1,conv(t,V)),(1,U7))
bad=next(((i,x) for i,x in enumerate(t) if x),None)
if bad: raise SystemExit('first residual %r'% (bad,))
out={'N':N,'U':U,'V':V,'U_sha256_csv':digest(U),'V_sha256_csv':digest(V),'relation':'Y^4 - X^3*Y^3 + (X - X^3 - 2*X^4)*Y^2 - (X^3 + X^4 + 2*X^5)*Y + X^7','relation_verified_through':N,'cleared_common_order':29,'cleared_verified_through':N+29,'U_first_80':U[:81],'V_first_80':V[:81]}
open('q6154/full.json','w').write(json.dumps(out,separators=(',',':')))
open('q6154/full-summary.json','w').write(json.dumps({k:v for k,v in out.items() if k not in ('U','V')},indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ('U','V')},indent=2),flush=True)
print('PASS FULL Q5000',flush=True)
