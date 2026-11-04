"""Anomaly and slope-zero check of the CICY 5421 candidate in the favourable 8-factor description (requires numpy, scipy).
Usage: python3 check_geometry.py"""
import numpy as np
from scipy.optimize import least_squares
from cicy_index import CICY
FAV=[[1,1,1,0,0,0,0,0,0],[0,0,1,0,1,0,0,0,0],[0,0,0,1,1,1,0,0,0],[0,1,0,1,0,0,0,0,0],[1,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,2],[0,0,0,0,0,1,1,0,0],[0,0,0,0,0,0,1,1,2]]
X=CICY([sum(r)-1 for r in FAV],FAV); d=np.array(X.d,float); c2=np.array([float(x) for x in X.c2J])
print("Euler number",X.euler,"(expected -48); c2.J =",c2.tolist())
SIG=(4,6,5,1,7)  # original ambient class r -> favourable class SIG[r] (matches c2 and all triple intersections)
MODELS={"A":[(-2,0,1,0,0),(0,-1,-1,0,1),(0,-1,0,1,0),(1,0,0,-2,0),(1,2,0,1,-1)],
        "B":[(-1,0,-2,1,1),(-1,1,1,0,0),(0,-1,0,1,0),(1,0,0,-2,0),(1,0,1,0,-1)]}
for name,M in MODELS.items():
    ks=[]
    for b in M:
        k=np.zeros(8)
        for r in range(5): k[SIG[r]]=b[r]
        ks.append(k)
    c1=sum(ks); c2V=np.array([-0.5*sum(np.einsum('st,s,t',d[r],k,k) for k in ks) for r in range(8)])
    print("model",name,"c1(V)=",c1.tolist(),"\n  c2(TX)-c2(V) =",(c2-c2V).tolist(),"(all > 0 => anomaly satisfiable)")
    S=[np.einsum('rst,r->st',d,k) for k in ks]; rng=np.random.default_rng(1); best=None
    for _ in range(300):
        u0=rng.uniform(-1.5,1.5,8)
        r=least_squares(lambda u:[np.exp(u)@s@np.exp(u)/np.exp(u).sum()**2 for s in S[:4]],u0,bounds=(-3.9,3.9),xtol=1e-14,ftol=1e-14,gtol=1e-14)
        if best is None or r.cost<best.cost: best=r
    t=np.exp(best.x); t/=t.max()
    J=np.array([2*s@t for s in S])
    print("  slope-zero residual %.1e at Kahler parameters t="%best.cost,np.round(t,3).tolist(),"(all >0), Jacobian rank",np.linalg.matrix_rank(J,1e-8))
