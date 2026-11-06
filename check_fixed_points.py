import sys,itertools as it,random,time
sys.path.insert(0,".")
import check_spectrum_v2 as e, check_smoothness as s
o=e.old; NF=e.NF; N=e.N
def fixed_comp_points(p,CW):
    # per factor: list of points lying in an eigenspace (parity 0 or 1)
    out=[]
    for f in range(NF):
        ws=CW[e.OFF[f]:e.OFF[f]+N[f]+1]; pts=set()
        for par in (0,1):
            idx=[j for j,w in enumerate(ws) if w==par]
            if not idx: continue
            n=len(idx)
            if n==1:
                v=[0]*(N[f]+1); v[idx[0]]=1; pts.add(tuple(v))
            else:
                for c in it.product(range(p),repeat=n):
                    if not any(c): continue
                    v=[0]*(N[f]+1)
                    for k,j in enumerate(idx): v[j]=c[k]
                    pts.add(s.normalize(v,p) if True else tuple(v))
        out.append(sorted(pts))
    return out
def run(p,act,seed):
    e.P=p; rng=random.Random(seed)
    F=[[(ex,c%p or 1) for ex,c in e.poly(i,o.CW,o.PCS[act],rng)] for i in range(4)]
    comps=fixed_comp_points(p,o.CW); hits=0
    for pt in it.product(*comps):
        if all(s.evalpoly(F[i],pt,p)==0 for i in range(4)): hits+=1
    return hits, [len(c) for c in comps]
if __name__=="__main__":
    for p in (31,1009,10007):
        t=time.time()
        for act in range(4):
            r=[run(p,act,sd)[0] for sd in range(10,30)]
            print("p",p,"action",act,"fixed points of g lying on X over 20 draws:",r,"total",sum(r),"%.0fs"%(time.time()-t),flush=True)
