"""Independent re-implementation (from scratch) of Z2-equivariant Koszul/Bott cohomology of line bundles on a CICY in a product of projective spaces."""
import itertools as it, random, sys
sys.path.insert(0,".")
import check_spectrum as old   # only for data: DIMS, DEGS, CW, PCS, A, B
P=998244353
NF=len(old.DIMS); N=[d for d in old.DIMS]
DEG=[[old.DEGS[f][i] for f in range(NF)] for i in range(4)]   # DEG[i][f]
OFF=[sum(n+1 for n in N[:f]) for f in range(NF)]
def comps(total,parts):
    if parts==1: yield (total,); return
    for a in range(total+1):
        for r in comps(total-a,parts-1): yield (a,)+r
def factor_basis(f,d):
    n=N[f]
    if d>=0: return 0,[tuple(c) for c in comps(d,n+1)]
    if d<=-(n+1): return n,[tuple(-x-1 for x in c) for c in comps(-d-(n+1),n+1)]  # exponents e=-a, a_i>=1  -> e_i=-(a_i)
    return None,[]
def wt_exp(CW,f,e): return sum(CW[OFF[f]+i]*e[i] for i in range(len(e)))%2
def term(L,S):
    ds=[L[f]-sum(DEG[i][f] for i in S) for f in range(NF)]
    ty=[];lists=[]
    for f in range(NF):
        t,b=factor_basis(f,ds[f])
        if t is None: return None
        ty.append(t); lists.append(b)
    return ty,lists
def poly(i,CW,pc,rng):
    terms=[]
    for ex in it.product(*[list(comps(DEG[i][f],N[f]+1)) for f in range(NF)]):
        if sum(wt_exp(CW,f,ex[f]) for f in range(NF))%2==pc[i]: terms.append((ex,rng.randrange(1,P)))
    return terms
def rank_mod(M):
    import numpy as np
    M=M.copy()%P; r=0; rows,cols=M.shape
    for c in range(cols):
        if r>=rows: break
        piv=None
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        piv=r+nz[0]; M[[r,piv]]=M[[piv,r]]
        inv=pow(int(M[r,c]),P-2,P); M[r]=(M[r]*inv)%P
        idx=np.nonzero(M[:,c])[0]
        for k in idx:
            if k!=r: M[k]=(M[k]-M[k,c]*M[r])%P
        r+=1
    return r
def run(L,CW,pc,seed=1):
    import numpy as np
    rng=random.Random(seed); F=[poly(i,CW,pc,rng) for i in range(4)]
    subsets={j:[S for S in it.combinations(range(4),j)] for j in range(5)}
    # basis per (j,q): list of (S, exps-tuple-of-tuples, weight)
    B={}; 
    for j in range(5):
        for S in subsets[j]:
            t=term(L,S)
            if t is None: continue
            ty,lists=t; q=sum(ty)
            gw=sum(pc[i] for i in S)%2
            for ex in it.product(*lists):
                w=(sum(wt_exp(CW,f,ex[f]) for f in range(NF))+gw)%2
                B.setdefault((j,q),[]).append((S,ex,w))
    idx={k:{(S,ex):n for n,(S,ex,w) in enumerate(v)} for k,v in B.items()}
    ranks={}
    for (j,q),cols in B.items():
        if j==0 or (j-1,q) not in B: continue
        rows=B[(j-1,q)]; ri=idx[(j-1,q)]
        M=np.zeros((len(rows),len(cols)),dtype=np.int64)
        for c,(S,ex,w) in enumerate(cols):
            for pos,i in enumerate(S):
                S2=tuple(s for s in S if s!=i); sign=-1 if pos%2 else 1
                for mon,co in F[i]:
                    new=[]; ok=True
                    for f in range(NF):
                        e=tuple(a+b for a,b in zip(ex[f],mon[f]))
                        # for negative-type factors all exps must stay <= -1
                        if ex[f] and any(a<0 for a in ex[f]) and any(a>=0 for a in e): ok=False;break
                        new.append(e)
                    if not ok: continue
                    r=ri.get((S2,tuple(new)))
                    if r is not None: M[r,c]=(M[r,c]+sign*co)%P
        for wv in (0,1):
            cc=[c for c,(S,ex,w) in enumerate(cols) if w==wv]; rr=[r for r,(S,ex,w) in enumerate(rows) if w==wv]
            ranks[(j,q,wv)]=rank_mod(M[np.ix_(rr,cc)]) if cc and rr else 0
    res={}
    for (j,q),v in B.items():
        for wv in (0,1):
            d=sum(1 for x in v if x[2]==wv)-ranks.get((j,q,wv),0)-ranks.get((j+1,q,wv),0)
            if d: res.setdefault(q-j,[0,0])[wv]+=d
    return {k:tuple(v) for k,v in sorted(res.items())}

if __name__=="__main__":
    # cross-check against check_spectrum.py: all 15 bundles x 2 models x 4 Z2 actions
    old.__dict__["PC"]=None
    bad=n=0
    for act in range(4):
        old.__dict__["PC"]=old.PCS[act]
        for name,M in (("A",old.A),("B",old.B)):
            for L in list(M)+[tuple(x+y for x,y in zip(M[a],M[b])) for a,b in it.combinations(range(5),2)]:
                n+=1
                if run(L,old.CW,old.PCS[act],3)!={k:tuple(v) for k,v in old.run(L,3).items()}: bad+=1; print("MISMATCH",act,name,L)
    print("compared",n,"cases; mismatches:",bad)
