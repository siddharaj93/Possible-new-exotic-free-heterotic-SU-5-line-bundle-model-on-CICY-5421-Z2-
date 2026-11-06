import sys,itertools as it,random,time
sys.path.insert(0,".")
import check_spectrum_v2 as e
o=e.old; NF=e.NF; N=e.N
def P1pts(p): return [(1,t) for t in range(p)]+[(0,1)]
def evalpoly(F,pt,p):
    s=0
    for ex,c in F:
        t=c
        for f in range(NF):
            for j,k in enumerate(ex[f]):
                if k: t=t*pow(pt[f][j],k,p)%p
        s+=t
    return s%p
def jac(F,pt,p):
    J=[[0]*12 for _ in range(4)]
    for i,Fi in enumerate(F):
        for ex,c in Fi:
            for f in range(NF):
                for j in range(N[f]+1):
                    k=ex[f][j]
                    if not k: continue
                    t=c*k
                    for f2 in range(NF):
                        for j2,k2 in enumerate(ex[f2]):
                            kk=k2-(1 if (f2==f and j2==j) else 0)
                            if kk: t=t*pow(pt[f2][j2],kk,p)%p
                    J[i][e.OFF[f]+j]=(J[i][e.OFF[f]+j]+t)%p
    return J
def rank(M,p):
    M=[r[:] for r in M]; r=0
    for c in range(len(M[0])):
        piv=next((x for x in range(r,len(M)) if M[x][c]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[v*inv%p for v in M[r]]
        for x in range(len(M)):
            if x!=r and M[x][c]%p:
                f=M[x][c]; M[x]=[(a-f*b)%p for a,b in zip(M[x],M[r])]
        r+=1
        if r==len(M): break
    return r
def nullspace(rows,n,p):
    # all vectors y in F_p^n with rows.y=0 ; return basis
    M=[r[:] for r in rows]; piv=[]; r=0
    for c in range(n):
        pv=next((x for x in range(r,len(M)) if M[x][c]%p),None)
        if pv is None: continue
        M[r],M[pv]=M[pv],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[v*inv%p for v in M[r]]
        for x in range(len(M)):
            if x!=r and M[x][c]%p:
                f=M[x][c]; M[x]=[(a-f*b)%p for a,b in zip(M[x],M[r])]
        piv.append(c); r+=1
    free=[c for c in range(n) if c not in piv]; basis=[]
    for fc in free:
        v=[0]*n; v[fc]=1
        for i,pc_ in enumerate(piv): v[pc_]=(-M[i][fc])%p
        basis.append(v)
    return basis
def normalize(v,p):
    for x in v:
        if x%p: inv=pow(x,p-2,p); return tuple(a*inv%p for a in v)
def run(p,act,seed,sym=True):
    e.P=p; rng=random.Random(seed)
    CW=o.CW if sym else [0]*len(o.CW); pc=o.PCS[act] if sym else [0]*4
    F=[[(ex,c%p or 1) for ex,c in e.poly(i,CW,pc,rng)] for i in range(4)]
    pts1=P1pts(p); npts=0; sing=[]; 
    for x1 in pts1:
      for x2 in pts1:
        for x4 in pts1:
            base=[None]*5; base[0]=x1; base[1]=x2; base[3]=x4
            # poly0 depends only on factors 0,1,3
            pt0=[x1,x2,(0,0),x4,(0,0,0,0)]
            if evalpoly(F[0],pt0,p)!=0: continue
            # poly1 linear in y given x1; poly2 linear in y given x2: coefficient vector by evaluating at unit vectors
            def lin(Fi,fi_idx,xv):
                out=[]
                for j in range(4):
                    y=[0]*4; y[j]=1
                    pt=[(0,0)]*5; pt[fi_idx]=xv; pt[4]=tuple(y)
                    out.append(evalpoly(Fi,pt,p))
                return out
            h1=lin(F[1],0,x1); h2=lin(F[2],1,x2)
            ns=nullspace([h1,h2],4,p)
            if not ns: continue
            # enumerate projective points of span(ns)
            ys=set()
            for co in it.product(range(p),repeat=len(ns)):
                if not any(co): continue
                v=[sum(c*b[i] for c,b in zip(co,ns))%p for i in range(4)]
                ys.add(normalize(v,p))
            for y in ys:
                for x3 in pts1:
                    pt=[x1,x2,x3,x4,y]
                    if evalpoly(F[3],pt,p)==0:
                        npts+=1
                        if rank(jac(F,pt,p),p)<4: sing.append(pt)
    return npts,sing
if __name__=="__main__":
    p=int(sys.argv[1]); sym=sys.argv[2]=="sym"; seeds=range(int(sys.argv[3]),int(sys.argv[4])); acts=[int(a) for a in sys.argv[5].split(",")]
    t=time.time()
    for act in acts:
        res=[]
        for s in seeds:
            n,sg=run(p,act,s,sym); res.append((n,len(sg)))
        print("p",p,"sym" if sym else "generic(no Z2)","action",act,"(points on X, singular):",res,"%.0fs"%(time.time()-t),flush=True)
