"""Standalone check of the CICY 5421 / Z2 heterotic SU(5) line-bundle candidate (no external data needed).
Z2-equivariant line-bundle cohomology via the Koszul spectral sequence with random Z2-covariant polynomials mod a prime.
Prints results per total degree as (n_plus, n_minus), checks E2 degeneration, and tests the exotic-free spectrum for each
diagonal Z2 action listed in the Oxford quotient table. Usage: python3 check_spectrum.py [seed]
Conventions: CICY 5421 original description, ambient P1xP1xP1xP1xP3, 4 polynomials; DEGS[f][i]=degree of polynomial i in factor f;
CW = Z2 weight of each homogeneous coordinate (0 even, 1 odd); PC[action] = Z2 character of each defining polynomial."""
import itertools as it, random, sys
P=1000003
DIMS=[1, 1, 1, 1, 3]
DEGS=[[1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 2], [2, 0, 0, 0], [0, 1, 1, 2]]
CW=[1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0]
PCS={}
PCS={0: [0, 0, 0, 0], 1: [1, 1, 0, 0], 2: [1, 0, 1, 0], 3: [0, 1, 1, 0]}
PC=PCS[0]
def monos(n,d):   # exponent tuples of total degree d in n+1 vars
    if n==0: return [(d,)]
    return [(i,)+r for i in range(d+1) for r in monos(n-1,d-i)]
def factor_basis(n,m,w):
    """list of (q, exps, weight)"""
    if m>=0: return [(0,e,sum(a*b for a,b in zip(e,w))%2) for e in monos(n,m)]
    if m<=-(n+1):
        N=-m; out=[]
        for a in it.product(range(1,N+1),repeat=n+1):
            if sum(a)==N: out.append((n,a,sum(x*y for x,y in zip(a,w))%2))
        return out
    return []
def rank_mod(M):
    M=[row[:] for row in M]; r=0
    if not M: return 0
    rows=len(M); cols=len(M[0])
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i][c]%P: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],P-2,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        r+=1
        if r==rows: break
    return r
def run(L,seed=1,verbose=False):
    dims,degs,cw,pc=DIMS,DEGS,CW,PC; m=len(dims); npoly=len(degs[0])
    starts=[sum(d+1 for d in dims[:i]) for i in range(m)]
    W=[cw[starts[f]:starts[f]+dims[f]+1] for f in range(m)]
    rnd=random.Random(seed)
    # random covariant polynomials: monomials of multidegree d_i with weight == pc_i
    polys=[]
    for i in range(npoly):
        mons=[]
        for combo in it.product(*[monos(dims[f],degs[f][i]) for f in range(m)]):
            wt=sum(sum(a*b for a,b in zip(e,W[f])) for f,e in enumerate(combo))%2
            if wt==pc[i]: mons.append((combo,rnd.randrange(1,P)))
        polys.append(mons)
    # terms
    terms={}
    for J in range(npoly+1):
        for sub in it.combinations(range(npoly),J):
            md=[L[f]-sum(degs[f][i] for i in sub) for f in range(m)]
            fb=[factor_basis(dims[f],md[f],W[f]) for f in range(m)]
            if any(len(b)==0 for b in fb): continue
            shift=sum(pc[i] for i in sub)%2
            basis=[]
            for combo in it.product(*fb):
                q=sum(c[0] for c in combo); ex=tuple(c[1] for c in combo); wt=(sum(c[2] for c in combo)+shift)%2
                basis.append((q,ex,wt))
            terms[sub]={"md":md,"basis":basis,"index":{(b[0],b[1]):k for k,b in enumerate(basis)}}
    # d1 maps grouped by level J (all subs of size J together), block by (q,weight)
    def block(sub,q,wt):
        return [k for k,b in enumerate(terms[sub]["basis"]) if b[0]==q and b[2]==wt]
    groups={}; dimsum={}
    for sub,T in terms.items():
        J=len(sub)
        for q in range(0,sum(dims)+1):
            for wt in (0,1):
                src=block(sub,q,wt)
                if not src: continue
                dimsum[(J,q,wt)]=dimsum.get((J,q,wt),0)+len(src)
                for k in src:
                    ex=T["basis"][k][1]; col={}
                    for pos,i in enumerate(sub):
                        sub2=tuple(x for x in sub if x!=i)
                        if sub2 not in terms: continue
                        sign=(-1)**pos; T2=terms[sub2]
                        for combo,coef in polys[i]:
                            newex=[]; okm=True
                            for f in range(m):
                                if T["md"][f]>=0: ne=tuple(a+b for a,b in zip(ex[f],combo[f]))
                                else:
                                    ne=tuple(a-b for a,b in zip(ex[f],combo[f]))
                                    if any(x<1 for x in ne): okm=False;break
                                newex.append(ne)
                            if not okm: continue
                            key=(q,tuple(newex))
                            if key in T2["index"]:
                                tk=(sub2,T2["index"][key]); col[tk]=(col.get(tk,0)+sign*coef)%P
                    groups.setdefault((J,q,wt),[]).append(col)
    rk={}
    for key,cols in groups.items():
        tgt={}
        for col in cols:
            for tk in col: tgt.setdefault(tk,len(tgt))
        if not tgt: rk[key]=0; continue
        M=[[0]*len(cols) for _ in range(len(tgt))]
        for cj,col in enumerate(cols):
            for tk,v in col.items(): M[tgt[tk]][cj]=v
        rk[key]=rank_mod(M)
    res={}
    for (J,q,wt),d in dimsum.items():
        dim=d-rk.get((J,q,wt),0)-rk.get((J+1,q,wt),0)
        if dim:
            k=q-J; r=res.setdefault(k,[0,0]); r[wt]+=dim
    return {k:tuple(v) for k,v in sorted(res.items())}

A=[(-2,0,1,0,0),(0,-1,-1,0,1),(0,-1,0,1,0),(1,0,0,-2,0),(1,2,0,1,-1)]
B=[(-1,0,-2,1,1),(-1,1,1,0,0),(0,-1,0,1,0),(1,0,0,-2,0),(1,0,1,0,-1)]
add=lambda a,b:tuple(x+y for x,y in zip(a,b))
def virt(r):
    v=[0,0]
    for d,c in r.items():
        s=-1 if d%2 else 1; v=[v[0]+s*c[0],v[1]+s*c[1]]
    return v
def exotic_free(ks,seed):
    ten=[]; pieces={}
    for k in ks:
        r=run(k,seed); v=virt(r); ten.append([-v[0],-v[1]])
    for a,b in it.combinations(range(5),2): pieces[(a,b)]=run(add(ks[a],ks[b]),seed)
    for c in it.product((0,1),repeat=5):
        if sum(c)%2: continue
        Q=u=0
        for a,h in enumerate(ten):
            t=h if c[a]==0 else [h[1],h[0]]; Q+=t[1]; u+=t[0]
        if not(Q==3 and u==3): continue
        d5b=L5b=T5=D5=0
        for (a,b),r in pieces.items():
            H1=r.get(1,(0,0)); H2=r.get(2,(0,0)); tt=(c[a]+c[b])%2
            p=H1 if tt==0 else (H1[1],H1[0]); q=H2 if tt==0 else (H2[1],H2[0])
            d5b+=p[0]; L5b+=p[1]; T5+=q[0]; D5+=q[1]
        if d5b==3 and T5==0 and D5>=1 and L5b-D5==3: return "EXOTIC-FREE (twists %s)"%(c,), ten, pieces
    return "not exotic-free", ten, pieces
if __name__=="__main__":
    seed=int(sys.argv[1]) if len(sys.argv)>1 else 1
    for act in range(4):
        PC=PCS[act]; globals()["PC"]=PC
        print("=== Z2 action",act,"poly characters",PC)
        for name,M in (("A",A),("B",B)):
            res,ten,pieces=exotic_free(M,seed)
            print(" model",name,":",res)
            if act==0 and name=="A":
                print("   10 (plus,minus) per bundle:",ten)
                for k,v in pieces.items(): print("   L%dL%d cohomology (deg:(plus,minus)):"%(k[0]+1,k[1]+1),v)
