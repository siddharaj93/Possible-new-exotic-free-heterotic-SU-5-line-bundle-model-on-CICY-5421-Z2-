# Candidate: heterotic SU(5) line-bundle model on CICY 5421 / Z2 (status: unverified candidate)

**Not a confirmed result.** Found by an independent computer search; all cohomology computed by home-made code (checked as listed below). Novelty has not been established.

## Model
- Manifold: CICY #5421 (h11=8, h21=32, Euler -48), freely acting Z2 (quotient h11=6, h21=18). Original 5-factor ambient description (non-favourable); favourable 8-factor description exists (favourcicylist.m).
- Bundle V = sum of 5 line bundles, ambient classes of the original description (conf [[1,1,1,0,0],[1,1,0,1,0],[1,0,0,0,2],[1,2,0,0,0],[3,0,1,1,2]]):
  - Model A: (-2,0,1,0,0), (0,-1,-1,0,1), (0,-1,0,1,0), (1,0,0,-2,0), (1,2,0,1,-1)
  - Model B: (-1,0,-2,1,1), (-1,1,1,0,0), (0,-1,0,1,0), (1,0,0,-2,0), (1,0,1,0,-1)
- Wilson line Z2: diag(1,1,1,-1,-1); equivariant twists (A, actions 0,1): (0,1,0,0,1); (A, actions 2,3): none.

## Properties computed
- c1=0, index -6 on the cover (3 families on the quotient), distinct bundles, 4-dim massless-U(1)-free (K full rank), no unbroken SU(2)'.
- Anomaly: c2(TX)-c2(V) > 0 in all 8 favourable directions (A: 20,14,18,12,14,16,10,24).
- Slope zero: solutions strictly inside the ample/positive orthant, Jacobian rank 4.
- Spectrum after Z2 with Wilson line: Q=u=e=3, no 10bar, d^c=3, no Higgs triplets, 1 Higgs pair.

## Checks and limits
- Cover cohomology dims: match pyCICY (user's run) and an independent 8-factor description.
- Z2 characters: stable over 4 primes x several seeds; Serre-duality consistent; exotic-free for all 4 diagonal Z2 actions (5th is non-diagonal, untested).
- Assumes the Koszul spectral sequence degenerates at E2. Z2 characters from one code base only.
- Literature: not found in papers reachable by search (comprehensive 1307.4787 scan covered only favourable h11<=6, Wilson lines deferred). Model databases not searched.

## Novelty check vs arXiv:2607.00078 (Yip, Mininno, Shiu; LB-Explorer)
- Full text read (user-attached). CICY 5421 does not appear anywhere: Table 20 lists the ~54 CICYs searched (h11 in [4,15], |Gamma| in [2,4]); 5421 is not one of them. Their database (Anderson-Gao-Gray-Lee) contains 5421, but it was not searched.
- Their search restricts to |k_a| <= 2 plus equivariance; our models A and B satisfy |k_a| <= 2 in the original basis.
- Other literature (earlier ML/GA scans, model databases) not checked in full text.

## E2 degeneration check (2026-11-04)
Explicit E2 tables (by Koszul level J and Bott degree q) for all nonzero bundles of model A and their duals: no two E2 entries sit at positions connected by any d_r with r>=2, so the Koszul spectral sequence degenerates at E2 for every piece (given the computed E2 table). The "assumes E2" caveat is removed for model A.
