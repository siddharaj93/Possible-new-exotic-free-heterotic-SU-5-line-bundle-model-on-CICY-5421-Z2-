# Possible new exotic-free heterotic SU(5) line-bundle model on CICY 5421 (Z2) — please try to break it

I'm an independent researcher (Siddharaj Pande). I ran a computer search for heterotic E8xE8 SU(5) line-bundle models on complete-intersection Calabi-Yau threefolds with free Z2 quotients and found a candidate on CICY #5421 (h11=8, h21=32, Euler -48), a manifold that the scans I could find (arXiv:1307.4787, 2003.13339, 2306.03147, 2607.00078 and others) do not appear to cover. I could not find it in those papers' text, but I could not search every database, so I am NOT claiming it is new, only that I did not find it.

## Claim
With V = sum of five line bundles (ambient classes of the original description, ambient P1xP1xP1xP1xP3):
- A: (-2,0,1,0,0), (0,-1,-1,0,1), (0,-1,0,1,0), (1,0,0,-2,0), (1,2,0,1,-1)
- B: (-1,0,-2,1,1), (-1,1,1,0,0), (0,-1,0,1,0), (1,0,0,-2,0), (1,0,1,0,-1)
the Z2 quotient with Wilson line diag(1,1,1,-1,-1) has 3 families, one Higgs pair, no Higgs triplets, no 10bar, no extra exotics, no extra massless U(1), no unbroken SU(2)'. c1=0; anomaly c2(TX)-c2(V) positive in all 8 Kahler directions (favourable description); slope-zero solutions inside the Kahler cone (Jacobian rank 4). Holds for all four diagonal Z2 actions in the Oxford quotient table.

## What was checked
- Cover cohomology dimensions: my code = pyCICY = an independent 8-factor description of the same manifold.
- Koszul spectral sequence degenerates at E2 (explicit E2 tables).
- Z2 characters: stable across primes/seeds, consistent with Serre duality.
- Reproduced by me on two machines (Linux, Windows).

## What I need
An independent recomputation of the Z2 characters (L1L2: H1,H2 both "plus", dimension 1; L2L5: H1 = 2 plus + 2 minus; L3L4: H1 = 1 + 1) with a different package, and anyone who knows whether this model is already known. Code and README: <link to repository>. Run `python check_spectrum.py` and `python check_geometry.py`.

## Limits
All Z2 characters come from my own Koszul code. The fifth (non-diagonal) Z2 action is untested. Novelty is unproven.
