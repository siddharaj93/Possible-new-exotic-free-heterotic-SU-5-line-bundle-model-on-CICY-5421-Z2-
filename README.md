# Candidate heterotic SU(5) line-bundle model on CICY 5421 / Z2 (unverified)

**Status: unverified candidate. Please try to break it.** This is not a claimed discovery.

On CICY #5421 (h11=8, h21=32, Euler number -48) with a freely acting Z2, a sum of five line bundles appears to give, after the Z2 quotient and a Z2 Wilson line, a spectrum with 3 families, one Higgs pair and no exotics.

## Bundles (ambient classes of the original description, ambient P1xP1xP1xP1xP3)
- A: (-2,0,1,0,0), (0,-1,-1,0,1), (0,-1,0,1,0), (1,0,0,-2,0), (1,2,0,1,-1)
- B: (-1,0,-2,1,1), (-1,1,1,0,0), (0,-1,0,1,0), (1,0,0,-2,0), (1,0,1,0,-1)

## Reproduce
```
pip install numpy scipy          # numpy/scipy only needed for check_geometry.py
python check_spectrum.py [seed]  # Z2-equivariant cohomology, exotic-free test for the 4 diagonal Z2 actions
python check_geometry.py         # anomaly and slope-zero check in the favourable 8-class description
```
(Run with plain `python`, not `python -I`.)

## What to verify independently
Z2 characters of H^q(L_aL_b): L1L2: H1 and H2 each 1-dimensional "plus"; L2L5: H1 = 2 plus + 2 minus; L3L4: H1 = 1 + 1; 10 multiplicities (1+1) for L1, L4, L5.

## Checks done
- Cover cohomology dimensions: this code = pyCICY = an independent 8-factor (favourable) description.
- E2 degeneration of the Koszul spectral sequence checked from explicit E2 tables.
- Z2 characters stable across primes/seeds and consistent with Serre duality; exotic-free for all four diagonal Z2 actions.
- Anomaly c2(TX)-c2(V) > 0 in all 8 Kahler directions; slope-zero solutions inside the Kahler cone, Jacobian rank 4.

## Limits
Z2 characters come from this repository's own code. The fifth (non-diagonal) Z2 action is untested. Novelty is not established. Earlier scans I checked do not mention CICY 5421 (see docs/candidate_note.md).

## Layout
- `check_spectrum.py`, `check_geometry.py`, `cicy_index.py`: code
- `docs/candidate_note.md`: detailed note; `docs/public_post_draft.md`: forum post draft

## Citation
See `CITATION.cff`. Add the Zenodo DOI once available.
