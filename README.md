# Candidate heterotic SU(5) line-bundle model on CICY 5421 / Z2 (unverified)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23993984.svg)](https://doi.org/10.5281/zenodo.23993984)

Archived at Zenodo: https://doi.org/10.5281/zenodo.23993984
Status: unverified numerical candidate, not a claimed discovery.

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
See `CITATION.cff`.

## Cross-check of the two engines
`python check_spectrum_v2.py` re-implements the cohomology calculation from scratch (different data structures, different prime) and compares it with `check_spectrum.py` on 15 bundles x 2 models x 4 Z2 actions (120 cases): 0 mismatches. This checks for implementation bugs; both use the same Koszul/Bott method.

## Freeness and smoothness (numerical evidence, not a proof)
- Freeness: for each of the four diagonal Z2 actions, a dimension count over all 32 fixed-point components of g in the ambient space shows that a generic covariant member misses them (the count is also run for the trivial action, where it correctly fails).
- `python check_fixed_points.py`: counts g-fixed points on X over F_p for random covariant equations. At p=31 such points occur by coincidence (about 14 hits in 20 draws per action); at p=1009 there are none in 80 draws (4 actions x 20), as expected for a free action. Edit the prime list in the file to rerun (p=10007 is slow).
- `python check_smoothness.py P sym|generic SEED0 SEED1 ACTIONS`, for example `python check_smoothness.py 31 sym 10 13 0,1,2,3`: enumerates all F_p points of X and counts those where the Jacobian drops rank. Small primes give many coincidental singular points (also for generic equations with no symmetry); non-fixed singular points per draw were about 1.0 (p=13), 2.5 (p=19), 0.29 (p=31), 0 (p=37, 10 draws) and 0 (p=43, 4 draws). This is consistent with generic smoothness, but it is evidence, not a proof; a Groebner-basis check of the singular locus over a large field would settle it.
