# CLE15 v1.0 modern Python adapter

This directory contains an auditable Python 3 compatibility candidate for the frozen official CLE15 v1.0 release. It is derived from the official Python file, with one branch corrected against the authoritative MATLAB implementation in the same official ZIP. The untouched source remains under `TCRoad/vendor/CLE15/10.4231_CZ4P-D448/v1.0/`.

## What is closed

All three official default entry cases run under modern Python, NumPy, SciPy, and Shapely and reproduce the vector curves in the three official PDFs to at most 0.001609 m s-1 at 38 checked anchors:

- `ER11E04_nondim_r0input`: `Vmax=50 m s-1`, `r0=900 km`, `f=5e-5 s-1`.
- `ER11E04_nondim_rmaxinput`: `Vmax=50 m s-1`, `rmax=25 km`, `f=5e-5 s-1`.
- `ER11E04_nondim_rfitinput`: `Vmax=50 m s-1`, `rfit=250 km`, `Vfit=12 m s-1`, `f=5e-5 s-1`.
- Shared official example environment: constant `Cd=1.5e-3`, `w_cool=0.002 m s-1`, constant `Ck/Cd=1`, and `eye_adj=0`.

The exact edits and the one MATLAB-authoritative correction are recorded in `TRANSFORM_LOG.md`.

One official-Python API inconsistency is deliberately preserved: `ER11E04_nondim_r0input` returns `(rr, VV, rmerge, Vmerge, rmax)`, whereas the MATLAB entry returns `(rr, VV, rmax, rmerge, Vmerge)`. Callers must use the Python order shown here; silently reordering it inside the mechanical adapter would break compatibility with the frozen Python file.

## What is not claimed

This is not yet proof of pointwise equivalence for every possible storm and parameter combination. The official Python file leaves its eye-adjustment implementation commented out, so `eye_adj=1` is not validated and must not be used through this adapter. The adapter also does not add translation asymmetry, surface reduction, environmental flow, rainfall physics, or any TCR coupling.

## Test

From the project root, in an environment satisfying `requirements.txt`:

```bash
python -m unittest discover -s TCRoad/adapters/c15_python3/tests -v
```

The test compares the adapter directly with anchors extracted from the official vector PDFs; it does not compare the adapter with itself.

## Verified runtimes

The committed tests passed in both of these existing local environments on 2026-08-13:

- `/opt/anaconda3/bin/python`: Python 3.12.2, NumPy 1.26.4, SciPy 1.16.3, Shapely 2.0.5.
- `/opt/anaconda3/envs/tcra/bin/python`: Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, Shapely 2.1.2.

The macOS system/Homebrew `python3` in this workspace does not currently provide SciPy or Shapely and was not used to claim numerical validation.
