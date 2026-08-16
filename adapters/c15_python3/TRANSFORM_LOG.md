# CLE15 v1.0 Python 3 adapter transform log

## Source boundary

- Frozen Python source: `../../vendor/CLE15/10.4231_CZ4P-D448/v1.0/original/CLE15_2020-06-23.py`.
- Authoritative same-release MATLAB source for the one corrected branch: `../../vendor/CLE15/10.4231_CZ4P-D448/v1.0/extracted/CLE15_windprofile_PUBLIC_2020-06-23/mfiles/ER11E04/ER11E04_nondim_rmaxinput.m`.
- Official numerical figures: `CLE15_plot_r0input.pdf`, `CLE15_plot_rmaxinput.pdf`, and `CLE15_plot_rfitinput.pdf` in the extracted official release.
- Neither `vendor/original` nor `vendor/extracted` was modified.

## Mechanical compatibility edits

The following edits were applied to the frozen official Python file and nothing else was intentionally changed:

1. Expanded 333 tab characters on 253 lines to eight-column spaces so Python 3 parses the original indentation consistently.
2. Changed four Python 2 `print` statements to function calls.
3. Wrapped six `zip` iterators in `list(...)`, preserving Python 2 `zip` materialisation for Shapely `LineString` construction.
4. Replaced ten removed `np.float(...)` aliases with builtin `float(...)` and one removed `np.int(...)` alias with builtin `int(...)`.
5. Replaced three WKT-string emptiness checks and one inverse check with Shapely 2 `is_empty`.
6. Replaced six WKT-prefix geometry-type checks with `geom_type` and three `MultiPoint[0]` operations with `MultiPoint.geoms[0]`.
7. Corrected the shebang and added a provenance docstring. No equation or constant is changed by this edit.

## MATLAB-authoritative correction to the separate Python file

The official Python `ER11E04_nondim_rmaxinput` contained a fallback that returned the raw, unmerged ER11 profile whenever the *last* bisection trial had no intersection. The same-release official MATLAB function has no such fallback: after bisection it always constructs the profile from the most recently stored ER11--E04 merge point.

The adapter disables only that Python-only fallback and executes the MATLAB branch. This is not a fitted correction and introduces no new equation, parameter, grid, threshold, or iteration. It is required for the default `rmax` case to reproduce the same-release official PDF.

## Explicit non-edits

- No physical equation, empirical coefficient, default parameter, grid spacing, convergence threshold, or loop bound was tuned.
- The official Python interpolation choices were not replaced wholesale by a fresh implementation.
- The official Python's accepted but inactive `eye_adj` block remains inactive. Therefore only the official default `eye_adj=0` is validated here.
- The official Python return order for `ER11E04_nondim_r0input` remains `(rr, VV, rmerge, Vmerge, rmax)`, even though the MATLAB entry uses `(rr, VV, rmax, rmerge, Vmerge)`.
- Runtime warnings caused by evaluating formulas at radius zero remain; the official code subsequently sets wind at radius zero to zero.

## Validation status

The three official default MATLAB scripts supplied all inputs. Vector blue curves in the three official PDFs supplied independent wind-speed anchors. Across 38 anchors, the largest absolute deviations are:

| Entry | Anchors | Maximum absolute deviation |
|---|---:|---:|
| `r0input` | 12 | 0.000127 m s-1 |
| `rmaxinput` | 13 | 0.001609 m s-1 |
| `rfitinput` | 13 | 0.000722 m s-1 |

The committed test tolerance is 0.003 m s-1, reflecting PDF/SVG coordinate quantisation rather than a model tolerance.

This closes the three *official default examples*. It does not prove pointwise equivalence over the full parameter domain, and it does not validate `eye_adj=1`.
