# Emanuel synthetic-event TCR input code audit

Status: **primary-code implementation found; Nature production equivalence unresolved**  
Audit date: 2026-08-13

## Decisive result

The official Emanuel v6.4 MATLAB package contains exact public code for the
warm-core moisture adapter, the synthetic-event shear proxy, storm translation,
time interpolation, and rainfall conversion. This is materially stronger than
inferring these steps from paper prose. However, Xi, Lin & Gori (2023) did not
identify this archive as their production hazard code, and the archive calls an
ER11-family wind profile rather than C15. Hence these equations are an
**official Emanuel public implementation**, not proof of the exact Nature 2023
implementation.

Frozen source and checksums:
`../vendor/Emanuel_TCR/scripts_ver6.4/`.

## 1. Warm-core saturation moisture (`qs900b.m`)

Call site: `raingen.m:141-145` constructs storm datetimes and calls
`qs900b(T600store, vstore)` before passing its output to TCR.

Despite the output variable name `q900`, the code sets `pref=950` hPa
(`qs900b.m:11`). Therefore this version computes saturation mixing ratio at
**950 hPa**, not 900 hPa. Inputs are environmental 600-hPa temperature `T600`
(K) and circular maximum wind `vmax` (knots).

Constants (`qs900b.m:13-22`):

```text
cp = 1005 J kg-1 K-1
Rv = 491 J kg-1 K-1
Rd = 287 J kg-1 K-1
Lv = 2.5e6 J kg-1
pref = 950 hPa
c3 = 1.6/100
vmax_ms = vmax_knots * 1852/3600
```

At 600 hPa (`qs900b.m:23-26`), with `Tc=max(T600-273.15,-50)`:

```text
es600 = 6.112 exp[17.67 Tc/(243.5+Tc)]       (hPa)
q600  = 0.622 es600/(600-es600)              (kg kg-1 mixing ratio)
```

The initial 950-hPa temperature guess is `T=T600+20 K`. Five Newton iterations
are applied (`qs900b.m:29-44`). At each iteration:

```text
es = 6.122 exp[17.67 (T-273.15)/(243.5+T-273.15)]
qs = 0.622 es/(950-es)

F(T) = cp ln(T/T600)
     + Lv (qs/T - q600/T600)
     - Rd ln(950/600)
     - (1.6/100) vmax_ms^2 = 0

F'(T) = [cp T + Lv qs (Lv/(Rv T)-1)]/T^2
T <- T - F/F'
```

The last term is the explicit intensity-dependent warm-core/eyewall entropy
correction. The code comment says 1.6 converts squared surface wind to squared
gradient wind and 100 K represents `Ts-Tt` (`qs900b.m:20-21`). Missing/fill
inputs with `T600<=100 K` return zero (`qs900b.m:35-45`). The output is a mixing
ratio, although code and documentation call it specific humidity. That
variable-definition difference must remain explicit when comparing the code
with papers that use the term specific humidity.

Version boundary: older `qs900.m:6-36` sets `pref=1000 hPa`, has no Vmax/warm-
core term, and uses a different iterative approximation. Zhu et al. (2013)
explicitly said warm-core effects were omitted; v6.4 `raingen.m` instead selects
the newer `qs900b`. These versions must not be conflated.

## 2. Synthetic-event shear proxy (`raingen.m`)

`utrans` supplies smoothed translation components in knots. If event-set fields
`u850store` and `v850store` exist, `raingen.m:56-63` computes:

```text
knotfac = 1852/3600
vdrift  = 1.5 m s-1 * sign(latitude_of_basin)

ush = 5 * knotfac * (ut - u850store)
vsh = 5 * knotfac * (vt - vdrift*cos(latitude) - v850store)
```

Thus this public synthetic adapter does **not** read a direct 200-minus-850 hPa
wind vector. It estimates the vector supplied to the TCR baroclinic term from
storm motion relative to 850-hPa environmental flow, removes a meridional
beta-drift contribution, then multiplies by five. If the 850-hPa fields are
absent, both components remain zero.

This does not establish that Nature 2023 used the same proxy. The Nature paper
does not publish this equation, and its deposited Zenodo code is limited to
statistical/probabilistic analysis and visualization, not hazard generation.

## 3. Translation and temporal interpolation

`utrans.m` assumes native two-hour track points. It calculates centered
differences (`utrans.m:30-40`), repairs Greenwich crossings (`:31-36`), linearly
extrapolates the first and last valid velocities (`:41-48`), and applies one
second-difference smoothing pass with `smfac=0.4` (`:50-53`). Output is knots.

`pointwshortnqdx.m:20-25` sets `nsteps=round(2/timeres)`. Lines 196-227 then
linearly interpolate every required storm variable from the native two-hour
steps: Vmax, RMW, secondary eyewall fields, latitude, longitude, datetime,
translation, shear, and moisture. The package defaults are `timeres=0.5 h` and
`timelength=96 h` (`params.m:44-45`). These are code defaults, not a scientific
requirement that every event must last 48 hours on either side.

### 3.1 RMW input boundary and the superseding fixed-r0 Lin closure

The public v6.4 rainfall scripts consume an RMW series but do not derive one for
the project's Lin v1.1 tracks; those files contain `v_trks` and `vmax_trks`, but
neither RMW nor outer radius. The subsequent Gori/Chavas source review therefore
supersedes the former time-step Knaff closure: following Gori et al. (2022), each
Lin event receives one lifecycle-fixed outer radius and official C15 `r0input`
derives time-varying RMW from `v_trks`, that `r0`, and `abs(f)`.

The source-constrained public catalogue matches Chavas et al. (2016) global
Table 1 `median/Q1/Q3 = 881.0/740.7/1054.4 km`, giving
`mu=6.78105762593618` and `sigma=0.261776756893889`; NumPy `PCG64`, seed
`20260810`, draws once in frozen `event_position=0..9999` order. It does not
claim Gori's private fit, seed, RNG, or event-to-draw mapping. Lin `v_trks`
supplies both C15 and the official `qs900b.m` warm-core calculation;
`vmax_trks` remains only an event-threshold/catalogue variable. Knaff,
Nederhoff, the project-fit outer-first model, and residual assignment are not
Lin production paths. The exact contract is
`CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md`.

## 4. TCR use of the inputs and disclosed numerical choices

`pointwshortnqdx.m:249-263` calls `windprofiles`, and `params.m:40` defaults to
profile 3, Emanuel & Rotunno (2011). It does not call C15.

The shear enters the baroclinic/entropy-slope component at
`pointwshortnqdx.m:413-416`. Translation and orographic effects are modified
radially at `:418-435`. Crucially, the public code then applies:

```text
w = min(w, 7)                       # pointwshortnqdx.m:436
wq = q_interp * max(w - wrad, 0)    # pointwshortnqdx.m:437
```

The defaults include `wrad=0.005 m s-1`, `eprecip=0.9`, `Htrop=4000 m`,
`deltar=2 km`, and `radcity=300 km` (`params.m:64-73,100`). Finally,
`raingen.m:147-149` converts vapor flux to rain rate and integrates:

```text
rainrate = eprecip * 1000 * 3600 * 0.00117 * wq   (mm h-1)
rain      = timeres * sum(rainrate over time)      (mm)
```

Therefore the earlier statement “no public positive-w cap exists” is false for
this official package. The accurate statement is: **the reviewed papers do not
disclose a positive-w cap, but official Emanuel v6.4 code caps w at 7 m s-1; it
remains unresolved whether Xi et al. (2023) retained this cap.**

The same source exposes additional numerical choices that cannot be inferred
from the papers alone:

- drag values at the radial finite-difference points are clipped to `[0, 0.005]`
  (`pointwshortnqdx.m:296-302`);
- both radial angular-momentum gradients are floored at `10` (`:314-320`;
  repeated for secondary eyewalls at `:376-382`);
- time-dependent radial motion is multiplied by the code's eye/time factors
  (`:321-330`), with an additional vanishing-secondary-eyewall guard at
  `:386-405`;
- the shear/baroclinic term has a concrete code form at `:413-416`;
- `radcity=300 km` produces a continuous 50-km-wide radial reduction of the
  **translation contribution** (`:418-424`), not a total-rainfall hard cutoff;
- the orographic gradients receive a continuous radial weight (`:426-435`),
  rather than the printed Lu `|V|<30 m s-1` switch;
- rainfall conversion uses `0.00117`, not the Lu/Xi paper value `0.0012`.

These are exact v6.4 implementation facts. They are candidates for author
comparison, not evidence that Nature retained them in its separate C15
production chain.

### 4.1 Bundled static fields and exact drag preprocessing

The official ZIP also bundles the actual arrays consumed by this public code:

| File | Variable and shape | Native grid | SHA-256 |
|---|---|---|---|
| `C_Drag500.mat` | `cd`, `1440 x 721`, float64 | 0.25-degree global | `2c9807ecf14fa786014021212461a1449d877f53e5e15e41c3fc5384f0c423d5` |
| `bathymetry.mat` | `bathy`, `1440 x 721`, float64 | 0.25-degree global | `c973d3601eb8078258918272cd3ce7dab8c896b932e57bba3b9c0b7cb7161ea6` |
| `bathymetry_high.mat` | `bathy`, `3600 x 1800`, float32 | 0.1-degree global | `66e2fe30da69806f1316b5a40782a60f539792c4d2c5a96ca85935f1a2123df4` |

`pointwshortnqdx.m:52-57` loads `C_Drag500.mat`, transforms it as
`cd=0.9*cd/(1+50*cd)`, and applies a `1e-3` floor. Lines 63-110 compute spatial
gradients on the 0.25-degree grid and bilinearly interpolate drag and its
gradients. Lines 113-119 then apply the code's coastal response adjustment.
`raingen.m:36-43` selects the 0.1-degree or 0.25-degree bathymetry/topography
array; lines 83-138 compute and interpolate terrain and gradients.

These bundled fields eliminate the need to guess an ERA-Interim parameter ID
when reproducing **Emanuel v6.4 itself**. They do not establish that Xi et al.
(2020) or Nature 2023 used the byte-identical arrays or every v6.4 clipping and
coastal rule; that identity remains an author/source question.

## 5. What remains unpublished / not proven

The following may not be filled in by assumption:

1. Whether Xi et al. (2023) used Emanuel v6.4 `qs900b` exactly, including its
   actual 950-hPa target, five Newton iterations, and constants.
2. Whether their production shear was the v6.4 motion-minus-850 proxy above,
   a direct 200--850 hPa shear, or another supplied downscaling-model field.
3. The exact C15 gradient/surface-wind conversion and C15-to-TCR radial profile
   adapter. The official v6.4 code uses ER11-family profiles.
4. Whether the Nature production run retained `w<=7`, density ratio `0.00117`,
   `dM/dr>=10`, drag floors/clips, time/eye factors, the translation-only
   `radcity=300 km` response, the orographic radial weighting, the 96-hour
   fixed-point window, and other v6.4 defaults.
5. Nature production handling outside C15 outer radius and at C15/TCR spatial
   and temporal derivative boundaries.
6. Gori's production lognormal fit, RNG, seed and event-to-draw sequence; none
   is published with the released Lin output.

These items remain source-identity boundaries, not execution gates. The single
public reconstruction frozen in `C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md`
proceeds without author correspondence and must not be relabelled as the
unpublished Nature production code.

## Primary links

- Official package index: https://texmex.mit.edu/pub/emanuel/ForAlex/
- Official ZIP: https://texmex.mit.edu/pub/emanuel/ForAlex/scripts_ver6.4.zip
- Official script README: https://texmex.mit.edu/pub/emanuel/ForAlex/Readme_matlab_scripts_ver6.4.pdf
- Xi, Lin & Gori (2023): https://doi.org/10.1038/s41558-023-01595-7
- Their deposited analysis code/data: https://doi.org/10.5281/zenodo.7407013
- Zhu et al. (2013): https://doi.org/10.1002/2013GL058284
- Xi & Lin (2022): https://doi.org/10.1029/2022GL099196
- Gori et al. (2022), event-fixed outer radius and C15 storm sets: https://doi.org/10.1038/s41558-021-01272-7
- Chavas et al. (2016), global observed outer-radius distribution: https://doi.org/10.1175/JCLI-D-15-0731.1
