# Nature synthetic-event environment: local and server code audit

Status: `SUPERSEDED_IN_PART_BY_OFFICIAL_EMANUEL_V6_4`  
Audit date: 2026-08-13  
Supersession note: this audit describes the local/server code state. The later discovery of official Emanuel v6.4 source closes a public warm-core/shear implementation candidate; see `EMANUEL_SYNTHETIC_TCR_INPUT_CODE_AUDIT.md`. The later Gori/Chavas source review also replaces this audit's former Knaff RMW closure with the event-fixed `r0` contract. The conclusion that the exact Nature-era production adapter remains unproven is unchanged.
Scope: Nature Climate Change (2023) synthetic-event TCR inputs only: storm-centred lower-tropospheric saturation humidity and environmental vertical wind shear.

## 1. Bottom line

The original local/server audit found two reproducible but distinct pieces of
code; a subsequent source search added a third, more authoritative public
candidate:

1. **Lin v1.1 contains the exact synthetic 250–850-hPa environmental-wind generator used by the downscaling model.** It samples a four-component wind vector from monthly means/covariances and saves hourly `u250`, `v250`, `u850`, and `v850` along every accepted track.
2. **The frozen pyTCR 1.2.3 wheel contains an explicit `T600 + Vmax -> q*` solver and an explicit rainfall-shear adapter.** However, pyTCR 1.2.3 was released in June 2025, after the 2023 Nature paper; the function called `calculate_qs900` actually solves at **950 hPa**; and its documented wind-speed units conflict across caller and callee. This is therefore a high-priority, source-coded candidate for author verification, **not evidence that the Nature production code used the same implementation**.
3. **The later-frozen official Emanuel v6.4 MATLAB package contains the upstream `qs900b`, synthetic shear, translation, interpolation and TCR code choices.** It supersedes pyTCR as the strongest public implementation candidate, but still does not prove Nature production equivalence because its default wind profile is ER11 rather than C15.

No reviewed code establishes the exact Nature-2023 production adapter for:

- environmental `T600` sampled along a synthetic event;
- the intensity-dependent warm-core conversion at the paper's stated near-900-hPa level; or
- the choice between direct Lin `250–850 hPa` shear and pyTCR's reconstructed `5 x (translation - 850-hPa wind)` shear.

Consequently, the old fixed-`q=0.01` event result remains outside the scientific production chain. This audit does **not** authorize a new synthetic-event rain run.

## 2. Frozen artifacts and identities

### 2.1 Project repository

- Server project root: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk`
- Project commit at audit: `1024c3bc8fe8c057631e27cde70ed3fa23ab4a68`
- Project worktree at audit: no tracked modifications reported by `git status --short`.

### 2.2 Lin downscaling source

- Frozen source: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/lin2023_tropical_cyclone_risk`
- Upstream origin: `https://github.com/linjonathan/tropical_cyclone_risk.git`
- Commit: `63a760f5ff95d84e7b4f388c40e83f11bcdc418f`
- Tag: `v1.1`
- Frozen source worktree: clean.
- Patched execution copy: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk/scratch/lin_runtime/v1.1-global-periodic-seed20260809-stream0`
- Runtime provenance: `runtime_manifest.json:150-153` records the same upstream commit and a clean source.
- The runtime patch changes deterministic RNG, global periodic longitude, calendar-month labelling, covariance `ddof`, interpolation and accepted-output buffering. It does not introduce a new humidity-to-rainfall adapter or replace the upstream 250–850-hPa shear definition.

### 2.3 pyTCR

- Frozen wheel: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/vendor/pytcr-1.2.3-py3-none-any.whl`
- SHA-256: `d09245629fbe704d76cbbc36729692804e6540fc6b2eb8ee557b0f0446602c69`
- Size: `60,575,952` bytes
- Frozen patched runtime: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/pytcr-runtimes/pyTCR-1.2.3-shared-cle-snapshot-v1`
- Runtime manifest: `runtime_manifest.json:15-21` records pyTCR `1.2.3`, the same wheel SHA/size, and patch scope `instantaneous rainfall primary radial profile only`.
- Official PyPI release JSON reports the same wheel SHA/size and upload time `2025-06-11T19:57:48.852243Z`: <https://pypi.org/pypi/pyTCR/1.2.3/json>.
- The public pyTCR repository was created in September 2024, and the JOSS software paper was published in June 2025: <https://github.com/levuvietphong/pyTCR>, <https://doi.org/10.21105/joss.08074>.

Therefore, this public pyTCR artifact postdates Xi, Lin and Gori (2023). Its code may encode the same scientific idea, but its existence cannot prove byte-for-byte or equation-for-equation identity with the Nature production implementation.

### 2.4 tcwindfields

- Frozen wheel: `/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/vendor/tcwindfields-0.2.0-py3-none-any.whl`
- SHA-256: `000564f8ed15f2e8d937ab41fdf76cce0a56fe181a3a2c85c39fdf5e0a4aa175`
- A scan of every Python member of this wheel found no `Emanuel 2017`, `Gori 2022`, `T600`, `q900`, warm-core, or environmental-shear implementation. It provides the wind-field/CLE layer, not the missing synthetic humidity adapter.

## 3. Exact Lin synthetic environmental-wind chain

All paths in this section are relative to:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/scratch/lin_runtime/v1.1-global-periodic-seed20260809-stream0`

The corresponding scientific definitions are inherited from frozen Lin v1.1; runtime-only differences are recorded above.

### 3.1 Variables and ordering

`namelist.py:65-71` freezes the steering levels as `[250, 850]` hPa.  
`track/env_wind.py:22-26` expands these into this ordered mean vector:

```text
[ua250_Mean, va250_Mean, ua850_Mean, va850_Mean]
```

The Lin README independently states the output meanings and units:

- `README.md:93`: `u250_trks`, 250-hPa environmental zonal wind, m s-1;
- `README.md:94`: `v250_trks`, 250-hPa environmental meridional wind, m s-1;
- `README.md:95`: `u850_trks`, 850-hPa environmental zonal wind, m s-1;
- `README.md:96`: `v850_trks`, 850-hPa environmental meridional wind, m s-1.

### 3.2 Monthly source statistics

`track/env_wind.py:173-217` implements the statistics:

- lines 184-188 choose `250/850 hPa` or `25000/85000 Pa` from the pressure-coordinate units;
- lines 190-199 reduce subdaily inputs to daily means when that branch is triggered;
- lines 201-208 select `ua250`, `va250`, `ua850`, and `va850`;
- lines 211-217 calculate complete-calendar-month means, variances and pairwise covariances (`ddof=1` in the patched runtime).

The accepted environment manifest confirms this production configuration:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/runs/lin_environment/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/manifest.json`

- lines 14-22: `p_midlevel=60000 Pa`, steering levels `250/850 hPa`;
- lines 37-40: each complete-calendar-month statistic is labelled at day 15;
- lines 106-134: accepted wind artifact, 240 months, 14 statistics, SHA-256 `e545a6d6b08b9a3bf5ebeec1c0e4621c4efab0f9795076ffe81b7f89119f7cd3`.

### 3.3 Synthetic time evolution along a track

`track/bam_track.py:18-30` constructs, independently for each wind component, a 15-term random Fourier series with spectral amplitude proportional to `n^-1.5`.  
`track/bam_track.py:103-106` creates that series for the four wind components.  
`track/bam_track.py:108-123` produces the environmental vector at a track point:

```text
env_winds = monthly_spatial_mean + cholesky(monthly_covariance) @ F(time)
```

The returned values have the four-component order listed in section 3.1 and units of m s-1.

### 3.4 Lin FAST shear

`track/env_wind.py:44-55` extracts the four components.  
`intensity/coupled_fast.py:113-122` defines the FAST environmental shear as:

```text
S_vector = [u250 - u850, v250 - v850]       (m s-1)
S = norm(S_vector)                           (m s-1)
```

`intensity/coupled_fast.py:124-131` multiplies the magnitude by the normalized midlevel saturation-entropy deficit for FAST ventilation. This shear is part of the Lin intensity integration; the code does not state that the same vector is passed to TCR rainfall.

### 3.5 Saved track fields

`util/compute.py:244-263` reconstructs the accepted event's environmental winds and commits them only after the event passes the Lin maximum-wind gate.  
`util/compute.py:312-320` writes the four hourly environmental-wind components, `v_trks`, and `vmax_trks` to the track NetCDF.

The accepted 100,000-track file is:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/data/lin/tracks/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/GLx5000peryear_stream0/output/tracks_GL_MPI-ESM1-2-LR_historical_r1i1p1f1_199501_201412.nc`

It contains `100000 x 361` values for each of `u250_trks`, `v250_trks`, `u850_trks`, and `v850_trks`; the time coordinate is 0 to 1,296,000 s at 3,600-s spacing. The NetCDF variables themselves have no attached unit attributes, so the units depend on the frozen Lin README/code contract above.

## 4. Exact `T600 + Vmax -> q*` code found in pyTCR 1.2.3

### 4.1 Function and units as written

File:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/pytcr-runtimes/pyTCR-1.2.3-shared-cle-snapshot-v1/tcr/terrain_boundary.py`

Function: `calculate_qs900(T600, vmax)` at lines `289-350`.

The function's own contract is:

- `T600`: temperature at 600 hPa, K (`lines 296-297`);
- `vmax`: maximum wind speed, knots (`lines 298-299`);
- first return: named `q900`, documented as saturation humidity at **950 hPa**, g g-1 (`lines 301-305`).

This is not a transcription error in this audit: `line 310` explicitly sets `pref = 950 hPa`. The public source therefore has a semantic mismatch between the function/variable name (`qs900`, `q900`) and its actual pressure (`950 hPa`).

### 4.2 Constants and exact equations

`terrain_boundary.py:309-315` uses:

```text
pref = 950 hPa
cp   = 1005 J kg-1 K-1
Rv   = 491 J kg-1 K-1
Rd   = 287 J kg-1 K-1
Lv   = 2.5e6 J kg-1
```

At lines `317-321`, for pressure `p` in hPa and temperature `T` in K:

```text
Tc = max(T - 273.15, -50 degC)
es = 6.112 * exp(17.67 * Tc / (243.5 + Tc))       [hPa]
q  = 0.622 * es / (p - es)                        [code-labelled g/g]
```

At lines `323-345`, the code converts `vmax` from knots to m s-1 and performs exactly five Newton updates, initialized by `T = T600 + 20 K`, to solve:

```text
cp * ln(T / T600)
+ Lv * (q950(T) / T - q600(T600) / T600)
- Rd * ln(950 / 600)
- 0.016 * vmax_ms^2
= 0
```

The intensity term is the source-coded warm-core contribution. Only elements with `T600 > 100 K` are updated (`lines 334-345`); other returned values remain zero (`lines 347-350`).

### 4.3 Caller and unresolved unit conflict

File:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/pytcr-runtimes/pyTCR-1.2.3-shared-cle-snapshot-v1/tcr/rainfall.py`

`generate_rainfall_point` is defined at lines `343-348`. It:

- documents its `velocity` input as m s-1 at lines `364-367`;
- documents `T600` in K at lines `383-384`;
- passes `velocity` unchanged to `calculate_qs900` at lines `459-463`.

The callee, however, documents `vmax` as knots and multiplies it by `1852/3600` internally. Thus the public caller and callee have incompatible documented velocity units. The frozen code contains no assertion or conversion at their boundary. This conflict must be resolved from the intended production interface or by author confirmation before the formula can be frozen for this study.

### 4.4 Limited call coverage

The `T600` branch is used only by `generate_rainfall_point` (`rainfall.py:459-468`). The two gridded APIs used by the old event work do not accept event-varying `T600`:

- `calculate_rainfall_rate`, `rainfall.py:15-20`, defaults to scalar `q900=0.01`;
- `calculate_etr_swath`, `rainfall.py:191-196`, also defaults to scalar `q900=0.01`.

Therefore, finding `calculate_qs900` does not close the production chain. A source-faithful gridded adapter, with verified pressure and velocity units, is still absent.

## 5. Exact rainfall-shear code found in pyTCR 1.2.3

### 5.1 Adapter

In `tcr/rainfall.py:430-440`, `generate_rainfall_point` constructs the two shear components as follows:

```text
vdrift_knots = 1.5 * hemisphere_sign
ush_ms = 5 * knots_to_ms * (utrans_knots - u850_knots)
vsh_ms = 5 * knots_to_ms * (
    vtrans_knots - vdrift_knots * cos(latitude) - v850_knots
)
```

The instantaneous gridded path contains the same relationship at `tcr/rainfall.py:157-163`; the event-total gridded path contains it at `tcr/rainfall.py:283-294`.

This is **not** a direct subtraction of the stored Lin `u250/v250` and `u850/v850`. It reconstructs a shear-like vector from storm translation, the 850-hPa flow, a beta-drift correction in the meridional component, and the factor 5.

### 5.2 Entry into the TCR vertical-motion equation

File:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/pytcr-runtimes/pyTCR-1.2.3-shared-cle-snapshot-v1/tcr/wind.py`

- `calculate_upward_velocity_field`, lines `708-750`, declares `us` and `vs` as vertical-shear components in m s-1;
- lines `803-816` broadcast those components to the TCR grid/time stencil;
- lines `833-847` use them in the baroclinic vertical-motion contribution;
- the time-series path applies the corresponding term at lines `1316-1337`.

This proves the formula's role in pyTCR 1.2.3. It does not prove that the 2023 Nature calculation used this later public adapter.

## 6. What the old server event runner actually did

File:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/code/run_shared_cle_pytcr_event.py`

The runner:

- reads only `u850_trks` and `v850_trks` from the Lin track at lines `526-535`; it does not read `u250_trks` or `v250_trks`;
- reproduces pyTCR's factor-5 translation/850-hPa adapter at lines `721-739`;
- multiplies the vertical velocity by the scalar `tcr_params.q900` at lines `748-754`;
- `tcr/parameters.py:14-18` fixes that default to `q900 = 0.01` and `Htrop = 4000 m`.

So the old `2929.97 mm` event result was generated by a fixed-humidity calculation, not by the `T600 + intensity-dependent warm core` chain. The factor-5 shear relationship was inherited from pyTCR 1.2.3 rather than invented inside the runner, but its post-2023 provenance still does not establish Nature-method identity.

The runner's `Gori (2022)` reference at `run_shared_cle_pytcr_event.py:312-367` and `998-1003` concerns the translation-vector rotation in the wind-field asymmetry. It does not provide the missing humidity formula or authenticate the TCR environmental-shear adapter.

## 7. T600 data availability is not an existing adapter

The accepted prepared GCM temperature file is:

`/mnt/sdb_test/tang/zengjun/TC_Road_Risk/data/cmip6/derived/lin/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/input/MPI-ESM1-2-LR_historical_r1i1p1f1_ta_1995-2014.nc`

Its `ta` variable is monthly mean air temperature in K on 19 pressure levels, including 600 hPa, for 240 months on a `96 x 192` grid.

Lin uses that field internally for thermodynamics:

- frozen source `thermo/calc_thermo.py:45-74` extracts midlevel temperature and humidity near `p_midlevel=60000 Pa`;
- `thermo/calc_thermo.py:107-116` saves only `vmax`, `chi`, and `rh_mid`.

It does **not** save event-centred `T600` to the synthetic track file. No searched project function samples monthly `ta600` along each accepted synthetic track and passes it to pyTCR's `calculate_qs900`. Having the source field on disk is therefore necessary but not sufficient evidence of the Nature adapter.

## 8. Search coverage and negative evidence

### 8.1 Local TCRoad

Before the Emanuel v6.4 archive was added to `vendor/`, the local project was searched with case-insensitive combinations of:

```text
Emanuel 2017; Gori 2022; q900; qs900; T600; warm core;
saturation; specific humidity; moist adiabat;
u250; u200; u850; v250; v200; v850; wind shear; shear
```

No executable humidity or shear adapter was found in that pre-v6.4 snapshot. The later official archive is audited separately in `EMANUEL_SYNTHETIC_TCR_INPUT_CODE_AUDIT.md`; this historical negative result must not be read as overriding that discovery. The earlier hits were method-boundary documents, especially:

- `methods/C15_TCRM_SOURCE_PARAMETER_MATRIX.md:50-59`;
- `methods/LU2018_TCR_EQUATION_CONTRACT.md:168-190`;
- `methods/NATURE_C15_TCR_REPRODUCTION_CONTRACT.md:46-53`;
- `methods/C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md`.

### 8.2 Server project and dependencies

The same terms were searched across:

- root tracked source under `/mnt/sdb_test/tang/zengjun/TC_Road_Risk`;
- every tracked text file in frozen Lin v1.1;
- the patched Lin runtime;
- every `.py` member of `pytcr-1.2.3-py3-none-any.whl`;
- every `.py` member of `tcwindfields-0.2.0-py3-none-any.whl`.

Within the original local/server search scope, the only executable `T600 -> q*` hit was pyTCR 1.2.3's `calculate_qs900`; the only complete synthetic-wind generator was Lin v1.1; and the only gridded rainfall-shear adapter was pyTCR's factor-5 translation/850-hPa relationship. The subsequently frozen Emanuel v6.4 archive supersedes that narrow code-candidate finding, but still contains no Nature-2023-labelled or Gori-2022-labelled C15 production adapter.

Representative reproducible read-only commands:

```bash
git -C /mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/lin2023_tropical_cyclone_risk \
  grep -n -i -E 'emanuel|gori|q.?900|T600|warm.?core|saturation|shear|u250|u200|u850'

unzip -p /mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/vendor/pytcr-1.2.3-py3-none-any.whl \
  tcr/terrain_boundary.py | nl -ba

unzip -p /mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/vendor/pytcr-1.2.3-py3-none-any.whl \
  tcr/rainfall.py | nl -ba
```

## 9. Decision for the paper workflow

### Reusable now

- The accepted hourly Lin environmental vector and its 250–850-hPa shear can be reproduced exactly from the frozen source.
- The pyTCR 1.2.3 `T600 + intensity` equation and factor-5 shear adapter can be reproduced exactly as a **2025 software implementation**.
- The official Emanuel v6.4 warm-core, shear, translation and numerical candidates can be reproduced exactly as that public code version.
- The prepared monthly GCM `ta600` source field exists with a frozen manifest and hash.
- Lin source semantics separate `v_trks` (circular/azimuthal-mean maximum wind) from `vmax_trks` (total maximum intensity). The frozen public reconstruction uses `v_trks` for official C15 `r0input` and the actual-950-hPa `qs900b.m` warm-core term; `vmax_trks` remains only an event-threshold/catalogue variable and does not enter the scale closure.
- Following the Gori et al. (2022) architecture, every Lin event receives one lifecycle-fixed `r0`. The public reconstruction matches Chavas et al. (2016) global Table 1 quartiles with `mu=6.78105762593618` and `sigma=0.261776756893889`, then draws once per frozen event position with NumPy `PCG64`, seed `20260810`. C15 `r0input` derives time-varying RMW. Knaff/Nederhoff and the project-fit outer-first model remain audit history only.

### Not yet authorized as Nature-2023-equivalent

- treating the pyTCR function's 950-hPa output as the paper's `q*900` without confirmation;
- choosing knots or m s-1 for the warm-core intensity term despite the caller/callee conflict;
- presenting the public reconstruction's `v_trks` warm-core choice as confirmed byte-identical Nature production behaviour;
- presenting the quantile-matched fixed-`r0` catalogue as Gori's private fit, seed, RNG or event-to-draw sequence;
- replacing the pyTCR factor-5 shear adapter with direct `u250-u850`, or vice versa, without confirming the production method;
- inventing a monthly-`T600` track sampler/interpolation rule and presenting it as the Nature method.

### Unpublished Nature identity points (not execution gates)

1. Was the target lower-tropospheric pressure 900, 925, or 950 hPa?
2. Is the source-coded warm-core equation above the equation used in Xi, Lin and Gori (2023), including all constants and five Newton iterations?
3. Whether Nature used the same circular-wind field and SI-unit conversion as the frozen public `qs900b.m` reconstruction.
4. Did the Nature TCR calculation use direct `250–850-hPa` synthetic shear or the factor-5 translation/850-hPa reconstruction, including the beta-drift correction?
5. How was monthly environmental `T600` sampled/interpolated along each hourly synthetic track?

Those points limit only a `source-identical Nature reproduction` claim. They do not block the single method-faithful public reconstruction already frozen in the C15–CLIMADA and Chavas fixed-`r0` contracts: the next computation is the fixed-`r0` rerun of the first Lin synthetic event, with no email dependency and no parallel physics branch.
