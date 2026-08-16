# Lu et al. (2018) TCR equation contract

Status: `SOURCE_CONTRACT_V2_PUBLIC_RECONSTRUCTION`  
Scope: equations and published parameters required by the Xi/Nature C15–TCR reproduction. This is not executable code.

Primary source: Lu et al. (2018), *Journal of the Atmospheric Sciences* 75, 2337–2358, <https://doi.org/10.1175/JAS-D-17-0264.1>.

## 1. Rainfall conversion

For positive total vertical velocity,

\[
P_{rate}=\epsilon_p\frac{\rho_{air}}{\rho_{liquid}}q_s w.
\]

Published values inherited by Xi et al. (2020):

- precipitation efficiency: \(\epsilon_p=0.9\);
- density ratio: \(\rho_{air}/\rho_{liquid}=0.0012\) in the Lu/Xi paper contract;
- when total \(w\le0\), set rainfall to zero.

The SI result is converted from m s\(^{-1}\) water depth to mm h\(^{-1}\) by the dimensional factor \(3.6\times10^6\). This is a unit conversion, not a calibration parameter. For the Xi historical branch, the public reconstruction explicitly uses the Lu/Xi paper density ratio `0.0012`; it must override the CLIMADA/Emanuel code constant `0.00117`. The positive vertical velocity is capped at `7 m s^-1`, as implemented independently in the official Emanuel v6.4 code and CLIMADA Petals v6.2. These are frozen public-reconstruction choices, not claims that Nature's unpublished source uses identical constants.

## 2. Angular momentum

The physically and dimensionally correct absolute angular momentum is

\[
M=rV+\frac12fr^2,
\]

where \(V\) is the axisymmetric tangential gradient wind. Lu's printed text contains \(\tfrac12fV^2\), which is dimensionally wrong; the implementation must use the standard \(\tfrac12fr^2\) expression and record this transparent correction.

The boundary-layer balance and continuity equations are

\[
u\frac{\partial M}{\partial r}\simeq-r\frac{\partial\tau_\theta}{\partial z},
\]

\[
\frac{\partial w}{\partial z}=-\frac1r\frac{\partial(ru)}{\partial r},
\]

leading to the frictional vertical-velocity structure used below.

## 3. Five vertical-velocity components

The total vertical velocity is

\[
w=w_f+w_h+w_t+w_s+w_r.
\]

### 3.1 Frictional convergence

\[
w_f=-\frac1r\frac{\partial}{\partial r}
\left[r^2\frac{\tau_{\theta s}}{\partial M/\partial r}\right],
\]

with

\[
\tau_{\theta s}=-C_d|\mathbf V|V.
\]

Here \(|\mathbf V|\) is the magnitude of the total horizontal wind, while the final \(V\) is the tangential gradient-wind component. Lu explicitly notes that the relevant stress is based on gradient-level wind rather than standard 10 m wind.

### 3.2 Topographic ascent

\[
w_h=\mathbf V\cdot\nabla h.
\]

Lu sets this component to zero where \(|\mathbf V|<V_{th}\), with

\[
V_{th}=30\ \mathrm{m\,s^{-1}}.
\]

This threshold applies only to the topographic term; it is not an event-acceptance threshold and not a total-vertical-velocity cap.

### 3.3 Vortex time dependence/stretching

\[
u=-\frac{\partial M/\partial t}{\partial M/\partial r},
\]

\[
w_t\simeq H_b\frac1r\frac{\partial}{\partial r}
\left[r\frac{\partial M/\partial t}{\partial M/\partial r}\right].
\]

Lu uses \(H_b=1\) km and notes that theory suggests a value closer to 3 km. Xi et al. (2020) changes it to 4 km. The Xi PDF prints “4000 km”, an obvious unit error; the reproducible setting is recorded as `4000 m` while retaining the printed-source discrepancy in provenance.

### 3.4 Environmental shear/baroclinicity

The full published expression is

\[
w_s\simeq
\frac{g}{c_p(T_s-T_t)(1-\epsilon_p)N^2}
V\left(f+\frac Vr+\frac{\partial V}{\partial r}\right)
(\Delta\mathbf V_e\cdot\mathbf j).
\]

Lu then uses approximate constants

\[
\epsilon_p\simeq0.5,\quad c_p\simeq1000,\quad g\simeq10,
\quad N^2\simeq4\times10^{-4}\ \mathrm{s^{-2}},
\quad T_s-T_t\simeq100\ \mathrm{K}
\]

to write the simplified form

\[
w_s\simeq0.5fV(\mathbf V_{g,200}-\mathbf V_{g,850})\cdot\mathbf j.
\]

For the method-faithful public reconstruction, this term follows the CLIMADA Petals v6.2 implementation in the frozen `tc_rainfield.py` rather than introducing another local derivation. The runner must record that implementation identity and its file SHA. This resolves the executable public route while preserving the evidence boundary: the unpublished Xi/Nature production branch is not claimed to be source-identical.

### 3.5 Radiative cooling

\[
w_r=-0.005\ \mathrm{m\,s^{-1}}.
\]

This is a published fixed setting, even though Lu notes that it lacks a particular theoretical justification.

## 4. Xi et al. (2020) numerical inheritance

Xi et al. (2020) states that parameters otherwise follow Lu, except for:

- \(H_b=4\) km;
- a spatially varying drag coefficient derived from ECMWF 0.25° surface roughness;
- horizontal TCR resolution `0.05° × 0.05°`.

Feldmann et al. (2019) gives

\[
C_D'=\left[\frac{\kappa}{\ln(500/z_0)}\right]^2,
\qquad \kappa=0.35,
\]

\[
C_D=\frac{0.9C_D'}{1+50C_D'}.
\]

The public papers do not specify the exact regridding, coastline and missing-value rules used by Xi's production files. This reconstruction therefore uses CLIMADA's own documented sampling of its package data, rather than fabricating an Xi/ECMWF identity: `c_drag_500.tif` (ERA5 roughness-derived, 0.25°) and `topography_land_360as.tif` (SRTM-derived, 0.1°). Their official URLs and SHA-256 hashes are frozen in [the open reconstruction contract](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md).

## 5. Published environmental branches

### Historical validation in Xi et al. (2020)

At each six-hourly storm centre, Xi et al. first extracts the spatially averaged environmental parameters,

\[
q_{input}=\overline{q_{925}}^{\,r\le200\,km},
\]

\[
\Delta\mathbf V_e=\overline{\mathbf V_{200}-\mathbf V_{850}}^{\,600\le r\le800\,km}.
\]

The resulting six-hourly environmental-parameter series are then linearly interpolated to one hour. This order follows the paper text and must not be replaced by first interpolating the full NCEP fields. These are historical-event validation inputs, not the Nature synthetic-storm environment.

### Synthetic events in Nature (2023)

Nature points to Emanuel (2017) for its synthetic environmental field. Xi & Lin (2022) further states that storm-centred saturation humidity near 900 hPa is derived from environmental 600-hPa temperature and an intensity-dependent storm warm core. Emanuel's official public MATLAB package `scripts_ver6.4` supplies one exact implementation: `qs900b.m` actually targets 950 hPa, includes an intensity term `(1.6/100)Vmax_ms^2`, and performs five Newton iterations; `raingen.m` constructs shear as five times storm motion relative to the 850-hPa flow, with a beta-drift correction. This closes the public Emanuel implementation, but the Nature paper does not identify v6.4 as its production code and v6.4 uses ER11 rather than C15. The source-identical Nature branch therefore remains unverified; the executable public branch is separately frozen in [the C15–CLIMADA contract](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md) and has no external-correspondence prerequisite.

The Lin v1.1 output contains neither outer radius nor RMW. Following Gori et al. (2022), each Lin event therefore receives one outer radius that remains fixed for its lifecycle, and official C15 `r0input` derives time-varying RMW from `v_trks`, that event-level `r0`, and `abs(f)`. The transparent public reconstruction matches the Chavas et al. (2016) global Table 1 quartiles with `mu=6.78105762593618` and `sigma=0.261776756893889`, then draws once per frozen `event_position` with NumPy `PCG64`, seed `20260810`. `v_trks` supplies C15 and `qs900b.m`; `vmax_trks` remains only an event-threshold/catalogue variable. Knaff/Nederhoff, the project outer-first model, and residual assignment do not enter production. This is a source-constrained public closure, not a claim that Gori's private fit or random sequence has been reproduced byte for byte; see [the fixed-r0 contract](CHAVAS2016_FIXED_R0_PUBLIC_RECONSTRUCTION_CONTRACT.md).

## 6. Frozen public numerical boundaries

Published or supported by the frozen public implementation and safe to encode within the cited branch:

- the equations above;
- C15 wind profile;
- \(\epsilon_p=0.9\) and \(H_b=4000\) m;
- Lu/Xi paper density ratio `0.0012` for the Xi 2020 historical-paper branch; the Nature synthetic branch must record whichever frozen public environmental implementation it uses rather than imply source-identical inheritance from Nature;
- \(V_{th}=30\) m s\(^{-1}\), \(w_r=-0.005\) m s\(^{-1}\);
- negative total \(w\) gives zero rainfall;
- event rainfall is the time integral of gridded rain rate;
- a 1-hour time step and `0.05°` grid;
- a moving `300 km` computation domain around each hourly TC centre;
- one lifecycle-fixed `r0` per Lin event, with time-varying RMW returned by official C15 `r0input`;
- `2 km` radial finite differences, `dM/dr>=10 m s^-1`, `Cd>=0.001` and `w<=7 m s^-1`;
- CLIMADA core v6.1.0 plus Petals v6.2.0, with the frozen `tc_rainfield.py` SHA.

Not source-identical to Xi/Nature and therefore not claims this project makes:

- Xi/Nature's unpublished first/last time derivative code;
- Xi's exact ECMWF static-file bytes and production interpolation;
- whether Nature used the official v6.4 warm-core humidity and shear formulas unchanged, including the actual 950-hPa target and input units.
- Gori's unpublished lognormal fit parameters, RNG, seed and event-to-draw sequence.

The frozen public route uses Petals' explicit numerical rules rather than treating them as tunable alternatives. The C15 provider mechanically appends the official definition point `(r0,V=0)` when needed; a query with `r>r0` raises an error. It does not extrapolate a tail, clip `r0`, or add a taper. The C15-to-`nocoriolis` companion preserves the official normalized absolute-angular-momentum shape as specified in [the open reconstruction contract](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md); it is not `Vd=V`.

No reviewed source authorizes tuning the cap, adding rain-rate clipping, imposing a fixed 48 h duration, inventing a C15 outer taper, disabling \(w_t\), or applying Lu's 20 km smoothing of noisy WRF wind profiles to analytic C15. The frozen `300 km` rule is not an arbitrary hard mask: Xi et al. (2023, JAMC, DOI `10.1175/JAMC-D-22-0131.1`) explicitly uses rainfall within 300 km of the moving TC centre, and Petals implements the same moving-domain limit.
