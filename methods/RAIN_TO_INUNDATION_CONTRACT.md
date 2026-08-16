# 降水转淹没合同（方法冻结，今夜不跑 SFINCS）

Status: **method frozen 2026-08-16**. No SFINCS production, no bathtub depth cube, no global flood layer tonight.

This note answers: we already invert **rain** on a moving 0.05° disk. How does that become **water depth on a road**, and what must not be faked?

---

## Decision (plain language)

Rain on a 5 km cell is **not** flood depth on a road.

**Intended next physical step:** a 2-D reduced-physics inundation model of the SFINCS class (Leijnse, van Ormondt, Nederhoff & van Dongeren, 2021, *Coastal Engineering*), forced by the rain fields we already store, plus topography, and — where the track is coastal — surge.

**Not tonight:** no SFINCS run, no LISFLOOD-FP production, no “rain × runoff coefficient = depth” layer, no 0.1° bathtub.

Until that model exists, road **asset** dollars can be computed, and wind can be overlaid, but **flood damage and flood-driven user delay are not to be published as if depth were known**.

---

## What the hazard worker actually keeps

Each historical event in `lin_road_domain_300km_v1` writes a compact footprint on the moving 300 km / 0.05° disk:

| Field | Meaning | Use for inundation |
|---|---|---|
| Event-total rain (mm) | TCR-integrated storm rainfall | Pluvial water **volume** available to pond / run off |
| Maximum 24 h rain (mm) | Peak daily total inside the event | Intensity control; design-storm analogue |
| Event-maximum near-surface wind | Model-native wind, not a 10-minute WMO claim | Wind damage; **not** a depth |

Hourly 3-D rain cubes were deleted on purpose (storage). Depth does **not** require keeping every hour of every event on disk. SFINCS (or any 2-D inundation model) needs:

1. a hyetograph **or** a conservative event-total / 24 h forcing,
2. a DEM,
3. roughness / drainage,
4. coastal water-level boundary if the event is on the shelf.

Those can be rebuilt from the compact fields plus the already-frozen track and TCR settings. Re-running TCR for a subset of landfalling events is cheaper than storing 100k hourly cubes.

Koks et al. (2019) tropical-cyclone layer was **GAR 2015 wind only**. They said rainfall flood was in the separate flood maps, not in the cyclone wind field. This study is not repeating that shortcut: rain is inverted, but rain is still not depth.

---

## Why a 2-D inundation model, and why SFINCS-class

Tropical-cyclone water on roads is **compound**:

- pluvial (TCR rain falling on the catchment),
- fluvial (rivers already high),
- coastal (surge and waves) on low coasts.

A bathtub (“every cell deeper than X is wet if it is below the rain-equivalent head”) ignores drainage, channels, and barriers. van Ginkel et al. (2021) showed that even in Europe, object-based road loss needs a real inundation map, not a land-use percentage. A bathtub on 0.05° (~5 km) would smear water across ridges and miss underpasses.

SFINCS is the intended engine because it is:

- open (Deltares),
- built for **compound** coastal / pluvial / fluvial flooding,
- reduced-physics, so a regional domain is computationally plausible after the 100k wind+rain batch,
- already used for tropical-cyclone compound flood studies (Leijnse et al. 2021 and subsequent coastal applications).

LISFLOOD-FP / CaMa-Flood remain acceptable **sensitivity** engines later. They are not the reason to delay choosing a class.

---

## How rain will enter SFINCS (when we implement it)

1. **Pilot first**, not 99k events. One basin we already have as an OSM inset (Gulf Coast **or** Pearl River Delta) × a handful of landfalling historical tracks.
2. Forcing: event-total rain distributed in time with a simple triangular or TCR-consistent hyetograph, **and** a max-24 h intensity cap so we do not dump a 3-day storm in one hour. If a subset of events is re-run, hourly TCR rain can replace that schematic.
3. DEM: the same public topography already used for TCR (`topography_land_360as`) is **too coarse for road depth**. The inundation pilot must use a ~30–90 m public DEM (FABDEM / Copernicus DEM), not the 360-arc-second TCR static.
4. Roads: sample **depth on the OSM way**, not on the 0.1° density cube. Bridges tagged in OSM are not inundated as pavement (van Ginkel).
5. Damage: depth → fraction of replacement cost with a road-specific curve (van Ginkel C1–C6 in Europe; Huizinga transport curve only as a named sensitivity). User delay: Pregnolato et al. (2017) depth–speed function on the assigned WorldOD flow, once that flow exists.

---

## Alternatives rejected as the *main* depth method

| Method | Why it is not the main case |
|---|---|
| Depth = event-total rain (mm) | Confuses rainfall height with flood depth. A 200 mm storm is not 20 cm of water on every road. |
| 0.1° bathtub from rain | Wrong geometry; ignores drainage and elevation barriers. |
| Koks-style “use GAR flood maps instead of our rain” | Breaks the event coupling we are building (same track → wind + rain + flood). |
| Full 3-D / Delft3D / SCHISM on 100k events | Physically richer, not finishable for this paper’s event set. |
| Starting SFINCS tonight on the whole historical set | Would fight the in-flight 192-way wind+rain workers and produce an unreviewed cube. |

---

## What is explicitly not claimed tonight

- No inundation NetCDF is produced.
- Compact rain fields remain **rain**.
- Wind+rain overlay on the 0.1° road-density cube remains **descriptive exposure**, not flood loss.
- Future SSP tracks, once finished, will reuse this same inundation recipe; they do not get a different depth shortcut.
