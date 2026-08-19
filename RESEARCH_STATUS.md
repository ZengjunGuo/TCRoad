# TC–road risk: what this study is, and where it stands (18 August 2026)

This is the public status note for the tropical-cyclone / urban-road work.
It is **not** a results paper. Hazard file counts below are the live
host inventory on `3090-2`. Crowther is hashed
locally and the GeoTIFF+manifest are on the host (1 157 937 187 bytes).

Large files (planet PBF, event NetCDFs, Lin tracks) stay on the server
`/mnt/sdb_test/tang/zengjun/TC_Road_Risk`. This GitHub tree holds methods,
code, unit-cost books, and tests.

---

## 1. One-paragraph research line

We estimate **direct reconstruction loss** and, later, **user (traffic) loss**
on the global motor-road network when tropical cyclones hit, first in the
current climate and then under the same GCM in four SSPs.

Hazard is a public reconstruction of the Nature / Gori / Xi line: Lin
synthetic tracks in MPI-ESM1-2-LR `r1i1p1f1`, official C15 wind, CLIMADA
TCR rain, compact event fields on a moving 300 km / 0.05° disk. Exposure is
OpenStreetMap motor roads frozen at **3 August 2026** (`planet-260803`).
Replacement cost is object-level 2025 USD. Flow will be WorldOD / GlODGen
commute demand assigned onto that graph, **not** observed AADT. Rain becomes
depth only through a 2-D inundation model (SFINCS-class), which has **not**
been started.

---

## 2. Live production (read from the host, not from an old README)

### 2.1 Historical wind + rain (`lin_road_domain_300km_v1`)

Live host inventory **2026-08-19 02:53 UTC** on `3090-2` (EasyConnect `status=4`, SSH OK):

| Item | 2026-08-19 02:53 UTC |
|---|---|
| Domain | 99,242 historical tracks within 300 km of a motor road (of 100,000) |
| Compact + road-overlap | **99,234** / 99,242 (`compact_nc=99234` `overlap_nc=99234`) |
| Closed leftover | **14472** compact+overlap present (mtime 02:39 UTC). `eq14472` session at shell prompt. |
| Leftover event ids | **8**: 11902, 11944, 12357, 50194, 62311, 68925, 72126, 86977 (all `METHOD_DOMAIN_PENDING`; no compact/overlap) |
| Frozen, not resampled | 12357, 68925, 72126 |
| In-flight hazard | **none**. Historical C15–TCR closed: 99,234 + 8 pending = 99,242. |
| Wind-asset ledger | **in-flight** tmux `windasset` / `score-historical` since 03:30 UTC. Inputs: 99,234 compact files (pending eight have no compact), 114 valued shards + extract lon/lat, Crowther sha256 `1812e5cbb1…`. Totals pending until `historical_wind_asset.summary.json` exists. Record: [`methods/HISTORICAL_WIND_ASSET_LEDGER.md`](methods/HISTORICAL_WIND_ASSET_LEDGER.md). |

This is the current-climate 1995–2014, 5,000 accepted tracks/year, stream 0.

### 2.2 Future Lin windows (same GCM, four SSPs × two 20-year slices)

All **eight** future environment files are published (plus historical):

- ssp126 / 245 / 370 / 585 × 2041–2060 and 2081–2100
- each `env_wnd_*.nc` is 495,477,475 bytes

Tracks (5,000/year, stream 0). Last host-answered pipeline log
(**2026-08-17 19:08 UTC**): **all eight future windows published**.
Not re-listed tonight (SSH down). No future C15–TCR production has been
started; that wait is method, not a missing env file.

| Window | Tracks | Notes |
|---|---|---|
| historical 1995–2014 | published (2.61 GB) | baseline |
| ssp126 2041–2060 | published 14 Aug 21:53 UTC | |
| ssp126 2081–2100 | published 15 Aug 09:55 UTC | |
| ssp245 2041–2060 | published 15 Aug 21:20 UTC | |
| ssp245 2081–2100 | published 16 Aug 08:43 UTC | |
| ssp370 2041–2060 | published 16 Aug 20:47 UTC | |
| ssp370 2081–2100 | published 17 Aug 05:14 UTC | |
| ssp585 2041–2060 | published 17 Aug 12:30 UTC | |
| ssp585 2081–2100 | published 17 Aug 19:08 UTC | |

Dead pane: `lin_s0_spatial_2500_5000` (historical tracks already published
10 Aug). Do not confuse it with live work.

### 2.3 OSM snapshot

- Source planet `planet-260803.osm.pbf` (88 GB, md5 `156085691b8f5cce296e36c35a6ba57b`)
- Filtered roads PBF 27 GB, sha256 `66b58129bfedfca870db95a946dc44e17979f5a5371eb753eef7dc61516c5dd0`
- Date **2026-08-03**. This is current OSM, frozen. We will not re-download
  and we will not rewind the map to 1995/2014.
- 113,198,514 accepted motor ways; 52.96 million km.

---

## 3. Replacement cost (implemented and tested)

Contract: [`methods/ROAD_ASSET_VALUATION_CONTRACT.md`](methods/ROAD_ASSET_VALUATION_CONTRACT.md)

**Replacement cost (2025 USD) = length (km) × unit cost (2025 USD/km)**
on each OSM way. The 0.1° density cube is a map, not a dollar grid.

- Developing countries: ROCKS 2018 **Actual new-build** medians by World
  Bank region, each row inflated with the US GDP deflator to 2025. Gravel
  resurfacing is not a rebuild price. One-row oddities (South Asia
  “2-lane” at $42k/km) are dropped.
- United States, Japan, Australia, Europe, China: **national / European
  Court of Auditors books**, not ROCKS Europe-and-Central-Asia.
- Missing lanes: country × class median, never a world average.
- Links default to 1 lane and are floored at half the mainline price.
- Main totals always include local motor roads **and** a no-local
  sensitivity. Unclassified share > 60% flags a country (GIRI rule).

Code: one script, `code/road_replacement_value.py`. Tests:
`code/tests/test_road_replacement_value_contract.py`.

Object-level 2025 USD layer is **done** on the host extract of
`planet-260803` (`data/valuation/global/global_replacement_value.summary.json`):
**110,822,264** accepted ways; **52.250 million km**;
**$49.230 T** all-motor / **$35.805 T** no-local. Eleven countries trip the
GIRI unclassified-share flag. This is not a 0.1° dollar grid.

---

## 4. Flow (method only tonight)

Contract: [`methods/ROAD_FLOW_ASSIGNMENT_CONTRACT.md`](methods/ROAD_FLOW_ASSIGNMENT_CONTRACT.md)

**Chosen:** WorldOD / GlODGen commute origin–destination demand, then
assignment onto the same OSM graph.

**Rejected as the main layer:** observed AADT.

WorldOD is a zone-to-zone commute matrix, not a traffic count. Indirect
loss needs rerouting; counts on a few rich-country highways cannot do that.
No global OD or assignment is in this repository yet.

---

## 5. Rain → inundation (method only tonight)

Contract: [`methods/RAIN_TO_INUNDATION_CONTRACT.md`](methods/RAIN_TO_INUNDATION_CONTRACT.md)

**Chosen:** 2-D reduced-physics inundation of the **SFINCS** class
(Leijnse et al. 2021), forced by stored event-total and max-24 h rain,
a ~30–90 m DEM, and surge where the track is coastal.

**Rejected as the main depth:** “depth = rainfall (mm)”, 0.1° bathtub,
or Koks’ use of a separate GAR flood map in place of our event rain.

Hourly rain cubes were deleted; they are not required to start SFINCS.
**SFINCS is not started tonight.** Compact rain remains rain. Do not
publish flood-dollar losses until depth exists.

---

## 6. Wind asset loss (method frozen 2026-08-18)

Contract: [`methods/WIND_ASSET_IMPACT_CONTRACT.md`](methods/WIND_ASSET_IMPACT_CONTRACT.md)

**Following Koks et al. (2019).** Operational rules are the public
`gmtra` implementation. Pavement is not rebuilt by wind.

- Cleanup on ways: Escobedo 2009 moderate volume × 2005 USD/m³, inflated
  to 2025, times Crowther density factor \(P=\min(N,10^{4})/10^{4}\),
  if converted gust ≥ 151 km/h (C15 cut = 91.0 km/h using Harper
  In-Land \(G_{3s/10min}=1.66\)).
- Bridges in the main text: full GIRI 2025 replacement if gust exceeds
  the `gmtra` class threshold **and** the event is rarer than the
  `gmtra` design return period.
- Kernel coded: `code/road_wind_asset_impact.py`. Object join:
  `code/road_wind_object_join.py` samples compact event-max C15 at the
  extract lon/lat (periodic, 0.05° cell) and Crowther \(N\), then calls
  the kernel. Pending eight IDs contribute no wind. Written record:
  [`methods/HISTORICAL_WIND_ASSET_LEDGER.md`](methods/HISTORICAL_WIND_ASSET_LEDGER.md).
- Crowther biome WGS84 GeoTIFF sha256
  `1812e5cbb17f91f3a1dfc3033e9cc402bc557ad6ed3827c84ce1fc3f8f05c338`
  (`data/impact/crowther.manifest.json`). Not a 0.1° dollar grid. Not
  Koks SI 5k–50k bands. Historical object-level apply is the live
  production step after the 99,234-event compact set closed.

---

## 7. What is deliberately not being done tonight

- No re-extract of `planet-260803`.
- No global traffic assignment, no MRIO.
- No SFINCS production.
- No Koks Supplementary Table 8 unit costs.
- No observed AADT as the flow layer.
- No future-window C15–TCR.

---

## 8. How to reproduce the dollar kernel locally

```bash
python3 -m unittest tests.test_road_replacement_value_contract
python3 road_replacement_value.py write-book
python3 road_replacement_value.py value \
  ../data/valuation/fixtures/representative_motor_roads.csv \
  --output ../data/valuation/fixtures/representative_motor_roads.valued.csv
```

Expect: US motorway priced from the national book (~$12 million/km),
Philippine primary from ROCKS East Asia, China motorway from the China
book, the bridge above pavement, the footway rejected, and
`replacement_usd_no_local` < `replacement_usd` because a residential
row is present.
