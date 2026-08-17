# TC–road risk: what this study is, and where it stands (16 August 2026)

This is the public status note for the tropical-cyclone / urban-road work.
It is **not** a results paper. Numbers below for live jobs were read from the
production host `tc-road-risk` (`3090-2`) at **2026-08-16 15:58 UTC**.

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

| Item | 16 Aug 2026 15:58 UTC |
|---|---|
| Domain | 99,242 historical tracks that come within 300 km of a motor road (of 100,000) |
| Compact hazard footprints | **83,536** |
| Road-overlap files | **83,536** |
| Event-side files | **83,687** |
| Locks held | ~83,869 |
| Workers still running | ~386 `run_lin_event_worker.py` processes |
| Remaining (approx.) | ~15,700 events |
| Progress | **~84%** of the road-domain sample |
| Audit queue | **3** events, all `METHOD_DOMAIN_PENDING` (`C15_FIXED_R0_BELOW_302KM_NUMERICAL_DOMAIN`): 12357 (r0=297.2 km), 68925 (278.4 km), 72126 (279.3 km). Not resampled. |
| Other live fail tally | Full 83k-attempt scan not finished in this session; successful attempts use `status=completed`. Antimeridian leftovers stay a post-batch audit. |
| In-flight worker change | **none** — do not touch the 192-way launcher |

Compact files are still being written (timestamps 15:57 UTC). This batch is
the current-climate 1995–2014, 5,000 accepted tracks/year, stream 0.

Antimeridian / method-domain failures remain a post-batch audit (the three
audit-queue IDs plus any lock-without-footprint remainder). They are **not**
being re-sampled mid-flight.

### 2.2 Future Lin windows (same GCM, four SSPs × two 20-year slices)

All **eight** future environment files are published (plus historical):

- ssp126 / 245 / 370 / 585 × 2041–2060 and 2081–2100
- each `env_wnd_*.nc` is 495,477,475 bytes

Tracks (5,000/year, stream 0), sequential in tmux `lin_future_windows`:

| Window | Tracks | Notes |
|---|---|---|
| historical 1995–2014 | published (2.61 GB) | baseline |
| ssp126 2041–2060 | published 14 Aug 21:53 UTC | |
| ssp126 2081–2100 | published 15 Aug 09:55 UTC | |
| ssp245 2041–2060 | published 15 Aug 21:20 UTC | |
| ssp245 2081–2100 | published 16 Aug 08:43 UTC | |
| **ssp370 2041–2060** | **running** (started 16 Aug 08:43 UTC) | directory exists; no final nc yet |
| ssp370 2081–2100 | not started | env ready |
| ssp585 2041–2060 | not started | env ready |
| ssp585 2081–2100 | not started | env ready |

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

A full 113-million-way dollar layer still needs a server extract of
`roads.osm.pbf`. That extract is **pending**; it is not required to freeze
the method, and it must not interrupt the wind+rain workers.

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

## 6. What is deliberately not being done tonight

- No change to the in-flight historical hazard launcher.
- No re-extract of `planet-260803`.
- No global traffic assignment, no MRIO.
- No SFINCS production.
- No Koks Supplementary Table 8 unit costs.
- No observed AADT as the flow layer.

---

## 7. How to reproduce the dollar kernel locally

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
