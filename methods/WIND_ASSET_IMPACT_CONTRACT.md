# Wind asset-impact contract (main case, frozen)

Status: **main case frozen 2026-08-18**.
This is how tropical-cyclone **wind** enters object-level road **asset** dollars.
It follows Koks et al. (2019). The operational definition is the public
`gmtra` implementation (`ElcoK/gmtra`, `damage.py` / `parallel.py`).

Do not invent a pavement wind MDR. Do not multiply ROCKS / national
replacement cost by a wind damage ratio. Do not use Koks Supplementary
Table 1 cleanup bands (5k–50k USD/km and the rest) as priced dollars.
Do not apply 150 or 151 km/h as a raw number on unlabeled C15 wind.
Do not write a 0.1° wind-loss grid.

Wind-driven **user** delay is not in this contract.

---

## 1. What we are computing

Two separate wind-asset ledgers, both in **2025 USD**, both on the OSM
object (way or bridge). They are **not** reconstruction of the carriageway.

1. **Road cleanup.** Tree debris removal and minor repair on motor-road
   ways, following Koks: wind above the tree-break threshold, scaled by
   Crowther tree density as in `gmtra`.
2. **Bridge collapse.** Full replacement of a bridge when the event gust
   exceeds the `gmtra` failure threshold **and** the event is rarer than
   the `gmtra` design return period.

Flood reconstruction (van Ginkel × the frozen 2025 unit-cost book) stays
on its own ledger after SFINCS depth exists. Do not mix the three.

Tunnels take no wind-cleanup and no wind-collapse dollars.

---

## 2. Method source

**Following Koks et al. (2019), *Nat. Commun.*** Tropical-cyclone wind
does not destroy the pavement. Road dollars are cleanup of fallen trees
and minor repair. Bridges can collapse at extreme gusts if they are also
beyond their design return period.

We apply that method with the rules in the public `gmtra` code that
produced their calculation.

Sources (for the contract, not for a methods polemic):

- Koks et al. (2019), DOI `10.1038/s41467-019-10442-3`, and SI Table 1.
- `https://github.com/ElcoK/gmtra` — `gmtra/damage.py`, `gmtra/parallel.py`.
- Crowther et al. (2015), *Nature* **525**, 201–205, DOI `10.1038/nature14967`.
- Escobedo et al. (2009), *Arboriculture & Urban Forestry* **35**, 100–106.
- Harper, Kepert & Ginger (2010), WMO/TD-No. 1555, Table 1.1.
- Virot et al. (2016), *Phys. Rev. E* **93**, 023001 (42 m/s = 151 km/h).
- GIRI bridge unit cost already frozen in
  [`ROAD_ASSET_VALUATION_CONTRACT.md`](ROAD_ASSET_VALUATION_CONTRACT.md).

---

## 3. Wind metric

Koks / GAR 2015 thresholds are **peak 3-s gusts** (km/h).

This study’s C15 field is stored as **model-native near-surface wind**
([`C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md`](C15_CLIMADA_OPEN_REPRODUCTION_CONTRACT.md)).
The Nature-family hazard is a 10-minute sustained wind. Roads and bridges
sit on land.

**Frozen conversion (threshold only; do not rename C15 files):**

Harper et al. (2010) Table 1.1, **In-Land**, \(T_0=600\) s, \(\tau=3\) s:

\[
G_{3s/10min}=1.66
\]

A Koks / Virot gust threshold \(V_{3s}\) is applied to C15 as
\(V^{\mathrm{C15}*}=V_{3s}/1.66\).

| Role | 3-s gust (Koks) | C15 threshold |
|---|---:|---:|
| Tree-break / cleanup | 151 km/h (42 m/s) | **91.0 km/h (25.3 m/s)** |

Sensitivity (SI only): \(G=1.52\) (Off-Land) and \(G=1.49\) (In-Land 3-s / 1-min).
Do not put the unconverted 151 km/h cut on C15 in the main text.

---

## 4. Road cleanup (main wind-road dollar)

For each motor-road way that is not a tunnel, and each event:

\[
L_{\mathrm{cleanup}}
=
L_{\mathrm{km}}
\times c_{2025}
\times \mathbf{1}\{V^{\mathrm{eq}}_{3s}\ge 151\,\mathrm{km/h}\}
\times P(N)
\]

with the `gmtra` density factor

\[
P(N)=\begin{cases}
0 & N\le 0\\
\min(N,\,10^{4})/10^{4} & N>0
\end{cases}
\]

- \(L_{\mathrm{km}}\): way length already used in valuation.
- \(N\): Crowther et al. (2015) tree density (stems km⁻², DBH ≥ 10 cm)
  at the way. Use the Yale EliScholar **biome** map, WGS84 GeoTIFF
  (`Revision_01` if that is the file we hash). Sample the pixel that
  contains the way’s representative point (same point as valuation).
- \(P(N)\): `damage.py` `regional_cyclone` — drop non-positive density,
  cap at 10 000 km⁻², divide by 10 000.
- \(c_{2025}\): Escobedo et al. (2009) **moderate** debris volume times
  their removal-and-disposal unit price, inflated 2005 → 2025 with the
  same World Bank US GDP deflator already used for valuation
  (`data/valuation/raw/wb_defl_usa.json`).

Escobedo arithmetic (frozen):

| Piece | Value |
|---|---|
| Volume, moderate | 3.40 m³ per 30.5 m of street |
| Volume per km | \(3.40\times 1000/30.5=111.475\) m³ km⁻¹ |
| Unit price | 28.25 USD m⁻³ (2005; abstract; body 28.11 is a footnote) |
| Deflator 2005 | 77.3945446986751 |
| Deflator 2025 | 122.361624831936 |
| Factor | 1.581014 |
| \(c_{2025}\) | **4,979 USD km⁻¹** |

SI band: Escobedo low / high volumes, same unit price and inflator
(864 and 25,583 USD km⁻¹). Do not vary \(c\) by OSM class. Spatial
variation is \(P(N)\), not Koks’s assumed 5k–50k class bands.

OSM class mapping for any Koks class label we still need (bridges)
is SI Table 4: `motorway` / `trunk` → primary, and so on.

---

## 5. Bridge collapse (main-text dollars)

For each OSM bridge on a motor road, and each event:

\[
L_{\mathrm{bridge}}
=
\begin{cases}
L_{\mathrm{br}}\times c_{\mathrm{br},2025}
&
\text{if }V^{\mathrm{eq}}_{3s}>V^{*}_{\mathrm{class}}
\text{ and }RP(V)>RP_{\mathrm{design}}\\
0
&
\text{otherwise.}
\end{cases}
\]

- \(c_{\mathrm{br},2025}\): GIRI 9.84 million USD/km (2023), inflated
  2023 → 2025 with the same US GDP deflator. Collapse is 100% loss, as
  in Koks / `gmtra`.
- \(V^{*}_{\mathrm{class}}\): `parallel.py` `wind_threshs` row mid-points
  (3-s gust), then divide by 1.66 for C15.

| `gmtra` class (SI Table 4) | 3-s \(V^{*}\) (km/h) | C15 \(V^{*}\) (km/h) |
|---|---:|---:|
| primary | 362.5 | 218.4 |
| secondary | 337.5 | 203.3 |
| other | 312.5 | 188.3 |

SI: the four Morris rows `[[400,375,350],[375,350,325],[350,325,300],[350,300,275]]`.

- \(RP_{\mathrm{design}}\): `parallel.py` `design_tables`, by World Bank
  income group and class (`HIC` / `UMC` / else):

| | primary | secondary | other |
|---|---:|---:|---:|
| HIC | 1/200 | 1/100 | 1/50 |
| UMC | 1/100 | 1/50 | 1/20 |
| LMC+LIC | 1/50 | 1/20 | 1/10 |

- \(RP(V)\): empirical exceedance of the converted gust at that bridge
  from the historical-window event set (same Lin / C15 catalogue).
  SI: move the design period one class step up and down.

Do not apply \(P(N)\) to bridges.

---

## 6. What is rejected

| Rejected | Why |
|---|---|
| \(D(v)\times\) 2025 replacement cost on pavement | Not Koks; no published global pavement wind MDR |
| Koks SI cleanup bands as 2018 or 2025 USD/km | SI labels them assumptions; Escobedo does not contain those numbers |
| Raw 151 km/h on C15 | Wrong averaging period |
| Hansen GFC or WorldCover as the main tree field | Different variable; Koks uses Crowther density |
| Bridge dollars only in SI | Main text includes them |
| Wind-flow / closure duration dollars | Separate contract, later |
| Starting this ledger before historical wind+rain is closed | Hazard first |

---

## 7. Implementation

Code: `code/road_wind_asset_impact.py` (dollar kernel) and
`code/road_wind_object_join.py` (compact wind + Crowther onto objects).
Tests: `code/tests/test_road_wind_asset_impact_contract.py`,
`code/tests/test_road_wind_object_join_contract.py`.
Plain-language winds / roads / how:
[`HISTORICAL_WIND_ASSET_LEDGER.md`](HISTORICAL_WIND_ASSET_LEDGER.md).

Valued ways keep `lon`/`lat` when those columns are present. Host valued
shards omitted them; the join fills lon/lat from the planet-260803 extract
by `way_id`. Impact `apply` reads a table that already has `v_c15_ms` and
`tree_dens_km2`. `score-event` / `score-historical` join compact
`event_maximum_near_surface_wind_speed` then call the kernel. Crowther
sampling is `sample_crowther_density` / `sample-trees` / the join batch
sampler; hash the GeoTIFF with `hash-trees` into
`data/impact/crowther.manifest.json`. The global GeoTIFF is not in git.

Historical 300 km wind+rain is closed (99,234 compact files; eight
`METHOD_DOMAIN_PENDING` IDs contribute no wind). Crowther is hashed on
the host. Object-level historical apply is the authorized next production
step.
