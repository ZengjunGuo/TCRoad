# 道路重建成本合同（主方案，已冻结）

Status: **main case frozen 2026-08-16**.
This is the replacement-cost method that will be implemented and published.
Do not switch back to a 0.1° dollar grid, Koks Table 8 numbers, CPI inflation,
world-average missing tags, or ROCKS prices for the United States / Japan /
Australia / Europe / China.

Traffic and flow are **not** in this contract. They are a later WorldOD problem.

---

## 1. What we are computing

For every motor-road segment in the frozen OpenStreetMap extract:

> **Replacement cost (2025 USD) = length (km) × unit cost (2025 USD/km)**

This is the cost of rebuilding that road to the same class and number of lanes
at 2025 construction prices. It is not market value, not land value, and not
traffic delay.

The 0.1° length-density cube remains a **descriptive map only**. Loss and
asset totals are computed on the road object (the OSM way, or a
class × country piece of that way). Hazard is sampled from the native
0.05° event field onto the way.

---

## 2. Which roads

**In (main case — all motor roads):**

| Group | OSM `highway` values |
|---|---|
| motorway | `motorway`, `motorway_link` |
| trunk | `trunk`, `trunk_link` |
| primary | `primary`, `primary_link` |
| secondary | `secondary`, `secondary_link` |
| tertiary | `tertiary`, `tertiary_link` |
| local | `residential`, `unclassified` |

**Out:** `footway`, `path`, `cycleway`, `steps`, `pedestrian`, `track`,
`service`, `living_street` (not motor), `construction`, `proposed`,
`bridleway`, `corridor`, `busway`, `raceway`, `road` (too ambiguous; counted
in QA only).

`*_link` ramps inherit the parent class but default to **1 lane**, so they
are not priced as a full 4-lane mainline.

Local roads are about two-thirds of global length. The main case **includes**
them. Every published table also reports a **trunk-and-above / no-local**
sensitivity so the dollar total is not silently dominated by streets.

If a country's `unclassified` length exceeds 60% of its motor-road length,
flag the country (GIRI quality rule) and do not treat that OSM extract as
complete.

---

## 3. OSM vintage — freeze this dump, do not pin to 1995

The road file we already have is `planet-260803`. That name is
**YYMMDD = 2026-08-03**. It is a snapshot of OpenStreetMap as it existed
on 3 August 2026, not a live feed and not a 1990s network.

**Decision (frozen):**

1. **Use this snapshot.** Do not download a newer planet for the paper.
   Re-downloading would change every length and every tag.
2. **Do not rewind OSM to 1995, 2014 or 2018.** Historical global OSM is
   incomplete. The 2017 completeness study found the map only ~83% done,
   with the missing part mostly low-order roads. A 1995 or 2014 “road
   year” cannot be reconstructed globally from OSM.
3. **Write the exposure year in the paper as 2026-08-03.** The climate
   baseline is still 1995–2014 (and the future windows are 2041–2060 /
   2081–2100). This is the standard static-exposure design: today’s network
   × current or future hazard. Future road building is a separate
   socio-economic scenario, not something we invent by using an older map.
4. Record the planet file SHA-256 in the valuation manifest. That hash
   *is* the year stamp.

So: yes, we are using **current** OSM. Yes, it must be **fixed to one
date**. That date is already 2026-08-03. We do not keep updating it.

---

## 4. How a segment gets a unit cost

Work through the segment in this order.

### 4.1 Bridge or tunnel

If OSM says `bridge=yes` (or a real bridge value) or `tunnel=yes`:

- bridge = GIRI 9.84 million USD/km, brought to 2025 USD
- tunnel = GIRI 19.80 million USD/km, brought to 2025 USD

Do not apply the pavement price. Do not apply the terrain multiplier
(those GIRI prices already mix settings).

### 4.2 Ordinary pavement — choose the construction type

| OSM class | Default construction type | Default lanes |
|---|---|---|
| motorway | new 4-lane expressway | 4 |
| trunk | new 4-lane if `lanes` ≥ 4, else new 2-lane highway | 2 |
| primary, secondary | new 2-lane highway | 2 |
| tertiary, paved | new 2-lane highway | 2 |
| tertiary, unpaved | new 1-lane road | 1 |
| local | new 1-lane road | 1 |

Never use ROCKS “Gravel Resurfacing” (~7,000 USD/km). That is a wearing
course, not a rebuild.

### 4.3 Lanes

- If `lanes` is a usable integer: use it.
- If missing: use the **median of tagged lanes for that class in that
  country**. If the country has no tagged sample: use the default in the
  table above.
- Never use a world average of tagged roads.
- Adjust price: **±25% per lane** relative to the default for that class.
  Floor at 0.5 × the default-lane price.

### 4.4 Paved or not

- If `surface` is clearly paved (`asphalt`, `concrete`, `paved`, …): paved.
- If clearly unpaved (`gravel`, `dirt`, `ground`, `unpaved`, `sand`, …): unpaved.
- If missing: use the country’s paved-road share, assigned from high class
  to low class (higher classes are assumed paved first). For a single
  untagged way, use the **expected** cost of that class in that country
  (`paved_fraction × paved price + (1 − paved_fraction) × unpaved price`)
  so national totals are not biased.

The global paved-road statistical series is old and incomplete. Where a
country has no figure, use an income-group default (high income 0.90,
upper-middle 0.65, lower-middle 0.40, low 0.20) and mark it.

### 4.5 Terrain

Mean slope along the way:

| Class | Slope | Multiplier |
|---|---:|---:|
| plain | 0–10° | 1.00 |
| rolling | 10–25° | 1.57 |
| steep | >25° | 1.94 |

Ratios follow GIRI Table 9 (highway plain / rolling / steep =
1.09 / 1.715 / 2.12 million USD/km). If slope has not yet been computed,
publish the flat price and a terrain-on sensitivity; do not pretend the
segment is flat forever.

### 4.6 Which price book (this is the important fork)

ROCKS 2018 is still the only public global file of *actual completed*
new-road costs. It is almost entirely development-bank projects. It has
**no North America**, almost no Japan / Australia / Western Europe, and
very few new-build rows (24 four-lane, 7 two-lane, 6 one-lane).

Therefore:

| Country group | Price book |
|---|---|
| Low- and middle-income countries, except China | ROCKS 2018 **Actual new-build** median for that World Bank region and work type, each row inflated from its project year to 2025 |
| United States, Canada | National / North American construction costs, 2025 USD |
| Europe high-income (EU, EFTA, UK) | van Ginkel / European Court of Auditors 2015 EUR costs, scaled by 2015 PPP GDP per capita relative to the EU-28 mean, then converted to 2025 USD |
| Japan, Korea, Australia, New Zealand, Singapore, Taiwan, Hong Kong | That country’s published construction costs, 2025 USD |
| China | Chinese published expressway / highway unit costs, 2025 USD |
| Other high-income (Gulf, Chile, …) | US class structure × (country PPP GDP per capita / US), clipped to 0.40–1.30 of the US price — **not** a ROCKS “rich-country average” |

If a region has fewer than 3 Actual new-build rows for a work type, do
not use a 1-row oddity (South Asia “new 2-lane” at 42,000 USD/km is
discarded). Fill from the global developing-country median of that work
type after outlier removal, and write `n` and `fill_rule` in the book.

### 4.7 Price year = 2025 construction dollars

- ROCKS values are US dollars in the **project completion year**.
  Inflate each row with the US GDP deflator from that year to 2025
  (complete, official, 1960–2025). This is a general-price inflator, not
  a highway bid index. Report a sensitivity that multiplies 2018–2025 by
  an extra 1.35 to reflect the published US highway-construction spike
  (NHCCI +71.5% from end-2020 to early 2024 versus a much smaller GDP
  deflator).
- Do **not** use CPI.
- Do **not** apply the US highway index to *cross-country levels*.
  African and East Asian observed ROCKS dollars already are the levels.
- European 2015 EUR: convert with the 2015 average EUR/USD, then apply
  the same 2015→2025 US GDP deflator (the costs were already
  EU-average real 2015 prices).
- GIRI bridge/tunnel: treat as 2023 USD and inflate 2023→2025.

---

## 5. What is published with the first value layer

1. Main case: all motor roads, 2025 USD replacement cost.
2. Same without local roads.
3. Developing countries on ROCKS, rich countries on national books —
   versus a rejected “everyone on ROCKS” check (for the appendix, to
   show the high-income error).
4. Low / median / high unit-cost band. Same road type can differ by
   a factor of three to four across countries after terrain is removed
   (Collier et al. 2016).
5. Country flags: unclassified share, missing lane tag share, price-book
   used, paved-share source.

---

## 6. Implementation order

1. Freeze this contract and the 2025 unit-cost book (done in-repo).
2. Value function + tests on synthetic segments (in-repo).
3. Extract every accepted way from `planet-260803` with tags, length,
   and a representative point.
4. Assign country (Natural Earth 50m admin-0).
5. Compute country × class lane medians; apply the book.
6. Add slope / terrain when the DEM pass is ready.
7. Only then overlay 0.05° wind and rain on the valued ways.

Code and numbers live in:

- `code/road_replacement_value.py`
- `code/build_road_unit_cost_book.py`
- `data/valuation/`

---

## 7. Explicitly rejected

- Koks Supplementary Table 8 as our unit costs (we recompute from ROCKS).
- Gravel resurfacing as a replacement cost.
- Pricing the 0.1° density cube.
- World-average lanes or world-average prices for untagged ways.
- Europe-and-Central-Asia ROCKS as a proxy for the US, Japan or Australia.
- CPI as a construction inflator.
- A new OSM download after `planet-260803`.
- Treating 1995–2014 climate as if the road network were also 1995–2014.
