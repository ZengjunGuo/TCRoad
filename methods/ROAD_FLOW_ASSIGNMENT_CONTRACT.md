# 道路流量赋值合同（方法冻结，今夜不生产全球 OD）

Status: **method frozen 2026-08-16**. Implementation of a global assigned-flow cube is **not** started tonight.

This note answers: how should this tropical-cyclone / road study put a *use* number on each motor-road segment? It is separate from replacement cost.

---

## Decision (plain language)

**We will not use observed traffic counts (AADT / loop detectors / TomTom / HERE).**

**We will use a commute origin–destination matrix of the WorldOD / GlODGen class, then assign those trips onto the frozen OSM motor-road graph.**

The number on a segment is therefore:

> **assigned daily commute flow = how many home-to-work trips are routed over this OSM way**

not “what a counter measured last Tuesday”.

Tonight we freeze that choice. We do **not** download WorldOD, do **not** run a global traffic assignment, and do **not** publish a flow cube.

---

## Why observed AADT is the wrong primary source here

1. **It does not exist globally on OSM ways.** High-income highway agencies publish AADT on numbered routes. Most of the 113 million motor ways in `planet-260803` have no count. Filling holes with “nearby AADT” or national averages would invent traffic on African and island networks that are exactly the tropical-cyclone core.
2. **AADT is not a network state.** A count on one link does not say where those vehicles came from or where they can divert. Indirect (user) loss from a typhoon is a rerouting problem. That needs an origin–destination matrix and a graph, which Koks et al. (2019) *Nature Communications* explicitly left to future work.
3. **Counts mix freight, through-traffic and commuting, with inconsistent years.** They cannot be aligned to the 2026-08-03 OSM vintage or to a single socio-economic year.
4. Commercial speed/volume products are not openly reproducible at planet scale.

Using AADT as the headline “flow” would make the paper look empirical where the map is actually empty.

---

## What WorldOD / GlODGen actually are

Two related Tsinghua products, both built so that **any city** can get a commute OD without local surveys:

- **WorldOD / WorldCommuting-OD** (Rong et al., 2025): an open global commute origin–destination dataset covering on the order of 1,625 cities in 179 countries. The object is *people moving from home zones to work zones on a typical day*, not a highway agency count.
- **GlODGen** (Rong et al., NeurIPS 2025; code `tsinghua-fib-lab/generate-od-pubtools`): a generator that takes **public satellite imagery + population** and produces an OD matrix for a city that was not in the training set. The authors report that public inputs recover about 98% of the skill of hard-to-get local urban features.

That is why this study named WorldOD in the first place: it is the only **open, global, zone-to-zone commute demand** that can sit on top of OSM.

They are **not** a ready-made AADT layer. They are demand. Assignment is a second step.

---

## How assignment will work (when we implement it)

1. Freeze zonal geography (GHSL / urban clusters consistent with WorldOD cities). Rural areas outside the 1,625-city set are either generated with GlODGen where population warrants it, or flagged as *no commute OD* — not silently filled with AADT.
2. Build a directed graph from the **same** `planet-260803` motor-road ways used for replacement cost (object-level, not 0.1°).
3. Assign each OD pair with a standard all-or-nothing or user-equilibrium router (travel time = free-flow time from OSM `maxspeed` / class defaults). Freight is **out of scope** for the first flow layer.
4. Store on each way: assigned commute trips / day, plus the list of OD pairs that use it (or a compressed edge-use table). That is what later gets multiplied by a depth–disruption function (Pregnolato et al. 2017) once inundation exists.
5. Publish two numbers, always: **asset loss** (replacement cost × damage ratio) and **user loss** (extra travel time × value of time). They must not be added as if they were the same dollar.

Pregnolato, Ford, Wilkinson & Dawson (2017), *Transportation Research Part D*: vehicle speed falls continuously with flood depth and roads are effectively impassable near 30 cm. That function needs **depth** and **who was going to use the link**. AADT cannot supply the second; assigned OD can.

---

## Alternatives considered and rejected

| Alternative | Why not as the main flow |
|---|---|
| Observed AADT | No global coverage on OSM ways; not a reroutable demand. |
| Population / road-density gravity with no OD | Invents traffic proportional to people nearby; misses corridors and bottlenecks. |
| Night-time lights or mobile-phone traces as the layer | Not openly reproducible at planet scale; licensing and year mismatch. |
| IRF / national vehicle-km totals smeared onto all roads | A national scalar, not a segment flow. Useful only as a country-level check. |
| Assuming every way carries the class-average flow | Same error as “untagged lanes = world average”. |

A country-level vehicle-km check against IRF World Road Statistics **is** allowed later, as a sanity bound, not as the segment value.

---

## What is explicitly not claimed tonight

- No global OD file is in this repository.
- No traffic assignment has been run.
- Future climate windows still use **today’s** network and **today’s** commute pattern unless a separate socio-economic scenario is added. That is the same static-exposure rule as the 2026-08-03 OSM snapshot.

Next implementation step (not tonight): download WorldOD for a pilot basin (Gulf Coast or Pearl River Delta — we already have those OSM insets), assign onto the inset graph, and only then scale the method.
