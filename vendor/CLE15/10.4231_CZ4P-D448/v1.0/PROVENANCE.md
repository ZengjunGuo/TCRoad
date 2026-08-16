# Official C15 source provenance

This directory freezes the official source distribution for the complete tropical-cyclone radial wind profile of Chavas, Lin and Emanuel (2015).

## Identity

- DOI: `10.4231/CZ4P-D448`
- Title: *Code for tropical cyclone wind profile model of Chavas et al. (2015, JAS)*
- Repository: Purdue University Research Repository (PURR)
- Repository version: `1.0`
- Issued: 2022
- PURR record: <https://purr.purdue.edu/publications/4066/1>
- DOI landing page: <https://doi.org/10.4231/CZ4P-D448>
- Rights in the PURR/DataCite record: `CC0-1.0`

The PURR record describes two files (Python and ZIP) and identifies them as the primary version first distributed by Daniel Chavas's Purdue laboratory on 2020-06-23.

## Acquisition

The live PURR attachment endpoint closed the TLS connection in the acquisition environment on 2026-08-13. The byte-identical dated primary files were therefore retrieved from the Internet Archive captures of their original Purdue University URLs:

- Original ZIP URL: <https://web.ics.purdue.edu/~dchavas/download/code/CLE15_windprofile_PUBLIC_2020-06-23.zip>
- Archived ZIP: <https://web.archive.org/web/20201006213832id_/https://web.ics.purdue.edu/~dchavas/download/code/CLE15_windprofile_PUBLIC_2020-06-23.zip>
- Original Python URL: <https://web.ics.purdue.edu/~dchavas/download/code/CLE15_2020-06-23.py>
- Archived Python: <https://web.archive.org/web/20201006213833id_/https://web.ics.purdue.edu/~dchavas/download/code/CLE15_2020-06-23.py>
- PURR bundle endpoint attempted: <https://purr.purdue.edu/publications/4066/serve/1?render=archive>

The authoritative DataCite metadata response is preserved as `original/datacite_10.4231_cz4p-d448.json`. File sizes and checksums are:

| File | Bytes | SHA-256 |
|---|---:|---|
| `original/CLE15_windprofile_PUBLIC_2020-06-23.zip` | 1,455,158 | `5ebb6e7c253e927653b4515278ace60cb43d17504c720e2b6bc03ab2f1519dc2` |
| `original/CLE15_2020-06-23.py` | 22,762 | `6f1306fae71d0e772f17dbf67d5c8cfd94fa543dd122cbb390c6d50161325113` |
| `original/datacite_10.4231_cz4p-d448.json` | 14,584 | `c103149466604d72058a79fb2a77acc4bec6bb236d9f2b8811e00a4c8ce902a9` |

## Source boundary

- `original/` is immutable source evidence.
- `extracted/` is a byte-preserving extraction of the official MATLAB bundle.
- Future compatibility code must live outside these directories and must never overwrite the originals.
- The MATLAB functions and bundled example PDFs are the numerical reference. The separate Python file uses Python 2 syntax, deprecated NumPy aliases and a commented eye-adjustment block; it is preserved as source evidence, not silently modernized in place.

The ZIP's older README contains a historical non-redistribution instruction, whereas the later PURR DOI record explicitly assigns CC0-1.0. Both records are retained without alteration; the formal repository rights metadata is reported as the governing release metadata.

## Official entry points

The MATLAB bundle supplies three physical entry modes:

- known outer radius `r0`: `ER11E04_nondim_r0input`
- known radius of maximum wind `rmax`: `ER11E04_nondim_rmaxinput`
- known intermediate wind radius `(rfit, Vfit)`: `ER11E04_nondim_rfitinput`

Each merges the ER11 inner solution with the E04 outer solution. No TCR rainfall code is included in this archive.
