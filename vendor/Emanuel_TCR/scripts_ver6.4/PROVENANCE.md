# Emanuel synthetic-TC MATLAB/TCR scripts v6.4: provenance

Retrieved: 2026-08-13 (America/Los_Angeles)

Official public directory maintained under Kerry Emanuel's MIT site:

- Index: https://texmex.mit.edu/pub/emanuel/ForAlex/
- Archive: https://texmex.mit.edu/pub/emanuel/ForAlex/scripts_ver6.4.zip
- README: https://texmex.mit.edu/pub/emanuel/ForAlex/Readme_matlab_scripts_ver6.4.pdf
- User guide: https://texmex.mit.edu/pub/emanuel/ForAlex/UsersGuide_ver6.4.pdf
- Quick-start guide: https://texmex.mit.edu/pub/emanuel/ForAlex/QuickStart_Guide_ver6.4.pdf

The directory index reports the ZIP timestamp as 2021-01-19 15:20 and size as
102,703,919 bytes. The locally frozen ZIP has that exact size and passed
`unzip -t` with "No errors detected in compressed data". Selected scripts were
extracted byte-for-byte into `extracted/`; no upstream file was edited.

## Authority and scope boundary

This is an official, publicly downloadable Emanuel/WindRiskTech MATLAB package
for processing synthetic tropical-cyclone event sets. It is primary evidence
for one public implementation of the synthetic-event TCR input adapter and
rainfall calculation. It is **not** identified by Xi et al. (2023) as the exact
production code used in the Nature Climate Change study. In particular, this
package defaults to the Emanuel--Rotunno (2011) radial wind profile, not C15.
Therefore it does not close the unpublished C15-to-TCR coupling boundary.

The frozen archive also contains the exact static arrays consumed by this public
implementation: `C_Drag500.mat` (`cd`, 1440 x 721; 0.25-degree global grid),
`bathymetry.mat` (`bathy`, 1440 x 721; 0.25 degree), and
`bathymetry_high.mat` (`bathy`, 3600 x 1800; 0.1 degree). They are preserved
byte-for-byte under `extracted/scripts/`. Their presence closes a source-exact
v6.4 run without guessing an ERA-Interim parameter ID, but the package metadata
does not prove that these were the exact static files used for Xi et al. (2020)
or Xi et al. (2023).

## Rights boundary

No package-level open-source license was found in the downloaded archive or
guides during this audit. Several source files carry WindRiskTech copyright
notices. The original site provides public download access, but redistribution,
modification, and relicensing rights are **UNKNOWN**. The frozen files are kept
unaltered for internal scholarly verification; do not publish or redistribute
them without resolving the rights status.
