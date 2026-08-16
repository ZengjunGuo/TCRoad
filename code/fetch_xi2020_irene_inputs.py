#!/usr/bin/env python3
"""Freeze the public IBTrACS and NCEP inputs for the Xi et al. (2020) Irene case.

This downloader implements only the source contract documented in
``methods/XI2020_HISTORICAL_CASE_INPUT_PLAN.md``.  It deliberately does not
select or fabricate the unresolved ERA-Interim topography/roughness fields.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
import urllib.request

import numpy as np
from netCDF4 import Dataset, chartostring, num2date


SID = "2011233N15301"
TIME_START = "2011-08-21 00:00:00"
TIME_END = "2011-08-30 00:00:00"
METHOD_CONTRACT_SHA256 = (
    "102f7a6f2462058948eb3e9997f9a0bbd3c3253c9e162123084c23023984c5b9"
)

OBJECTS = {
    "ibtracs_na": {
        "filename": "IBTrACS.NA.v04r01.nc",
        "url": (
            "https://www.ncei.noaa.gov/data/international-best-track-archive-for-"
            "climate-stewardship-ibtracs/v04r01/access/netcdf/"
            "IBTrACS.NA.v04r01.nc"
        ),
        "kind": "ibtracs",
    },
    "ncep_r1_shum925": {
        "filename": (
            "ncep_r1_shum925_irene2011_20110821T00_20110830T00_"
            "bbox12p5_55_277p5_305.nc"
        ),
        "url": (
            "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/"
            "pressure/shum.2011.nc?var=shum&north=55&west=277.5&east=305&"
            "south=12.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&"
            "time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=925&"
            "accept=netcdf4"
        ),
        "kind": "ncep",
        "variable": "shum",
        "level": 925.0,
        "lat": [12.5, 55.0, 18],
        "lon": [277.5, 305.0, 12],
        "units": "kg/kg",
    },
    "ncep_r1_uwnd200": {
        "filename": (
            "ncep_r1_uwnd200_irene2011_20110821T00_20110830T00_"
            "bbox7p5_62p5_272p5_312p5.nc"
        ),
        "url": (
            "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/"
            "pressure/uwnd.2011.nc?var=uwnd&north=62.5&west=272.5&east=312.5&"
            "south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&"
            "time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=200&"
            "accept=netcdf4"
        ),
        "kind": "ncep",
        "variable": "uwnd",
        "level": 200.0,
        "lat": [7.5, 62.5, 23],
        "lon": [272.5, 312.5, 17],
        "units": "m/s",
    },
    "ncep_r1_uwnd850": {
        "filename": (
            "ncep_r1_uwnd850_irene2011_20110821T00_20110830T00_"
            "bbox7p5_62p5_272p5_312p5.nc"
        ),
        "url": (
            "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/"
            "pressure/uwnd.2011.nc?var=uwnd&north=62.5&west=272.5&east=312.5&"
            "south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&"
            "time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=850&"
            "accept=netcdf4"
        ),
        "kind": "ncep",
        "variable": "uwnd",
        "level": 850.0,
        "lat": [7.5, 62.5, 23],
        "lon": [272.5, 312.5, 17],
        "units": "m/s",
    },
    "ncep_r1_vwnd200": {
        "filename": (
            "ncep_r1_vwnd200_irene2011_20110821T00_20110830T00_"
            "bbox7p5_62p5_272p5_312p5.nc"
        ),
        "url": (
            "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/"
            "pressure/vwnd.2011.nc?var=vwnd&north=62.5&west=272.5&east=312.5&"
            "south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&"
            "time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=200&"
            "accept=netcdf4"
        ),
        "kind": "ncep",
        "variable": "vwnd",
        "level": 200.0,
        "lat": [7.5, 62.5, 23],
        "lon": [272.5, 312.5, 17],
        "units": "m/s",
    },
    "ncep_r1_vwnd850": {
        "filename": (
            "ncep_r1_vwnd850_irene2011_20110821T00_20110830T00_"
            "bbox7p5_62p5_272p5_312p5.nc"
        ),
        "url": (
            "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/"
            "pressure/vwnd.2011.nc?var=vwnd&north=62.5&west=272.5&east=312.5&"
            "south=7.5&horizStride=1&time_start=2011-08-21T00%3A00%3A00Z&"
            "time_end=2011-08-30T00%3A00%3A00Z&timeStride=1&vertCoord=850&"
            "accept=netcdf4"
        ),
        "kind": "ncep",
        "variable": "vwnd",
        "level": 850.0,
        "lat": [7.5, 62.5, 23],
        "lon": [272.5, 312.5, 17],
        "units": "m/s",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decoded_variables_sha256(data: Dataset, names: set[str]) -> str:
    """Hash decoded variable identity and values, independent of HDF5 packing."""
    digest = hashlib.sha256()
    for name in sorted(names):
        variable = data.variables[name]
        values = np.ma.asarray(variable[:])
        mask = np.ma.getmaskarray(values)
        payload = np.ascontiguousarray(np.ma.filled(values, 0))
        identity = {
            "name": name,
            "dimensions": list(variable.dimensions),
            "shape": list(variable.shape),
            "dtype": str(variable.dtype),
            "units": getattr(variable, "units", None),
        }
        digest.update(json.dumps(identity, sort_keys=True).encode("utf-8"))
        digest.update(np.ascontiguousarray(mask, dtype=np.uint8).tobytes())
        digest.update(payload.tobytes())
    return digest.hexdigest()


def download(url: str, target: Path) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "TCRoad/1.0"})
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".part", dir=target.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=180) as response, temporary.open(
            "wb"
        ) as output:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                output.write(block)
            headers = {
                key: response.headers.get(key)
                for key in (
                    "Content-Type",
                    "Content-Length",
                    "ETag",
                    "Last-Modified",
                )
                if response.headers.get(key) is not None
            }
            final_url = response.geturl()
            status = response.status
        with temporary.open("rb") as handle:
            magic = handle.read(8)
        if not (magic.startswith(b"CDF") or magic == b"\x89HDF\r\n\x1a\n"):
            raise ValueError(f"download is not NetCDF/HDF5: {magic!r}")
        os.replace(temporary, target)
        return {"http_status": status, "final_url": final_url, "headers": headers}
    finally:
        if temporary.exists():
            temporary.unlink()


def validate_ncep(path: Path, contract: dict[str, object]) -> dict[str, object]:
    variable = str(contract["variable"])
    with Dataset(path) as data:
        required = {"time", "level", "lat", "lon", variable}
        if not required.issubset(data.variables):
            raise ValueError(f"missing NCEP variables: {required - set(data.variables)}")
        time = data.variables["time"]
        dates = num2date(time[:], units=time.units, calendar=getattr(time, "calendar", "standard"))
        if len(dates) != 37:
            raise ValueError(f"expected 37 six-hourly times, got {len(dates)}")
        first = dates[0].strftime("%Y-%m-%d %H:%M:%S")
        last = dates[-1].strftime("%Y-%m-%d %H:%M:%S")
        if (first, last) != (TIME_START, TIME_END):
            raise ValueError(f"unexpected time range: {(first, last)}")
        level = np.asarray(data.variables["level"][:], dtype=float)
        if level.shape != (1,) or not np.isclose(level[0], contract["level"]):
            raise ValueError(f"unexpected pressure level: {level}")
        lat = np.asarray(data.variables["lat"][:], dtype=float)
        lon = np.asarray(data.variables["lon"][:], dtype=float)
        lat_min, lat_max, nlat = contract["lat"]
        lon_min, lon_max, nlon = contract["lon"]
        if len(lat) != nlat or not np.allclose([lat.min(), lat.max()], [lat_min, lat_max]):
            raise ValueError(f"unexpected latitude axis: {lat}")
        if len(lon) != nlon or not np.allclose([lon.min(), lon.max()], [lon_min, lon_max]):
            raise ValueError(f"unexpected longitude axis: {lon}")
        values = np.ma.asarray(data.variables[variable][:])
        if values.shape != (37, 1, nlat, nlon):
            raise ValueError(f"unexpected {variable} shape: {values.shape}")
        if np.ma.count_masked(values) or not np.all(np.isfinite(values.filled(np.nan))):
            raise ValueError(f"{variable} contains missing or non-finite values")
        units = getattr(data.variables[variable], "units", "")
        if units != contract["units"]:
            raise ValueError(f"unexpected {variable} units: {units!r}")
        return {
            "dimensions": {name: len(dim) for name, dim in data.dimensions.items()},
            "variable": variable,
            "units": units,
            "pressure_hpa": float(level[0]),
            "time_start": first + " UTC",
            "time_end": last + " UTC",
            "time_count": len(dates),
            "lat_range": [float(lat.min()), float(lat.max())],
            "lon_range": [float(lon.min()), float(lon.max())],
            "value_range": [float(values.min()), float(values.max())],
            "decoded_variables_sha256": decoded_variables_sha256(data, required),
        }


def validate_ibtracs(path: Path) -> dict[str, object]:
    with Dataset(path) as data:
        required = {
            "sid",
            "iso_time",
            "iflag",
            "usa_atcf_id",
            "usa_lat",
            "usa_lon",
            "usa_wind",
            "usa_rmw",
        }
        if not required.issubset(data.variables):
            raise ValueError(f"missing IBTrACS variables: {required - set(data.variables)}")
        sid = np.asarray(chartostring(data.variables["sid"][:])).astype(str)
        positions = np.flatnonzero(sid == SID)
        if len(positions) != 1:
            raise ValueError(f"expected one {SID}, found {len(positions)}")
        index = int(positions[0])
        times = np.asarray(chartostring(data.variables["iso_time"][index])).astype(str)
        standard = []
        for time_index, text in enumerate(times):
            if not text:
                continue
            stamp = datetime.strptime(text, "%Y-%m-%d %H:%M:%S")
            if stamp.minute == 0 and stamp.second == 0 and stamp.hour in (0, 6, 12, 18):
                if datetime(2011, 8, 21) <= stamp <= datetime(2011, 8, 30):
                    standard.append((time_index, text))
        if len(standard) != 37 or standard[0][1] != TIME_START or standard[-1][1] != TIME_END:
            raise ValueError(f"unexpected Irene six-hourly sequence: {standard}")
        indices = np.asarray([item[0] for item in standard], dtype=int)
        wind = np.ma.asarray(data.variables["usa_wind"][index, indices])
        rmw = np.ma.asarray(data.variables["usa_rmw"][index, indices])
        return {
            "storm_dimension": len(data.dimensions["storm"]),
            "sid": SID,
            "storm_index": index,
            "six_hourly_time_start": standard[0][1] + " UTC",
            "six_hourly_time_end": standard[-1][1] + " UTC",
            "six_hourly_time_count": len(standard),
            "usa_wind_available_count": int(wind.count()),
            "usa_rmw_available_count": int(rmw.count()),
            "usa_wind_units": getattr(data.variables["usa_wind"], "units", ""),
            "usa_rmw_units": getattr(data.variables["usa_rmw"], "units", ""),
            "decoded_variables_sha256": decoded_variables_sha256(data, required),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    raw_dir = output_dir / "raw"
    manifest_path = output_dir / "irene2011_inputs.manifest.json"
    if manifest_path.exists() and not args.overwrite:
        raise FileExistsError(f"manifest already exists: {manifest_path}")

    records: dict[str, object] = {}
    for key, contract in OBJECTS.items():
        path = raw_dir / str(contract["filename"])
        if path.exists():
            if not args.overwrite:
                raise FileExistsError(f"input already exists: {path}")
            path.unlink()
        transfer = download(str(contract["url"]), path)
        validation = (
            validate_ibtracs(path)
            if contract["kind"] == "ibtracs"
            else validate_ncep(path, contract)
        )
        records[key] = {
            "relative_path": str(path.relative_to(output_dir)),
            "source_url": contract["url"],
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "transfer": transfer,
            "validation": validation,
        }

    manifest = {
        "schema_version": "1.0",
        "status": "pass",
        "purpose": "Xi et al. (2020) Irene historical C15-TCRM input freeze",
        "case": {"name": "IRENE", "sid": SID, "atcf_id": "AL092011"},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "downloader_sha256": sha256(Path(__file__).resolve()),
        "method_contract": "TCRoad/methods/XI2020_HISTORICAL_CASE_INPUT_PLAN.md",
        "method_contract_sha256": METHOD_CONTRACT_SHA256,
        "scope_note": (
            "This package intentionally excludes unresolved ERA-Interim "
            "topography and roughness inputs and is not yet a runnable Xi2020 case."
        ),
        "artifacts": records,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{manifest_path.name}.", suffix=".tmp", dir=output_dir
    )
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, manifest_path)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(json.dumps({"status": "pass", "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
