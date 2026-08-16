#!/usr/bin/env python3
"""Build the immutable event-level outer-radius catalogue for the Lin 10k sample.

The catalogue implements the public reconstruction fixed for this study:

* one outer radius of vanishing wind (``r0``) is drawn per synthetic storm;
* the draw follows a global lognormal distribution matched to the median and
  interquartile range in Chavas et al. (2016), Table 1; and
* NumPy's PCG64 bit generator is initialized once; a standard-normal
  sequence is drawn in frozen 10,000-event sample-position order; and
  ``r0 = exp(mu + sigma * Z)`` is stored to the nearest millimetre so the
  catalogue is bitwise identical across libm/NumPy versions.

No draw is truncated, rejected, resampled, clipped, or otherwise altered.
The frozen sample SHA-256 and unique event identifiers are checked before any
output is written.  The NetCDF is published first by same-directory rename and
the JSON manifest is published last as the completion marker.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
from typing import Any, Sequence
import uuid

from netCDF4 import Dataset, __version__ as netcdf4_version, chartostring
import numpy as np


SCRIPT_VERSION = "1.1.0"
SCHEMA_VERSION = "tcr-fixed-r0-catalogue-v1"

EXPECTED_SAMPLE_COUNT = 10_000
FROZEN_SAMPLE_SHA256 = (
    "856ed368466cf4f8a1f0b8e351bcc8f44eae32d9d60c55623c6ad2217275d1af"
)

CHAVAS2016_TITLE = "Observed Tropical Cyclone Size Revisited"
CHAVAS2016_AUTHORS = (
    "Daniel R. Chavas; Ning Lin; Wenhao Dong; Yanluan Lin"
)
CHAVAS2016_JOURNAL = "Journal of Climate"
CHAVAS2016_DOI = "10.1175/JCLI-D-15-0731.1"
CHAVAS2016_URL = (
    "https://journals.ametsoc.org/view/journals/clim/29/8/"
    "jcli-d-15-0731.1.xml"
)
CHAVAS2016_GLOBAL_N = 578
CHAVAS2016_R0_Q1_KM = 740.7
CHAVAS2016_R0_MEDIAN_KM = 881.0
CHAVAS2016_R0_Q3_KM = 1054.4

STANDARD_NORMAL_Q75 = 0.6744897501960817
LOGNORMAL_MU_LN_KM = math.log(CHAVAS2016_R0_MEDIAN_KM)
LOGNORMAL_SIGMA_LN_KM = (
    math.log(CHAVAS2016_R0_Q3_KM) - math.log(CHAVAS2016_R0_Q1_KM)
) / (2.0 * STANDARD_NORMAL_Q75)

RNG_BIT_GENERATOR = "PCG64"
RNG_SEED = 20_260_810
LOW_R0_DIAGNOSTIC_THRESHOLD_KM = 302.0
# Storage quantum only.  It does not truncate, clip, or resample the draw.
R0_STORAGE_QUANTUM_M = 0.001
STANDARD_NORMAL_SEQUENCE_SHA256 = (
    "c503c52f7f84a93f90b3d85fc8dddf5eb2e244652cfdbf6b2f99b83229abdf10"
)
QUANTIZED_OUTER_RADIUS_M_SEQUENCE_SHA256 = (
    "be647f18bb05dafd0bf54e1343ad1e52a4e6fa7f39921faced5a22869a8eef71"
)


def sha256(path: Path) -> str:
    """Return the SHA-256 of one file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _nc_strings(variable: Any) -> np.ndarray:
    values = variable[:]
    if values.dtype.kind in {"U", "O"}:
        return np.asarray(values, dtype=str)
    if values.dtype.kind == "S" and values.ndim == 2:
        return np.asarray(chartostring(values), dtype=str)
    return np.asarray(values).astype(str)


def canonical_event_id_order_sha256(event_ids: Sequence[str]) -> str:
    """Hash a precisely defined, encoding-independent event-ID sequence."""

    encoded = json.dumps(
        list(event_ids), ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _float64_le_sha256(values: np.ndarray) -> str:
    """Hash a numeric sequence as contiguous IEEE-754 little-endian float64."""

    canonical = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    return hashlib.sha256(canonical.tobytes(order="C")).hexdigest()


def load_frozen_event_ids(
    sample_path: Path,
    *,
    expected_sample_sha256: str = FROZEN_SAMPLE_SHA256,
    expected_count: int = EXPECTED_SAMPLE_COUNT,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Validate the frozen sample identity and return IDs in sample order."""

    sample_path = sample_path.resolve()
    if not sample_path.is_file():
        raise FileNotFoundError(sample_path)
    actual_sha = sha256(sample_path)
    if actual_sha != expected_sample_sha256:
        raise ValueError(
            "frozen sample SHA-256 mismatch: "
            f"expected {expected_sample_sha256}, got {actual_sha}"
        )

    with Dataset(sample_path) as sample:
        if "event" not in sample.dimensions:
            raise ValueError("frozen sample lacks the event dimension")
        if "event_id" not in sample.variables:
            raise ValueError("frozen sample lacks the event_id variable")
        count = len(sample.dimensions["event"])
        if count != expected_count:
            raise ValueError(f"expected {expected_count} events, found {count}")
        event_ids = _nc_strings(sample.variables["event_id"])

    if event_ids.shape != (expected_count,):
        raise ValueError(
            "event_id must be a one-dimensional vector over the event dimension"
        )
    event_ids = np.asarray([str(value) for value in event_ids], dtype=object)
    empty_positions = np.flatnonzero(
        np.asarray([not value.strip() for value in event_ids], dtype=bool)
    )
    if empty_positions.size:
        raise ValueError(
            "blank event_id values at sample positions "
            f"{empty_positions.astype(int).tolist()}"
        )

    unique_ids, counts = np.unique(event_ids.astype(str), return_counts=True)
    duplicates = unique_ids[counts > 1]
    if duplicates.size:
        preview = duplicates[:10].astype(str).tolist()
        raise ValueError(
            f"event_id values are not unique; duplicate examples: {preview}"
        )

    return event_ids, {
        "path": str(sample_path),
        "bytes": sample_path.stat().st_size,
        "sha256": actual_sha,
        "expected_sha256": expected_sample_sha256,
        "event_count": expected_count,
        "event_id_unique_count": int(unique_ids.size),
        "event_id_order_sha256": canonical_event_id_order_sha256(
            event_ids.astype(str).tolist()
        ),
        "event_id_order_sha256_definition": (
            "SHA-256 of the UTF-8 canonical JSON array of event IDs, with "
            "ensure_ascii=false and separators=(',', ':')"
        ),
    }


def draw_standard_normal(count: int) -> np.ndarray:
    """Draw the frozen standard-normal sequence used to construct r0."""

    if count < 1:
        raise ValueError("count must be positive")
    generator = np.random.Generator(np.random.PCG64(RNG_SEED))
    values = np.asarray(generator.standard_normal(count), dtype=np.float64)
    if values.shape != (count,) or not np.all(np.isfinite(values)):
        raise RuntimeError("PCG64 standard-normal draw did not return finite values")
    return values


def r0_millimetres(r0_km: np.ndarray) -> np.ndarray:
    """Return integer-valued millimetres on the published storage grid."""

    return np.rint(np.asarray(r0_km, dtype=np.float64) * 1_000_000.0)


def quantize_r0_km(r0_km: np.ndarray) -> np.ndarray:
    """Store r0 to the published millimetre grid without altering the draw."""

    return r0_millimetres(r0_km) / 1_000_000.0


def r0_km_to_metres(r0_km: np.ndarray) -> np.ndarray:
    """Convert quantized kilometres to metres through integer millimetres."""

    return r0_millimetres(r0_km) / 1000.0


def draw_fixed_r0_km(count: int) -> np.ndarray:
    """Draw the unmodified event-level r0 sequence under the frozen RNG.

    The scientific draw is the explicit lognormal transform
    ``exp(mu + sigma * Z)`` of a PCG64 standard-normal sequence.  NumPy's
    ``Generator.lognormal`` is deliberately not used: that helper is not
    bitwise stable across NumPy versions.  The published catalogue then
    stores each draw to the nearest millimetre so ``exp`` libm rounding
    cannot change the frozen bytes.
    """

    z_values = draw_standard_normal(count)
    raw = np.exp(LOGNORMAL_MU_LN_KM + LOGNORMAL_SIGMA_LN_KM * z_values)
    if not np.all(np.isfinite(raw)) or np.any(raw <= 0.0):
        raise RuntimeError("explicit lognormal transform returned a non-positive r0")
    values = quantize_r0_km(raw)
    if count == EXPECTED_SAMPLE_COUNT:
        z_sha = _float64_le_sha256(z_values)
        r0_sha = _float64_le_sha256(r0_km_to_metres(values))
        if z_sha != STANDARD_NORMAL_SEQUENCE_SHA256:
            raise RuntimeError(
                "PCG64 standard-normal sequence drifted from the frozen contract"
            )
        if r0_sha != QUANTIZED_OUTER_RADIUS_M_SEQUENCE_SHA256:
            raise RuntimeError(
                "quantized r0 sequence drifted from the frozen millimetre contract"
            )
    return values


def _diagnostics(
    event_ids: np.ndarray, r0_km: np.ndarray
) -> dict[str, Any]:
    low_positions = np.flatnonzero(r0_km < LOW_R0_DIAGNOSTIC_THRESHOLD_KM)
    quantiles = np.quantile(r0_km, [0.25, 0.5, 0.75])
    return {
        "event_count": int(r0_km.size),
        "minimum_r0_km": float(np.min(r0_km)),
        "maximum_r0_km": float(np.max(r0_km)),
        "mean_r0_km": float(np.mean(r0_km)),
        "sample_standard_deviation_r0_km": float(np.std(r0_km, ddof=1)),
        "sample_q1_r0_km": float(quantiles[0]),
        "sample_median_r0_km": float(quantiles[1]),
        "sample_q3_r0_km": float(quantiles[2]),
        "low_r0_check": {
            "criterion": "r0_km < threshold_km",
            "threshold_km": LOW_R0_DIAGNOSTIC_THRESHOLD_KM,
            "count": int(low_positions.size),
            "fraction": float(low_positions.size / r0_km.size),
            "event_positions": low_positions.astype(int).tolist(),
            "event_ids": event_ids[low_positions].astype(str).tolist(),
            "r0_km": r0_km[low_positions].astype(float).tolist(),
        },
    }


def _distribution_contract() -> dict[str, Any]:
    return {
        "scientific_source": {
            "citation": (
                "Chavas, D. R., Lin, N., Dong, W., and Lin, Y. (2016), "
                "Observed Tropical Cyclone Size Revisited, Journal of "
                "Climate, 29, 2923-2939."
            ),
            "title": CHAVAS2016_TITLE,
            "authors": CHAVAS2016_AUTHORS,
            "year": 2016,
            "journal": CHAVAS2016_JOURNAL,
            "doi": CHAVAS2016_DOI,
            "url": CHAVAS2016_URL,
            "table": "Table 1, global r0 distribution",
            "published_global_sample_count": CHAVAS2016_GLOBAL_N,
            "published_r0_q1_km": CHAVAS2016_R0_Q1_KM,
            "published_r0_median_km": CHAVAS2016_R0_MEDIAN_KM,
            "published_r0_q3_km": CHAVAS2016_R0_Q3_KM,
        },
        "distribution": {
            "family": "lognormal",
            "radius_definition": "outer radius of vanishing wind (r0)",
            "radius_unit": "km",
            "logarithm": "natural",
            "fit": (
                "quantile matched to the published global median and "
                "interquartile range"
            ),
            "mu_ln_km": LOGNORMAL_MU_LN_KM,
            "sigma_ln_km": LOGNORMAL_SIGMA_LN_KM,
            "standard_normal_q75": STANDARD_NORMAL_Q75,
            "formula_mu": "ln(881.0 km)",
            "formula_sigma": (
                "[ln(1054.4 km) - ln(740.7 km)] / "
                "[2 * Phi^-1(0.75)]"
            ),
            "draw_formula": "r0_km = exp(mu_ln_km + sigma_ln_km * Z)",
            "standard_normal_source": (
                "numpy.random.Generator(PCG64(seed)).standard_normal"
            ),
            "storage_quantization": {
                "applied": True,
                "role": "cross-platform byte identity only",
                "quantum_m": R0_STORAGE_QUANTUM_M,
                "rule": "round-to-nearest millimetre, ties to even",
                "does_not_truncate_or_resample": True,
            },
            "truncation_applied": False,
            "rejection_or_resampling_applied": False,
            "clipping_applied": False,
        },
    }


def _method_record() -> dict[str, Any]:
    contract = _distribution_contract()
    return {
        **contract,
        "distribution_contract_sha256": _canonical_json_sha256(contract),
        "distribution_contract_sha256_definition": (
            "SHA-256 of canonical UTF-8 JSON for scientific_source and "
            "distribution, with sorted keys and compact separators"
        ),
        "rng": {
            "library": "numpy",
            "numpy_version": np.__version__,
            "bit_generator": RNG_BIT_GENERATOR,
            "seed_decimal": RNG_SEED,
            "generator_call": (
                "numpy.random.Generator(PCG64(seed)).standard_normal; "
                "r0_km = exp(mu + sigma * Z); store nearest millimetre"
            ),
            "standard_normal_sequence_sha256": STANDARD_NORMAL_SEQUENCE_SHA256,
            "quantized_outer_radius_m_sequence_sha256": (
                QUANTIZED_OUTER_RADIUS_M_SEQUENCE_SHA256
            ),
            "draw_count_per_event": 1,
            "ordering": (
                "one sequential draw in frozen sample event-position order, "
                "positions 0 through 9999"
            ),
            "stream_reinitialized_per_event": False,
        },
    }


def _write_catalogue_netcdf(
    path: Path,
    *,
    event_ids: np.ndarray,
    r0_km: np.ndarray,
    source_record: dict[str, Any],
    diagnostics: dict[str, Any],
) -> None:
    method = _method_record()
    outer_radius_m = r0_km_to_metres(r0_km)
    draw_sha = _float64_le_sha256(outer_radius_m)
    binding_sha = _canonical_json_sha256(
        [
            [str(event_id), float(radius_m).hex()]
            for event_id, radius_m in zip(event_ids, outer_radius_m, strict=True)
        ]
    )
    with Dataset(path, "w", format="NETCDF4") as output:
        output.createDimension("event", r0_km.size)
        positions = output.createVariable("event_position", "i4", ("event",))
        identifiers = output.createVariable("event_id", str, ("event",))
        radius = output.createVariable(
            "outer_radius_m",
            "f8",
            ("event",),
            zlib=True,
            complevel=4,
            shuffle=True,
            chunksizes=(min(4096, r0_km.size),),
        )

        positions.long_name = "zero-based position in the frozen 10k sample"
        identifiers.long_name = "frozen synthetic-storm event identifier"
        radius.long_name = "event-level outer radius of vanishing wind"
        radius.units = "m"
        radius.draw_role = "one fixed value per storm"
        radius.source_distribution_unit = "km"

        positions[:] = np.arange(r0_km.size, dtype=np.int32)
        identifiers[:] = event_ids.astype(object)
        radius[:] = outer_radius_m

        output.title = "Frozen event-level r0 catalogue for the Lin 10k sample"
        output.schema_version = SCHEMA_VERSION
        output.script_version = SCRIPT_VERSION
        output.status = "FROZEN_IMMUTABLE"
        output.source_sample_sha256 = source_record["sha256"]
        output.source_event_id_order_sha256 = source_record[
            "event_id_order_sha256"
        ]
        output.scientific_source_doi = CHAVAS2016_DOI
        output.scientific_source_table = "Table 1, global r0 distribution"
        output.distribution = "lognormal in natural-log kilometres"
        output.lognormal_mu_ln_km = LOGNORMAL_MU_LN_KM
        output.lognormal_sigma_ln_km = LOGNORMAL_SIGMA_LN_KM
        output.distribution_contract_sha256 = method[
            "distribution_contract_sha256"
        ]
        output.rng_bit_generator = RNG_BIT_GENERATOR
        output.rng_seed_decimal = np.int64(RNG_SEED)
        output.rng_numpy_version = method["rng"]["numpy_version"]
        output.draw_order = "frozen sample event positions 0..9999"
        output.draws_per_event = np.int32(1)
        output.truncation_applied = np.int8(0)
        output.rejection_or_resampling_applied = np.int8(0)
        output.clipping_applied = np.int8(0)
        output.outer_radius_m_sequence_sha256 = draw_sha
        output.event_outer_radius_binding_sha256 = binding_sha
        output.low_r0_diagnostic_threshold_km = (
            LOW_R0_DIAGNOSTIC_THRESHOLD_KM
        )
        output.count_r0_less_than_302_km = np.int32(
            diagnostics["low_r0_check"]["count"]
        )


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def _temporary_path(target: Path) -> Path:
    return target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")


def build_catalogue(
    sample_path: Path,
    output_nc: Path,
    manifest_path: Path | None = None,
    *,
    expected_sample_sha256: str = FROZEN_SAMPLE_SHA256,
    expected_count: int = EXPECTED_SAMPLE_COUNT,
    parent_draw_count: int | None = None,
) -> dict[str, Any]:
    """Build and atomically publish the immutable catalogue and manifest."""

    sample_path = sample_path.resolve()
    output_nc = output_nc.resolve()
    manifest_path = (
        manifest_path.resolve()
        if manifest_path is not None
        else output_nc.with_suffix(".manifest.json")
    )
    if output_nc == manifest_path:
        raise ValueError("NetCDF and manifest paths must differ")
    existing = [path for path in (output_nc, manifest_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "immutable catalogue output already exists: "
            + ", ".join(str(path) for path in existing)
        )

    event_ids, source_record = load_frozen_event_ids(
        sample_path,
        expected_sample_sha256=expected_sample_sha256,
        expected_count=expected_count,
    )
    if parent_draw_count is None:
        r0_km = draw_fixed_r0_km(expected_count)
    else:
        if parent_draw_count < expected_count:
            raise ValueError("parent draw count is smaller than the sample")
        parent = draw_fixed_r0_km(parent_draw_count)
        with Dataset(sample_path) as sample:
            if "source_track_index" not in sample.variables:
                raise ValueError("sample lacks source_track_index for parent-draw binding")
            track_index = np.asarray(
                sample.variables["source_track_index"][:], dtype=np.int64
            )
        if track_index.size != expected_count:
            raise ValueError("source_track_index length differs from the sample")
        if np.any(track_index < 0) or np.any(track_index >= parent_draw_count):
            raise ValueError("source_track_index is outside the parent draw sequence")
        r0_km = parent[track_index]
    diagnostics = _diagnostics(event_ids, r0_km)
    outer_radius_m = r0_km_to_metres(r0_km)
    exact_sequence = {
        "outer_radius_m_sequence_sha256": _float64_le_sha256(outer_radius_m),
        "outer_radius_m_sequence_sha256_definition": (
            "SHA-256 over event-position-ordered contiguous IEEE-754 "
            "little-endian float64 outer_radius_m bytes"
        ),
        "event_outer_radius_binding_sha256": _canonical_json_sha256(
            [
                [str(event_id), float(radius_m).hex()]
                for event_id, radius_m in zip(
                    event_ids, outer_radius_m, strict=True
                )
            ]
        ),
        "event_outer_radius_binding_sha256_definition": (
            "SHA-256 of canonical UTF-8 JSON containing ordered "
            "[event_id, float.hex(outer_radius_m)] pairs"
        ),
    }

    output_nc.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_nc = _temporary_path(output_nc)
    temporary_manifest = _temporary_path(manifest_path)
    published_nc = False
    try:
        _write_catalogue_netcdf(
            temporary_nc,
            event_ids=event_ids,
            r0_km=r0_km,
            source_record=source_record,
            diagnostics=diagnostics,
        )
        with temporary_nc.open("rb") as stream:
            os.fsync(stream.fileno())
        catalogue_record = {
            "path": str(output_nc),
            "bytes": temporary_nc.stat().st_size,
            "sha256": sha256(temporary_nc),
        }
        script_path = Path(__file__).resolve()
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "status": "FROZEN_IMMUTABLE",
            "publication_contract": (
                "The NetCDF is atomically renamed first; this manifest is "
                "atomically renamed last and is the completed-publication marker."
            ),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_sample": source_record,
            "method": _method_record(),
            "exact_catalogue_sequence": exact_sequence,
            "diagnostics": diagnostics,
            "artifacts": {"fixed_r0_catalogue_netcdf": catalogue_record},
            "software": {
                "builder": {
                    "path": str(script_path),
                    "version": SCRIPT_VERSION,
                    "sha256": sha256(script_path),
                },
                "python_version": platform.python_version(),
                "python_implementation": platform.python_implementation(),
                "numpy_version": np.__version__,
                "netCDF4_python_version": netcdf4_version,
            },
        }
        _write_json(temporary_manifest, manifest)

        os.replace(temporary_nc, output_nc)
        published_nc = True
        os.replace(temporary_manifest, manifest_path)
        published_nc = False
        return manifest
    except Exception:
        # Roll back only the just-published NetCDF if the final completion
        # marker could not be published.  Pre-existing files are rejected above.
        if published_nc and output_nc.exists() and not manifest_path.exists():
            output_nc.unlink()
        raise
    finally:
        for temporary in (temporary_nc, temporary_manifest):
            if temporary.exists():
                temporary.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        type=Path,
        required=True,
        help="frozen 10,000-event sample NetCDF (exact SHA-256 required)",
    )
    parser.add_argument(
        "--output-nc",
        type=Path,
        required=True,
        help="new immutable fixed-r0 catalogue NetCDF",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="JSON manifest (default: output filename with .manifest.json)",
    )
    parser.add_argument(
        "--expected-count",
        type=int,
        default=EXPECTED_SAMPLE_COUNT,
        help="required sample length; use the road-domain count for production",
    )
    parser.add_argument(
        "--expected-sample-sha256",
        type=str,
        default=FROZEN_SAMPLE_SHA256,
        help="required SHA-256 of the sample NetCDF",
    )
    parser.add_argument(
        "--parent-draw-count",
        type=int,
        default=None,
        help="draw this many r0 values in source_track_index order, then subset",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = build_catalogue(
        args.sample,
        args.output_nc,
        args.manifest,
        expected_sample_sha256=args.expected_sample_sha256,
        expected_count=args.expected_count,
        parent_draw_count=args.parent_draw_count,
    )
    result = {
        "status": manifest["status"],
        "catalogue": manifest["artifacts"]["fixed_r0_catalogue_netcdf"],
        "manifest": str(
            args.manifest.resolve()
            if args.manifest is not None
            else args.output_nc.resolve().with_suffix(".manifest.json")
        ),
        "count_r0_less_than_302_km": manifest["diagnostics"][
            "low_r0_check"
        ]["count"],
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
