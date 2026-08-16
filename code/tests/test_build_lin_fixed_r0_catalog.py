"""Narrow contract tests for the frozen event-level r0 catalogue."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from netCDF4 import Dataset
import numpy as np


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from build_lin_fixed_r0_catalog import (  # noqa: E402
    CHAVAS2016_R0_MEDIAN_KM,
    LOGNORMAL_MU_LN_KM,
    LOGNORMAL_SIGMA_LN_KM,
    QUANTIZED_OUTER_RADIUS_M_SEQUENCE_SHA256,
    RNG_SEED,
    STANDARD_NORMAL_SEQUENCE_SHA256,
    build_catalogue,
    draw_fixed_r0_km,
    draw_standard_normal,
    load_frozen_event_ids,
    r0_km_to_metres,
    sha256,
)


def _sample_fixture(path: Path, event_ids: list[str]) -> str:
    with Dataset(path, "w", format="NETCDF4") as sample:
        sample.createDimension("event", len(event_ids))
        variable = sample.createVariable("event_id", str, ("event",))
        variable[:] = np.asarray(event_ids, dtype=object)
    return sha256(path)


class FixedR0CatalogueContractTest(unittest.TestCase):
    def test_frozen_parameters_and_10k_draw_are_exact_and_unmodified(self):
        self.assertEqual(RNG_SEED, 20_260_810)
        self.assertAlmostEqual(LOGNORMAL_MU_LN_KM, np.log(881.0), places=15)
        self.assertAlmostEqual(
            LOGNORMAL_SIGMA_LN_KM,
            (np.log(1054.4) - np.log(740.7))
            / (2.0 * 0.6744897501960817),
            places=15,
        )
        self.assertEqual(CHAVAS2016_R0_MEDIAN_KM, 881.0)

        z_values = draw_standard_normal(10_000)
        values = draw_fixed_r0_km(10_000)
        expected_first5_mm = np.array(
            [1000708405, 601141571, 847777465, 660960756, 786728288],
            dtype=np.float64,
        )
        np.testing.assert_array_equal(values[:5], expected_first5_mm / 1_000_000.0)
        self.assertEqual(int(np.rint(values.min() * 1_000_000.0)), 308736422)
        self.assertEqual(int(np.count_nonzero(values < 302.0)), 0)
        z_canonical = np.ascontiguousarray(z_values, dtype="<f8")
        self.assertEqual(
            hashlib.sha256(z_canonical.tobytes(order="C")).hexdigest(),
            STANDARD_NORMAL_SEQUENCE_SHA256,
        )
        canonical_m = np.ascontiguousarray(r0_km_to_metres(values), dtype="<f8")
        self.assertEqual(
            hashlib.sha256(canonical_m.tobytes(order="C")).hexdigest(),
            QUANTIZED_OUTER_RADIUS_M_SEQUENCE_SHA256,
        )
        np.testing.assert_allclose(
            values,
            np.exp(LOGNORMAL_MU_LN_KM + LOGNORMAL_SIGMA_LN_KM * z_values),
            rtol=0.0,
            atol=5e-7,
        )

    def test_sample_sha_and_unique_event_ids_are_hard_requirements(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.nc"
            actual_sha = _sample_fixture(sample, ["e0", "e1", "e2"])
            event_ids, record = load_frozen_event_ids(
                sample,
                expected_sample_sha256=actual_sha,
                expected_count=3,
            )
            self.assertEqual(event_ids.tolist(), ["e0", "e1", "e2"])
            self.assertEqual(record["sha256"], actual_sha)

            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                load_frozen_event_ids(
                    sample,
                    expected_sample_sha256="0" * 64,
                    expected_count=3,
                )

            duplicate = root / "duplicate.nc"
            duplicate_sha = _sample_fixture(duplicate, ["same", "same"])
            with self.assertRaisesRegex(ValueError, "not unique"):
                load_frozen_event_ids(
                    duplicate,
                    expected_sample_sha256=duplicate_sha,
                    expected_count=2,
                )

    def test_outputs_are_hash_linked_atomic_and_immutable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.nc"
            sample_sha = _sample_fixture(sample, ["storm-a", "storm-b", "storm-c"])
            output = root / "catalogue" / "fixed_r0.nc"
            manifest_path = root / "catalogue" / "fixed_r0.manifest.json"

            manifest = build_catalogue(
                sample,
                output,
                manifest_path,
                expected_sample_sha256=sample_sha,
                expected_count=3,
            )
            on_disk_manifest = json.loads(manifest_path.read_text())

            self.assertEqual(manifest, on_disk_manifest)
            self.assertEqual(manifest["status"], "FROZEN_IMMUTABLE")
            self.assertEqual(
                manifest["artifacts"]["fixed_r0_catalogue_netcdf"]["sha256"],
                sha256(output),
            )
            self.assertFalse(
                manifest["method"]["distribution"]["truncation_applied"]
            )
            self.assertFalse(
                manifest["method"]["distribution"][
                    "rejection_or_resampling_applied"
                ]
            )
            self.assertFalse(
                manifest["method"]["distribution"]["clipping_applied"]
            )
            with Dataset(output) as catalogue:
                self.assertEqual(
                    list(catalogue.variables["event_id"][:]),
                    ["storm-a", "storm-b", "storm-c"],
                )
                np.testing.assert_array_equal(
                    catalogue.variables["event_position"][:], [0, 1, 2]
                )
                np.testing.assert_allclose(
                    catalogue.variables["outer_radius_m"][:],
                    r0_km_to_metres(draw_fixed_r0_km(3)),
                    rtol=0.0,
                    atol=0.0,
                )
                self.assertEqual(
                    catalogue.variables["outer_radius_m"].units, "m"
                )
                self.assertEqual(
                    catalogue.outer_radius_m_sequence_sha256,
                    manifest["exact_catalogue_sequence"][
                        "outer_radius_m_sequence_sha256"
                    ],
                )
                self.assertEqual(
                    catalogue.event_outer_radius_binding_sha256,
                    manifest["exact_catalogue_sequence"][
                        "event_outer_radius_binding_sha256"
                    ],
                )
                self.assertEqual(
                    catalogue.distribution_contract_sha256,
                    manifest["method"]["distribution_contract_sha256"],
                )
                self.assertEqual(catalogue.truncation_applied, 0)
                self.assertEqual(catalogue.rejection_or_resampling_applied, 0)
                self.assertEqual(catalogue.clipping_applied, 0)

            before = (sha256(output), sha256(manifest_path))
            with self.assertRaises(FileExistsError):
                build_catalogue(
                    sample,
                    output,
                    manifest_path,
                    expected_sample_sha256=sample_sha,
                    expected_count=3,
                )
            self.assertEqual(before, (sha256(output), sha256(manifest_path)))
            self.assertEqual(list(output.parent.glob(".*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
