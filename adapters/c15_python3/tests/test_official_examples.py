"""Numerical anchors extracted from the three official CLE15 example PDFs.

The PDFs are vector figures in the frozen official ZIP.  Their blue model paths
were converted to SVG with Poppler ``pdftocairo`` and mapped from the plotted
axes (0--1000 km, 0--55 m s-1).  The 0.003 m s-1 tolerance is wider than the
six-decimal SVG coordinate quantisation but far narrower than a scientifically
meaningful wind-field difference.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import sys
import unittest
import warnings

import numpy as np


ADAPTER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ADAPTER_DIR))
import c15  # noqa: E402


PDF_ANCHORS = {
    "r0input": {
        "radius_km": [0, 25, 50, 100, 200, 250, 300, 400, 500, 600, 700, 800],
        "wind_ms": [
            0,
            49.978462,
            38.766216,
            22.251600,
            12.989852,
            10.969427,
            9.528959,
            7.509536,
            6.040860,
            4.804634,
            3.627265,
            2.323625,
        ],
    },
    "rmaxinput": {
        "radius_km": [0, 25, 50, 100, 200, 250, 300, 400, 500, 600, 700, 800, 900],
        "wind_ms": [
            0,
            50.003650,
            39.497859,
            22.777607,
            13.293275,
            11.229014,
            9.760671,
            7.708981,
            6.225819,
            4.987582,
            3.823262,
            2.567866,
            0.676808,
        ],
    },
    "rfitinput": {
        "radius_km": [0, 25, 50, 100, 200, 250, 300, 400, 500, 600, 700, 800, 900],
        "wind_ms": [
            0,
            49.839924,
            41.519627,
            24.339012,
            14.190605,
            11.996426,
            10.443343,
            8.293378,
            6.764135,
            5.514011,
            4.374450,
            3.214876,
            1.807023,
        ],
    },
}


def _official_default(case: str):
    common = (5e-5, 0, 1.5e-3, 2 / 1000, 0, 1.0, 0, 0.15)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if case == "r0input":
            return c15.ER11E04_nondim_r0input(50.0, 900e3, *common)
        if case == "rmaxinput":
            return c15.ER11E04_nondim_rmaxinput(50.0, 25e3, *common)
        if case == "rfitinput":
            return c15.ER11E04_nondim_rfitinput(50.0, 250e3, 12.0, *common)
    raise AssertionError(case)


class OfficialExampleCurvesTest(unittest.TestCase):
    def test_three_official_pdf_curves(self):
        for case, reference in PDF_ANCHORS.items():
            with self.subTest(case=case):
                result = _official_default(case)
                radius_m = result[0]
                wind_ms = result[1]
                self.assertTrue(np.all(np.diff(radius_m) > 0))
                self.assertTrue(np.all(np.isfinite(wind_ms)))
                query_m = np.asarray(reference["radius_km"], dtype=float) * 1000
                expected = np.asarray(reference["wind_ms"], dtype=float)
                actual = np.interp(query_m, radius_m, wind_ms)
                np.testing.assert_allclose(actual, expected, rtol=0, atol=0.003)

    def test_official_default_constraints_are_recovered(self):
        r0_result = _official_default("r0input")
        rmax_result = _official_default("rmaxinput")
        rfit_result = _official_default("rfitinput")

        self.assertAlmostEqual(r0_result[0][np.argmax(r0_result[1])], r0_result[4], delta=1e-8)
        self.assertAlmostEqual(rmax_result[0][np.argmax(rmax_result[1])], 25e3, delta=1e-8)
        self.assertAlmostEqual(
            np.interp(250e3, rfit_result[0], rfit_result[1]),
            12.0,
            delta=0.01,
        )


if __name__ == "__main__":
    unittest.main()
