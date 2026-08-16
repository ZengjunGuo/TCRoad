"""Public-code-derived bridge from official C15 to CLIMADA Petals TCR.

This module has one deliberately narrow responsibility: provide the radial
wind quantities required by the public CLIMADA Petals 6.2.0 TCR numerical
implementation from the frozen official C15 solver.  Two official input modes
are exposed: ``rmaxinput`` for the observed-RMW Irene benchmark and
``r0input`` for a fixed event outer radius.  It does not modify CLIMADA, add a
radial taper, extrapolate a wind tail, or substitute the ER11 profile used by
CLIMADA's default TCR configuration.

The bridge is not claimed to be source-code-identical to the unpublished
Xi/Lin production implementation.  Its two wind quantities follow the public
C15 and TCR angular-momentum definitions:

* ``gradient``: the official C15 gradient-wind profile ``V(r)``;
* ``cyclostrophic``/``nocoriolis`` companion:

  ``mu = (r*V + 0.5*f*r**2) / (rm*Vmax + 0.5*f*rm**2)``

  ``Vd = rm*Vmax*mu/r`` (and ``Vd=0`` at ``r=0``).

The latter is the normalized-angular-momentum companion used by the public TCR
formulation; it is not ``M/r`` and is not a square-root conversion.
"""

from __future__ import annotations

from contextlib import contextmanager, redirect_stdout
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import importlib.metadata
import importlib.util
import inspect
import io
from pathlib import Path
import threading
from types import ModuleType
from typing import Any, Callable, Iterator, Mapping, Protocol
import warnings

import numpy as np


# Frozen source identity -----------------------------------------------------

C15_SOURCE_DOI = "10.4231/CZ4P-D448"
C15_SOURCE_VERSION = "v1.0"
C15_SOURCE_RELEASE = "2020-06-23"
C15_OFFICIAL_PYTHON_SHA256 = (
    "6f1306fae71d0e772f17dbf67d5c8cfd94fa543dd122cbb390c6d50161325113"
)
C15_PYTHON3_ADAPTER_SHA256 = (
    "823aeca59f1faa4ea118ea1f39d137d640a2611be053f77194749789da19f73b"
)

CLIMADA_PETALS_TARGET_VERSION = "6.2.0"
CLIMADA_PETALS_TARGET_COMMIT = "6ecd7af096f126df2da1023fbc5013765566d5e9"
CLIMADA_TCRAINFIELD_SHA256 = (
    "6f0bd30dc5532d907401a862f9d8b560c3feca6356c79ac7de444a52b315e062"
)

XI2020_H_TROP_M = 4000.0
XI2020_RHO_A_OVER_RHO_L = 0.0012

# Defaults are the values in the official C15 rmaxinput example.  They are
# explicit here rather than inherited from an unrelated wind/rain package.
C15_RMAXINPUT_DEFAULTS: Mapping[str, float | int] = {
    "Cdvary": 0,
    "C_d": 1.5e-3,
    "w_cool": 2.0e-3,
    "CkCdvary": 0,
    "CkCd": 1.0,
    "eye_adj": 0,
    "alpha_eye": 0.15,
}

# The official r0input and rmaxinput entry points expose the same physical
# parameters.  Keep a separately named immutable interface so metadata cannot
# silently describe one solver mode as the other.
C15_R0INPUT_DEFAULTS: Mapping[str, float | int] = dict(C15_RMAXINPUT_DEFAULTS)

C15_WIND_MODEL_NAME = "C15"
C15_WIND_MODEL_SENTINEL = 1515


class C15ProfileDomainError(ValueError):
    """Raised when an active TCR query lies beyond C15's finite outer radius."""


class C15ProfileInputError(ValueError):
    """Raised when a track row cannot define an official C15 profile."""


@dataclass(frozen=True)
class C15RadialProfile:
    """One cached official C15 radial-profile solution in SI units."""

    radius_m: np.ndarray
    gradient_wind_ms: np.ndarray
    outer_radius_m: float
    radius_max_wind_m: float
    maximum_wind_ms: float
    coriolis_s: float


class WindProfileProvider(Protocol):
    """Callable seam consumed by :func:`patched_c15_tcr`."""

    def __call__(
        self,
        si_track: object,
        radius_m: np.ndarray,
        active_mask: np.ndarray,
        *,
        cyclostrophic: bool,
    ) -> np.ndarray:
        """Return an array with shape ``(time, point)`` in m s-1."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


@lru_cache(maxsize=1)
def _load_frozen_c15_adapter() -> ModuleType:
    """Load the local frozen adapter without mutating ``sys.path``."""

    adapter_path = _project_root() / "adapters" / "c15_python3" / "c15.py"
    actual_hash = _sha256(adapter_path)
    if actual_hash != C15_PYTHON3_ADAPTER_SHA256:
        raise RuntimeError(
            "Frozen C15 adapter hash mismatch: "
            f"expected {C15_PYTHON3_ADAPTER_SHA256}, got {actual_hash}"
        )
    spec = importlib.util.spec_from_file_location("tcroad_frozen_c15_v1", adapter_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load frozen C15 adapter at {adapter_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _official_rmaxinput_solver(
    vmax_ms: float,
    rmax_m: float,
    coriolis_s: float,
    parameters: Mapping[str, float | int],
) -> tuple[np.ndarray, np.ndarray, float]:
    """Run the frozen official solver and return ``(r, V, r0)``."""

    c15 = _load_frozen_c15_adapter()
    # The official implementation prints convergence counters.  They are not
    # scientific output and are suppressed without changing warning policy or
    # numerical execution.
    try:
        with redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            radius_m, wind_ms, r0_m, _, _ = c15.ER11E04_nondim_rmaxinput(
                vmax_ms,
                rmax_m,
                abs(coriolis_s),
                parameters["Cdvary"],
                parameters["C_d"],
                parameters["w_cool"],
                parameters["CkCdvary"],
                parameters["CkCd"],
                parameters["eye_adj"],
                parameters["alpha_eye"],
            )
    except SystemExit as error:
        raise RuntimeError(
            "Official C15 rmaxinput aborted via SystemExit for "
            f"Vmax={vmax_ms:.9g} m s-1, Rmax={rmax_m:.9g} m, "
            f"|f|={abs(coriolis_s):.9g} s-1; exit={error.code!r}"
        ) from error
    return np.asarray(radius_m, dtype=float), np.asarray(wind_ms, dtype=float), float(r0_m)


def _official_r0input_solver(
    vmax_ms: float,
    r0_m: float,
    coriolis_s: float,
    parameters: Mapping[str, float | int],
) -> tuple[np.ndarray, np.ndarray, float]:
    """Run official r0input and return ``(r, V, rmax)``.

    The frozen Python adapter's exact return order is
    ``rr, VV, rmerge, Vmerge, rmax``.  The explicit unpacking below prevents
    the merge-radius fields from being mistaken for the inferred RMW.
    """

    c15 = _load_frozen_c15_adapter()
    try:
        with redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            radius_m, wind_ms, _, _, rmax_m = c15.ER11E04_nondim_r0input(
                vmax_ms,
                r0_m,
                abs(coriolis_s),
                parameters["Cdvary"],
                parameters["C_d"],
                parameters["w_cool"],
                parameters["CkCdvary"],
                parameters["CkCd"],
                parameters["eye_adj"],
                parameters["alpha_eye"],
            )
    except SystemExit as error:
        raise RuntimeError(
            "Official C15 r0input aborted via SystemExit for "
            f"Vmax={vmax_ms:.9g} m s-1, r0={r0_m:.9g} m, "
            f"|f|={abs(coriolis_s):.9g} s-1; exit={error.code!r}"
        ) from error
    return (
        np.asarray(radius_m, dtype=float),
        np.asarray(wind_ms, dtype=float),
        float(rmax_m),
    )


class C15WindProfileProvider:
    """Cached official-C15 provider for CLIMADA Petals TCR radial queries.

    Parameters
    ----------
    parameters
        Explicit official ``rmaxinput`` parameter values.  ``eye_adj`` must
        remain zero because the frozen official Python release does not expose
        a validated active eye-adjustment branch.
    solver
        Test seam only.  Production callers leave this as ``None`` to run the
        hash-verified frozen official C15 adapter.
    """

    def __init__(
        self,
        parameters: Mapping[str, float | int] | None = None,
        *,
        solver: Callable[
            [float, float, float, Mapping[str, float | int]],
            tuple[np.ndarray, np.ndarray, float],
        ]
        | None = None,
    ) -> None:
        merged = dict(C15_RMAXINPUT_DEFAULTS)
        if parameters is not None:
            unknown = set(parameters).difference(merged)
            if unknown:
                raise C15ProfileInputError(
                    f"Unknown official C15 rmaxinput parameters: {sorted(unknown)}"
                )
            merged.update(parameters)
        if merged["eye_adj"] != 0:
            raise C15ProfileInputError(
                "eye_adj=1 is not validated in the frozen official Python adapter"
            )
        self.parameters = merged
        self._solver = solver or _official_rmaxinput_solver
        self._profiles: dict[tuple[float, float, float], C15RadialProfile] = {}
        self._cache_lock = threading.RLock()
        self._cache_hits = 0
        self._cache_misses = 0

    @property
    def cache_info(self) -> tuple[int, int, int]:
        """Return ``(hits, misses, current_size)`` for the profile cache."""

        with self._cache_lock:
            return self._cache_hits, self._cache_misses, len(self._profiles)

    def profile_for(
        self, vmax_ms: float, rmax_m: float, coriolis_s: float
    ) -> C15RadialProfile:
        """Return one official profile, solving it only once per exact SI tuple."""

        vmax_ms = float(vmax_ms)
        rmax_m = float(rmax_m)
        coriolis_s = abs(float(coriolis_s))
        values = np.asarray([vmax_ms, rmax_m, coriolis_s])
        if not np.all(np.isfinite(values)):
            raise C15ProfileInputError("C15 inputs must be finite")
        if vmax_ms <= 0 or rmax_m <= 0 or coriolis_s <= 0:
            raise C15ProfileInputError(
                "C15 requires Vmax>0, Rmax>0, and a non-zero Coriolis magnitude"
            )

        key = (vmax_ms, rmax_m, coriolis_s)
        with self._cache_lock:
            cached = self._profiles.get(key)
            if cached is not None:
                self._cache_hits += 1
                return cached

            radius_m, wind_ms, r0_m = self._solver(
                vmax_ms, rmax_m, coriolis_s, self.parameters
            )
            if (
                radius_m.ndim != 1
                or wind_ms.shape != radius_m.shape
                or radius_m.size < 2
                or not np.all(np.isfinite(radius_m))
                or not np.all(np.isfinite(wind_ms))
                or not np.all(np.diff(radius_m) > 0)
                or not np.isfinite(r0_m)
                or r0_m < radius_m[-1]
            ):
                raise RuntimeError("Official C15 rmaxinput returned an invalid radial profile")

            # r0 is the official finite outer radius and V(r0)=0 by definition.
            # The solver's regular radial grid can stop just inside r0; adding
            # that exact defining point is neither a tail nor a taper.
            if r0_m > radius_m[-1]:
                radius_m = np.append(radius_m, r0_m)
                wind_ms = np.append(wind_ms, 0.0)

            profile = C15RadialProfile(
                radius_m=radius_m,
                gradient_wind_ms=wind_ms,
                outer_radius_m=r0_m,
                radius_max_wind_m=rmax_m,
                maximum_wind_ms=vmax_ms,
                coriolis_s=coriolis_s,
            )
            self._profiles[key] = profile
            self._cache_misses += 1
            return profile

    def __call__(
        self,
        si_track: object,
        radius_m: np.ndarray,
        active_mask: np.ndarray,
        *,
        cyclostrophic: bool,
    ) -> np.ndarray:
        """Evaluate official C15 only at active TCR points.

        The first row is explicitly zeroed to preserve the behavior of
        CLIMADA's ``compute_angular_windspeeds`` wrapper, including calls on
        time-shifted track slices made by TCR's finite differences.
        """

        radius_m = np.asarray(radius_m, dtype=float)
        active_mask = np.asarray(active_mask, dtype=bool)
        if radius_m.ndim != 2 or active_mask.shape != radius_m.shape:
            raise ValueError("radius_m and active_mask must share shape (time, point)")

        ntime = radius_m.shape[0]
        for name in ("vmax", "rad", "cp"):
            if name not in si_track:
                raise C15ProfileInputError(f"SI track is missing required variable {name!r}")
            if np.asarray(si_track[name].values).shape != (ntime,):
                raise C15ProfileInputError(
                    f"SI track variable {name!r} must have shape ({ntime},)"
                )

        result = np.zeros_like(radius_m, dtype=float)
        # CLIMADA intentionally zeros the first row for every wind-profile
        # call.  Start at 1 so no C15 solve/query occurs for that row.
        for time_index in range(1, ntime):
            row_mask = active_mask[time_index]
            if not np.any(row_mask):
                continue
            query_m = radius_m[time_index, row_mask]
            if not np.all(np.isfinite(query_m)) or np.any(query_m < 0):
                raise C15ProfileInputError("Active C15 radii must be finite and non-negative")

            vmax_ms = float(si_track["vmax"].values[time_index])
            rmax_m = float(si_track["rad"].values[time_index])
            coriolis_s = abs(float(si_track["cp"].values[time_index]))
            profile = self.profile_for(vmax_ms, rmax_m, coriolis_s)
            outside = query_m > profile.outer_radius_m
            if np.any(outside):
                maximum_query = float(np.max(query_m[outside]))
                raise C15ProfileDomainError(
                    "Active TCR radius exceeds the finite official C15 domain: "
                    f"query={maximum_query:.6f} m, r0={profile.outer_radius_m:.6f} m, "
                    f"time_index={time_index}"
                )

            gradient_ms = np.interp(
                query_m, profile.radius_m, profile.gradient_wind_ms
            )
            if not cyclostrophic:
                result[time_index, row_mask] = gradient_ms
                continue

            companion_ms = np.zeros_like(query_m)
            nonzero = query_m > 0
            radius_nonzero = query_m[nonzero]
            absolute_momentum = (
                radius_nonzero * gradient_ms[nonzero]
                + 0.5 * coriolis_s * radius_nonzero**2
            )
            maximum_momentum = (
                rmax_m * vmax_ms + 0.5 * coriolis_s * rmax_m**2
            )
            mu = absolute_momentum / maximum_momentum
            companion_ms[nonzero] = rmax_m * vmax_ms * mu / radius_nonzero
            result[time_index, row_mask] = companion_ms

        return result


class C15FixedR0WindProfileProvider:
    """Official-C15 ``r0input`` provider with one outer radius per event.

    ``outer_radius_m`` is immutable for the provider's lifetime.  C15 infers a
    potentially different RMW from each track row's circular maximum wind and
    Coriolis magnitude.  Active queries beyond the finite event r0 are rejected
    rather than padded, tapered, or extrapolated.

    Parameters
    ----------
    outer_radius_m
        Fixed event outer radius (C15 ``r0``) in metres.
    parameters
        Explicit official ``r0input`` parameter values.  ``eye_adj`` must
        remain zero because the frozen official Python release does not expose
        a validated active eye-adjustment branch.
    solver
        Test seam only.  Production callers leave this as ``None`` to run the
        hash-verified frozen official C15 adapter.
    """

    def __init__(
        self,
        outer_radius_m: float,
        parameters: Mapping[str, float | int] | None = None,
        *,
        solver: Callable[
            [float, float, float, Mapping[str, float | int]],
            tuple[np.ndarray, np.ndarray, float],
        ]
        | None = None,
    ) -> None:
        outer_radius_m = float(outer_radius_m)
        if not np.isfinite(outer_radius_m) or outer_radius_m <= 0:
            raise C15ProfileInputError(
                "C15 fixed event outer_radius_m must be finite and positive"
            )
        merged = dict(C15_R0INPUT_DEFAULTS)
        if parameters is not None:
            unknown = set(parameters).difference(merged)
            if unknown:
                raise C15ProfileInputError(
                    f"Unknown official C15 r0input parameters: {sorted(unknown)}"
                )
            merged.update(parameters)
        if merged["eye_adj"] != 0:
            raise C15ProfileInputError(
                "eye_adj=1 is not validated in the frozen official Python adapter"
            )
        self.outer_radius_m = outer_radius_m
        self.parameters = merged
        self._solver = solver or _official_r0input_solver
        self._profiles: dict[tuple[float, float], C15RadialProfile] = {}
        self._cache_lock = threading.RLock()
        self._cache_hits = 0
        self._cache_misses = 0

    @property
    def cache_info(self) -> tuple[int, int, int]:
        """Return ``(hits, misses, current_size)`` for the profile cache."""

        with self._cache_lock:
            return self._cache_hits, self._cache_misses, len(self._profiles)

    def profile_for(
        self, vmax_ms: float, coriolis_s: float
    ) -> C15RadialProfile:
        """Return one r0input profile with C15-inferred, row-specific RMW."""

        vmax_ms = float(vmax_ms)
        coriolis_s = abs(float(coriolis_s))
        values = np.asarray([vmax_ms, self.outer_radius_m, coriolis_s])
        if not np.all(np.isfinite(values)):
            raise C15ProfileInputError("C15 inputs must be finite")
        if vmax_ms <= 0 or coriolis_s <= 0:
            raise C15ProfileInputError(
                "C15 r0input requires Vmax>0, r0>0, and a non-zero Coriolis magnitude"
            )

        key = (vmax_ms, coriolis_s)
        with self._cache_lock:
            cached = self._profiles.get(key)
            if cached is not None:
                self._cache_hits += 1
                return cached

            radius_m, wind_ms, rmax_m = self._solver(
                vmax_ms, self.outer_radius_m, coriolis_s, self.parameters
            )
            if (
                radius_m.ndim != 1
                or wind_ms.shape != radius_m.shape
                or radius_m.size < 2
                or not np.all(np.isfinite(radius_m))
                or not np.all(np.isfinite(wind_ms))
                or not np.all(np.diff(radius_m) > 0)
                or not np.isfinite(rmax_m)
                or rmax_m <= 0
                or rmax_m >= self.outer_radius_m
                or radius_m[-1] > self.outer_radius_m
            ):
                raise RuntimeError("Official C15 r0input returned an invalid radial profile")

            # C15 defines V(r0)=0.  Its regular r/rmax grid normally stops just
            # inside r0, so add the exact defining endpoint.  This is neither a
            # tail nor a taper and no value beyond r0 is manufactured.
            if radius_m[-1] < self.outer_radius_m:
                radius_m = np.append(radius_m, self.outer_radius_m)
                wind_ms = np.append(wind_ms, 0.0)
            else:
                radius_m = radius_m.copy()
                wind_ms = wind_ms.copy()
                radius_m[-1] = self.outer_radius_m
                wind_ms[-1] = 0.0

            profile = C15RadialProfile(
                radius_m=radius_m,
                gradient_wind_ms=wind_ms,
                outer_radius_m=self.outer_radius_m,
                radius_max_wind_m=rmax_m,
                maximum_wind_ms=vmax_ms,
                coriolis_s=coriolis_s,
            )
            self._profiles[key] = profile
            self._cache_misses += 1
            return profile

    def __call__(
        self,
        si_track: object,
        radius_m: np.ndarray,
        active_mask: np.ndarray,
        *,
        cyclostrophic: bool,
    ) -> np.ndarray:
        """Evaluate fixed-r0 C15 profiles only at active TCR points."""

        radius_m = np.asarray(radius_m, dtype=float)
        active_mask = np.asarray(active_mask, dtype=bool)
        if radius_m.ndim != 2 or active_mask.shape != radius_m.shape:
            raise ValueError("radius_m and active_mask must share shape (time, point)")

        ntime = radius_m.shape[0]
        for name in ("vmax", "cp"):
            if name not in si_track:
                raise C15ProfileInputError(f"SI track is missing required variable {name!r}")
            if np.asarray(si_track[name].values).shape != (ntime,):
                raise C15ProfileInputError(
                    f"SI track variable {name!r} must have shape ({ntime},)"
                )

        result = np.zeros_like(radius_m, dtype=float)
        # Preserve CLIMADA compute_angular_windspeeds semantics, including its
        # finite-difference calls on shifted track slices.
        for time_index in range(1, ntime):
            row_mask = active_mask[time_index]
            if not np.any(row_mask):
                continue
            query_m = radius_m[time_index, row_mask]
            if not np.all(np.isfinite(query_m)) or np.any(query_m < 0):
                raise C15ProfileInputError("Active C15 radii must be finite and non-negative")
            outside = query_m > self.outer_radius_m
            if np.any(outside):
                maximum_query = float(np.max(query_m[outside]))
                raise C15ProfileDomainError(
                    "Active TCR radius exceeds the fixed finite official C15 domain: "
                    f"query={maximum_query:.6f} m, r0={self.outer_radius_m:.6f} m, "
                    f"time_index={time_index}"
                )

            vmax_ms = float(si_track["vmax"].values[time_index])
            coriolis_s = abs(float(si_track["cp"].values[time_index]))
            profile = self.profile_for(vmax_ms, coriolis_s)
            gradient_ms = np.interp(
                query_m, profile.radius_m, profile.gradient_wind_ms
            )
            if not cyclostrophic:
                result[time_index, row_mask] = gradient_ms
                continue

            rmax_m = profile.radius_max_wind_m
            companion_ms = np.zeros_like(query_m)
            nonzero = query_m > 0
            radius_nonzero = query_m[nonzero]
            absolute_momentum = (
                radius_nonzero * gradient_ms[nonzero]
                + 0.5 * coriolis_s * radius_nonzero**2
            )
            maximum_momentum = (
                rmax_m * vmax_ms + 0.5 * coriolis_s * rmax_m**2
            )
            mu = absolute_momentum / maximum_momentum
            companion_ms[nonzero] = rmax_m * vmax_ms * mu / radius_nonzero
            result[time_index, row_mask] = companion_ms

        return result


_PATCH_LOCK = threading.RLock()


def _validate_tc_rainfield_module(tc_rainfield: ModuleType) -> None:
    module_path = Path(inspect.getsourcefile(tc_rainfield) or "")
    if not module_path.is_file():
        raise RuntimeError("Cannot identify CLIMADA Petals tc_rainfield.py")
    actual_hash = _sha256(module_path)
    if actual_hash != CLIMADA_TCRAINFIELD_SHA256:
        try:
            installed = importlib.metadata.version("climada_petals")
        except importlib.metadata.PackageNotFoundError:
            installed = "unknown"
        raise RuntimeError(
            "Unsupported CLIMADA Petals private TCR seam: "
            f"installed={installed}, tc_rainfield.py sha256={actual_hash}; "
            f"expected target {CLIMADA_PETALS_TARGET_VERSION} sha256="
            f"{CLIMADA_TCRAINFIELD_SHA256}"
        )
    expected_parameters = (
        "si_track",
        "d_centr",
        "mask_centr_close",
        "model",
        "cyclostrophic",
        "matlab_ref_mode",
    )
    actual_parameters = tuple(inspect.signature(tc_rainfield._windprofile).parameters)
    if actual_parameters != expected_parameters:
        raise RuntimeError(
            "CLIMADA Petals _windprofile signature changed: "
            f"expected {expected_parameters}, got {actual_parameters}"
        )


def assert_environmental_pressure_schema_only() -> None:
    """Assert that the selected C15--TCR path never consumes pressure deficit.

    CLIMADA's generic track-to-SI converter requires central and environmental
    pressures and derives ``pdelta``.  The public TCR implementation selected
    here must not read any of them after conversion.  This static assertion is
    intentionally hash-pinned to Petals 6.2.0 and avoids a second science run.
    """

    from climada_petals.hazard import tc_rainfield

    _validate_tc_rainfield_module(tc_rainfield)
    functions = (
        tc_rainfield.compute_rain,
        tc_rainfield._tcr,
        tc_rainfield._compute_vertical_velocity,
        tc_rainfield._horizontal_winds,
        tc_rainfield._windprofile,
        tc_rainfield._w_frict_stretch,
        tc_rainfield._w_topo,
        tc_rainfield._w_shear,
    )
    forbidden = ("pdelta", "environmental_pressure", '["env"]', "['env']")
    for function in functions:
        source = inspect.getsource(function)
        hits = [token for token in forbidden if token in source]
        if hits:
            raise RuntimeError(
                f"CLIMADA TCR path {function.__name__} reads schema-only "
                f"environmental pressure tokens: {hits}"
            )


@contextmanager
def _patched_xi2020_physical_constants(
    *, lower_troposphere_height_m: float, rho_air_over_rho_liquid: float
) -> Iterator[None]:
    """Temporarily apply the two explicit Xi/Lu TCR physical constants."""

    from climada_petals.hazard import tc_rainfield

    _validate_tc_rainfield_module(tc_rainfield)
    if lower_troposphere_height_m != XI2020_H_TROP_M:
        raise ValueError(
            f"this reconstruction requires H_TROP={XI2020_H_TROP_M:g} m"
        )
    if rho_air_over_rho_liquid != XI2020_RHO_A_OVER_RHO_L:
        raise ValueError(
            "this reconstruction requires rho_air/rho_liquid="
            f"{XI2020_RHO_A_OVER_RHO_L:g}"
        )
    if float(tc_rainfield.H_TROP) != XI2020_H_TROP_M:
        raise RuntimeError(
            f"unexpected CLIMADA H_TROP={tc_rainfield.H_TROP!r}; "
            f"expected {XI2020_H_TROP_M:g} m"
        )

    original_h_trop = tc_rainfield.H_TROP
    original_density_ratio = tc_rainfield.RHO_A_OVER_RHO_L
    try:
        tc_rainfield.H_TROP = lower_troposphere_height_m
        tc_rainfield.RHO_A_OVER_RHO_L = rho_air_over_rho_liquid
        yield
    finally:
        tc_rainfield.H_TROP = original_h_trop
        tc_rainfield.RHO_A_OVER_RHO_L = original_density_ratio


def run_tcr_public_reconstruction(
    *,
    track: object,
    centroid_lat: np.ndarray,
    centroid_lon: np.ndarray,
    elevation_tif: Path | str,
    c_drag_tif: Path | str,
    e_precip: float,
    lower_troposphere_height_m: float,
    rho_air_over_rho_liquid: float,
    max_w_foreground: float,
    res_radial_m: float,
    min_c_drag: float,
    max_dist_eye_km: float,
    provider: WindProfileProvider | None = None,
) -> dict[str, Any]:
    """Run one C15-driven public CLIMADA-Petals TCR reconstruction.

    The wrapper intentionally calls Petals' public numerical functions once,
    while retaining the supplied union-grid ordering.  It returns the raw
    one-hour rain-rate matrix; accumulation semantics remain the runner's
    explicit responsibility.  The default ``rmaxinput`` provider is retained
    solely for the observed-RMW Irene benchmark.  Synthetic-event callers must
    inject their explicitly constructed provider, such as
    :class:`C15FixedR0WindProfileProvider`.
    """

    from climada.hazard import Centroids
    from climada_petals.hazard import tc_rainfield

    _validate_tc_rainfield_module(tc_rainfield)
    assert_environmental_pressure_schema_only()
    elevation_tif = Path(elevation_tif).resolve()
    c_drag_tif = Path(c_drag_tif).resolve()
    if not elevation_tif.is_file() or not c_drag_tif.is_file():
        raise FileNotFoundError("explicit elevation and drag GeoTIFFs are required")
    if max_dist_eye_km != 300.0:
        raise ValueError("this published reconstruction requires a 300-km TCR domain")
    latitude = np.asarray(centroid_lat, dtype=float)
    longitude = np.asarray(centroid_lon, dtype=float)
    if latitude.ndim != 1 or longitude.shape != latitude.shape:
        raise ValueError("centroid_lat and centroid_lon must share one-dimensional shape")
    centroids = Centroids.from_lat_lon(latitude, longitude)
    all_indices = np.arange(latitude.size, dtype=np.int64)
    model_kwargs = {
        "wind_model": C15_WIND_MODEL_NAME,
        "elevation_tif": elevation_tif,
        "c_drag_tif": c_drag_tif,
        "e_precip": float(e_precip),
        "max_w_foreground": float(max_w_foreground),
        "res_radial_m": float(res_radial_m),
        "min_c_drag": float(min_c_drag),
    }

    provider_was_injected = provider is not None
    active_provider: WindProfileProvider = provider or C15WindProfileProvider()
    with _PATCH_LOCK:
        with _patched_xi2020_physical_constants(
            lower_troposphere_height_m=lower_troposphere_height_m,
            rho_air_over_rho_liquid=rho_air_over_rho_liquid,
        ):
            with patched_c15_tcr(active_provider):
                _, rainrates_sparse = tc_rainfield._compute_rain_sparse(
                    track=track,
                    centroids=centroids,
                    idx_centr_filter=all_indices,
                    model="TCR",
                    model_kwargs=model_kwargs,
                    store_rainrates=True,
                    metric="geosphere",
                    intensity_thres=0.0,
                    max_dist_eye_km=max_dist_eye_km,
                )
    if rainrates_sparse is None:
        raise RuntimeError("CLIMADA did not return requested rain rates")
    rainrate = np.asarray(rainrates_sparse.toarray(), dtype=float)
    cache_info = getattr(active_provider, "cache_info", (0, 0, 0))
    hits, misses, size = cache_info
    provider_class = (
        f"{type(active_provider).__module__}.{type(active_provider).__qualname__}"
    )
    if isinstance(active_provider, C15FixedR0WindProfileProvider):
        profile_metadata: dict[str, Any] = {
            "wind_profile": "official C15 r0input with fixed event outer radius",
            "c15_input_mode": "r0input",
            "c15_outer_radius_m": active_provider.outer_radius_m,
            "c15_r0input_parameters": dict(active_provider.parameters),
            "c15_rmax_inference": "official r0input; varies with track-row Vmax and |f|",
            "c15_radial_domain": "0 <= r <= fixed event r0; no tail or taper",
        }
    elif isinstance(active_provider, C15WindProfileProvider):
        profile_metadata = {
            "wind_profile": "official C15 rmaxinput",
            "c15_input_mode": "rmaxinput",
            "c15_rmaxinput_parameters": dict(active_provider.parameters),
            "c15_radial_domain": "0 <= r <= row-specific inferred r0; no tail or taper",
        }
    else:
        profile_metadata = {
            "wind_profile": f"injected provider {provider_class}",
            "c15_input_mode": "injected_provider",
        }
    return {
        "rainfall_rate_mm_h": rainrate,
        "metadata": {
            **profile_metadata,
            "wind_profile_provider_class": provider_class,
            "wind_profile_provider_injected": provider_was_injected,
            "legacy_rmaxinput_default_for_irene_compatibility": (
                not provider_was_injected
                and isinstance(active_provider, C15WindProfileProvider)
            ),
            "bridge_module_sha256": _sha256(Path(__file__).resolve()),
            "c15_source_doi": C15_SOURCE_DOI,
            "c15_source_version": C15_SOURCE_VERSION,
            "c15_official_python_sha256": C15_OFFICIAL_PYTHON_SHA256,
            "c15_python3_adapter_sha256": C15_PYTHON3_ADAPTER_SHA256,
            "climada_petals_version": CLIMADA_PETALS_TARGET_VERSION,
            "climada_petals_commit": CLIMADA_PETALS_TARGET_COMMIT,
            "climada_tc_rainfield_sha256": CLIMADA_TCRAINFIELD_SHA256,
            "metric": "geosphere",
            "h_trop_m": lower_troposphere_height_m,
            "rho_air_over_rho_liquid": rho_air_over_rho_liquid,
            "profile_cache_hits": hits,
            "profile_cache_misses": misses,
            "profile_cache_size": size,
            "environmental_pressure_schema_only_asserted": True,
        },
    }
@contextmanager
def patched_c15_tcr(
    provider: WindProfileProvider | None = None,
) -> Iterator[WindProfileProvider]:
    """Temporarily inject a ``C15`` wind-profile dispatch into Petals TCR.

    Only the hash-pinned private ``_windprofile`` seam is replaced, and only
    calls dispatched with ``wind_model='C15'`` are intercepted.  ER11 and other
    CLIMADA wind models retain their original implementation.  The function and
    model mapping are restored in ``finally``; no site-package file is edited.

    This is a process-global monkeypatch and must surround a single TCR run, not
    concurrent TCR calls in multiple threads.
    """

    from climada_petals.hazard import tc_rainfield

    _validate_tc_rainfield_module(tc_rainfield)
    active_provider: WindProfileProvider = provider or C15WindProfileProvider()

    with _PATCH_LOCK:
        original_windprofile = tc_rainfield._windprofile
        model_mapping = tc_rainfield.MODEL_VANG
        had_c15 = C15_WIND_MODEL_NAME in model_mapping
        previous_c15 = model_mapping.get(C15_WIND_MODEL_NAME)
        if had_c15 and previous_c15 != C15_WIND_MODEL_SENTINEL:
            raise RuntimeError(
                "CLIMADA already defines a different C15 wind-model dispatch; "
                "refusing to shadow it"
            )

        def _patched_windprofile(
            si_track: object,
            d_centr: np.ndarray,
            mask_centr_close: np.ndarray,
            model: int,
            cyclostrophic: bool = False,
            matlab_ref_mode: bool = False,
        ) -> np.ndarray:
            if model != C15_WIND_MODEL_SENTINEL:
                return original_windprofile(
                    si_track,
                    d_centr,
                    mask_centr_close,
                    model,
                    cyclostrophic=cyclostrophic,
                    matlab_ref_mode=matlab_ref_mode,
                )
            if matlab_ref_mode:
                raise ValueError(
                    "matlab_ref_mode is an ER11 reference option and is not part of "
                    "the public C15 bridge"
                )
            result = np.asarray(
                active_provider(
                    si_track,
                    d_centr,
                    mask_centr_close,
                    cyclostrophic=cyclostrophic,
                ),
                dtype=float,
            )
            if result.shape != np.asarray(d_centr).shape:
                raise RuntimeError(
                    "C15 provider returned shape "
                    f"{result.shape}, expected {np.asarray(d_centr).shape}"
                )
            result = result.copy()
            if result.shape[0]:
                result[0, :] = 0.0
            return result

        try:
            model_mapping[C15_WIND_MODEL_NAME] = C15_WIND_MODEL_SENTINEL
            tc_rainfield._windprofile = _patched_windprofile
            yield active_provider
        finally:
            tc_rainfield._windprofile = original_windprofile
            if had_c15:
                model_mapping[C15_WIND_MODEL_NAME] = previous_c15
            else:
                model_mapping.pop(C15_WIND_MODEL_NAME, None)


__all__ = [
    "C15ProfileDomainError",
    "C15ProfileInputError",
    "C15FixedR0WindProfileProvider",
    "C15RadialProfile",
    "C15WindProfileProvider",
    "C15_OFFICIAL_PYTHON_SHA256",
    "C15_PYTHON3_ADAPTER_SHA256",
    "C15_R0INPUT_DEFAULTS",
    "C15_RMAXINPUT_DEFAULTS",
    "C15_SOURCE_DOI",
    "C15_SOURCE_RELEASE",
    "C15_SOURCE_VERSION",
    "C15_WIND_MODEL_NAME",
    "CLIMADA_PETALS_TARGET_COMMIT",
    "CLIMADA_PETALS_TARGET_VERSION",
    "CLIMADA_TCRAINFIELD_SHA256",
    "XI2020_H_TROP_M",
    "XI2020_RHO_A_OVER_RHO_L",
    "assert_environmental_pressure_schema_only",
    "patched_c15_tcr",
    "run_tcr_public_reconstruction",
]
