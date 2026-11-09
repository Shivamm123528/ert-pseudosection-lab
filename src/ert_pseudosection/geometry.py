"""Survey geometry and half-space calculations for four common ERT arrays.

All electrodes are collinear on a flat surface, at x = 0, d, 2d, ... .
The negative electrode index -1 means an ideal, infinitely remote current
return electrode (pole-dipole only). None of the depths are inversion results.

Median-depth factors: M. H. Loke (1999), *Electrical imaging surveys for
 environmental and engineering studies*, Table 2 (after Edwards, 1977).
https://www.geometrics.com/wp-content/uploads/2018/10/Lokenote.pdf
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi
from typing import Literal

import numpy as np
import pandas as pd

ArrayName = Literal[
    "wenner", "wenner_schlumberger", "dipole_dipole", "pole_dipole"
]
PoleDirection = Literal["forward", "reverse", "both"]
ARRAYS: tuple[str, ...] = (
    "wenner", "wenner_schlumberger", "dipole_dipole", "pole_dipole"
)
ARRAY_LABELS = {
    "wenner": "Wenner-alpha",
    "wenner_schlumberger": "Wenner–Schlumberger",
    "dipole_dipole": "Dipole–Dipole",
    "pole_dipole": "Pole–Dipole",
}

# Published *rounded* ze/a coefficients for 1 <= n <= 6. Values at n = 1
# for Wenner-Schlumberger and pole-dipole are rounded to 0.52 in Loke's table.
# Wenner-alpha is recorded separately at its published precision (0.519).
MEDIAN_DEPTH_FACTORS = {
    "wenner": {1: 0.519},
    "wenner_schlumberger": {
        1: 0.520, 2: 0.930, 3: 1.320, 4: 1.710, 5: 2.090, 6: 2.480
    },
    "dipole_dipole": {
        1: 0.416, 2: 0.697, 3: 0.962, 4: 1.220, 5: 1.476, 6: 1.730
    },
    "pole_dipole": {
        1: 0.520, 2: 0.930, 3: 1.320, 4: 1.710, 5: 2.090, 6: 2.480
    },
}

COLUMNS = [
    "measurement_id", "array", "direction", "a_multiplier", "n",
    "A_idx", "B_idx", "M_idx", "N_idx",
    "x_A_m", "x_B_m", "x_M_m", "x_N_m",
    "a_m", "array_span_m", "plot_x_m", "pseudodepth_m",
    "k_signed_m", "k_magnitude_m", "current_a", "voltage_v",
    "apparent_resistivity_ohm_m",
]


@dataclass(frozen=True)
class SurveyConfig:
    """The survey has equally spaced electrodes on a horizontal straight line.

    max_a_multiplier: 'a' = a_multiplier * spacing_m. For Wenner, this
        controls its equal electrode spacing. For other arrays, it controls
        their dipole/potential-pair length; n controls the separate gap.
    max_n: maximum separation factor for DD, WS, PD. Loke's tabulated
        median-depth factors extend to n=6 in the 1999 table used here.
    """

    n_electrodes: int = 32
    spacing_m: float = 5.0
    max_a_multiplier: int = 2
    max_n: int = 6
    pole_direction: PoleDirection = "forward"

    def __post_init__(self) -> None:
        for name in ("n_electrodes", "max_a_multiplier", "max_n"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} must be an integer")
        if not (4 <= self.n_electrodes <= 256):
            raise ValueError("n_electrodes must be between 4 and 256")
        if not (isfinite(self.spacing_m) and self.spacing_m > 0):
            raise ValueError("spacing_m must be finite and positive")
        if not (1 <= self.max_a_multiplier <= 32):
            raise ValueError("max_a_multiplier must be between 1 and 32")
        if not (1 <= self.max_n <= 6):
            raise ValueError("max_n must be between 1 and 6 (published coefficients)")
        if self.pole_direction not in ("forward", "reverse", "both"):
            raise ValueError("pole_direction must be 'forward', 'reverse', or 'both'")

    @property
    def length_m(self) -> float:
        return (self.n_electrodes - 1) * self.spacing_m

    @property
    def electrode_positions_m(self) -> np.ndarray:
        return np.arange(self.n_electrodes, dtype=float) * self.spacing_m


def geometric_factor_m(
    x_a: float, x_b: float, x_m: float, x_n: float
) -> float:
    """Signed K (m) for a flat, homogeneous, isotropic half-space.

    x_b = +infinity denotes the remote return electrode. Keep the *sign*:
    rho_a = K * (V_M - V_N) / I. A-B-M-N dipole-dipole has negative K with
    this electrode numbering, and its expected measured delta-V is negative.
    """
    coordinates = (x_a, x_b, x_m, x_n)
    if any(np.isnan(x) for x in coordinates):
        raise ValueError("Electrode positions cannot be NaN")
    if not all(np.isfinite(x) or (i == 1 and x == float("inf"))
               for i, x in enumerate(coordinates)):
        raise ValueError("Only B may be an ideal remote electrode (+infinity)")

    def reciprocal_distance(a: float, b: float) -> float:
        r = abs(a - b)
        if r == 0:
            raise ValueError("Current and potential electrodes cannot coincide")
        return 0.0 if not np.isfinite(r) else 1.0 / r

    denom = (
        reciprocal_distance(x_a, x_m)
        - reciprocal_distance(x_a, x_n)
        - reciprocal_distance(x_b, x_m)
        + reciprocal_distance(x_b, x_n)
    )
    if abs(denom) < 1e-14:
        raise ValueError("Degenerate array: zero geometrical voltage")
    return 2.0 * pi / denom


def median_pseudodepth_m(array: ArrayName, n: int, a_m: float) -> float:
    """Homogeneous-earth median-depth plotting convention (NOT measured DOI)."""
    if array not in MEDIAN_DEPTH_FACTORS:
        raise ValueError(f"Unknown ERT array: {array}")
    if n not in MEDIAN_DEPTH_FACTORS[array]:
        raise ValueError(f"No published depth coefficient for {array}, n={n}")
    if not (isfinite(a_m) and a_m > 0):
        raise ValueError("a_m must be finite and positive")
    return a_m * MEDIAN_DEPTH_FACTORS[array][n]


def _new_row(
    config: SurveyConfig, array: ArrayName, a_idx: int, b_idx: int,
    m_idx: int, n_idx: int, a_multiplier: int, n: int,
    direction: str,
) -> dict:
    ids = [a_idx, m_idx, n_idx] + ([] if b_idx == -1 else [b_idx])
    if any(index < 0 or index >= config.n_electrodes for index in ids):
        raise ValueError("Electrode index out of range")
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicated electrode in a measurement")
    x_a, x_m, x_n = (index * config.spacing_m for index in (a_idx, m_idx, n_idx))
    x_b = (b_idx * config.spacing_m) if b_idx >= 0 else float("inf")
    active = [x_a, x_m, x_n] + ([] if b_idx < 0 else [x_b])
    a_m = a_multiplier * config.spacing_m
    signed_k = geometric_factor_m(x_a, x_b, x_m, x_n)
    return {
        "array": array, "direction": direction,
        "a_multiplier": a_multiplier, "n": n,
        "A_idx": a_idx, "B_idx": b_idx, "M_idx": m_idx, "N_idx": n_idx,
        "x_A_m": x_a, "x_B_m": np.nan if b_idx == -1 else x_b,
        "x_M_m": x_m, "x_N_m": x_n,
        "a_m": a_m,
        "array_span_m": max(active) - min(active),
        "plot_x_m": (max(active) + min(active)) / 2.0,
        "pseudodepth_m": median_pseudodepth_m(array, n, a_m),
        "k_signed_m": signed_k, "k_magnitude_m": abs(signed_k),
        "current_a": np.nan, "voltage_v": np.nan,
        "apparent_resistivity_ohm_m": np.nan,
    }


def generate_survey(config: SurveyConfig, array: ArrayName) -> pd.DataFrame:
    """Generate all configurations that fit, deduplicated by array parameters.

    Wenner: A --a-- M --a-- N --a-- B.
    Wenner-Schlumberger: A --na-- M --a-- N --na-- B.
    Dipole-Dipole: A --a-- B --na-- M --a-- N.
    Pole-Dipole forward: A --na-- M --a-- N; B is remote.
    Pole-Dipole reverse: N --a-- M --na-- A; B is remote.
    """
    if array not in ARRAYS:
        raise ValueError(f"Unknown ERT array: {array}")
    records: list[dict] = []
    size = config.n_electrodes
    for step in range(1, config.max_a_multiplier + 1):
        n_options = [1] if array == "wenner" else range(1, config.max_n + 1)
        for n in n_options:
            if array == "wenner":
                width = 3 * step
            elif array == "wenner_schlumberger":
                width = (2 * n + 1) * step
            elif array == "dipole_dipole":
                width = (n + 2) * step
            else:
                width = (n + 1) * step

            if width >= size:
                continue
            for start in range(size - width):
                if array == "wenner":
                    specs = [(start, start + 3*step, start + step,
                              start + 2*step, "symmetric")]
                elif array == "wenner_schlumberger":
                    specs = [(start, start + width, start + n*step,
                              start + (n+1)*step, "symmetric")]
                elif array == "dipole_dipole":
                    specs = [(start, start + step, start + (n+1)*step,
                              start + (n+2)*step, "forward")]
                else:
                    specs = []
                    if config.pole_direction in ("forward", "both"):
                        specs.append((start, -1, start + n*step,
                                      start + width, "forward"))
                    if config.pole_direction in ("reverse", "both"):
                        specs.append((start + width, -1, start + step,
                                      start, "reverse"))
                for a_idx, b_idx, m_idx, n_idx, direction in specs:
                    records.append(
                        _new_row(config, array, a_idx, b_idx, m_idx, n_idx,
                                 step, n, direction)
                    )

    if not records:
        return pd.DataFrame(columns=COLUMNS)
    frame = pd.DataFrame(records)
    frame.insert(0, "measurement_id", np.arange(1, len(frame) + 1, dtype=int))
    return frame[COLUMNS]


def generate_all(config: SurveyConfig) -> dict[str, pd.DataFrame]:
    """One deterministic DataFrame per array, suitable for comparisons."""
    return {array: generate_survey(config, array) for array in ARRAYS}


def apply_homogeneous_example(
    frame: pd.DataFrame, resistivity_ohm_m: float = 100.0,
    current_a: float = 0.1,
) -> pd.DataFrame:
    """Generate ideal *analytical* readings for a uniform half-space.

    Not a layered-earth simulator, full ERT forward modelling or inversion.
    The sign of V_M-V_N agrees with signed geometric factor K.
    """
    if not (isfinite(resistivity_ohm_m) and resistivity_ohm_m > 0):
        raise ValueError("resistivity_ohm_m must be positive")
    if not (isfinite(current_a) and current_a > 0):
        raise ValueError("current_a must be positive")
    out = frame.copy()
    out["current_a"] = current_a
    out["voltage_v"] = resistivity_ohm_m * current_a / out["k_signed_m"]
    out["apparent_resistivity_ohm_m"] = resistivity_ohm_m
    return out


def attach_measurements(template: pd.DataFrame, measurements: pd.DataFrame) -> pd.DataFrame:
    """Attach user readings only if the CSV matches the exported geometry.

    Template contains IDs and electrode indices for *all* survey rows. A row
    may contain either apparent_resistivity_ohm_m or BOTH voltage_v,current_a.
    Blank rows are kept blank and omitted from apparent-resistivity charts.
    Existing values are not silently preferred if both choices are provided.
    """
    keys = ["measurement_id", "A_idx", "B_idx", "M_idx", "N_idx"]
    missing = [c for c in keys if c not in measurements]
    if missing:
        raise ValueError(f"Missing CSV columns: {', '.join(missing)}")
    if len(measurements) != len(template):
        raise ValueError("CSV must contain every exported row (unused readings can be blank)")
    actual = measurements.copy()
    for key in keys:
        numeric = pd.to_numeric(actual[key], errors="coerce")
        if numeric.isna().any() or not np.isfinite(numeric).all() or (numeric % 1 != 0).any():
            raise ValueError(f"CSV column {key} must contain valid integer identifiers")
        actual[key] = numeric.astype(int)
    if actual["measurement_id"].duplicated().any():
        raise ValueError("CSV contains repeated measurement_id values")
    actual = actual.sort_values("measurement_id").reset_index(drop=True)
    expected = template.sort_values("measurement_id").reset_index(drop=True)
    if not actual[keys].equals(expected[keys]):
        raise ValueError("CSV electrode IDs do not match the current array configuration")

    def get_measurement_col(name: str) -> pd.Series:
        if name not in actual:
            return pd.Series(np.nan, index=actual.index, dtype=float)
        values = actual[name]
        parsed = pd.to_numeric(values, errors="coerce")
        # Distinguish genuinely blank cells from typo or nonnumerical text.
        not_blank = values.notna() & (values.astype(str).str.strip() != "")
        if (not_blank & parsed.isna()).any() or ((parsed.notna()) & ~np.isfinite(parsed)).any():
            raise ValueError(f"Invalid numeric entry in {name}")
        return parsed.astype(float)

    rho = get_measurement_col("apparent_resistivity_ohm_m")
    volts = get_measurement_col("voltage_v")
    amps = get_measurement_col("current_a")
    rho_set = rho.notna()
    volts_set = volts.notna()
    amps_set = amps.notna()
    if (rho_set & (volts_set | amps_set)).any():
        raise ValueError("For each row, supply either rho_a OR voltage/current, not both")
    if (volts_set != amps_set).any():
        raise ValueError("Each measured row needs both voltage_v and current_a")
    if (rho_set & (rho <= 0)).any() or (amps_set & (amps <= 0)).any():
        raise ValueError("Apparent resistivity and injected current must be positive")
    inferred = expected["k_signed_m"] * volts / amps
    if (volts_set & ((inferred <= 0) | ~np.isfinite(inferred))).any():
        raise ValueError(
            "Calculated rho_a is nonpositive. Check voltage sign: "
            "voltage_v must be V_M − V_N using the exported electrode ordering."
        )
    resolved = rho.copy()
    resolved.loc[volts_set] = inferred.loc[volts_set]
    expected["apparent_resistivity_ohm_m"] = resolved
    expected["current_a"] = amps
    expected["voltage_v"] = volts
    return expected[COLUMNS]
