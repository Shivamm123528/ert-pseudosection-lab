"""Regression tests for array positions, signed K, depth and CSV integrity."""
from math import pi

import numpy as np
import pandas as pd
import pytest

from ert_pseudosection.geometry import (
    ARRAYS, MEDIAN_DEPTH_FACTORS, SurveyConfig, apply_homogeneous_example,
    attach_measurements, generate_all, generate_survey, geometric_factor_m,
    median_pseudodepth_m,
)


def test_wenner_12_electrodes_has_18_positions() -> None:
    frame = generate_survey(SurveyConfig(12, 5, 3, 6), "wenner")
    assert len(frame) == 9 + 6 + 3
    first = frame.iloc[0]
    assert [first[k] for k in ("A_idx", "M_idx", "N_idx", "B_idx")] == [0, 1, 2, 3]
    assert first["plot_x_m"] == 7.5
    assert first["pseudodepth_m"] == pytest.approx(0.519 * 5)
    assert first["k_signed_m"] == pytest.approx(2 * pi * 5)
    assert set(frame["n"]) == {1}


@pytest.mark.parametrize("array, expected", [
    ("wenner", 2 * pi * 5),
    ("wenner_schlumberger", 2 * pi * 5),
    ("dipole_dipole", -6 * pi * 5),
    ("pole_dipole", 4 * pi * 5),
])
def test_known_n1_geometric_factors(array: str, expected: float) -> None:
    frame = generate_survey(SurveyConfig(24, 5, 1, 1), array)
    assert frame.iloc[0]["k_signed_m"] == pytest.approx(expected)


@pytest.mark.parametrize("array,n,expected", [
    ("wenner_schlumberger", 2, pi * 2 * 3 * 5),
    ("dipole_dipole", 2, -pi * 2 * 3 * 4 * 5),
    ("pole_dipole", 2, 2 * pi * 2 * 3 * 5),
])
def test_n2_geometric_factors(array: str, n: int, expected: float) -> None:
    frame = generate_survey(SurveyConfig(24, 5, 1, 2), array)
    row = frame.loc[(frame["n"] == n) & (frame["a_multiplier"] == 1)].iloc[0]
    assert row["k_signed_m"] == pytest.approx(expected)


def test_median_depth_table_and_rejected_extrapolation() -> None:
    assert median_pseudodepth_m("dipole_dipole", 3, 10) == pytest.approx(9.62)
    assert median_pseudodepth_m("wenner_schlumberger", 6, 5) == pytest.approx(12.4)
    assert MEDIAN_DEPTH_FACTORS["wenner"][1] == 0.519
    with pytest.raises(ValueError):
        median_pseudodepth_m("dipole_dipole", 7, 5)


def test_valid_indices_no_duplicates_no_off_line_x() -> None:
    config = SurveyConfig(24, 2.5, 3, 6, "both")
    for array, frame in generate_all(config).items():
        assert not frame.empty
        assert frame["plot_x_m"].between(0, config.length_m).all()
        assert (frame["pseudodepth_m"] > 0).all()
        assert frame["measurement_id"].is_unique
        for row in frame.itertuples():
            idxs = [row.A_idx, row.M_idx, row.N_idx]
            if row.B_idx >= 0:
                idxs.append(row.B_idx)
            assert len(set(idxs)) == len(idxs)
            assert all(0 <= i < 24 for i in idxs)
            if array == "pole_dipole":
                assert row.B_idx == -1
                assert np.isnan(row.x_B_m)
            else:
                assert row.B_idx >= 0


def test_forward_reverse_pole_dipole_are_mirror_images() -> None:
    fwd = generate_survey(SurveyConfig(12, 5, 1, 2, "forward"), "pole_dipole")
    rev = generate_survey(SurveyConfig(12, 5, 1, 2, "reverse"), "pole_dipole")
    assert len(fwd) == len(rev)
    assert np.allclose(fwd["k_signed_m"], rev["k_signed_m"])
    assert np.allclose(fwd["pseudodepth_m"], rev["pseudodepth_m"])
    assert np.allclose(fwd["plot_x_m"], rev["plot_x_m"])
    assert set(fwd["direction"]) == {"forward"}
    assert set(rev["direction"]) == {"reverse"}
    both = generate_survey(SurveyConfig(12, 5, 1, 2, "both"), "pole_dipole")
    assert len(both) == len(fwd) * 2


def test_homogeneous_example_signed_voltage_produces_positive_rho() -> None:
    for array, frame in generate_all(SurveyConfig(18, 5, 2, 4)).items():
        synthetic = apply_homogeneous_example(frame, 150, 0.08)
        assert (synthetic["apparent_resistivity_ohm_m"] == 150).all()
        inferred = synthetic["k_signed_m"] * synthetic["voltage_v"] / synthetic["current_a"]
        assert np.allclose(inferred, 150)
        if array == "dipole_dipole":
            assert (synthetic["voltage_v"] < 0).all()
            assert (synthetic["k_signed_m"] < 0).all()


def test_attach_uploaded_voltage_data() -> None:
    base = generate_survey(SurveyConfig(12, 5, 1, 3), "dipole_dipole")
    real = apply_homogeneous_example(base, 150, 0.08)
    upload = real.drop(columns=["apparent_resistivity_ohm_m"])
    result = attach_measurements(base, upload)
    assert np.allclose(result["apparent_resistivity_ohm_m"], 150)
    assert np.allclose(result["voltage_v"], real["voltage_v"])


def test_attach_uploaded_rho_and_blank_rows() -> None:
    base = generate_survey(SurveyConfig(12, 5, 1, 1), "wenner")
    upload = base.copy()
    upload.loc[0, "apparent_resistivity_ohm_m"] = 87
    result = attach_measurements(base, upload)
    assert result.loc[0, "apparent_resistivity_ohm_m"] == 87
    assert result["apparent_resistivity_ohm_m"].notna().sum() == 1


@pytest.mark.parametrize("mutate", [
    lambda df: df.assign(A_idx=df["A_idx"] + 1),
    lambda df: df.iloc[:-1],
    lambda df: df.assign(measurement_id=1),
])
def test_reject_misaligned_csv(mutate) -> None:
    base = generate_survey(SurveyConfig(12, 5, 1, 1), "wenner")
    with pytest.raises(ValueError):
        attach_measurements(base, mutate(base.copy()))


def test_reject_voltage_sign_and_both_input_types() -> None:
    base = generate_survey(SurveyConfig(12, 5, 1, 2), "dipole_dipole")
    filled = base.copy()
    filled.loc[0, "voltage_v"] = +0.05  # incorrectly signed for DD
    filled.loc[0, "current_a"] = 0.1
    with pytest.raises(ValueError, match="voltage sign"):
        attach_measurements(base, filled)
    filled.loc[0, "voltage_v"] = -0.05
    filled.loc[0, "apparent_resistivity_ohm_m"] = 100
    with pytest.raises(ValueError, match="either rho_a OR"):
        attach_measurements(base, filled)


def test_reject_invalid_config_and_degenerate_positions() -> None:
    with pytest.raises(ValueError):
        SurveyConfig(n_electrodes=3)
    with pytest.raises(ValueError):
        SurveyConfig(spacing_m=-1)
    with pytest.raises(ValueError):
        SurveyConfig(max_n=7)
    with pytest.raises(ValueError):
        geometric_factor_m(0, 10, 0, 5)
    with pytest.raises(ValueError):
        generate_survey(SurveyConfig(), "unknown")


def test_direct_signed_geometric_factor_remote() -> None:
    assert geometric_factor_m(0, float("inf"), 5, 10) == pytest.approx(20*pi)


def test_all_arrays_empty_possible_with_four_electrodes() -> None:
    # Some n / a combinations cannot physically fit, but there is at least
    # one valid array geometry for 4 electrodes and the code must not crash.
    schemes = generate_all(SurveyConfig(4, 5, 3, 6))
    assert all(len(data) >= 1 for data in schemes.values())
