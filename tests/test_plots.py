"""Smoke tests for each visualization and the example CSV export/import loop."""
from io import StringIO

import pandas as pd

from ert_pseudosection.geometry import (
    SurveyConfig, apply_homogeneous_example, attach_measurements, generate_all,
)
from ert_pseudosection.plots import (
    comparison_chart, comparison_png, coverage_chart, density_chart,
    electrode_diagram, png_bytes,
)


def test_interactive_charts_and_png_exports() -> None:
    schemes = generate_all(SurveyConfig(n_electrodes=12, spacing_m=5, max_a_multiplier=2, max_n=3))
    for array, frame in schemes.items():
        assert "plotly" in coverage_chart(frame, array).to_html()
        assert "plotly" in density_chart(frame, 5).to_html()
        assert "plotly" in electrode_diagram(frame.iloc[0]).to_html()
        assert png_bytes(frame, array).startswith(b"\x89PNG")
        filled = apply_homogeneous_example(frame, 100, 0.1)
        assert png_bytes(filled, array, show_resistivity=True).startswith(b"\x89PNG")
        assert "plotly" in coverage_chart(filled, array, show_resistivity=True).to_html()
    assert len(comparison_chart(schemes).data) == 4
    assert comparison_png(schemes).startswith(b"\x89PNG")


def test_user_export_import_csv_round_trip() -> None:
    scheme = generate_all(SurveyConfig(n_electrodes=12))["dipole_dipole"]
    simulated = apply_homogeneous_example(scheme, 92, 0.1)
    text = simulated.drop(columns=["apparent_resistivity_ohm_m"]).to_csv(index=False)
    uploaded = pd.read_csv(StringIO(text))
    validated = attach_measurements(scheme, uploaded)
    assert validated["apparent_resistivity_ohm_m"].round(7).eq(92).all()
