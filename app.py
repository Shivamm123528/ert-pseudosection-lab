"""Streamlit dashboard for geometry-first ERT pseudosection survey design."""
from __future__ import annotations

from io import BytesIO

import pandas as pd
import streamlit as st

from ert_pseudosection.geometry import (
    ARRAYS, ARRAY_LABELS, SurveyConfig, apply_homogeneous_example,
    attach_measurements, generate_all,
)
from ert_pseudosection.plots import (
    comparison_chart, comparison_png, coverage_chart, density_chart,
    electrode_diagram, png_bytes,
)

st.set_page_config(
    page_title="ERT Pseudosection Lab",
    page_icon="📡",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("ERT Pseudosection Lab")
st.caption("Electrical resistivity survey design · Four electrode arrays · Reproducible Python calculations")
st.info(
    "**Interpretation boundary:** The plotted depth is the published **homogeneous-earth "
    "median-depth convention** (Loke, 1999, after Edwards, 1977). It is **not** a "
    "measured geological depth, a sensitivity map, or a resistivity inversion result.",
    icon="ℹ️",
)

with st.sidebar:
    st.header("Survey controls")
    n_electrodes = st.number_input("Number of electrodes", min_value=4, max_value=256,
                                    value=32, step=1)
    spacing = st.number_input("Electrode interval (m)", min_value=0.1,
                              max_value=1000.0, value=5.0, step=0.5, format="%.2f")
    a_multiplier = st.slider("Maximum a multiplier", min_value=1, max_value=8, value=2,
                              help="a = multiplier × electrode interval. All multipliers from 1 to this value are used.")
    n_factor = st.slider("Maximum separation factor n", min_value=1, max_value=6, value=6,
                         help="Applied to Dipole–Dipole, Wenner–Schlumberger and Pole–Dipole. Published depth factors stop at n=6.")
    pole_direction = st.selectbox("Pole–Dipole direction", options=["forward", "reverse", "both"],
                                  help="Forward: A–M–N. Reverse: N–M–A. B is an ideal remote electrode.")
    st.caption("Horizontal line, flat topography, equal electrode intervals.")
    st.divider()
    st.markdown("**Electrodes:** A and B inject/return current; M and N measure potential.")
    st.markdown("**Array n:** a relative gap. **a:** physical spacing in metres.")

config = SurveyConfig(
    n_electrodes=int(n_electrodes), spacing_m=float(spacing),
    max_a_multiplier=int(a_multiplier), max_n=int(n_factor),
    pole_direction=pole_direction,
)

schemes = generate_all(config)
page1, page2, page3 = st.tabs(["Survey designer", "Compare arrays", "Theory, sources & limitations"])

with page1:
    array = st.selectbox("Choose electrode arrangement", options=list(ARRAYS),
                         format_func=lambda key: ARRAY_LABELS[key])
    base = schemes[array]
    if base.empty:
        st.warning("This electrode count cannot fit the selected geometry. Increase the number of electrodes.")
    else:
        cols = st.columns(4)
        cols[0].metric("Survey length", f"{config.length_m:,.1f} m")
        cols[1].metric("Configurations", f"{len(base):,}")
        cols[2].metric("Deepest median pseudodepth", f"{base['pseudodepth_m'].max():.2f} m")
        cols[3].metric("Maximum |K|", f"{base['k_magnitude_m'].max():,.1f} m")

        st.subheader("Choose data mode")
        mode = st.radio("Measurement input", options=[
            "Geometry only", "Ideal homogeneous-earth example", "Upload readings CSV"
        ], horizontal=True, label_visibility="collapsed")
        frame = base
        if mode == "Ideal homogeneous-earth example":
            c1, c2 = st.columns(2)
            with c1:
                rho = st.number_input("Uniform half-space resistivity (Ω·m)",
                                      min_value=0.1, value=100.0, step=10.0)
            with c2:
                current = st.number_input("Injected current (A)",
                                          min_value=0.00001, value=0.1, format="%.4f")
            frame = apply_homogeneous_example(base, float(rho), float(current))
            st.caption("Analytical homogeneous half-space: V = ρI/K. "
                       "This is not a layered-earth forward solver or field data.")
        elif mode == "Upload readings CSV":
            st.markdown(
                "Download the **blank measurement template**, then populate either "
                "`apparent_resistivity_ohm_m`, or **both** `current_a` and "
                "`voltage_v` (where `voltage_v = V_M − V_N`). "
                "Keep all IDs and electrode columns unchanged; blank rows are allowed."
            )
            st.download_button("Download blank CSV template", base.to_csv(index=False).encode("utf-8"),
                               file_name=f"{array}_blank_template.csv", mime="text/csv")
            uploaded = st.file_uploader("Upload the completed template", type=["csv"])
            if uploaded is not None:
                try:
                    frame = attach_measurements(base, pd.read_csv(uploaded))
                    count = int(frame["apparent_resistivity_ohm_m"].notna().sum())
                    st.success(f"Validated {count:,} populated resistivity readings against {len(base):,} geometry rows.")
                except (ValueError, pd.errors.ParserError, UnicodeDecodeError) as exc:
                    st.error(f"CSV validation failed: {exc}")
                    frame = base
        st.subheader("Measurement pseudosection")
        colored = frame["apparent_resistivity_ohm_m"].notna().any()
        plot = coverage_chart(frame, array, show_resistivity=bool(colored))
        st.plotly_chart(plot, use_container_width=True)
        if not colored:
            st.caption("Points show measurement positions only; they do not represent inferred geological layers.")

        d1, d2, d3 = st.columns(3)
        d1.download_button(
            "Download full survey CSV", frame.to_csv(index=False).encode("utf-8"),
            file_name=f"{array}_measurements.csv", mime="text/csv")
        d2.download_button(
            "Download graph PNG",
            png_bytes(frame, f"{ARRAY_LABELS[array]} — ERT pseudosection", show_resistivity=bool(colored)),
            file_name=f"{array}_pseudosection.png", mime="image/png")
        d3.download_button(
            "Download interactive HTML",
            plot.to_html(full_html=True, include_plotlyjs=True).encode("utf-8"),
            file_name=f"{array}_pseudosection.html", mime="text/html")

        with st.expander("Point density map — counts, not sensitivity", expanded=False):
            st.plotly_chart(density_chart(base, config.spacing_m), use_container_width=True)
            st.caption("Each cell counts assigned pseudosection points. The map does not quantify "
                       "subsurface measurement sensitivity or resolution.")

        st.subheader("Inspect one electrode configuration")
        selected_id = st.slider("Measurement ID", min_value=1, max_value=len(base), value=1)
        selected = frame.iloc[selected_id - 1]
        st.plotly_chart(electrode_diagram(selected), use_container_width=True)
        st.write(
            f"**Spacing a:** {selected['a_m']:g} m · **n:** {int(selected['n'])} · "
            f"**K signed:** {selected['k_signed_m']:.4f} m · "
            f"**Plotting point:** ({selected['plot_x_m']:.2f}, {selected['pseudodepth_m']:.2f}) m"
        )
        if int(selected["B_idx"]) == -1:
            st.caption("B = −1 in CSV means an ideal remote return electrode; "
                       "the remote B is not part of the surface-line span.")
        with st.expander("View the complete electrode table"):
            st.dataframe(frame, hide_index=True, use_container_width=True)

with page2:
    st.subheader("Compare the four arrays with identical survey settings")
    st.plotly_chart(comparison_chart(schemes), use_container_width=True)
    summary = pd.DataFrame([
        {
            "Array": ARRAY_LABELS[name],
            "Measurements": len(data),
            "Maximum median pseudodepth (m)": round(float(data["pseudodepth_m"].max()), 3)
                if len(data) else 0,
            "Maximum active span (m)": round(float(data["array_span_m"].max()), 3)
                if len(data) else 0,
            "Maximum |K| (m)": round(float(data["k_magnitude_m"].max()), 3)
                if len(data) else 0,
        } for name, data in schemes.items()
    ])
    st.dataframe(summary, hide_index=True, use_container_width=True)
    st.download_button("Download comparison PNG", comparison_png(schemes),
                       file_name="ert_four_array_comparison.png", mime="image/png")
    st.download_button("Download summary CSV", summary.to_csv(index=False).encode("utf-8"),
                       file_name="ert_array_summary.csv", mime="text/csv")
    st.caption("A deeper median pseudodepth or more assigned data points does not by itself "
               "mean better reconstruction: uncertainty, noise, sensitivity and inversion also matter.")

with page3:
    st.subheader("Electrode geometry")
    st.markdown(r"""
    | Array | Surface electrode order | Physical lengths |
    |---|---|---|
    | Wenner-alpha | A – M – N – B | AM = MN = NB = a |
    | Wenner–Schlumberger | A – M – N – B | AM = NB = na; MN = a |
    | Dipole–Dipole | A – B – M – N | AB = MN = a; BM = na |
    | Pole–Dipole | A – M – N (or reversed) | AM = na; MN = a; B is ideally remote |

    **a = multiplier × fixed electrode interval.** For Wenner, increasing the multiplier
    increases all three equal distances. For other arrays, independently increasing the
    separation factor **n** moves the current and potential electrodes farther apart.
    """)
    st.subheader("Theory")
    st.latex(r"\rho_a = K \frac{\Delta V}{I}\,,\qquad \Delta V = V_M - V_N")
    st.latex(r"K = \frac{2\pi}{\frac{1}{AM}-\frac{1}{AN}-\frac{1}{BM}+\frac{1}{BN}}")
    st.markdown(
        "The program retains the **signed** geometric factor. For a Dipole–Dipole "
        "configuration labelled A–B–M–N, both K and the ideal voltage V_M−V_N "
        "are negative, giving a **positive** apparent resistivity. "
        "For Pole–Dipole, terms containing the remote B are zero (B → ∞)."
    )
    st.latex(r"z_e = c_{\mathrm{array},n}\,a")
    coefficients = pd.DataFrame({
        "n": [1, 2, 3, 4, 5, 6],
        "Wenner-alpha": [0.519, None, None, None, None, None],
        "Wenner–Schlumberger": [0.520, 0.930, 1.320, 1.710, 2.090, 2.480],
        "Dipole–Dipole": [0.416, 0.697, 0.962, 1.220, 1.476, 1.730],
        "Pole–Dipole": [0.520, 0.930, 1.320, 1.710, 2.090, 2.480],
    })
    st.dataframe(coefficients, hide_index=True, use_container_width=True)
    st.caption("Coefficients rounded as in M. H. Loke (1999), Table 2, "
               "citing Edwards (1977). Wenner-alpha uses 0.519 for every a multiplier.")
    st.subheader("Limitations")
    st.markdown(
        "- Geometric plotting locations are **not** actual measurement locations in rock. "
        "A datum responds to a broad subsurface volume.\n"
        "- Median depths assume a **flat, homogeneous, isotropic half-space**. "
        "Significant resistivity contrasts and topography change practical DOI.\n"
        "- The pole return electrode is assumed infinitely far away. "
        "A real finite-distance return needs its true position and an adjusted K.\n"
        "- Point-count density is **not** an electrical sensitivity calculation.\n"
        "- Homogeneous readings are **analytical teaching examples**. The tool "
        "does not perform 2-D forward modelling or resistivity inversion.\n"
        "- Data coverage describes geometric sampling; it is not a guarantee of resolution."
    )
    st.subheader("References")
    st.markdown(
        "1. M. H. Loke (1999), *Electrical imaging surveys for environmental and "
        "engineering studies*, Table 2 and Section 2.3 — "
        "[read original tutorial](https://www.geometrics.com/wp-content/uploads/2018/10/Lokenote.pdf).\n"
        "2. L. S. Edwards (1977), *A modified pseudosection for resistivity and IP*, "
        "*Geophysics*, DOI: [10.1190/1.1440762](https://doi.org/10.1190/1.1440762).\n"
        "3. US EPA, [Electrical Resistivity](https://www.epa.gov/environmental-geophysics/electrical-resistivity).\n"
        "4. pyGIMLi, [ERT array-scheme implementation](https://www.pygimli.org/_modules/pygimli/physics/ert/ertScheme/)."
    )

st.divider()
st.caption("ERT Pseudosection Lab · Reproducible geometry and theoretical median-depth plotting · No inversion")
