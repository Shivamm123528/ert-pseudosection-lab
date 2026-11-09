"""Plotting helpers: Plotly for Streamlit, Matplotlib for PNG exports."""
from __future__ import annotations

from io import BytesIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .geometry import ARRAY_LABELS


def coverage_chart(frame: pd.DataFrame, array: str, *, show_resistivity: bool = False) -> go.Figure:
    """Interactive *point-placement* pseudosection, not subsurface imaging."""
    has_rho = (
        show_resistivity
        and "apparent_resistivity_ohm_m" in frame
        and frame["apparent_resistivity_ohm_m"].notna().any()
    )
    displayed = frame.dropna(subset=["apparent_resistivity_ohm_m"]) if has_rho else frame
    colors = (
        np.log10(displayed["apparent_resistivity_ohm_m"].astype(float))
        if has_rho else displayed["a_multiplier"].astype(float)
    )
    ids = np.stack([
        displayed["measurement_id"].astype(int),
        displayed["A_idx"].astype(int), displayed["B_idx"].astype(int),
        displayed["M_idx"].astype(int), displayed["N_idx"].astype(int),
        displayed["n"].astype(int), displayed["a_m"].astype(float),
        displayed["k_signed_m"].astype(float),
    ], axis=1) if len(displayed) else np.empty((0, 8))
    plot_title = ("Apparent-resistivity pseudosection" if has_rho
                  else "Geometry-only data coverage")
    hover = (
        "Measurement %{customdata[0]:.0f}<br>"
        "A / B / M / N: %{customdata[1]:.0f} / %{customdata[2]:.0f} / "
        "%{customdata[3]:.0f} / %{customdata[4]:.0f}<br>"
        "n: %{customdata[5]:.0f}; a: %{customdata[6]:.2f} m<br>"
        "K (signed): %{customdata[7]:.2f} m<br>"
        "x: %{x:.2f} m; pseudodepth: %{y:.2f} m"
    )
    if has_rho:
        hover += "<br>ρₐ: %{text} Ω·m"
    hover += "<extra></extra>"
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=displayed["plot_x_m"], y=displayed["pseudodepth_m"],
        mode="markers", customdata=ids,
        text=(displayed["apparent_resistivity_ohm_m"].map(lambda x: f"{x:.3g}")
              if has_rho else None),
        marker=dict(
            size=10, opacity=0.88, color=colors, colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="log₁₀(ρₐ)" if has_rho else "a multiplier"),
            line=dict(width=0.35, color="rgba(30,35,55,0.6)"),
        ),
        hovertemplate=hover,
    ))
    fig.update_layout(
        title=f"{ARRAY_LABELS[array]} — {plot_title}",
        xaxis_title="Horizontal position (m)",
        yaxis_title="Median pseudodepth (m; positive downward)",
        height=490, margin=dict(l=30, r=20, t=60, b=35),
        template="plotly_white",
    )
    fig.update_yaxes(autorange="reversed", rangemode="tozero")
    return fig


def density_chart(frame: pd.DataFrame, spacing_m: float) -> go.Figure:
    """2-D histogram of assigned points, NOT physical sensitivity or resolution."""
    if frame.empty:
        return go.Figure()
    dx = max(spacing_m / 2, 0.05)
    dz = max(spacing_m / 4, 0.025)
    xs = frame["plot_x_m"].to_numpy(dtype=float)
    zs = frame["pseudodepth_m"].to_numpy(dtype=float)
    x_edges = np.arange(0, xs.max() + 2 * dx, dx)
    z_edges = np.arange(0, zs.max() + 2 * dz, dz)
    counts, xb, zb = np.histogram2d(xs, zs, bins=[x_edges, z_edges])
    centers_x = 0.5 * (xb[:-1] + xb[1:])
    centers_z = 0.5 * (zb[:-1] + zb[1:])
    fig = go.Figure(go.Heatmap(
        x=centers_x, y=centers_z, z=counts.T, colorscale="YlGnBu",
        colorbar=dict(title="Point count"), hovertemplate=(
            "x: %{x:.1f} m<br>pseudodepth: %{y:.2f} m"
            "<br>Points in bin: %{z:.0f}<extra></extra>"
        ),
    ))
    fig.update_layout(
        title="Assigned measurement-point density (geometry only)",
        xaxis_title="Horizontal position (m)",
        yaxis_title="Median pseudodepth (m)",
        height=450, template="plotly_white", margin=dict(l=30, r=15, t=60, b=35),
    )
    fig.update_yaxes(autorange="reversed")
    return fig


def comparison_chart(schemes: dict[str, pd.DataFrame]) -> go.Figure:
    """Four comparable plots using the same horizontal and depth axes."""
    orders = list(ARRAY_LABELS)
    fig = make_subplots(rows=2, cols=2, subplot_titles=[ARRAY_LABELS[a] for a in orders],
                        horizontal_spacing=0.1, vertical_spacing=0.16)
    all_x, all_z = [], []
    for number, array in enumerate(orders):
        row, col = divmod(number, 2)
        frame = schemes[array]
        fig.add_trace(go.Scattergl(
            x=frame["plot_x_m"], y=frame["pseudodepth_m"],
            mode="markers", marker=dict(size=6, color="#2668ae", opacity=0.7),
            name=ARRAY_LABELS[array], showlegend=False,
            hovertemplate="x: %{x:.1f} m<br>pseudodepth: %{y:.2f} m<extra></extra>",
        ), row=row + 1, col=col + 1)
        if len(frame):
            all_x.extend(frame["plot_x_m"].tolist())
            all_z.extend(frame["pseudodepth_m"].tolist())
    if all_x and all_z:
        bound_x = [0, max(all_x) * 1.04]
        bound_z = [max(all_z) * 1.07, 0]
        fig.update_xaxes(range=bound_x)
        fig.update_yaxes(range=bound_z)
    fig.update_layout(
        title="Array comparison — positions in median-depth convention",
        height=740, template="plotly_white", margin=dict(l=35, r=20, t=90, b=50),
    )
    for row in (1, 2):
        for col in (1, 2):
            fig.update_xaxes(title_text="Horizontal position (m)", row=row, col=col)
            fig.update_yaxes(title_text="Pseudodepth (m)", row=row, col=col)
    return fig


def electrode_diagram(record: pd.Series) -> go.Figure:
    """One selected electrode configuration and its assigned plotting point."""
    active = []
    for label in ("A", "B", "M", "N"):
        index = int(record[f"{label}_idx"])
        if index >= 0:
            active.append((label, float(record[f"x_{label}_m"]), index))
    active.sort(key=lambda pair: pair[1])
    fig = go.Figure()
    left, right = active[0][1], active[-1][1]
    margin = max((right-left) * 0.16, float(record["a_m"]) / 2)
    fig.add_hline(y=0, line_color="#64748b", line_width=2)
    for label, x, index in active:
        color = "#ef8a3d" if label in ("A", "B") else "#227db5"
        fig.add_trace(go.Scatter(
            x=[x], y=[0], mode="markers+text", text=[f"{label} (E{index+1})"],
            textposition="top center", marker=dict(size=14, color=color),
            showlegend=False, hovertemplate=f"{label}, electrode {index+1}, x={x:g} m<extra></extra>",
        ))
    fig.add_trace(go.Scatter(
        x=[float(record["plot_x_m"])], y=[-float(record["pseudodepth_m"])],
        mode="markers+text", text=["Assigned plotting point"],
        textposition="bottom center", marker=dict(symbol="diamond", color="#22a879", size=12),
        showlegend=False,
        hovertemplate="x: %{x:.2f} m<br>z: %{y:.2f} m (negative plotting axis)<extra></extra>",
    ))
    if int(record["B_idx"]) == -1:
        fig.add_annotation(x=right, y=0.32, showarrow=False, text="B: ideal remote electrode (∞)")
    fig.update_layout(
        height=275, template="plotly_white", showlegend=False,
        xaxis=dict(title="Electrode position (m)", range=[left-margin, right+margin]),
        yaxis=dict(title="Diagram z (m)", range=[-float(record["pseudodepth_m"])*1.65-0.2, 1.4]),
        margin=dict(l=30, r=20, t=20, b=35),
    )
    return fig


def png_bytes(frame: pd.DataFrame, title: str, show_resistivity: bool = False) -> bytes:
    """Static PNG export without browser dependencies like Kaleido."""
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    rho = show_resistivity and frame["apparent_resistivity_ohm_m"].notna().any()
    points = frame.dropna(subset=["apparent_resistivity_ohm_m"]) if rho else frame
    color = (np.log10(points["apparent_resistivity_ohm_m"].to_numpy(dtype=float))
             if rho else points["a_multiplier"].to_numpy(dtype=float))
    plot = ax.scatter(points["plot_x_m"], points["pseudodepth_m"],
                      c=color, s=30, cmap="viridis", edgecolors="black", linewidths=0.2)
    fig.colorbar(plot, ax=ax, label="log10(ρa / Ω·m)" if rho else "a multiplier")
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("Horizontal position (m)")
    ax.set_ylabel("Median pseudodepth (m)")
    ax.grid(alpha=0.20)
    buf = BytesIO()
    fig.savefig(buf, dpi=180, format="png", facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def comparison_png(schemes: dict[str, pd.DataFrame]) -> bytes:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (array, frame) in zip(axes.flat, schemes.items()):
        ax.scatter(frame["plot_x_m"], frame["pseudodepth_m"],
                   color="#2668ae", s=15, alpha=0.8)
        ax.set_title(f"{ARRAY_LABELS[array]} ({len(frame)} points)")
        ax.set_xlabel("Horizontal position (m)")
        ax.set_ylabel("Median pseudodepth (m)")
        ax.grid(alpha=0.18)
    axes[0, 0].invert_yaxis()
    fig.suptitle("ERT pseudosection coverage — geometry only", fontsize=14)
    buf = BytesIO()
    fig.savefig(buf, dpi=180, format="png", facecolor="white")
    plt.close(fig)
    return buf.getvalue()
