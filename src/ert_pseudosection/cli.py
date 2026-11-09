"""Generate four ERT array files and images without launching Streamlit."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .geometry import SurveyConfig, apply_homogeneous_example, generate_all
from .plots import comparison_png, png_bytes


def main() -> None:
    parser = argparse.ArgumentParser(description="ERT pseudosection geometry generator")
    parser.add_argument("--electrodes", type=int, default=32)
    parser.add_argument("--spacing", type=float, default=5.0, help="Metres")
    parser.add_argument("--max-a", type=int, default=2, dest="max_a")
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument("--pole-direction", choices=["forward", "reverse", "both"], default="forward")
    parser.add_argument("--resistivity", type=float, default=None,
                        help="Optional homogeneous half-space rho in ohm-m")
    parser.add_argument("--current", type=float, default=0.1,
                        help="Example current in amperes (used with --resistivity)")
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    config = SurveyConfig(n_electrodes=args.electrodes, spacing_m=args.spacing,
                          max_a_multiplier=args.max_a, max_n=args.max_n,
                          pole_direction=args.pole_direction)
    args.output.mkdir(parents=True, exist_ok=True)
    schemes = generate_all(config)
    summary = []
    for array, geometry in schemes.items():
        frame = (apply_homogeneous_example(geometry, args.resistivity, args.current)
                 if args.resistivity is not None and len(geometry) else geometry)
        frame.to_csv(args.output / f"{array}.csv", index=False, float_format="%.10g")
        if len(frame):
            chart = png_bytes(frame, f"{array} — ERT pseudosection",
                              show_resistivity=args.resistivity is not None)
            (args.output / f"{array}.png").write_bytes(chart)
        summary.append({"array": array, "measurements": len(frame),
                        "max_median_pseudodepth_m": (
                            float(frame["pseudodepth_m"].max()) if len(frame) else 0.0)})
    pd.DataFrame(summary).to_csv(args.output / "summary.csv", index=False)
    (args.output / "comparison.png").write_bytes(comparison_png(schemes))
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"Saved project outputs to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
