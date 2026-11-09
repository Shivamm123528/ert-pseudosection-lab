"""ERT Pseudosection Lab: transparent survey-design and plotting utilities."""
from .geometry import (
    ARRAYS, ARRAY_LABELS, MEDIAN_DEPTH_FACTORS, SurveyConfig,
    apply_homogeneous_example, attach_measurements, generate_all,
    generate_survey, geometric_factor_m, median_pseudodepth_m,
)

__all__ = [
    "ARRAYS", "ARRAY_LABELS", "MEDIAN_DEPTH_FACTORS", "SurveyConfig",
    "apply_homogeneous_example", "attach_measurements", "generate_all",
    "generate_survey", "geometric_factor_m", "median_pseudodepth_m",
]
