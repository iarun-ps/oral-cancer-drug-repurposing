"""Analyze how binding-score weights affect drug rankings."""

import numpy as np
import pandas as pd

WEIGHTS = (0.0, 0.25, 0.50, 0.75, 1.0)


def sensitivity_analysis(shortlist, top_k=10):
    required = {
        "drug_name",
        "reversal_rank",
        "binding_rank",
        "combined_rank",
    }

    if not required.issubset(shortlist.columns):
        raise ValueError("Missing required ranking columns")

    if shortlist.empty:
        raise ValueError("Shortlist is empty")

    df = shortlist.copy()

    if df["drug_name"].isna().any() or not df["drug_name"].is_unique:
        raise ValueError("Drug names must be present and unique")

    for col in ("reversal_rank", "binding_rank", "combined_rank"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

        if not np.isfinite(df[col].to_numpy(dtype=float)).all():
            raise ValueError(f"Invalid values in {col}")

        if (df[col] <= 0).any():
            raise ValueError(f"{col} must be positive")

    expected = df["reversal_rank"] * df["binding_rank"]

    if not np.allclose(df["combined_rank"], expected):
        raise ValueError("Combined rank must equal rank product")

    baseline = df.sort_values(
        ["combined_rank", "reversal_rank", "binding_rank", "drug_name"]
    ).reset_index(drop=True)

    baseline_order = {
        name: rank
        for rank, name in enumerate(baseline["drug_name"], start=1)
    }

    k = min(top_k, len(df))
    baseline_top = set(baseline["drug_name"].head(k))

    summaries = []
    rankings = []

    for weight in WEIGHTS:
        current = df.copy()

        current["weighted_score"] = (
            (1 - weight) * current["reversal_rank"]
            + weight * current["binding_rank"]
        )

        current = current.sort_values(
            ["weighted_score", "reversal_rank", "binding_rank", "drug_name"]
        ).reset_index(drop=True)

        current["weighted_rank"] = np.arange(1, len(current) + 1)
        current["binding_weight"] = weight

        current["baseline_rank"] = current["drug_name"].map(
            baseline_order
        )

        correlation = (
            float(np.corrcoef(
                current["baseline_rank"],
                current["weighted_rank"]
            )[0, 1])
            if len(current) > 1 else 1.0
        )

        overlap = len(
            baseline_top.intersection(current["drug_name"].head(k))
        )

        summaries.append({
            "binding_weight": weight,
            "top_k": k,
            "top_k_overlap": overlap,
            "top_k_overlap_fraction": overlap / k,
            "spearman_correlation": correlation,
        })

        rankings.append(current[
            [
                "binding_weight",
                "drug_name",
                "weighted_score",
                "weighted_rank",
                "baseline_rank",
            ]
        ])

    return (
        pd.DataFrame(summaries),
        pd.concat(rankings, ignore_index=True),
    )
