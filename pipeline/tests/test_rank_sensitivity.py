"""Tests for drug ranking sensitivity analysis."""

import pandas as pd
import pytest

from pipeline.rank_sensitivity import sensitivity_analysis


@pytest.fixture
def shortlist():
    return pd.DataFrame({
        "drug_name": ["Drug_A", "Drug_B", "Drug_C"],
        "reversal_rank": [1, 2, 3],
        "binding_rank": [2, 1, 3],
        "combined_rank": [2, 2, 9],
    })


def test_five_binding_weights(shortlist):
    summary, _ = sensitivity_analysis(shortlist)
    assert summary["binding_weight"].tolist() == [
        0.0, 0.25, 0.5, 0.75, 1.0
    ]


def test_total_ranking_records(shortlist):
    _, rankings = sensitivity_analysis(shortlist)
    assert len(rankings) == 15


def test_reversal_only_ranking(shortlist):
    _, rankings = sensitivity_analysis(shortlist)
    result = rankings[rankings["binding_weight"] == 0.0]
    assert result["drug_name"].tolist() == [
        "Drug_A", "Drug_B", "Drug_C"
    ]


def test_binding_only_ranking(shortlist):
    _, rankings = sensitivity_analysis(shortlist)
    result = rankings[rankings["binding_weight"] == 1.0]
    assert result["drug_name"].tolist() == [
        "Drug_B", "Drug_A", "Drug_C"
    ]


def test_top_one_overlap(shortlist):
    summary, _ = sensitivity_analysis(shortlist, top_k=1)
    assert summary["top_k_overlap"].iloc[0] == 1
    assert summary["top_k_overlap"].iloc[-1] == 0


def test_invalid_rank_product(shortlist):
    shortlist.loc[0, "combined_rank"] = 99
    with pytest.raises(ValueError, match="rank product"):
        sensitivity_analysis(shortlist)


def test_duplicate_drugs(shortlist):
    shortlist.loc[1, "drug_name"] = "Drug_A"
    with pytest.raises(ValueError, match="unique"):
        sensitivity_analysis(shortlist)


def test_empty_shortlist():
    with pytest.raises(ValueError, match="empty"):
        sensitivity_analysis(pd.DataFrame(columns=[
            "drug_name", "reversal_rank",
            "binding_rank", "combined_rank"
        ]))


def test_missing_columns(shortlist):
    shortlist = shortlist.drop(columns=["binding_rank"])
    with pytest.raises(ValueError, match="Missing"):
        sensitivity_analysis(shortlist)


def test_negative_rank(shortlist):
    shortlist.loc[0, "binding_rank"] = -1
    with pytest.raises(ValueError, match="positive"):
        sensitivity_analysis(shortlist)
