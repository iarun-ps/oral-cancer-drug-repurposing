"""Tests for the drug repurposing integration pipeline."""

from pathlib import Path

import pandas as pd
import pytest

from pipeline.run_pipeline import run_pipeline

MOCK = Path("pipeline/mock")


def execute(tmp_path, gold_path=None, mock=True):
    return run_pipeline(
        MOCK / "oral_signature.csv",
        MOCK / "candidate_drugs.csv",
        MOCK / "binding_scores.csv",
        tmp_path / "shortlist.csv",
        gold_path=gold_path,
        mock=mock,
    )


def test_pipeline_runs(tmp_path):
    result = execute(tmp_path)

    assert len(result) == 3
    assert (tmp_path / "shortlist.csv").exists()


def test_rank_product(tmp_path):
    result = execute(tmp_path)

    assert result["combined_rank"].tolist() == [2, 2, 9]


def test_known_targets_excluded(tmp_path):
    result = execute(tmp_path)

    targets = dict(zip(result["drug_name"], result["best_target"]))

    assert targets["MockDrug_A"] == "TP53"
    assert targets["MockDrug_B"] == "AKT1"
    assert targets["MockDrug_C"] == "MTOR"


def test_mock_results_labelled(tmp_path):
    result = execute(tmp_path)

    assert result["evidence_note"].eq("SYNTHETIC TEST DATA").all()


def test_real_pipeline_requires_gold_list(tmp_path):
    with pytest.raises(ValueError, match="require --gold"):
        execute(tmp_path, mock=False)


def test_gold_list_annotation(tmp_path):
    gold = tmp_path / "known_drugs.csv"

    pd.DataFrame({
        "drug_name": ["MockDrug_A"]
    }).to_csv(gold, index=False)

    result = execute(tmp_path, gold_path=gold, mock=False)

    flags = dict(zip(
        result["drug_name"],
        result["known_oral_cancer_drug"]
    ))

    assert flags["MockDrug_A"] == "yes"
    assert flags["MockDrug_B"] == "no"


def test_missing_input_rejected(tmp_path):
    with pytest.raises(ValueError, match="File not found"):
        run_pipeline(
            tmp_path / "missing.csv",
            MOCK / "candidate_drugs.csv",
            MOCK / "binding_scores.csv",
            tmp_path / "shortlist.csv",
            mock=True,
        )


def test_deterministic_output(tmp_path):
    first = execute(tmp_path / "run1")
    second = execute(tmp_path / "run2")

    pd.testing.assert_frame_equal(first, second)
