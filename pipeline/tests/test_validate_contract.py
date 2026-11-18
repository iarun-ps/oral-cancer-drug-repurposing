"""Automated tests for the pipeline data contract."""

from pathlib import Path

import pandas as pd
import pytest

from pipeline.validate_contract import (
    SIGNATURE_COLUMNS,
    load_csv,
    validate_binding,
    validate_drugs,
    validate_signature,
)

MOCK = Path("pipeline/mock")
REAL = Path("signature/oral_cancer/results/oral_signature.csv")


@pytest.fixture
def signature():
    return pd.read_csv(MOCK / "oral_signature.csv", dtype=str)


@pytest.fixture
def drugs():
    return pd.read_csv(MOCK / "candidate_drugs.csv", dtype=str)


@pytest.fixture
def binding():
    return pd.read_csv(MOCK / "binding_scores.csv", dtype=str)


def test_mock_signature(signature):
    validate_signature(signature)


def test_real_signature():
    df = load_csv(REAL, SIGNATURE_COLUMNS)
    validate_signature(df)


def test_mock_drugs(drugs):
    validate_drugs(drugs)


def test_mock_binding(binding, drugs):
    validate_binding(binding, drugs)


def test_missing_gene_symbol(signature):
    signature.loc[0, "gene_symbol"] = ""

    with pytest.raises(ValueError, match="missing gene_symbol"):
        validate_signature(signature)


def test_invalid_direction(signature):
    signature.loc[0, "direction"] = "invalid"

    with pytest.raises(ValueError):
        validate_signature(signature)


def test_invalid_padj(signature):
    signature.loc[0, "padj"] = "1.5"

    with pytest.raises(ValueError):
        validate_signature(signature)


def test_duplicate_gene(signature):
    signature.loc[1, "gene_symbol"] = signature.loc[0, "gene_symbol"]

    with pytest.raises(ValueError, match="duplicate gene"):
        validate_signature(signature)


def test_missing_smiles(drugs):
    drugs.loc[0, "smiles"] = ""

    with pytest.raises(ValueError, match="missing smiles"):
        validate_drugs(drugs)


def test_invalid_smiles(drugs):
    drugs.loc[0, "smiles"] = "NOT_A_SMILES"

    with pytest.raises(ValueError, match="invalid SMILES"):
        validate_drugs(drugs)


def test_reversal_score_direction(drugs):
    drugs.loc[0, "reversal_score"] = "0.1"

    with pytest.raises(ValueError, match="reversal_score"):
        validate_drugs(drugs)


def test_invalid_drug_rank(drugs):
    drugs.loc[0, "rank"] = "-1"

    with pytest.raises(ValueError):
        validate_drugs(drugs)


def test_training_data_leakage(binding, drugs):
    binding.loc[0, "in_training_data"] = "yes"

    with pytest.raises(ValueError, match="in_training_data"):
        validate_binding(binding, drugs)


def test_invalid_confidence(binding, drugs):
    binding.loc[0, "confidence"] = "1.5"

    with pytest.raises(ValueError):
        validate_binding(binding, drugs)


def test_binding_smiles_mismatch(binding, drugs):
    binding.loc[0, "smiles"] = "CCCC"

    with pytest.raises(ValueError, match="SMILES mismatch"):
        validate_binding(binding, drugs)


def test_unknown_candidate(binding, drugs):
    binding.loc[0, "drug_name"] = "UnknownDrug"

    with pytest.raises(ValueError):
        validate_binding(binding, drugs)


def test_known_target_mismatch(binding, drugs):
    binding.loc[0, "is_known_target"] = "no"

    with pytest.raises(ValueError, match="known-target mismatch"):
        validate_binding(binding, drugs)


def test_missing_columns(tmp_path):
    path = tmp_path / "invalid.csv"
    path.write_text("gene_symbol,rank\nMMP1,1\n")

    with pytest.raises(ValueError, match="incorrect columns"):
        load_csv(path, SIGNATURE_COLUMNS)


def test_invalid_numeric_value(signature):
    signature.loc[0, "log2fc"] = "invalid"

    with pytest.raises(ValueError):
        validate_signature(signature)


def test_real_binding_accepts_documented_extra_columns():
    from pipeline.validate_contract import BINDING_COLUMNS

    df = load_csv(
        Path("binding/results/binding_scores.csv"),
        BINDING_COLUMNS,
    )

    assert len(df) == 784
    assert "tree_sd" in df.columns
    assert "max_similarity_to_train" in df.columns


def test_categorical_confidence(binding, drugs):
    binding["confidence"] = "low"
    validate_binding(binding, drugs)


def test_random_forest_method(binding, drugs):
    binding["method"] = "RandomForest_Morgan2048_pChEMBL"
    validate_binding(binding, drugs)


def test_additional_known_target(binding, drugs):
    drug_name = binding.loc[0, "drug_name"]
    drugs.loc[drugs["drug_name"] == drug_name, "known_targets"] = ""
    binding.loc[0, "is_known_target"] = "yes"

    validate_binding(binding, drugs)


def test_invalid_similarity(binding, drugs):
    binding["tree_sd"] = "0.5"
    binding["max_similarity_to_train"] = "1.5"

    with pytest.raises(ValueError):
        validate_binding(binding, drugs)


def test_duplicate_drug_target_pair(binding, drugs):
    duplicate = binding.iloc[[0]].copy()

    binding = pd.concat(
        [binding, duplicate],
        ignore_index=True
    )

    with pytest.raises(ValueError, match="duplicate drug-target"):
        validate_binding(binding, drugs)
