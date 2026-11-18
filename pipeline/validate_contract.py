"""Validate CSV files used by the drug repurposing pipeline."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

SIGNATURE_COLUMNS = (
    "gene_symbol", "entrez_id", "log2fc", "padj",
    "direction", "rank", "n_cohorts"
)

DRUG_COLUMNS = (
    "rank", "drug_name", "pubchem_cid", "smiles",
    "reversal_score", "raw_score", "cell_line",
    "source_tool", "known_targets", "approval_status",
    "used_in_oral_cancer"
)

BINDING_COLUMNS = (
    "drug_name", "smiles", "target_gene", "pdb_id",
    "predicted_affinity", "confidence", "method",
    "is_known_target", "in_training_data"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_csv(path, expected):
    require(path.is_file(), f"File not found: {path}")

    df = pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig"
    )

    actual_columns = list(df.columns)
    accepted_columns = list(expected)

    if tuple(expected) == BINDING_COLUMNS:
        for optional in ("tree_sd", "max_similarity_to_train"):
            if optional in actual_columns:
                accepted_columns.append(optional)

    require(
        actual_columns == accepted_columns,
        f"{path}: incorrect columns. Expected {accepted_columns}; "
        f"found {actual_columns}"
    )

    require(not df.empty, f"{path}: CSV is empty")
    return df


def require_text(df, columns, label):
    for col in columns:
        values = df[col].str.strip()

        require(
            values.ne("").all(),
            f"{label}: missing {col} at rows "
            f"{(df.index[values.eq('')] + 2).tolist()}"
        )

        require(
            values.eq(df[col]).all(),
            f"{label}: extra whitespace in {col}"
        )


def number(
    df, column, label, *,
    integer=False, minimum=None, maximum=None
):
    values = pd.to_numeric(df[column], errors="coerce")

    require(
        values.notna().all()
        and np.isfinite(values.to_numpy(dtype=float)).all(),
        f"{label}: {column} must contain valid finite numbers"
    )

    if integer:
        require(
            (values % 1 == 0).all(),
            f"{label}: {column} must contain integers"
        )

    if minimum is not None:
        require(
            (values >= minimum).all(),
            f"{label}: {column} must be >= {minimum}"
        )

    if maximum is not None:
        require(
            (values <= maximum).all(),
            f"{label}: {column} must be <= {maximum}"
        )

    return values


def allowed(df, column, choices, label):
    require(
        df[column].isin(choices).all(),
        f"{label}: invalid {column}. Allowed: {sorted(choices)}"
    )


def validate_signature(df):
    label = "signature"

    require_text(df, ["gene_symbol", "direction"], label)

    require(
        df["gene_symbol"].is_unique,
        "signature: duplicate gene symbols"
    )

    number(df, "entrez_id", label, integer=True, minimum=1)
    log2fc = number(df, "log2fc", label)
    number(df, "padj", label, minimum=0, maximum=1)
    ranks = number(df, "rank", label, integer=True, minimum=1)
    number(df, "n_cohorts", label, integer=True, minimum=1, maximum=2)

    allowed(df, "direction", {"up", "down"}, label)

    correct_direction = (
        ((df["direction"] == "up") & (log2fc > 0))
        | ((df["direction"] == "down") & (log2fc < 0))
    )

    require(
        correct_direction.all(),
        "signature: direction disagrees with log2fc"
    )

    require(
        not df.assign(_rank=ranks)
        .duplicated(["direction", "_rank"])
        .any(),
        "signature: duplicate rank within direction"
    )


def validate_drugs(df):
    label = "candidates"

    require_text(
        df,
        [
            "drug_name", "smiles", "cell_line",
            "source_tool", "approval_status",
            "used_in_oral_cancer"
        ],
        label
    )

    require(
        df["drug_name"].is_unique,
        "candidates: duplicate drug names"
    )

    ranks = number(df, "rank", label, integer=True, minimum=1)

    require(
        ranks.is_unique,
        "candidates: duplicate ranks"
    )

    number(df, "pubchem_cid", label, integer=True, minimum=1)
    scores = number(df, "reversal_score", label)
    number(df, "raw_score", label)

    allowed(
        df, "approval_status",
        {"approved", "investigational", "unknown"},
        label
    )

    allowed(
        df, "used_in_oral_cancer",
        {"yes", "no", "unknown"},
        label
    )

    ordered_scores = scores.iloc[
        np.argsort(ranks.to_numpy())
    ].to_numpy()

    require(
        (np.diff(ordered_scores) <= 1e-12).all(),
        "candidates: higher reversal_score must have a better rank"
    )

    for drug, smiles in zip(df["drug_name"], df["smiles"]):
        require(
            Chem.MolFromSmiles(smiles) is not None,
            f"candidates: invalid SMILES for {drug}"
        )


def validate_binding(df, drugs):
    label = "binding"

    require(
        not df.duplicated(["drug_name", "target_gene"]).any(),
        "binding: duplicate drug-target pairs detected"
    )

    require_text(
        df,
        [
            "drug_name", "smiles", "target_gene",
            "pdb_id", "method", "is_known_target",
            "in_training_data"
        ],
        label
    )

    number(df, "predicted_affinity", label)

    confidence_levels = {"low", "medium", "high"}
    if not df["confidence"].isin(confidence_levels).all():
        number(df, "confidence", label, minimum=0, maximum=1)

    if "tree_sd" in df.columns:
        number(df, "tree_sd", label, minimum=0)

    if "max_similarity_to_train" in df.columns:
        number(
            df, "max_similarity_to_train", label,
            minimum=0, maximum=1
        )

    allowed(
        df, "method",
        {
            "docking", "qsar", "gnn",
            "RandomForest_Morgan2048_pChEMBL"
        },
        label
    )

    allowed(
        df, "is_known_target",
        {"yes", "no"},
        label
    )

    require(
        df["in_training_data"].eq("no").all(),
        "binding: in_training_data must be no"
    )

    require(
        set(df["drug_name"]) == set(drugs["drug_name"]),
        "binding: drug names must match the candidate list"
    )

    known = drugs.set_index("drug_name")

    for row in df.itertuples(index=False):
        drug = known.loc[row.drug_name]

        require(
            row.smiles == drug["smiles"],
            f"binding: SMILES mismatch for {row.drug_name}"
        )

        targets = {
            target.strip()
            for target in drug["known_targets"].split(";")
            if target.strip()
        }

        require(
            row.target_gene not in targets
            or row.is_known_target == "yes",
            f"binding: known-target mismatch for "
            f"{row.drug_name}/{row.target_gene}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Validate pipeline CSV contracts"
    )

    parser.add_argument(
        "--signature",
        type=Path,
        default=Path("pipeline/mock/oral_signature.csv")
    )

    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("pipeline/mock/candidate_drugs.csv")
    )

    parser.add_argument(
        "--binding",
        type=Path,
        default=Path("pipeline/mock/binding_scores.csv")
    )

    args = parser.parse_args()

    try:
        signature = load_csv(args.signature, SIGNATURE_COLUMNS)
        drugs = load_csv(args.candidates, DRUG_COLUMNS)
        binding = load_csv(args.binding, BINDING_COLUMNS)

        validate_signature(signature)
        print(f"PASS: Disease signature ({len(signature)} rows)")

        validate_drugs(drugs)
        print(f"PASS: Drug candidates ({len(drugs)} rows)")

        validate_binding(binding, drugs)
        print(f"PASS: Binding scores ({len(binding)} rows)")

        print("\nAll CSV contract checks passed.")
        print(
            "NOTE: Training-data and scaffold leakage still "
            "require an independent audit by Person C."
        )

    except (ValueError, OSError, pd.errors.ParserError) as error:
        parser.exit(1, f"FAIL: {error}\n")


if __name__ == "__main__":
    main()
