"""Integrate oral cancer drug repurposing and binding results."""

import argparse
from pathlib import Path

import pandas as pd

from pipeline.rank_sensitivity import sensitivity_analysis
from pipeline.validate_contract import (
    BINDING_COLUMNS,
    DRUG_COLUMNS,
    SIGNATURE_COLUMNS,
    load_csv,
    validate_binding,
    validate_drugs,
    validate_signature,
)


def run_pipeline(signature_path, drugs_path, binding_path, output_path,
                 gold_path=None, mock=False):

    # Load and validate inputs
    signature = load_csv(signature_path, SIGNATURE_COLUMNS)
    drugs = load_csv(drugs_path, DRUG_COLUMNS)
    binding = load_csv(binding_path, BINDING_COLUMNS)

    validate_signature(signature)
    validate_drugs(drugs)
    validate_binding(binding, drugs)

    print(f"PASS: Signature ({len(signature)} genes)")
    print(f"PASS: Candidates ({len(drugs)} drugs)")
    print(f"PASS: Binding ({len(binding)} records)")

    # Exclude already-known drug-target interactions
    eligible = binding[
        binding["is_known_target"] == "no"
    ].copy()

    if eligible.empty:
        raise ValueError("No eligible binding predictions found.")

    # Convert scores to numbers
    eligible["predicted_affinity"] = pd.to_numeric(
        eligible["predicted_affinity"]
    )

    eligible["confidence"] = pd.to_numeric(
        eligible["confidence"]
    )

    # Select strongest eligible binding prediction per drug
    eligible = eligible.sort_values(
        ["predicted_affinity", "confidence", "target_gene", "pdb_id"],
        ascending=[False, False, True, True]
    )

    best_binding = eligible.drop_duplicates(
        subset=["drug_name"],
        keep="first"
    ).copy()

    # Assign binding ranks
    best_binding = best_binding.sort_values(
        ["predicted_affinity", "confidence", "drug_name"],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    best_binding["binding_rank"] = range(
        1, len(best_binding) + 1
    )

    # Merge drug ranking with binding results
    merged = drugs.merge(
        best_binding[
            [
                "drug_name", "target_gene",
                "predicted_affinity", "confidence",
                "binding_rank"
            ]
        ],
        on="drug_name",
        how="inner",
        validate="one_to_one"
    )

    missing = sorted(
        set(drugs["drug_name"]) - set(merged["drug_name"])
    )

    if missing:
        print("WARNING: Drugs without eligible targets:", missing)

    # Calculate rank product
    merged["reversal_rank"] = pd.to_numeric(
        merged["rank"]
    ).astype(int)

    merged["reversal_score"] = pd.to_numeric(
        merged["reversal_score"]
    )

    merged["combined_rank"] = (
        merged["reversal_rank"] * merged["binding_rank"]
    )

    merged = merged.rename(
        columns={"target_gene": "best_target"}
    )

    # Add known oral cancer drug status
    if mock:
        merged["known_oral_cancer_drug"] = "no"
        merged["evidence_note"] = "SYNTHETIC TEST DATA"
    else:
        if gold_path is None or not gold_path.is_file():
            raise ValueError(
                "Real runs require --gold with Person E's "
                "frozen known-drug CSV."
            )

        gold = pd.read_csv(gold_path, dtype=str)

        if "drug_name" not in gold.columns:
            raise ValueError(
                "Gold-list CSV requires a drug_name column."
            )

        known_drugs = set(gold["drug_name"].dropna())

        merged["known_oral_cancer_drug"] = (
            merged["drug_name"]
            .isin(known_drugs)
            .map({True: "yes", False: "no"})
        )

        merged["evidence_note"] = "PENDING_PERSON_E_REVIEW"

    # Sort final results
    merged = merged.sort_values(
        ["combined_rank", "reversal_rank", "binding_rank", "drug_name"]
    )

    columns = [
        "drug_name",
        "reversal_rank",
        "reversal_score",
        "best_target",
        "predicted_affinity",
        "confidence",
        "binding_rank",
        "combined_rank",
        "known_oral_cancer_drug",
        "evidence_note"
    ]

    final = merged[columns]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Generate ranking sensitivity analysis
    summary, rankings = sensitivity_analysis(final)

    final.to_csv(output_path, index=False)

    summary_path = output_path.with_name(
        f"{output_path.stem}_sensitivity_summary.csv"
    )
    rankings_path = output_path.with_name(
        f"{output_path.stem}_sensitivity_rankings.csv"
    )

    summary.to_csv(summary_path, index=False)
    rankings.to_csv(rankings_path, index=False)

    print("PASS: Ranking sensitivity analysis")
    print(f"Saved: {summary_path}")
    print(f"Saved: {rankings_path}")

    print(f"\nPASS: Generated {len(final)} shortlisted drugs")
    print(f"Saved: {output_path}")

    return final


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mock", action="store_true")
    parser.add_argument("--gold", type=Path)
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    if args.mock:
        signature = Path("pipeline/mock/oral_signature.csv")
        drugs = Path("pipeline/mock/candidate_drugs.csv")
        binding = Path("pipeline/mock/binding_scores.csv")
        output = Path("pipeline/results/mock_shortlist.csv")
    else:
        signature = Path(
            "signature/oral_cancer/results/oral_signature.csv"
        )
        drugs = Path("repurposing/results/candidate_drugs.csv")
        binding = Path("binding/results/binding_scores.csv")
        output = Path("pipeline/results/final_shortlist.csv")

    if args.output:
        output = args.output

    if args.mock and output.name == "final_shortlist.csv":
        parser.error(
            "Mock data cannot be saved as final_shortlist.csv."
        )

    run_pipeline(
        signature,
        drugs,
        binding,
        output,
        gold_path=args.gold,
        mock=args.mock
    )


if __name__ == "__main__":
    main()
