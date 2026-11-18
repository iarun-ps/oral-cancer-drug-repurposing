# Oral Cancer Drug Repurposing — Data Contract

Status: DRAFT — awaiting confirmation from Persons A, B and C.

## General Rules

- All files must use CSV format, UTF-8 encoding and comma separators.
- Column names must match exactly.
- Required values must not be missing.
- Gene symbols must follow HGNC naming.
- Drug names must match across files.
- SMILES must be valid molecular structures.
- Do not change column names without team agreement.

## Stage 1: Disease Signature (Person A)

File: signature/oral_cancer/results/oral_signature.csv

| Column | Type | Description |
|---|---|---|
| gene_symbol | string | HGNC gene symbol |
| entrez_id | integer | NCBI gene identifier |
| log2fc | float | Tumour versus normal fold change |
| padj | float | Adjusted p-value |
| direction | string | up or down |
| rank | integer | Rank within each direction |
| n_cohorts | integer | Number of agreeing cohorts (1 or 2) |

Validation:
- gene_symbol must not be empty.
- direction must be up or down.
- Positive log2fc means upregulated.
- Negative log2fc means downregulated.
- padj must be between 0 and 1.
- rank must be a positive integer.

Additional files:
- oral_signature_up.txt
- oral_signature_down.txt

Each text file contains one gene symbol per line.

## Stage 2: Drug Repurposing (Person B)

File: repurposing/results/candidate_drugs.csv

| Column | Type | Description |
|---|---|---|
| rank | integer | Drug reversal rank |
| drug_name | string | Drug name |
| pubchem_cid | integer | PubChem compound ID |
| smiles | string | Canonical SMILES |
| reversal_score | float | Higher means stronger reversal |
| raw_score | float | Original tool score |
| cell_line | string | Cell line or all |
| source_tool | string | LINCS analysis tool |
| known_targets | string | Semicolon-separated gene symbols |
| approval_status | string | approved, investigational or unknown |
| used_in_oral_cancer | string | yes, no or unknown |

Validation:
- drug_name must not be empty.
- smiles must not be empty and must parse using RDKit.
- rank must be a positive integer.
- Higher reversal_score must represent stronger reversal.
- Person B must document the raw_score conversion.
- approval_status and used_in_oral_cancer must use the listed values.

## Stage 3: Binding Validation (Person C)

File: binding/results/binding_scores.csv

| Column | Type | Description |
|---|---|---|
| drug_name | string | Must match candidate_drugs.csv |
| smiles | string | Must match the candidate SMILES |
| target_gene | string | Target protein gene symbol |
| pdb_id | string | Protein structure identifier |
| predicted_affinity | float | Higher means stronger binding |
| confidence | float | Prediction confidence |
| method | string | docking, qsar or gnn |
| is_known_target | string | yes or no |
| in_training_data | string | Must be no |

Validation:
- Every candidate must have a binding result.
- Drug names and SMILES must match Stage 2.
- confidence must be between 0 and 1.
- in_training_data must always be no.
- Binding affinities must use a documented comparable scale.
- Known-target results are sanity checks, not novel findings.
- Candidate and gold-list drugs and matching core scaffolds must be excluded from model training.
- Person C must provide evidence of the training-data leakage checks.

## Stage 4: Final Shortlist (Person D)

File: pipeline/results/final_shortlist.csv

Required columns:

- drug_name
- reversal_rank
- reversal_score
- best_target
- predicted_affinity
- confidence
- binding_rank
- combined_rank
- known_oral_cancer_drug
- evidence_note

Ranking:
1. Use the drug reversal rank from Stage 2.
2. Select the strongest predicted binding result among targets not already known for that drug.
3. Rank eligible drugs by predicted affinity.
4. Calculate combined_rank = reversal_rank * binding_rank.
5. Lower combined_rank is better.
6. Sort deterministically when scores are tied.
7. Do not invent binding scores for drugs without eligible targets.

Sensitivity analysis:
- Compare binding weights of 0%, 25%, 50%, 75% and 100%.
- Report top-10 overlap and rank correlation.

The known_oral_cancer_drug field comes from the frozen gold list
prepared by Person E in docs/known_oral_cancer_drugs.csv.

The evidence_note field is supplied or reviewed by Person E.

## Integration Requirements

- Validate every input before merging.
- Reject missing required columns and invalid values.
- Reject unmatched drug names or SMILES.
- Reject binding results marked in_training_data=yes.
- Do not use known-target predictions as novel binding evidence.
- Do not silently substitute missing data with zero.
- Mock outputs must be clearly identified as synthetic.
- Final results must be reproducible from the documented pipeline.

## Contract Approval

- Person A: Pending
- Person B: Pending
- Person C: Pending
- Person D: Draft prepared
