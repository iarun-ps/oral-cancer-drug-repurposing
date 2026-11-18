# Omics-Driven Drug Repurposing for Indian Oral Cancer

## Project Overview

This BioAI capstone project investigates potential drug repurposing candidates for Indian oral squamous cell carcinoma (OSCC).

The project integrates:
1. RNA-seq differential gene expression analysis.
2. Drug signature reversal using LINCS L1000.
3. Computational drug-target binding validation.
4. Drug ranking, evaluation and interpretation.

The objective is to produce a reproducible, evidence-based shortlist of potential drug candidates for further research.

Predicted candidates are research hypotheses, not clinically validated treatments.

## Project Workflow

```text
Indian Oral Cancer RNA-seq Data
             |
             v
1. Disease Signature Analysis
   Identify upregulated and downregulated genes
             |
             v
2. Drug Repurposing
   Find drugs that reverse the disease signature
             |
             v
3. Binding Validation
   Predict interactions with cancer-relevant proteins
             |
             v
4. Pipeline Integration
   Validate inputs and combine drug rankings
             |
             v
5. Evaluation and Reporting
   Evaluate and explain shortlisted candidates
```

## Datasets

| Dataset | Purpose |
|---|---|
| GSE101547 | Indian OSCC RNA-seq, 11 matched pairs |
| GSE213862 | Indian OSCC RNA-seq, 46 tumour-normal pairs |
| TCGA-HNSC | Independent signature replication analysis |
| LINCS L1000 | Drug-induced gene expression signatures |
| BindingDB | Experimental drug-target binding measurements |
| RCSB PDB | Protein structures |
| DepMap | Supporting target dependency evidence |

## Team Responsibilities

| Member | Responsibility |
|---|---|
| Person A - Nivedh | Disease signature analysis |
| Person B - Laniya | Drug repurposing |
| Person C - Archana | Binding validation |
| Person D - Arun | Pipeline integration, testing and GitHub management |
| Person E - Savin | Evaluation, explainability and reporting |

## Repository Structure

```text
oral-cancer-drug-repurposing/
|
|-- data/                 Data sources and provenance
|-- signature/            Disease signature analysis
|-- repurposing/          Drug candidate identification
|-- binding/              Drug-target binding validation
|-- pipeline/             Integration and automated testing
|   |-- mock/             Synthetic test datasets
|   |-- tests/            Automated pytest tests
|   |-- results/          Generated shortlist and analysis
|   |-- validate_contract.py
|   |-- run_pipeline.py
|   |-- rank_sensitivity.py
|
|-- analysis/             Evaluation and explainability
|-- report/               Figures, report and slides
|-- docs/                 Data contracts and documentation
|-- scripts/              Analysis scripts
|-- .github/workflows/    GitHub Actions CI
|-- requirements.txt
|-- requirements-signature.txt
```

## Installation

Recommended Python version: 3.13.

Clone the repository:

```bash
git clone https://github.com/iarun-ps/oral-cancer-drug-repurposing.git
cd oral-cancer-drug-repurposing
```

Create a virtual environment:

```bash
py -3.13 -m venv .venv
```

Activate it in Windows Git Bash:

```bash
source .venv/Scripts/activate
```

Install core dependencies:

```bash
python -m pip install -r requirements.txt
```

Additional dependencies for disease signature analysis:

```bash
python -m pip install -r requirements-signature.txt
```

## Running the Pipeline

### 1. Validate Mock Data

```bash
python pipeline/validate_contract.py
```

### 2. Run Mock Integration Pipeline

```bash
python -m pipeline.run_pipeline --mock
```

This generates:

- pipeline/results/mock_shortlist.csv
- pipeline/results/mock_shortlist_sensitivity_summary.csv
- pipeline/results/mock_shortlist_sensitivity_rankings.csv

Mock data is synthetic and must not be interpreted as scientific evidence.

### 3. Run Automated Tests

```bash
python -m pytest pipeline/tests/ -q
```

### 4. Check Code Quality

```bash
python -m ruff check pipeline/
```

### 5. Run the Pipeline on Real Results

Required inputs:

- signature/oral_cancer/results/oral_signature.csv
- repurposing/results/candidate_drugs.csv
- binding/results/binding_scores.csv
- docs/known_oral_cancer_drugs.csv

Once all required inputs have been finalized and validated:

```bash
python -m pipeline.run_pipeline --gold docs/known_oral_cancer_drugs.csv
```

The integration script reads the outputs of the scientific stages.
It does not independently rerun RNA-seq analysis, LINCS queries or binding-model training.

## Drug Ranking Method

The pipeline combines two rankings:

- Reversal rank: How strongly a drug reverses the disease signature.
- Binding rank: How strongly a drug is predicted to bind an eligible target.

The combined score is:

```text
Rank Product = Reversal Rank × Binding Rank
```

Lower rank products indicate higher-priority candidates.

Known drug-target interactions are excluded when selecting novel binding hypotheses.

## Ranking Sensitivity Analysis

The pipeline also compares binding weights of:

- 0%
- 25%
- 50%
- 75%
- 100%

It calculates ranking correlation and top-10 overlap to assess ranking stability.

## Data Validation

The validator checks:

- Required CSV columns and data types.
- Missing gene symbols.
- Invalid molecular SMILES.
- Drug reversal-score ordering.
- Binding-score confidence.
- Drug-name and SMILES consistency.
- Declared training-data membership.
- Known versus predicted drug-target interactions.

The validator cannot independently prove that training-data or scaffold leakage has been prevented. Binding-model training requires a separate leakage audit.

CSV specifications are documented in:

`docs/data-contract.md`

The contract requires confirmation from the responsible team members.

## Automated Testing

GitHub Actions runs automated checks on pushes and pull requests.

These include:

- Python dependency installation.
- Ruff code-quality checks.
- Mock CSV validation.
- Disease signature validation.
- Mock pipeline execution.
- Pytest tests.

## Current Integration Status

The disease signature, integration pipeline and mock-data tests are available.

Before the final real-data run:

- Person B must resolve missing PubChem IDs for MBCQ and AZD-4547, or an approved contract revision must define how missing IDs are handled.
- Person C must provide validated binding results.
- Person E must provide the frozen known-drug reference list.
- The team must approve the data contract.
- The final integrated results must be evaluated and reviewed.

## Reproducibility

All pipeline inputs must follow the agreed data contract.

The project uses version-controlled code, automated tests and documented data sources.

Large raw datasets and sensitive information should not be committed to GitHub.

## Limitations

- Computational predictions do not demonstrate clinical effectiveness.
- Binding predictions require experimental validation.
- Results depend on dataset quality and model assumptions.
- The pipeline is intended for research and educational purposes.

## License

No license has been selected yet.
