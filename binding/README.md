# Binding and target validation (Part C)

Predicts how strongly each reversal candidate (from `repurposing/results/candidate_drugs.csv`) binds
8 cancer-relevant proteins, and flags toxicity. This is **computational evidence to support a shortlist,
not proof that any drug works.**

## 1. Targets (8)
| Gene | UniProt | PDB (for optional docking) |
|---|---|---|
| EGFR | P00533 | 5HG8 |
| PIK3CA | P42336 | 8EXL |
| CDK6 | Q00534 | 5L2S |
| MET | P08581 | 9SXJ (.cif only) |
| FGFR1 | P11362 | 5EW8 |
| KDR | P35968 | 3WZD |
| MAP2K1 | Q02750 | 7B7R (verify) |
| HDAC2 | Q92769 | 4LY1 (verify) |

Targets were chosen from the oral cancer literature and checked in DepMap (release not recorded in our notes (to be added)) CRISPR gene
effect for head and neck lines: CDK6 and EGFR show strong dependency, PIK3CA moderate, the rest none.
A target with no dependency is still scored, but is weaker evidence of relevance.

## 2. Method
1. **Training data:** measured activities from **ChEMBL** (not BindingDB as first planned; a deliberate change).
   Metric: pChEMBL (9 = about 1 nM, 6 = about 1 uM). Assay filter: none in the download. The raw 31,415 rows are 93% binding assays, 7% functional cell assays, plus a few ADME/toxicity records; measurement types are 88% IC50, 8% Ki, 4% Kd and a few EC50/Potency/ED50 values. A binding-only rerun (assay_type = B, IC50/Ki/Kd) is a possible robustness check.
   Download was capped at 5,000 rows per protein, so for several proteins the data is **not a random sample**.
   After cleaning: 23,281 target-molecule records.
2. **Leakage control:** all 98 candidates were removed from training (same molecule by InChIKey first block,
   largest fragment only), and also every molecule sharing a Bemis-Murcko scaffold with a candidate (Level 2).
   12 of 98 candidates had ChEMBL data. Level 2 removed at most 1.5% of any protein's data, so Level 2 is the
   main training set (`binding/data/train_level2.csv`); Level 1 is `train_level1.csv`. The 12 measured candidates
   are kept as a separate test (`heldout_candidates_measured.csv`).
3. **Model:** one random forest (300 trees) per protein on Morgan fingerprints (radius 2, 2048 bits).
   `confidence` combines tree disagreement and similarity to the nearest training molecule.
4. **Toxicity:** random forest per assay on Tox21 (12 assays), class-balanced.

## 3. Results
**Scaffold-split test (whole scaffold families held out; 80/20):**

| Target | n train | n test | RMSE (model) | RMSE (predict-the-mean) | Spearman |
|---|---|---|---|---|---|
| CDK6 | 668 | 167 | 0.81 | 0.97 | 0.59 |
| EGFR | 2500 | 625 | 0.83 | 1.26 | 0.75 |
| FGFR1 | 2947 | 737 | 0.86 | 1.29 | 0.73 |
| HDAC2 | 2260 | 554 | 0.88 | 1.08 | 0.53 |
| KDR | 2916 | 729 | 0.88 | 1.21 | 0.64 |
| MAP2K1 | 865 | 217 | 0.92 | 1.45 | 0.66 |
| MET | 2901 | 726 | 0.96 | 1.33 | 0.67 |
| PIK3CA | 3431 | 858 | 0.82 | 1.21 | 0.72 |

The model beats the predict-the-mean baseline for every protein. A typical error is about 0.8 to 0.96 pChEMBL
units (roughly 6- to 9-fold in binding strength).

**Held-out known drugs (39 measured pairs from 12 candidates, never in training):**
Spearman 0.45, RMSE 1.07, mean error (predicted - measured) -0.44. About 33 distinct pairs
(some repeat because ChEMBL lists several entries of one molecule). Small and noisy: a sanity check only.

**Toxicity (Tox21, scaffold split, AUC; 0.5 = coin flip):**

| Assay | AUC |
|---|---|
| NR-AR | 0.69 |
| NR-AR-LBD | 0.82 |
| NR-AhR | 0.81 |
| NR-Aromatase | 0.71 |
| NR-ER | 0.64 |
| NR-ER-LBD | 0.68 |
| NR-PPAR-gamma | 0.77 |
| SR-ARE | 0.70 |
| SR-ATAD5 | 0.76 |
| SR-HSE | 0.66 |
| SR-MMP | 0.81 |
| SR-p53 | 0.74 |

Mean AUC 0.73.

## 4. Limitations (please read before using the scores)
- **Use `predicted_affinity` as a coarse signal** (about 7 or higher: plausible binder; 6 or lower: probably weak),
  not as an exact potency. The model **under-predicts strong binders** (shrinks toward the mean) and
  **cannot tell that some kinase inhibitors do not bind PIK3CA** (a lipid kinase).
- **Training labels mix assay types** (about 7% are functional, not binding, assays; IC50, Ki and Kd are pooled as pChEMBL). This adds noise.
- **738 of 784 pairs are low confidence**: most candidates are chemically unlike the training molecules.
- `is_known_target` = "yes" if ChEMBL lists the gene as a target OR a measured pChEMBL of 6 or more exists.
  "no" means "not listed", not "does not bind".
- Antibody candidates (for example cetuximab) cannot go through this small-molecule pipeline.
- Toxicity predictions are weak warning flags from cell assays (mean AUC 0.73). Only 5 of 98 candidates
  are in Tox21 (those are partly lookups). No flag does **not** mean safe.
- Some candidates (for example menadione, actinomycin-D) probably rank high because they are generally cytotoxic,
  not because they are specific to oral cancer.
- Docking and a GNN were not done. No wet-lab validation.
- Rerun Step 5 (leakage filter) when `docs/known_oral_cancer_drugs.csv` is frozen, so gold-list drugs are also excluded.

## 5. Files
| File | What |
|---|---|
| `results/binding_scores.csv` | Data-contract output: 98 drugs x 8 targets = 784 rows (+ extra columns `tree_sd`, `max_similarity_to_train`) |
| `results/binding_summary.csv` | One row per drug: best target, confidence, toxicity flags |
| `results/baseline_metrics.csv`, `heldout_check.csv`, `tox_metrics.csv`, `tox_predictions.csv`, `leakage_filter_summary.csv` | Evidence tables |
| `train_models.py` | Rebuilds `models/rf_<GENE>.joblib` (about 10 min, fixed seed). Models are not in Git (about 570 MB); ask Laniya for the zip |
| `featurize.py` | Rebuilds the exact fingerprints for the models |
| `data/` | Cleaned ChEMBL data and the leakage-filtered training sets |
| `figures/` | Plots used in the report |
| `targets.md` | Target choice and notes |
