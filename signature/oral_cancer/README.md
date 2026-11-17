# Disease signature: Indian oral squamous cell carcinoma (Person A)

Turns two public Indian RNA-seq cohorts into one ranked list of up- and down-regulated genes
(tumour vs normal, paired by patient). This is the input for drug matching (Person B).

## Outputs

| File | What it is |
|---|---|
| `results/oral_signature.csv` | Final consensus signature in the data-contract format (Section 8) |
| `results/oral_signature_up.txt` / `_down.txt` | Top 150 gene symbols per direction, one per line |
| `results/preliminary/gse101547_preliminary_signature.csv` | Cohort-1-only signature, **preliminary**, for B's rough run |
| `results/consensus_details/oral_signature_annotated.csv` | Every consensus gene with its class flag and whether it was excluded |
| `results/consensus_details/replication_summary.txt` | Replication numbers (Section 11) |
| `results/cohort_data/` | Per-cohort sample tables, full DE results, sample QC tables |
| `results/figures/` | PCA and volcano plots for both cohorts |

`oral_signature.csv` columns: `gene_symbol, entrez_id, log2fc, padj, direction, rank, n_cohorts`.
A positive `log2fc` means higher in tumour. `rank` is 1 = strongest, **within each direction**.

## Data

| Cohort | Samples used | Notes |
|---|---|---|
| GSE101547 | 11 tumour + 11 matched normal | Patient 23 dropped: tumour sample missing from the counts table |
| GSE213862 | 46 tumour + 46 matched normal | 43 margin samples (title ending in M) ignored |

Links and accessions are in `data_sources.md`. Raw data are not stored in this repo.

## Methods

1. **Sample tables.** Condition and patient were taken from the sample titles (GSE101547: last letter D = tumour,
   otherwise normal; GSE213862: last letter T / N / M). Pairs were checked in code: every patient has exactly one
   tumour and one normal, 11 and 46 pairs respectively.
2. **Counts.** Raw counts were mapped to HGNC symbols with the NCBI GRCh38.p13 annotation table; duplicate symbols
   keep the highest-expressed row. [TODO: confirm whether the GSE213862 counts file holds integers; the script rounds
   it to integers.] Low-count filter: a gene is kept only if it has at least 10 reads in at least 10% of samples
   (minimum 3). This kept 17,724 genes in GSE101547 and 18,333 in GSE213862.
3. **Differential expression.** PyDESeq2, paired design `~ patient + condition`, Wald test, contrast tumour vs
   normal (normal = reference), default multiple-testing correction. Each cohort is analysed separately, never pooled.
4. **QC.** PCA on variance-stabilised counts with matched pairs joined by lines, volcano plots, and a per-sample
   muscle/fat marker score (MYH2, ACTA1, DES, MB, CASQ1, CKM, ADIPOQ, PLIN1, CIDEC, FABP4, mean z-score of VST values).
5. **Consensus.** Genes tested in both cohorts with baseMean >= 10 in both. A gene is kept if padj < 0.05 and
   |log2FC| >= 1 in both cohorts and the direction agrees.
   - `log2fc` = mean of the two cohorts' log2 fold changes.
   - `padj` = the **weaker (larger)** of the two cohorts' adjusted p-values.
   - **Ranking** within each direction by the combined Wald statistic |stat_101547 + stat_213862| / sqrt(2),
     so that precisely estimated genes outrank noisy low-count genes with large fold changes.
6. **Class filter.** Immunoglobulin genes, a hand-curated list of skeletal-muscle and adipose genes, and non-coding
   genes (LINC, LOC, MIR, -AS and similar) are left out of `oral_signature.csv` (kept, flagged, in the annotated file).
   Reason: immunoglobulin genes reflect immune infiltration, muscle/fat genes reflect tissue composition of the
   normal biopsies, and L1000 covers few non-coding genes. Genes without an Entrez ID are dropped.
7. **Stability (leave-one-cohort-out).** A signature of the top 150 genes per direction was built from one cohort and
   tested in the other.

## Results

| Item | Value |
|---|---|
| Shared tested genes (baseMean >= 10 in both) | 15,272 |
| Significant in GSE101547 / GSE213862 | 2,995 / 3,580 |
| Significant in both | 2,215 |
| Same direction in both | 2,202 (99.41%) |
| Spearman r of log2FC, all shared genes | 0.778 |
| Spearman r of log2FC, consensus genes | 0.883 |
| Genes in `oral_signature.csv` | 2,013 (657 up, 1,356 down) |
| Excluded by class | immunoglobulin 52, muscle/adipose 62, non-coding 75 |
| LOCO GSE101547 -> GSE213862 | 100% same direction, 98.7% also significant, Spearman 0.922 |
| LOCO GSE213862 -> GSE101547 | 100% same direction, 99.0% also significant, Spearman 0.942 |

Known cancer markers behave as expected (log2FC, GSE101547 / GSE213862): KRT4 -9.25 / -5.43, KRT13 -6.32 / -4.61
(down); KRT17 3.83 / 4.14, MMP1 7.51 / 6.58, MMP3 3.40 / 4.63, MMP10 5.69 / 5.19 (up).

## Limitations

- **Tissue composition.** Several "normal" biopsies are rich in skeletal muscle or fat (GSE213862: RP71N, RP86N,
  RP78N, RP74N, RP57N and others have the highest marker scores and form the far-right cluster on PC1 of the PCA).
  Down-regulated genes therefore partly reflect the type of normal tissue sampled, not only tumour-intrinsic change.
  The class filter removes the best-known muscle/fat genes but the curated list is not exhaustive.
- GSE101547 has only 11 pairs, so many genes are significant in GSE213862 only. The two cohorts come from different
  laboratories and were not batch-corrected against each other.
- Tobacco status is only partly labelled, so the tobacco YES/NO comparison was not run.
- Results are computational hypotheses for the downstream steps; nothing here is medical advice.

## Rerun

From the repo root, with the raw files listed in `data_sources.md` placed in `data/raw/oral_cancer/`:

```
pip install pydeseq2 pandas numpy scipy scikit-learn matplotlib adjustText
python scripts/run_de.py all
python scripts/signature2.py
```

Settings (thresholds, `TOP_N`, excluded classes) are at the top of `signature2.py`.

## Software versions

Python [TODO], pydeseq2 [TODO], pandas [TODO], numpy [TODO], scipy [TODO], scikit-learn [TODO], matplotlib [TODO].
Run date: 17 Nov 2026.
