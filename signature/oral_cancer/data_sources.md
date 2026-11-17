# Data sources (Part A)

| Dataset | Accession page | Files used |
|---|---|---|
| GSE101547 (RNA-seq, Indian OSCC) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101547 | `GSE101547_series_matrix.txt.gz`, `GSE101547_raw_counts_htseq.txt.gz` |
| GSE213862 (RNA-seq, Indian OSCC) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213862 | `GSE213862_series_matrix.txt.gz`, `GSE213862_PPM1_AC1_rawdata.txt.gz` |
| Gene annotation (GRCh38.p13) | linked from the GEO RNA-seq counts page | `Human.GRCh38.p13.annot.tsv.gz` |

For each dataset, copy the following from its GEO page before submission (needed for `docs/ethics-and-data.md`):

- Publication / PubMed ID: [TODO]
- Consent / ethics statement given by the original study: [TODO]
- Terms of use / citation requirement: [TODO]

Only public, de-identified sample IDs and group labels from GEO are used. No raw data are stored in this repository.
