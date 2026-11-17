import gzip
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ==========================================
# STEP 1: Build sample metadata & clean counts
# ==========================================

# 1. Read GSM IDs and sample titles from the series matrix
gsm_ids, titles = [], []
with gzip.open("GSE101547_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsm_ids = [x.strip('"') for x in line.strip().split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.strip().split("\t")[1:]]

# 2. Load raw counts and detect whether columns use titles (2D) or GSM IDs
counts_raw = pd.read_csv("GSE101547_raw_counts_htseq.txt.gz", sep="\t")
counts_raw.columns = [c.strip() for c in counts_raw.columns]
counts_raw = counts_raw.rename(columns={counts_raw.columns[0]: "gene_id"})

id_col = titles if any(t in counts_raw.columns for t in titles) else gsm_ids

meta = pd.DataFrame({
    "sample_id": id_col,
    "patient": [t[:-1] for t in titles],
    "condition": ["tumour" if t.endswith("D") else "normal" for t in titles],
})

# 3. Drop patient 23 and assert every remaining patient has 1 tumour & 1 normal
meta = meta[meta["patient"] != "23"].copy()

pair_counts = meta.groupby(["patient", "condition"]).size().unstack(fill_value=0)
assert (pair_counts["tumour"] == 1).all(), "Missing or duplicate tumour sample!"
assert (pair_counts["normal"] == 1).all(), "Missing or duplicate normal sample!"
assert len(meta) == 22, f"Expected 22 samples, got {len(meta)}"

meta.to_csv("gse101547_samples.csv", index=False)
meta = meta.set_index("sample_id")
print("Saved gse101547_samples.csv (22 samples, 11 matched pairs).")

# 4. Clean gene IDs, map to symbols, and deduplicate
# Drop HTSeq summary rows (e.g., __no_feature) if present
counts_raw = counts_raw[~counts_raw["gene_id"].astype(str).str.startswith("__")].copy()

# Strip Ensembl version suffix (.12) if present
counts_raw["gene_id"] = counts_raw["gene_id"].astype(str).str.split(".").str[0]

annot = pd.read_csv("Human.GRCh38.p13.annot.tsv.gz", sep="\t", low_memory=False)
annot_key = "EnsemblGeneID" if counts_raw["gene_id"].str.startswith("ENSG").any() else "GeneID"
annot[annot_key] = annot[annot_key].astype(str).str.split(".").str[0]

merged = counts_raw.merge(
    annot[[annot_key, "Symbol"]].dropna().drop_duplicates(subset=annot_key),
    left_on="gene_id",
    right_on=annot_key,
)
merged = merged[merged["Symbol"].str.strip() != ""]

# Keep only the 22 valid samples & resolve duplicate symbols (keep highest-expressed row)
valid_samples = meta.index.tolist()
merged["total_expr"] = merged[valid_samples].sum(axis=1)
clean_genes = (
    merged.sort_values("total_expr", ascending=False)
    .drop_duplicates(subset="Symbol", keep="first")
    .set_index("Symbol")[valid_samples]
)

# Transpose to (samples x genes) with integer counts as required by PyDESeq2
counts = clean_genes.T.astype(int)
print(f"Cleaned counts matrix shape (samples x genes): {counts.shape}")

# ==========================================
# STEP 2: Differential Expression (PyDESeq2)
# ==========================================
meta["patient"] = meta["patient"].astype(str)
meta["condition"] = pd.Categorical(meta["condition"], categories=["normal", "tumour"])

dds = DeseqDataSet(
    counts=counts,
    metadata=meta,
    design="~ patient + condition",
    refit_cooks=True,
)
dds.deseq2()

stat = DeseqStats(dds, contrast=["condition", "tumour", "normal"])
stat.summary()
res = stat.results_df
res.to_csv("gse101547_de_results.csv")
print("Done! Saved gse101547_de_results.csv")