import gzip
from pathlib import Path
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

RAW_DIR = Path("../data/raw/oral_cancer")
OUT_DIR = Path("../signature/oral_cancer/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Metadata
gsm_ids, titles = [], []
with gzip.open(RAW_DIR / "GSE101547_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsm_ids = [x.strip('"') for x in line.strip().split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.strip().split("\t")[1:]]

counts_raw = pd.read_csv(RAW_DIR / "GSE101547_raw_counts_htseq.txt.gz", sep="\t")
counts_raw.columns = [c.strip() for c in counts_raw.columns]
counts_raw = counts_raw.rename(columns={counts_raw.columns[0]: "gene_id"})

id_col = titles if any(t in counts_raw.columns for t in titles) else gsm_ids
meta = pd.DataFrame({
    "sample_id": id_col,
    "patient": [t[:-1] for t in titles],
    "condition": ["tumour" if t.endswith("D") else "normal" for t in titles],
})
meta = meta[meta["patient"] != "23"].copy()

pair_counts = meta.groupby(["patient", "condition"]).size().unstack(fill_value=0)
assert (pair_counts["tumour"] == 1).all() and (pair_counts["normal"] == 1).all()
meta.to_csv(OUT_DIR / "gse101547_samples.csv", index=False)
meta = meta.set_index("sample_id")

# 2. Clean counts & map symbols
counts_raw = counts_raw[~counts_raw["gene_id"].astype(str).str.startswith("__")].copy()
counts_raw["gene_id"] = counts_raw["gene_id"].astype(str).str.split(".").str[0]

annot = pd.read_csv(RAW_DIR / "Human.GRCh38.p13.annot.tsv.gz", sep="\t", low_memory=False)
annot_key = "EnsemblGeneID" if counts_raw["gene_id"].str.startswith("ENSG").any() else "GeneID"
annot[annot_key] = annot[annot_key].astype(str).str.split(".").str[0]

merged = counts_raw.merge(
    annot[[annot_key, "Symbol"]].dropna().drop_duplicates(subset=annot_key),
    left_on="gene_id", right_on=annot_key
)
merged = merged[merged["Symbol"].str.strip() != ""]

valid_samples = meta.index.tolist()
merged["total_expr"] = merged[valid_samples].sum(axis=1)
clean_genes = (
    merged.sort_values("total_expr", ascending=False)
    .drop_duplicates(subset="Symbol", keep="first")
    .set_index("Symbol")[valid_samples]
)

# Transpose to samples x genes & enforce alignment + >=10 count pre-filter
counts = clean_genes.T.astype(int)
assert list(counts.index) == list(meta.index), "Sample alignment mismatch!"
counts = counts.loc[:, counts.sum(axis=0) >= 10]

# 3. Run PyDESeq2
meta["patient"] = meta["patient"].astype(str)
meta["condition"] = pd.Categorical(meta["condition"], categories=["normal", "tumour"])

dds = DeseqDataSet(counts=counts, metadata=meta, design="~ patient + condition", refit_cooks=True)
dds.deseq2()

stat = DeseqStats(dds, contrast=["condition", "tumour", "normal"])
stat.summary()
res = stat.results_df
res.to_csv(OUT_DIR / "gse101547_de_results.csv")

# 4. Step 3: Preliminary signature for Person B (top 150 up & top 150 down)
sig = res.dropna(subset=["padj"]).query("padj < 0.05 and abs(log2FoldChange) >= 1").copy()
top_up = sig[sig["log2FoldChange"] > 0].sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(150)
top_down = sig[sig["log2FoldChange"] < 0].sort_values(["padj", "log2FoldChange"], ascending=[True, True]).head(150)

prelim = pd.concat([top_up.assign(direction="up"), top_down.assign(direction="down")])
prelim.to_csv(OUT_DIR / "gse101547_preliminary_signature.csv")
print(f"Saved preliminary signature: {len(top_up)} up, {len(top_down)} down.")