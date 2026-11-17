"""Paired tumour-vs-normal differential expression + QC plots for one cohort.

Usage (from any folder):
    python run_de.py GSE101547
    python run_de.py GSE213862
    python run_de.py all

Outputs (in signature/oral_cancer/results/):
    cohort_data/<gse>_samples.csv, <gse>_de_results.csv, <gse>_sample_qc.csv
    figures/<gse>_pca.png, <gse>_volcano.png
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

try:
    from adjustText import adjust_text  # pip install adjustText
except ImportError:
    adjust_text = None

# repo root = the nearest parent folder that contains data/raw (so the script can live in any subfolder)
_here = Path(__file__).resolve()
ROOT = next((p for p in _here.parents if (p / "data" / "raw").is_dir()), _here.parent.parent)
RAW_DIR = ROOT / "data" / "raw" / "oral_cancer"
RES_DIR = ROOT / "signature" / "oral_cancer" / "results"
DATA_DIR = RES_DIR / "cohort_data"   # per-cohort tables (samples, DE results, sample QC)
FIG_DIR = RES_DIR / "figures"        # PCA + volcano plots
for d in (RES_DIR, DATA_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Low-count filter: keep a gene only if it has >= MIN_COUNT reads in
# >= max(3, MIN_SAMPLE_FRAC * n_samples) samples (was: >=10 total reads, too weak)
MIN_COUNT = 10
MIN_SAMPLE_FRAC = 0.10

# Marker genes used only to flag muscle/fat-rich samples in the QC table
TISSUE_MARKERS = ["MYH2", "ACTA1", "DES", "MB", "CASQ1", "CKM", "ADIPOQ", "PLIN1", "CIDEC", "FABP4"]


# ---------------------------------------------------------------- metadata
def meta_gse101547(sample_ids, titles):
    meta = pd.DataFrame({
        "sample_id": sample_ids,
        "patient": [t[:-1] for t in titles],
        "condition": ["tumour" if t.endswith("D") else "normal" for t in titles],
    })
    return meta[meta["patient"] != "23"].copy()  # patient 23: tumour missing from counts


def meta_gse213862(sample_ids, titles):
    label_map = {"T": "tumour", "N": "normal", "M": "margin"}
    meta = pd.DataFrame({
        "sample_id": sample_ids,
        "patient": [t[:-1] for t in titles],
        "condition": [label_map.get(t[-1], "unknown") for t in titles],
    })
    print(f"Dropping {(meta['condition'] == 'margin').sum()} margin ('M') samples...")
    return meta[meta["condition"] != "margin"].copy()


COHORTS = {
    "GSE101547": dict(matrix="GSE101547_series_matrix.txt.gz",
                      counts="GSE101547_raw_counts_htseq.txt.gz",
                      build_meta=meta_gse101547, n_pairs=11, round_counts=False),
    "GSE213862": dict(matrix="GSE213862_series_matrix.txt.gz",
                      counts="GSE213862_PPM1_AC1_rawdata.txt.gz",
                      build_meta=meta_gse213862, n_pairs=46, round_counts=True),
}


# ---------------------------------------------------------------- loading
def load_inputs(cfg):
    gsm_ids, titles = [], []
    with gzip.open(RAW_DIR / cfg["matrix"], "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsm_ids = [x.strip('"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"').strip() for x in line.strip().split("\t")[1:]]

    counts_raw = pd.read_csv(RAW_DIR / cfg["counts"], sep="\t")
    counts_raw.columns = [c.strip() for c in counts_raw.columns]
    counts_raw = counts_raw.rename(columns={counts_raw.columns[0]: "gene_id"})

    id_col = titles if any(t in counts_raw.columns for t in titles) else gsm_ids
    meta = cfg["build_meta"](id_col, titles)

    assert set(meta["condition"]).issubset({"tumour", "normal"}), meta["condition"].unique()
    pairs = meta.groupby(["patient", "condition"]).size().unstack(fill_value=0)
    bad = pairs[(pairs["tumour"] != 1) | (pairs["normal"] != 1)]
    if not bad.empty:
        print("Incomplete pairs:\n", bad)
    assert bad.empty, "Not every patient has exactly 1 tumour and 1 normal!"
    assert len(pairs) == cfg["n_pairs"], f"Expected {cfg['n_pairs']} pairs, found {len(pairs)}"
    assert all(s in counts_raw.columns for s in meta["sample_id"]), "Sample IDs missing from counts table!"
    return meta, counts_raw


def counts_by_symbol(counts_raw, samples, round_counts):
    counts_raw = counts_raw[~counts_raw["gene_id"].astype(str).str.startswith("__")].copy()
    counts_raw["gene_id"] = counts_raw["gene_id"].astype(str).str.strip().str.split(".").str[0]

    annot = pd.read_csv(RAW_DIR / "Human.GRCh38.p13.annot.tsv.gz", sep="\t", low_memory=False)
    if counts_raw["gene_id"].str.startswith("ENSG").any():
        key = "EnsemblGeneID"
    elif counts_raw["gene_id"].str.isdigit().mean() > 0.5:
        key = "GeneID"
    else:
        key = None

    if key:
        annot[key] = annot[key].astype(str).str.split(".").str[0]
        merged = counts_raw.merge(
            annot[[key, "Symbol"]].dropna().drop_duplicates(subset=key),
            left_on="gene_id", right_on=key)
    else:
        merged = counts_raw.rename(columns={"gene_id": "Symbol"}).dropna(subset=["Symbol"])

    merged = merged[merged["Symbol"].astype(str).str.strip() != ""]
    merged["total_expr"] = merged[samples].sum(axis=1)
    clean = (merged.sort_values("total_expr", ascending=False)
             .drop_duplicates(subset="Symbol", keep="first")
             .set_index("Symbol")[samples])
    counts = clean.T.round().astype(int) if round_counts else clean.T.astype(int)
    return counts


# ---------------------------------------------------------------- plots
def pca_plot(coords, var, meta, out_png, title):
    df = meta.copy()
    df["PC1"], df["PC2"] = coords[:, 0], coords[:, 1]
    fig, ax = plt.subplots(figsize=(9, 8))
    for _, g in df.groupby("patient"):
        ax.plot(g["PC1"], g["PC2"], color="lightgray", linestyle="--", linewidth=0.7, zorder=1)
    colors = {"normal": "#1f77b4", "tumour": "#d62728"}
    texts = []
    for cond, g in df.groupby("condition", observed=False):
        ax.scatter(g["PC1"], g["PC2"], label=cond, color=colors[cond], s=45, alpha=0.85, zorder=2)
        for _, r in g.iterrows():
            texts.append(ax.text(r["PC1"], r["PC2"], r["patient"], fontsize=6))
    if adjust_text:
        adjust_text(texts, ax=ax)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    return df


def volcano_plot(res, out_png, title):
    vol = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    vol["nlp"] = -np.log10(vol["padj"].clip(lower=1e-300))
    up = (vol["padj"] < 0.05) & (vol["log2FoldChange"] >= 1)
    down = (vol["padj"] < 0.05) & (vol["log2FoldChange"] <= -1)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    ns = ~(up | down)
    ax.scatter(vol.loc[ns, "log2FoldChange"], vol.loc[ns, "nlp"], color="lightgray", alpha=0.5, s=12, label="Not significant")
    ax.scatter(vol.loc[up, "log2FoldChange"], vol.loc[up, "nlp"], color="#d62728", alpha=0.75, s=16, label="Up (padj<0.05, log2FC>=1)")
    ax.scatter(vol.loc[down, "log2FoldChange"], vol.loc[down, "nlp"], color="#1f77b4", alpha=0.75, s=16, label="Down (padj<0.05, log2FC<=-1)")
    ax.axvline(1, color="black", linestyle="--", linewidth=0.9)
    ax.axvline(-1, color="black", linestyle="--", linewidth=0.9)
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.9)

    # clip x-axis to the significant genes so a few near-zero-count genes don't stretch it
    xmax = np.ceil(vol.loc[up | down, "log2FoldChange"].abs().max()) + 1
    ax.set_xlim(-xmax, xmax)

    top = pd.concat([vol[up].sort_values("padj").head(8), vol[down].sort_values("padj").head(8)])
    texts = [ax.text(r["log2FoldChange"], r["nlp"], sym, fontsize=8) for sym, r in top.iterrows()]
    if adjust_text:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    ax.set_xlabel("log2 Fold Change (Tumour vs Normal)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# ---------------------------------------------------------------- main
def run_cohort(gse):
    cfg = COHORTS[gse]
    tag = gse.lower()
    print(f"\n===== {gse} =====")
    meta, counts_raw = load_inputs(cfg)
    meta.to_csv(DATA_DIR / f"{tag}_samples.csv", index=False)
    meta = meta.set_index("sample_id")
    samples = meta.index.tolist()

    counts = counts_by_symbol(counts_raw, samples, cfg["round_counts"])
    assert list(counts.index) == list(meta.index), "Sample alignment mismatch!"

    min_samples = max(3, int(np.ceil(MIN_SAMPLE_FRAC * counts.shape[0])))
    keep = (counts >= MIN_COUNT).sum(axis=0) >= min_samples
    print(f"Low-count filter (>= {MIN_COUNT} reads in >= {min_samples} samples): "
          f"kept {int(keep.sum())} of {counts.shape[1]} genes")
    counts = counts.loc[:, keep]

    meta["patient"] = meta["patient"].astype(str)
    meta["condition"] = pd.Categorical(meta["condition"], categories=["normal", "tumour"])

    dds = DeseqDataSet(counts=counts, metadata=meta, design="~ patient + condition", refit_cooks=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "tumour", "normal"])
    stat.summary()
    res = stat.results_df
    res.index.name = "Symbol"
    res.to_csv(DATA_DIR / f"{tag}_de_results.csv")

    markers = ["MMP1", "MMP3", "MMP10", "KRT17", "KRT4", "KRT13"]
    print("Marker sanity check:")
    print(res.loc[res.index.intersection(markers), ["log2FoldChange", "padj"]])

    # QC: VST -> PCA + tissue-composition marker score per sample
    if hasattr(dds, "vst"):
        dds.vst()
        vst_mat = np.asarray(dds.layers["vst_counts"])
    else:
        vst_mat = np.log2(np.asarray(dds.layers["normed_counts"]) + 1)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(vst_mat)
    qc = pca_plot(coords, pca.explained_variance_ratio_, meta,
                  FIG_DIR / f"{tag}_pca.png",
                  f"{gse} PCA (VST counts, {cfg['n_pairs']} matched pairs)")

    vst_df = pd.DataFrame(vst_mat, index=meta.index, columns=counts.columns)
    mk = [g for g in TISSUE_MARKERS if g in vst_df.columns]
    z = (vst_df[mk] - vst_df[mk].mean()) / vst_df[mk].std()
    qc["tissue_marker_score"] = z.mean(axis=1)  # high = muscle/fat-rich sample
    qc.reset_index().to_csv(DATA_DIR / f"{tag}_sample_qc.csv", index=False)
    print("Samples with highest muscle/fat marker score (check these against the PCA outliers):")
    print(qc.sort_values("tissue_marker_score", ascending=False).head(8)[["patient", "condition", "PC1", "PC2", "tissue_marker_score"]])

    volcano_plot(res, FIG_DIR / f"{tag}_volcano.png", f"{gse} Differential Expression Volcano Plot")
    print(f"{gse} complete.")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    for g in (COHORTS if arg == "all" else [arg]):
        run_cohort(g)
