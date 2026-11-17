"""Build the preliminary and consensus oral cancer signatures.

Run AFTER run_de.py has produced both *_de_results.csv files.

Outputs (signature/oral_cancer/results/):
    oral_signature.csv, oral_signature_up.txt, oral_signature_down.txt   (the hand-off files)
    preliminary/gse101547_preliminary_signature.csv   (cohort 1 only)
    consensus_details/oral_signature_annotated.csv    (all genes incl. excluded, with flags)
    consensus_details/replication_summary.txt         (numbers for the README)
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# repo root = the nearest parent folder that contains data/raw (so the script can live in any subfolder)
_here = Path(__file__).resolve()
ROOT = next((p for p in _here.parents if (p / "data" / "raw").is_dir()), _here.parent.parent)
RAW_DIR = ROOT / "data" / "raw" / "oral_cancer"
RES_DIR = ROOT / "signature" / "oral_cancer" / "results"
DATA_DIR = RES_DIR / "cohort_data"           # inputs from run_de.py
PRE_DIR = RES_DIR / "preliminary"            # cohort-1-only signature for Person B's rough run
DET_DIR = RES_DIR / "consensus_details"      # annotated table + replication numbers
for d in (RES_DIR, PRE_DIR, DET_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ----- settings (record these in the README) -----
PADJ_THRESH = 0.05
LFC_THRESH = 1.0
BASEMEAN_MIN = 10        # gene must have baseMean >= this in BOTH cohorts
TOP_N = 150              # per direction, for the .txt files (agree N with Person B)
# Genes of these classes are left out of oral_signature.csv (kept in the annotated file).
# Set to set() to keep everything.
EXCLUDE_CLASSES = {"immunoglobulin", "muscle_adipose", "noncoding"}

# ----- gene classes (heuristic, edit freely) -----
IG_RE = re.compile(r"^IG[HKL][VDJ]\d|^IGH[GAMED]\d*$|^IGKC$|^IGLC\d$")
NC_RE = re.compile(r"^(LINC\d|LOC\d|MIR\d|SNOR[AD]|RNU\d|RPSAP|HCG\d)|-AS\d*$|-IT\d+$|HG$|^(LNCAROD|LASTR|EPIC1)$")
MUSCLE_ADIPOSE = {
    # skeletal / cardiac muscle
    "MYH1", "MYH2", "MYH7", "MYH8", "MYL1", "MYL2", "MYL3", "MYLPF", "MYBPC1", "MYBPC2", "TNNC1", "TNNC2",
    "TNNI1", "TNNI2", "TNNT1", "TNNT3", "ACTA1", "DES", "MB", "CKM", "CASQ1", "CSRP3", "LDB3", "NRAP",
    "XIRP2", "ABRA", "SMPX", "PYGM", "TRDN", "LMOD2", "LMOD3", "MYOZ1", "KLHL40", "KLHL41", "SMYD1",
    "COX6A2", "CAV3", "CACNA1S", "SGCG", "PPP1R3A", "ATP1A2", "TTN", "NEB", "RYR1", "ATP2A1",
    "MIR133A1HG", "MYOT", "ACTN2", "ACTN3", "CASQ2", "HSPB6",
    # adipose
    "ADIPOQ", "PLIN1", "PLIN4", "CIDEC", "FABP4", "LEP", "LPL", "CFD", "GPD1", "THRSP", "PCK1", "RBP4",
    "ACADL", "FBP2", "GYS2",
}


def gene_class(sym):
    if IG_RE.search(sym):
        return "immunoglobulin"
    if sym in MUSCLE_ADIPOSE:
        return "muscle_adipose"
    if NC_RE.search(sym):
        return "noncoding"
    return "other"


def load_de(path):
    df = pd.read_csv(path, index_col=0)
    df.index.name = "gene_symbol"
    df = df.dropna(subset=["log2FoldChange", "padj", "stat"])
    return df[df["baseMean"] >= BASEMEAN_MIN]


def entrez_map():
    a = pd.read_csv(RAW_DIR / "Human.GRCh38.p13.annot.tsv.gz", sep="\t",
                    usecols=["Symbol", "GeneID"], low_memory=False).dropna()
    a["entrez_id"] = pd.to_numeric(a["GeneID"], errors="coerce")
    a = a.dropna(subset=["entrez_id"]).drop_duplicates(subset="Symbol")
    a["entrez_id"] = a["entrez_id"].astype(int)
    return a[["Symbol", "entrez_id"]].rename(columns={"Symbol": "gene_symbol"})


def finalise(df, entrez, label):
    """Filter classes, attach entrez IDs, rank within direction by |z| (1 = strongest)."""
    n0 = len(df)
    df = df[~df["gene_class"].isin(EXCLUDE_CLASSES)]
    n1 = len(df)
    df = df.merge(entrez, on="gene_symbol", how="inner")
    print(f"[{label}] {n0} genes -> {n1} after class filter -> {len(df)} with an Entrez ID")
    df["rank"] = df.groupby("direction")["z"].rank(ascending=False, method="first").astype(int)
    df = df.sort_values(["direction", "rank"], ascending=[False, True])  # up first
    return df


CONTRACT_COLS = ["gene_symbol", "entrez_id", "log2fc", "padj", "direction", "rank", "n_cohorts"]

entrez = entrez_map()
r1 = load_de(DATA_DIR / "gse101547_de_results.csv")
r2 = load_de(DATA_DIR / "gse213862_de_results.csv")

# =========================================================
# 1. Preliminary signature (GSE101547 only) for Person B
# =========================================================
p = r1[(r1["padj"] < PADJ_THRESH) & (r1["log2FoldChange"].abs() >= LFC_THRESH)].copy()
p = p.reset_index()
p["log2fc"] = p["log2FoldChange"]
p["direction"] = np.where(p["log2fc"] > 0, "up", "down")
p["z"] = p["stat"].abs()
p["n_cohorts"] = 1
p["gene_class"] = p["gene_symbol"].map(gene_class)
p = finalise(p, entrez, "preliminary")
p = p[p["rank"] <= TOP_N]
p[CONTRACT_COLS].to_csv(PRE_DIR / "gse101547_preliminary_signature.csv", index=False)

# =========================================================
# 2. Consensus signature
# =========================================================
cols = ["baseMean", "log2FoldChange", "stat", "pvalue", "padj"]
j = r1[cols].join(r2[cols], lsuffix="_a", rsuffix="_b", how="inner")  # a = GSE101547, b = GSE213862
j["gene_class"] = [gene_class(s) for s in j.index]

sig_a = (j["padj_a"] < PADJ_THRESH) & (j["log2FoldChange_a"].abs() >= LFC_THRESH)
sig_b = (j["padj_b"] < PADJ_THRESH) & (j["log2FoldChange_b"].abs() >= LFC_THRESH)
same = np.sign(j["log2FoldChange_a"]) == np.sign(j["log2FoldChange_b"])
both = j[sig_a & sig_b]
cons = j[sig_a & sig_b & same].copy().reset_index()

cons["log2fc"] = (cons["log2FoldChange_a"] + cons["log2FoldChange_b"]) / 2
cons["padj"] = np.maximum(cons["padj_a"], cons["padj_b"])        # weaker of the two cohorts
cons["z"] = (cons["stat_a"] + cons["stat_b"]).abs() / np.sqrt(2)  # combined evidence used for ranking
cons["direction"] = np.where(cons["log2fc"] > 0, "up", "down")
cons["n_cohorts"] = 2

# annotated file: every consensus gene with its class and whether it was excluded
ann = cons.copy()
ann["excluded"] = ann["gene_class"].isin(EXCLUDE_CLASSES)
ann.sort_values("z", ascending=False).to_csv(DET_DIR / "oral_signature_annotated.csv", index=False)

final = finalise(cons, entrez, "consensus")
final[CONTRACT_COLS].to_csv(RES_DIR / "oral_signature.csv", index=False)
up = final[(final["direction"] == "up") & (final["rank"] <= TOP_N)]["gene_symbol"]
dn = final[(final["direction"] == "down") & (final["rank"] <= TOP_N)]["gene_symbol"]
up.to_csv(RES_DIR / "oral_signature_up.txt", index=False, header=False)
dn.to_csv(RES_DIR / "oral_signature_down.txt", index=False, header=False)
top_final = set(up) | set(dn)

# =========================================================
# 3. Replication numbers + leave-one-cohort-out (LOCO)
# =========================================================
lines = []
rho_all, _ = spearmanr(j["log2FoldChange_a"], j["log2FoldChange_b"])
rho_cons, _ = spearmanr(cons["log2FoldChange_a"], cons["log2FoldChange_b"])
lines += [
    "=== REPLICATION SUMMARY ===",
    f"Thresholds: padj<{PADJ_THRESH}, |log2FC|>={LFC_THRESH}, baseMean>={BASEMEAN_MIN} in both cohorts",
    f"Shared tested genes:               {len(j)}",
    f"Significant in GSE101547:          {int(sig_a.sum())}",
    f"Significant in GSE213862:          {int(sig_b.sum())}",
    f"Significant in both:               {len(both)}",
    f"Same direction in both:            {len(cons)} / {len(both)} ({len(cons) / len(both) * 100:.2f}%)",
    f"Spearman r, all shared genes:      {rho_all:.3f}",
    f"Spearman r, consensus genes:       {rho_cons:.3f}",
    f"Consensus genes (all classes):     {len(cons)}",
    f"Genes in oral_signature.csv:       {len(final)} ({int((final['direction'] == 'up').sum())} up, {int((final['direction'] == 'down').sum())} down)",
    f"Excluded by class: " + ", ".join(f"{c}={int((cons['gene_class'] == c).sum())}" for c in sorted(EXCLUDE_CLASSES)),
    "",
    f"=== LEAVE-ONE-COHORT-OUT (signature built from one cohort, tested in the other; top {TOP_N} per direction) ===",
]

jf = j[~j["gene_class"].isin(EXCLUDE_CLASSES)]
for src, dst, name in [("a", "b", "GSE101547 -> GSE213862"), ("b", "a", "GSE213862 -> GSE101547")]:
    s = jf[(jf[f"padj_{src}"] < PADJ_THRESH) & (jf[f"log2FoldChange_{src}"].abs() >= LFC_THRESH)].copy()
    s["dir"] = np.sign(s[f"log2FoldChange_{src}"])
    s["zabs"] = s[f"stat_{src}"].abs()
    pick = pd.concat([s[s["dir"] > 0].nlargest(TOP_N, "zabs"), s[s["dir"] < 0].nlargest(TOP_N, "zabs")])
    agree = np.sign(pick[f"log2FoldChange_{dst}"]) == pick["dir"]
    sig_dst = agree & (pick[f"padj_{dst}"] < PADJ_THRESH)
    rho, _ = spearmanr(pick[f"log2FoldChange_{src}"], pick[f"log2FoldChange_{dst}"])
    lines += [
        f"{name}: {len(pick)} genes | same direction in held-out cohort: {agree.mean() * 100:.1f}% | "
        f"significant + same direction: {sig_dst.mean() * 100:.1f}% | Spearman r of log2FC: {rho:.3f} | "
        f"overlap with final top list: {len(set(pick.index) & top_final)}"
    ]

text = "\n".join(lines)
print(text)
(DET_DIR / "replication_summary.txt").write_text(text + "\n")
print(f"\nSaved outputs to {RES_DIR}")

# =========================================================
# 4. Self-check against the data contract (Section 8)
# =========================================================
chk = pd.read_csv(RES_DIR / "oral_signature.csv")
assert list(chk.columns) == CONTRACT_COLS, "column names/order differ from the contract"
assert chk.notna().all().all(), "missing values"
assert chk["gene_symbol"].str.strip().ne("").all() and not chk["gene_symbol"].duplicated().any()
assert pd.api.types.is_integer_dtype(chk["entrez_id"]) and pd.api.types.is_integer_dtype(chk["rank"])
assert set(chk["direction"]) <= {"up", "down"}
assert ((chk["log2fc"] > 0) == (chk["direction"] == "up")).all(), "direction does not match log2fc sign"
assert chk["n_cohorts"].isin([1, 2]).all()
for d, g in chk.groupby("direction"):
    assert sorted(g["rank"]) == list(range(1, len(g) + 1)), f"ranks for {d} are not 1..n"
assert len(pd.read_csv(RES_DIR / "oral_signature_up.txt", header=None)) <= TOP_N
print("Contract self-check passed (still ask Person D to run the official validator).")
