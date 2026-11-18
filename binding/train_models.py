"""Rebuild the binding models (same settings and seed as used for binding_scores.csv).
Run from the repo root:  python binding/train_models.py     (about 10 minutes)
Writes binding/models/rf_<GENE>.joblib  (not stored in Git, files are large)."""
import os, sys, joblib
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
sys.path.insert(0, "binding")
from featurize import smiles_to_fp

train = pd.read_csv("binding/data/train_level2.csv")
os.makedirs("binding/models", exist_ok=True)
for gene, g in train.groupby("target"):
    X = smiles_to_fp(g["smiles"].tolist()); y = g["pchembl"].values
    rf = RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=42).fit(X, y)
    joblib.dump(rf, f"binding/models/rf_{gene}.joblib", compress=3)
    print(gene, "saved,", len(g), "training molecules")
