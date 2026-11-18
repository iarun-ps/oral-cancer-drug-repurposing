"""Rebuild the exact features used by the binding models (Morgan fingerprint, radius 2, 2048 bits).
Usage:  import joblib; from featurize import smiles_to_fp
        rf = joblib.load("binding/models/rf_EGFR.joblib"); x = smiles_to_fp(["CCO"]); rf.predict(x)
Model files: binding/models/rf_<GENE>.joblib. Training data: binding/data/train_level2.csv
(columns: target, smiles, pchembl). Output unit: pChEMBL (9 = ~1 nM, 6 = ~1 uM)."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
def smiles_to_fp(smiles_list):
    return np.vstack([_gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles_list]).astype(np.uint8)
