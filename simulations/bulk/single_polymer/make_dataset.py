import os
import json
import numpy as np
import pandas as pd

base_dir = "results"



def extract_features(angles, probs):
    trans = sum(p for a, p in zip(angles, probs) if abs(a) > 150)
    gauche = sum(p for a, p in zip(angles, probs) if 30 < abs(a) < 90)
    cis = sum(p for a, p in zip(angles, probs) if abs(a) < 30)

    mean = sum(a * p for a, p in zip(angles, probs))
    var = sum(p * (a - mean)**2 for a, p in zip(angles, probs))

    angles_rad = np.radians(angles)

    C = sum(p * np.cos(a) for a, p in zip(angles_rad, probs))
    S = sum(p * np.sin(a) for a, p in zip(angles_rad, probs))

    R = np.sqrt(C**2 + S**2)

    circular_var = 1 - R

    return trans, gauche, cis, mean, var, circular_var


rows = []

for polymer_id in os.listdir(base_dir):
    polymer_path = os.path.join(base_dir, polymer_id)

    if not os.path.isdir(polymer_path):
        continue

    dih_file = os.path.join(polymer_path, "MD", "dihedral.json")
    rot_file = os.path.join(polymer_path, "rotacf.json")
    res_file = os.path.join(polymer_path, "results.json")

    if not os.path.exists(dih_file):
        continue

    # ----------- READ DIHEDRAL -----------
    with open(dih_file) as f:
        dih = json.load(f)

    angles = dih["angle"]
    probs = dih["probability"]

    trans, gauche, cis, mean, var, circular_var = extract_features(angles, probs)

    # ----------- READ TAU -----------
    tau = None
    if os.path.exists(rot_file):
        with open(rot_file) as f:
            tau = json.load(f).get("tau")

    # ----------- READ RG -----------
    rg = None
    if os.path.exists(res_file):
        with open(res_file) as f:
            rg = json.load(f).get("radius_of_gyration", {}).get("mean")

    rows.append({
        "polymer_id": polymer_id,
        "Rg": rg,
        "tau": tau,
        "trans": trans,
        "gauche": gauche,
        "cis": cis,
        "mean_angle": mean,
        "var_dihedral": var,
        "circular_var": circular_var
    })

df = pd.DataFrame(rows)
df.to_csv("final_dataset.csv", index=False)

print("✅ Dataset creato:", df.shape)
