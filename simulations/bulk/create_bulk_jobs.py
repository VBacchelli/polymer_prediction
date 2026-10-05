import os
import subprocess
import argparse
from dotenv import load_dotenv
from pymongo import MongoClient


def is_bad_smiles(smiles: str) -> bool:
    cleaned = smiles.replace("[*]", "")
    return "()" in cleaned


# ----------------------
# CLI args
# ----------------------
parser = argparse.ArgumentParser()
parser.add_argument("--temp", type=str, default=None, help="Optional temperature")
args = parser.parse_args()
temp = args.temp


# ----------------------
# Mongo setup
# ----------------------
load_dotenv("mongo.env")

MONGO_URI = os.getenv("MONGO_URI")

client = MongoClient(MONGO_URI)
db = client["PolymerPrediction"]
col = db["biceranoPolymers"]


# ----------------------
# Loop over polymers
# ----------------------
for job in col.find():

    polymer_id = str(job["_id"])
    smiles = job["smiles"]
    name = job.get("name")

    if name:
        print(f"Submitting {name} - {polymer_id}")
    else:
        print(f"Submitting {polymer_id}")

    if is_bad_smiles(smiles):
        print(f"[SKIP] bad SMILES (creates ()): {name}")
        continue

    if temp and temp != "300":
        job_id = f"{temp}_{polymer_id}"
    else:
        job_id = polymer_id

    subprocess.run([
        "sbatch",
        "--job-name", f"poly_{polymer_id}",
        "--output", f"logs/{temp or '300'}_%x.out",
        "--error", f"logs/{temp or '300'}_%x.err",
        "bulk_job_template.sbatch",
        job_id,
        smiles
    ])
