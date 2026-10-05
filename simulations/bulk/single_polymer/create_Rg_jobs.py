from pymongo import MongoClient
import os
from dotenv import load_dotenv
import subprocess

load_dotenv("mongo.env")
MONGO_URI = os.getenv("MONGO_URI")
client = MongoClient(MONGO_URI)
db = client["PolymerPrediction"]
col = db["extendedPolymers"]

jobs = col.find()

for job in jobs:
        polymer_id = str(job["_id"])
        smiles = job["smiles"]

        name = job.get("name")

        if name:
            print(f"Submitting {name} - {polymer_id}")
        else:
            print(f"Submitting {polymer_id}")

        result = subprocess.run(
            [
                "sbatch",
                "Rg_job_template.sbatch",
                polymer_id,
                smiles
            ],)
