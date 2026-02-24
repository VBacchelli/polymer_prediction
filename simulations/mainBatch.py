import hashlib
import json
import os
import pandas as pd
import re
import traceback

from simulate import setup_biosushy, run_simulation

def main():
    (
        WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    ) = setup_biosushy()

    from pymongo import MongoClient
    from dotenv import load_dotenv

    load_dotenv()
    MONGO_URI = os.getenv("MONGO_URI")

    client = MongoClient(MONGO_URI)
    db=client["PolymerPrediction"]
    collection=db["biceranoPolymers"]

    docs = list(collection.find({}))
    df = pd.DataFrame(docs)
    os.makedirs("simulation_resultsHf", exist_ok=True)
    os.chdir("simulation_resultsHf")

    results = []

    for _, row in df.iterrows():
        polymer_id = str(row._id) # mantengo l'id che ho sul dataset per poterli riallineare
        print(f"🚀 Running simulation for {polymer_id}")
        smiles=row.smiles
        try:
            results_path = run_simulation(
                WorkflowManager,
                polymer_builder,
                ff_builder,
                md_engine,
                analysis,
                polymer_name=polymer_id,
                smiles=smiles,
                n_monomers=3,
                stoichiometry=0.3,
                temperature=300,
                steps=40000,
                timestep_ps=0.001,
            )

            if not results_path or not os.path.exists(results_path):
                continue

            with open(results_path) as f:
                data = json.load(f)

            results.append({
                "polymer_id": polymer_id,
                "smiles": row.smiles,
                "Rg_mean": data.get("radius_of_gyration", {}).get("mean"),
            })

        except Exception as e:
            print(f"❌ Simulation failed for {polymer_id}: {e}")

    results_df = pd.DataFrame(results)
    print(results_df)
    results_df.to_csv("results_biceranoHf.csv", index=False)

    print("✅ All simulations completed.")


if __name__ == "__main__":
    main()
