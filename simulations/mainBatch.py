import hashlib
import json
import os
import pandas as pd

from simulate import setup_biosushy, run_simulation


def polymer_id_from_row(row):
    if hasattr(row, "polymer_name") and pd.notna(row.polymer_name):
        return str(row.polymer_name)

    h = hashlib.sha1(row.smiles.encode("utf-8")).hexdigest()[:8]
    return f"poly_{h}"


def main():
    (
        WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    ) = setup_biosushy()

    df = pd.DataFrame([
        {
            "polymer_name": "polyethylene", # funzionante
            "smiles": "[*]CC[*]",
        },
        {
            "polymer_name": "polystyrene", # funzionante
            "smiles": "[*]CC(c1ccccc1)[*]",
        },
        {
            "polymer_name": "polyethylene terephthalate", #  Parameters have not been assigned to all angles
            "smiles": "[*]OC(=O)c1ccc(cc1)C(=O)OCC[*]",
        },
        {
            "polymer_name": "nylon 6", # Found no types for atom numbered 21 which is atomic number 6
            "smiles": "[*]NCCCCCC(=O)[*]",
        }, 
        {
            "polymer_name": "Poly(acrylic acid)", # Pre-condition Violation - bond already exists
            "smiles": "*CC(*)C(=O)O",
        }
    ])

    os.makedirs("simulation_results", exist_ok=True)
    os.chdir("simulation_results")

    results = []

    for row in df.itertuples(index=False):
        polymer_id = polymer_id_from_row(row)
        print(f"🚀 Running simulation for {polymer_id}")

        try:
            results_path = run_simulation(
                WorkflowManager,
                polymer_builder,
                ff_builder,
                md_engine,
                analysis,
                polymer_name=polymer_id,
                smiles=row.smiles,
                n_monomers=8,
                stoichiometry=0.3,
                temperature=300,
                steps=20000,
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
    results_df.to_csv("simulation_results.csv", index=False)

    print("✅ All simulations completed.")


if __name__ == "__main__":
    main()
