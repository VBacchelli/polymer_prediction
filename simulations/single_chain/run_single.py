import json
import os
from pathlib import Path
from simulate import setup_biosushy, run_simulation


def main():
    output_dir = Path(__file__).resolve().parent / "simulation_results"
    original_dir = Path.cwd()

    (
        WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    ) = setup_biosushy()

    print("🚀 Running simulation...")

    output_dir.mkdir(exist_ok=True)
    os.chdir(output_dir)
    try:
        results_path=run_simulation(
            WorkflowManager,
            polymer_builder,
            ff_builder,
            md_engine,
            analysis,
            polymer_name="polyethylene",
            smiles="[*]CC[*]",
            n_monomers=8,
            stoichiometry=0.3,
            temperature=300,
            steps=20000,
            timestep_ps = 0.001
        )
    finally:
        os.chdir(original_dir)

    if not os.path.exists(results_path):
        return None

    with open(results_path) as f:
        data = json.load(f)

    print("Rg mean=", data.get("radius_of_gyration", {}).get("mean"))


if __name__ == "__main__":
    main()
