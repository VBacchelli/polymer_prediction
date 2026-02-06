import json
import os
from simulate import setup_biosushy, run_simulation


def main():
    (
        WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    ) = setup_biosushy()

    os.makedirs("simulations_test", exist_ok=True)

    print("🚀 Running simulation...")

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

    if not os.path.exists(results_path):
        return None

    with open(results_path) as f:
        data = json.load(f)

    print("Rg mean=", data.get("radius_of_gyration", {}).get("mean"))


if __name__ == "__main__":
    main()
