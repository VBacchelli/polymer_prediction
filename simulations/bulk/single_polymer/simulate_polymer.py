#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import importlib
import traceback


def setup_biosushy(repo_name="BIO-SUSHY-tutorials"):
    repo_url = "https://github.com/daimoners/BIO-SUSHY-tutorials.git"

    base_dir = os.path.dirname(os.path.abspath(__file__))
    repo_path = os.path.join(base_dir, "..", repo_name)
    repo_path = os.path.abspath(repo_path)

    #if os.path.exists(repo_path):
    #   subprocess.run(["git", "-C", repo_path, "pull"], check=False)
    #else:
    #    parent_dir = os.path.dirname(repo_path)
    #    subprocess.run(["git", "clone", repo_url, repo_path], check=True)

    if repo_path not in sys.path:
        sys.path.insert(0, repo_path)

    import workflow
    import polymer_builder
    import ff_builder
    import md_engine
    import analysis

    importlib.reload(workflow)
    importlib.reload(polymer_builder)
    importlib.reload(ff_builder)
    importlib.reload(md_engine)
    importlib.reload(analysis)

    return (
        workflow.WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    )


def run_single_polymer(name, smiles):

    # --- parameters ---
    stoichiometry = 1.0
    n_monomers = 10
    temperature = 300
    steps = 200000
    timestep_ps = 0.001
    # ----------------------------------------
    results_root = os.path.join(os.getcwd(), "results")
    os.makedirs(results_root, exist_ok=True)
    os.chdir(results_root)

    WorkflowManager, polymer_builder, ff_builder, md_engine, analysis = setup_biosushy()

    wm = WorkflowManager(name)

    wm.add_path("results_path", os.path.join(wm.base_dir, "results.json"))

    wm.update_state("parameters", {
        "polymer_type": "homopolymer",
        "smiles_A": smiles,
        "stoichiometry": stoichiometry,
        "target_degree_polymerization": n_monomers,
        "temperature_k": temperature,
        "steps": steps,
        "timestep_ps": timestep_ps,
        "ensemble": "NVT_vacuum"
    })

    # --- BUILD POLYMER ---
    polymer_folder = os.path.join(wm.base_dir, "polymer")
    os.makedirs(polymer_folder, exist_ok=True)

    output_name_path = os.path.join(polymer_folder, name)

    raw_pdb = polymer_builder.build_polymer(
        smiles_A=smiles,
        polymer_type="homopolymer",
        n_monomers=n_monomers,
        name=output_name_path,
        ratio=stoichiometry,
    )

    wm.add_path("raw_pdb", raw_pdb)

    # --- FORCE FIELD ---
    force_field_dir = os.path.join(wm.base_dir, "force_field")

    try:
        gro_file, top_file = ff_builder.build_opls_system(
            raw_pdb,
            output_name=force_field_dir
        )

        wm.add_path("gro_file", gro_file)
        wm.add_path("top_file", top_file)

    except Exception as e:
        print(f"❌ FF error: {e}")
        traceback.print_exc()
        return

    # --- MD RUN ---
    md_dir = os.path.join(wm.base_dir, "MD")
    os.makedirs(md_dir, exist_ok=True)

    try:
        os.chdir(md_dir)

        traj, pdb = md_engine.run_vacuum_simulation(
            gro_file=wm.get_path("gro_file"),
            top_file=wm.get_path("top_file"),
            output_prefix="polymer_vac",
            temp_k=temperature,
            n_steps=steps
        )

        wm.add_path("trajectory", os.path.abspath(traj))
        wm.add_path("final_pdb", os.path.abspath(pdb))

    except Exception as e:
        print(f"❌ MD error: {e}")
        traceback.print_exc()
        return

    # --- ANALYSIS ---
    try:
        analysis.analyze_trajectory(wm)
    except Exception as e:
        print(f"⚠️ Analysis failed: {e}")

    print(f"✅ DONE: {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--smiles", required=True)

    args = parser.parse_args()

    run_single_polymer(args.name, args.smiles)
