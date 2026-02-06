import os
import sys
import subprocess
import importlib
import tempfile


def setup_biosushy(repo_name="BIO-SUSHY-tutorials"):
    """
    Clone/update BIO-SUSHY and return loaded modules.
    """

    repo_url = "https://github.com/daimoners/BIO-SUSHY-tutorials.git"

    if os.path.exists(repo_name):
        print(f"🔄 Updating {repo_name}...")
        subprocess.run(["git", "-C", repo_name, "pull"], check=False)
    else:
        print(f"⬇️ Cloning {repo_name}...")
        subprocess.run(["git", "clone", repo_url], check=True)

    repo_path = os.path.abspath(repo_name)
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

    print("✅ BIO-SUSHY loaded")

    return (
        workflow.WorkflowManager,
        polymer_builder,
        ff_builder,
        md_engine,
        analysis,
    )


def run_simulation(
    WorkflowManager,
    polymer_builder,
    ff_builder,
    md_engine,
    analysis,
    polymer_name,
    smiles,
    stoichiometry,
    n_monomers,
    temperature,
    steps,
    timestep_ps
):
    """
    Run a single simulation.
    """
    wm = WorkflowManager(polymer_name)
    
    results_path = os.path.join(wm.base_dir, "results.json")
    wm.add_path("results_path", results_path)

    print("📂 Workflow base dir:", wm.base_dir)
    
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

    print(f"✅ Settings loaded for: {polymer_name}")
    print(f"   📂 Master Directory: {wm.base_dir}")
    print(f"   📝 State Initialized: {wm.state_file}")

    base_dir = wm.base_dir
    polymer_folder = os.path.join(base_dir, "polymer")
    os.makedirs(polymer_folder, exist_ok=True)

    output_name_path = os.path.join(polymer_folder, polymer_name)

    raw_pdb = polymer_builder.build_polymer(
        smiles_A=smiles,
        polymer_type="homopolymer",
        n_monomers=n_monomers,
        name=output_name_path,
        ratio=stoichiometry,
    )


    wm.add_path("raw_pdb", raw_pdb)
    wm.update_state("build_status", {"raw_build": "success"})

    polymer_pdb = wm.get_path("raw_pdb")
    force_field_dir = os.path.join(wm.base_dir, "force_field")

    try:
        if not polymer_pdb or not os.path.exists(polymer_pdb):
            print(f"❌ Error: Raw PDB not found in state. Run Step 5.")
        else:
            print(f"⚙️ Generating OPLS parameters...")

            # Run Builder
            gro_file, top_file = ff_builder.build_opls_system(polymer_pdb, output_name=force_field_dir)

            # --- WORKFLOW: UPDATE STATE ---
            wm.add_path("gro_file", gro_file)
            wm.add_path("top_file", top_file)
            wm.update_state("build_status", {"parameterization": "success"})
            # ------------------------------

            print(f"\n✅ Force Field Generated Successfully!")
            print(f"   📜 Topology: {os.path.basename(top_file)}")

    except Exception as e:
        wm.update_state("build_status", {"parameterization": "failed", "error": str(e)})
        print(f"\n❌ Error during parameterization: {e}")

    gro_path = wm.get_path("gro_file")
    top_path = wm.get_path("top_file")

    # Outputs
    md_dir = os.path.join(wm.base_dir, "MD")
    output_prefix = "polymer_vac"

    try:
        if not gro_path: raise FileNotFoundError("GRO file missing in state.")

        os.makedirs(md_dir, exist_ok=True)
        print(f"🚀 Starting Simulation in: {md_dir}")

        original_dir = os.getcwd()
        os.chdir(md_dir)

        try:
            # Run Simulation
            traj, pdb = md_engine.run_vacuum_simulation(
                gro_file=gro_path,
                top_file=top_path,
                output_prefix=output_prefix,
                temp_k=temperature,
                n_steps=steps
                # Assuming md_engine uses the standard 1fs timestep internally
            )

            # --- WORKFLOW: SAVE RESULTS ---
            wm.add_path("trajectory", os.path.abspath(traj))
            wm.add_path("final_pdb", os.path.abspath(pdb))

            # Save MD Details to results.json
            total_time_ps = steps * timestep_ps

            wm.save_result("temperature", temperature, units="K")
            wm.save_result("simulation_steps", steps, units="steps")
            wm.save_result("timestep", timestep_ps, units="ps")
            wm.save_result("simulated_time", total_time_ps, units="ps")
            # ------------------------------

        finally:
            os.chdir(original_dir)

        print(f"✅ Simulation Complete.")

    except Exception as e:
        wm.update_state("simulation_status", {"status": "failed", "error": str(e)})
        print(f"\n❌ Error: {e}")

    cwd = os.getcwd()
    os.chdir(md_dir)

    analysis.analyze_trajectory(wm)

    return wm.get_path("results_path")