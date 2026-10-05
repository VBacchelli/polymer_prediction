import json
import os
import subprocess

from lib.gromacs_utils import run_gromacs
from lib.bulk_analyzer import (
    evaluate_tg,
    evaluate_thermal_expansion,
    measure_bulk_modulus,
    measure_young_poisson,
    measure_bulk_properties,
)


def create_supercell(input_gro, input_top, output_gro, output_top, n=2):
    print(f"      [Mechanics] Creating {n}x{n}x{n} supercell...")
    subprocess.run(
        ["gmx", "genconf", "-f", input_gro, "-o", output_gro, "-nbox", str(n), str(n), str(n)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        check=True,
        text=True,
    )

    with open(input_top, "r", encoding="utf-8") as fin, open(output_top, "w", encoding="utf-8") as fout:
        in_mols = False
        for line in fin:
            stripped = line.strip()

            if stripped.lower() == "[ molecules ]":
                in_mols = True
                fout.write(line)
                continue

            if in_mols and stripped and not stripped.startswith(";"):
                if stripped.startswith("["):
                    in_mols = False
                    fout.write(line)
                    continue

                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        fout.write(f"{parts[0]:<10} {int(parts[1]) * n**3}\n")
                    except ValueError:
                        fout.write(line)
                else:
                    fout.write(line)
            else:
                fout.write(line)


def _to_number(value):
    """Convert numeric strings to floats, while leaving non-numeric values unchanged."""
    try:
        return float(value)
    except Exception:
        return value


def _add_property(results, key, value, unit):
    """
    Store only the useful screening output in the main JSON.

    Detailed methods, source XVGs, fit slopes, fit ranges, etc. are written to
    polymer_properties_diagnostics.json instead of being duplicated inside results.json.
    """
    results[key] = {
        "value": _to_number(value),
        "unit": unit,
    }


def _as_float(result_entry):
    try:
        return float(result_entry.get("value"))
    except Exception:
        return None


def _build_quality_flags(properties, diagnostics):
    """Internal quality checks used to generate compact warnings."""
    flags = {}

    density = _as_float(properties.get("density", {}))
    if density is not None:
        flags["density_positive"] = density > 0.0
        flags["density_reasonable_polymer_range"] = 500.0 <= density <= 2500.0

    bulk = _as_float(properties.get("bulk_modulus", {}))
    if bulk is not None:
        flags["bulk_modulus_positive"] = bulk > 0.0

    young = _as_float(properties.get("young_modulus", {}))
    if young is not None:
        flags["young_modulus_positive"] = young > 0.0

    poisson = _as_float(properties.get("poisson_ratio", {}))
    if poisson is not None:
        flags["poisson_ratio_physical"] = -0.1 <= poisson <= 0.6

    tg = _as_float(properties.get("Tg", {}))
    if tg is not None:
        tg_diag = diagnostics.get("Tg", {})
        flags["tg_positive"] = tg > 0.0
        if "tg_inside_temperature_range" in tg_diag:
            flags["tg_inside_temperature_range"] = bool(tg_diag["tg_inside_temperature_range"])

    alpha = _as_float(properties.get("volumetric_thermal_expansion_coefficient", {}))
    if alpha is not None:
        flags["thermal_expansion_positive"] = alpha > 0.0

    bm_r2 = diagnostics.get("bulk_modulus", {}).get("fit_r2")
    if bm_r2 is not None:
        flags["bulk_modulus_fit_r2_ge_0p8"] = bm_r2 >= 0.8

    ym_r2 = diagnostics.get("young_modulus", {}).get("young_modulus_fit", {}).get("fit_r2")
    if ym_r2 is not None:
        flags["young_modulus_fit_r2_ge_0p8"] = ym_r2 >= 0.8

    pr_r2 = diagnostics.get("young_modulus", {}).get("poisson_ratio_fit", {}).get("fit_r2")
    if pr_r2 is not None:
        flags["poisson_ratio_fit_r2_ge_0p8"] = pr_r2 >= 0.8

    return flags


def _collect_warnings(quality_flags):
    """Convert internal quality flags into short user-facing warnings."""
    warnings = []

    if quality_flags.get("density_positive") is False:
        warnings.append("density is not positive")
    if quality_flags.get("density_reasonable_polymer_range") is False:
        warnings.append("density is outside the expected polymer range")

    if quality_flags.get("bulk_modulus_positive") is False:
        warnings.append("bulk_modulus is not positive")
    if quality_flags.get("young_modulus_positive") is False:
        warnings.append("young_modulus is not positive")
    if quality_flags.get("poisson_ratio_physical") is False:
        warnings.append("poisson_ratio is outside the expected physical range")

    if quality_flags.get("tg_positive") is False:
        warnings.append("Tg was not detected")
    if quality_flags.get("tg_inside_temperature_range") is False:
        warnings.append("Tg is outside the simulated temperature range")

    if quality_flags.get("thermal_expansion_positive") is False:
        warnings.append("thermal expansion coefficient is not positive")

    if quality_flags.get("bulk_modulus_fit_r2_ge_0p8") is False:
        warnings.append("bulk_modulus fit quality is low")
    if quality_flags.get("young_modulus_fit_r2_ge_0p8") is False:
        warnings.append("young_modulus fit quality is low")
    if quality_flags.get("poisson_ratio_fit_r2_ge_0p8") is False:
        warnings.append("poisson_ratio fit quality is low")

    return warnings


def _extract_energy_series(edr_file, work_dir, output_name, selection):
    """Extract one GROMACS energy series to an XVG file and return its path."""
    os.makedirs(work_dir, exist_ok=True)
    out_xvg = os.path.join(work_dir, output_name)
    subprocess.run(
        ["gmx", "energy", "-f", edr_file, "-o", out_xvg],
        input=f"{selection}\n",
        text=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return out_xvg


def _write_debug_json(sim_dir, diagnostics, quality_flags):
    """Write detailed diagnostics outside the main results JSON."""
    debug_path = os.path.join(sim_dir, "polymer_properties_diagnostics.json")
    debug_payload = {
        "quality_flags": quality_flags,
        "fit_diagnostics": diagnostics,
    }
    with open(debug_path, "w", encoding="utf-8") as f:
        json.dump(debug_payload, f, indent=4)
    return debug_path


def run(json_path: str, mode="full"):
    # added
    FAST_MODE = (mode == "fast")
    # ============================================================
    # 1. LOAD WORKFLOW STATE AND RESOLVE REQUIRED INPUTS
    # ============================================================
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    state = data["state"]
    state_paths = state.get("paths", {})

    root_dir = data["root_dir"]
    sim_dir = data["sim_dir"]

    os.makedirs(sim_dir, exist_ok=True)

    mdp_dir = os.path.join(root_dir, "inputs", "mdp")

    bulk_gro = state_paths.get("bulk_structure")
    bulk_top = state_paths.get("bulk_top")
    bulk_tpr = state_paths.get("bulk_tpr")
    bulk_xtc = state_paths.get("bulk_xtc")
    bulk_edr = state_paths.get("bulk_edr")

    if not bulk_gro or not os.path.exists(bulk_gro):
        raise FileNotFoundError("polymer_properties: Missing equilibrated bulk GRO from state['paths'].")
    if not bulk_top or not os.path.exists(bulk_top):
        raise FileNotFoundError("polymer_properties: Missing bulk topology from state['paths'].")
    if not bulk_tpr or not os.path.exists(bulk_tpr):
        raise FileNotFoundError("polymer_properties: Missing equilibrated bulk TPR from state['paths'].")
    if not bulk_xtc or not os.path.exists(bulk_xtc):
        raise FileNotFoundError("polymer_properties: Missing equilibrated bulk XTC from state['paths'].")
    if not bulk_edr or not os.path.exists(bulk_edr):
        raise FileNotFoundError("polymer_properties: Missing equilibrated bulk EDR from state['paths'].")

    # ============================================================
    # 2. RUN MODULE WORKFLOW
    # ============================================================
    properties = {}
    fit_diagnostics = {}
    plot_paths = {}

    print("-> Calculating Polymer Properties from equilibrated bulk...")

    general_dir = os.path.join(sim_dir, "general_props")
    vals, units, general_diag = measure_bulk_properties(
        general_dir,
        bulk_edr,
        bulk_xtc,
        bulk_tpr,
        start=10,
        return_diagnostics=True,
    )
    fit_diagnostics["general_bulk_properties"] = general_diag

    for k, v in vals.items():
        _add_property(properties, k, v, units.get(k, ""))

    # -----------------------------
    # Tg
    # -----------------------------
    tg_mdp = os.path.join(mdp_dir, "bulk_tg.mdp")
    tg_dir = os.path.join(sim_dir, "Tg")
    if (not FAST_MODE) and os.path.exists(tg_mdp):
        print("   -> Running Tg Simulation...")
        _, _, tg_edr, _ = run_gromacs(tg_mdp, bulk_gro, bulk_top, "tg_sim", tg_dir)

        vol_xvg = _extract_energy_series(tg_edr, tg_dir, "vol.xvg", "Volume")
        temp_xvg = _extract_energy_series(tg_edr, tg_dir, "temp.xvg", "Temperature")
        tg_plot = os.path.join(tg_dir, "tg.png")

        tg_val, tg_diag = evaluate_tg(vol_xvg, temp_xvg, tg_plot, return_diagnostics=True)
        fit_diagnostics["Tg"] = tg_diag
        plot_paths["Tg"] = tg_plot

        if tg_val > 0:
            _add_property(properties, "Tg", float(f"{tg_val:.1f}"), "K")

    # -----------------------------
    # Optional thermal expansion
    # -----------------------------
    thermal_mdp = os.path.join(mdp_dir, "bulk_thermal_exp.mdp")
    thermal_dir = os.path.join(sim_dir, "thermal_expansion")
    if (not FAST_MODE) and os.path.exists(thermal_mdp):
        print("   -> Running Thermal Expansion Simulation...")
        _, _, thermal_edr, _ = run_gromacs(
            thermal_mdp,
            bulk_gro,
            bulk_top,
            "thermal_expansion",
            thermal_dir,
        )

        vol_xvg = _extract_energy_series(thermal_edr, thermal_dir, "vol.xvg", "Volume")
        temp_xvg = _extract_energy_series(thermal_edr, thermal_dir, "temp.xvg", "Temperature")
        thermal_plot = os.path.join(thermal_dir, "thermal_expansion.png")

        th_vals, th_units, th_diag = evaluate_thermal_expansion(
            vol_xvg,
            temp_xvg,
            thermal_plot,
            return_diagnostics=True,
        )
        fit_diagnostics["thermal_expansion"] = th_diag
        plot_paths["thermal_expansion"] = thermal_plot

        for k, v in th_vals.items():
            _add_property(properties, k, v, th_units.get(k, ""))

    # -----------------------------
    # Mechanical supercell
    # -----------------------------
    mech_gro = os.path.join(sim_dir, "supercell.gro")
    mech_top = os.path.join(sim_dir, "supercell.top")
    if (not FAST_MODE):
        create_supercell(bulk_gro, bulk_top, mech_gro, mech_top, n=2)

    # -----------------------------
    # Bulk modulus
    # -----------------------------
    bm_mdp = os.path.join(mdp_dir, "bulk_modulus.mdp")
    bm_dir = os.path.join(sim_dir, "bulk_mod")
    if (not FAST_MODE) and os.path.exists(bm_mdp):
        print("   -> Running Bulk Modulus Workflow...")
        _, _, bm_edr, _ = run_gromacs(bm_mdp, mech_gro, mech_top, "bulk_modulus", bm_dir)

        bm_vals, bm_units, bm_plot, bm_diag = measure_bulk_modulus(
            bm_dir,
            bm_edr,
            start=10,
            return_diagnostics=True,
        )
        fit_diagnostics["bulk_modulus"] = bm_diag
        plot_paths["bulk_modulus"] = bm_plot

        for k, v in bm_vals.items():
            _add_property(properties, k, v, bm_units.get(k, ""))

        # Isothermal compressibility is intentionally not reported: for coating
        # screening, bulk_modulus already contains the same information in a
        # more directly useful form.

    # -----------------------------
    # Young's modulus and Poisson ratio
    # -----------------------------
    pre_ym_mdp = os.path.join(mdp_dir, "bulk_pre_young.mdp")
    ym_mdp = os.path.join(mdp_dir, "bulk_young.mdp")

    if not os.path.exists(pre_ym_mdp):
        pre_ym_mdp = os.path.join(mdp_dir, "pre_young.mdp")
    if not os.path.exists(ym_mdp):
        ym_mdp = os.path.join(mdp_dir, "young_modulus.mdp")

    ym_dir = os.path.join(sim_dir, "young_mod")
    if (not FAST_MODE) and os.path.exists(pre_ym_mdp) and os.path.exists(ym_mdp):
        print("   -> Running Young's Modulus Workflow...")

        print("      1. Pre-Equilibration...")
        pre_ym_gro, _, _, _ = run_gromacs(
            pre_ym_mdp,
            mech_gro,
            mech_top,
            "pre_young_modulus",
            ym_dir,
        )

        print("      2. Deformation...")
        _, _, ym_edr, _ = run_gromacs(
            ym_mdp,
            pre_ym_gro,
            mech_top,
            "young_modulus",
            ym_dir,
        )

        ym_vals, ym_units, young_plot, poisson_plot, ym_diag = measure_young_poisson(
            ym_dir,
            ym_edr,
            start=100,
            return_diagnostics=True,
        )
        fit_diagnostics["young_modulus"] = ym_diag

        if young_plot:
            plot_paths["young_modulus"] = young_plot if os.path.isabs(young_plot) else os.path.join(ym_dir, young_plot)
        if poisson_plot:
            plot_paths["poisson_ratio"] = poisson_plot if os.path.isabs(poisson_plot) else os.path.join(ym_dir, poisson_plot)

        for k, v in ym_vals.items():
            _add_property(properties, k, v, ym_units.get(k, ""))

    quality_flags = _build_quality_flags(properties, fit_diagnostics)
    warnings = _collect_warnings(quality_flags)
    debug_json = _write_debug_json(sim_dir, fit_diagnostics, quality_flags)

    print("\n   [Summary]")
    for k, v in properties.items():
        print(f"   - {k}: {v['value']} {v.get('unit', '')}")

    if warnings:
        print("\n   [Warnings]")
        for warning in warnings:
            print(f"   - {warning}")

    # ============================================================
    # 3. RETURN STANDARDIZED MODULE OUTPUT
    # ============================================================
    paths = {
        "polymer_properties_dir": sim_dir,
        "polymer_general_properties_dir": general_dir,
        "polymer_supercell_structure": mech_gro,
        "polymer_supercell_top": mech_top,
        "polymer_properties_diagnostics": debug_json,
    }

    # Keep plot paths also in state['paths'] so they can be found later.
    for plot_name, plot_path in plot_paths.items():
        paths[f"{plot_name}_plot"] = plot_path

    if os.path.exists(tg_dir):
        paths["polymer_tg_dir"] = tg_dir
    if os.path.exists(thermal_dir):
        paths["polymer_thermal_expansion_dir"] = thermal_dir
    if os.path.exists(bm_dir):
        paths["polymer_bulk_modulus_dir"] = bm_dir
    if os.path.exists(ym_dir):
        paths["polymer_young_modulus_dir"] = ym_dir

    # Main JSON remains compact: no duplicated diagnostics, no source XVGs,
    # no methods, and no fit internals.
    #
    # IMPORTANT FOR FINAL AGGREGATION:
    # WorkflowManager.build_final_results_json() extracts the top-level
    # "results" block from each module JSON. Therefore "results" should contain
    # only the physical screening quantities that should appear in final
    # results.json. Plots, warnings, and diagnostics are kept outside "results"
    # so they remain available in polymer_properties.json but are not mixed into
    # the final aggregated results.
    return {
        "paths": paths,
        "metadata": {
            "module_used": "polymer_properties",
            "source_model": "equilibrated_bulk",
            "warnings": warnings,
        },
        "artifacts": {
            "plots": plot_paths,
            "diagnostics_json": debug_json,
        },
        "results": properties,
    }