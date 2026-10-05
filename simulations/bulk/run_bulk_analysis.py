import json
import os
import sys

# =============================
# PATH SETUP
# =============================

ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "BIO-SUSHY-workflows",
    "FOOD_v2"
)
sys.path.append(ROOT)

from modules import polymer_properties


# =============================
# MAIN PIPELINE
# =============================

def run(pid):

    res_file = os.path.join(ROOT, "results", pid, "results.json")

    if not os.path.exists(res_file):
        raise FileNotFoundError(f"Missing results.json for {pid}")

    with open(res_file) as f:
        data = json.load(f)

    try:
        paths = data["_modules"]["bulk_builder"]["paths"]
    except KeyError:
        raise RuntimeError(f"No bulk_builder data for {pid}")

    # =============================
    # RUN polymer_properties (FAST)
    # =============================

    sim_dir = os.path.join(ROOT, "simulations", pid, "analysis")

    input_data = {
        "state": {"paths": paths},
        "root_dir": ROOT,
        "sim_dir": sim_dir
    }

    tmp_json = f"/tmp/{pid}_analysis.json"

    with open(tmp_json, "w") as f:
        json.dump(input_data, f)

    print("[INFO] Running polymer_properties (fast)...")

    props_out = polymer_properties.run(tmp_json, mode="fast")

    # =============================
    # RESULTS
    # =============================

    results = props_out["results"]

    clean_results = {
        "density": results.get("density", {}).get("value"),
        "Rg": results.get("radius_of_gyration", {}).get("value"),
        "end_to_end": results.get("end_to_end_dist", {}).get("value"),
    }

    print("\n[RESULT]")
    for k, v in clean_results.items():
        print(f" - {k}: {v}")

    # =============================
    # SAVE RESULTS to analysis.json
    # =============================

    analysis_file = os.path.join(ROOT, "results", pid, "analysis.json")

    output = {
        "polymer_id": pid,
        "analysis": clean_results
    }

    tmp_out = analysis_file + ".tmp"

    with open(tmp_out, "w") as f:
        json.dump(output, f, indent=2)

    os.replace(tmp_out, analysis_file)
    print(f"[OK] Saved to {analysis_file}")


# =============================
# ENTRY POINT
# =============================

if __name__ == "__main__":
    run(sys.argv[1])
