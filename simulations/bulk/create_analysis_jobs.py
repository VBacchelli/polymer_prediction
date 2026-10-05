import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--temp", type=str, required=False, help="Temperature (optional)")
args = parser.parse_args()

temperature = args.temp

ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "BIO-SUSHY-workflows",
    "FOOD_v2"
)

RESULTS_DIR = os.path.join(ROOT, "results")

for folder_name in os.listdir(RESULTS_DIR):

    poly_dir = os.path.join(RESULTS_DIR, folder_name)

    # directories only
    if not os.path.isdir(poly_dir):
        continue

    if temperature:
        # expect format: <temp>_<polymer_id>
        parts = folder_name.split("_", 1)

        if len(parts) != 2:
            print(f"[SKIP] {folder_name} (invalid format)")
            continue

        folder_temp, polymer_id = parts

        if folder_temp != temperature:
            continue

        job_name_id = polymer_id
        sbatch_args = [folder_name, temperature]

    else:
        # only accept pure polymer_id (no temp prefix)
        if "_" in folder_name:
            continue

        polymer_id = folder_name
        job_name_id = polymer_id
        sbatch_args = [polymer_id]

    # check if result file exists
    res_file = os.path.join(poly_dir, "results.json")

    if not os.path.exists(res_file):
        print(f"[SKIP] {polymer_id} (no results.json)")
        continue

    analysis_file = os.path.join(poly_dir, "analysis.json")

    print(f"[SUBMIT] {polymer_id}")

    # log dir
    log_dir = os.path.join("logs", "analysis")
    os.makedirs(log_dir, exist_ok=True)

    # log naming
    if temperature:
        out_log = f"{log_dir}/{temperature}_%x_%j.out"
        err_log = f"{log_dir}/{temperature}_%x_%j.err"
    else:
        out_log = f"{log_dir}/%x_%j.out"
        err_log = f"{log_dir}/%x_%j.err"

    subprocess.run([
        "sbatch",
        "--job-name", f"analysis_{job_name_id}",
        "--output", out_log,
        "--error", err_log,
        "analysis_job_template.sbatch",
        *sbatch_args
    ])
