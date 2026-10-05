import os
import json
import csv
import argparse

# =============================
# ARGUMENTS
# =============================

parser = argparse.ArgumentParser(description="Summarize analysis.json files")
parser.add_argument("--temp", help="Temperature prefix (e.g. 300)")
args = parser.parse_args()

TEMP = args.temp
USE_PREFIX = TEMP is not None

# =============================
# PATHS
# =============================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

ROOT = os.path.join(
    BASE_DIR,
    "BIO-SUSHY-workflows",
    "FOOD_v2"
)

RESULTS_DIR = os.path.join(ROOT, "results")

# Output file name adapts
suffix = TEMP if TEMP else "base"
OUTPUT_CSV = os.path.join(BASE_DIR, f"analysis_summary_{suffix}.csv")

rows = []
all_keys = set()

# =============================
# COLLECT DATA
# =============================

for dirname in os.listdir(RESULTS_DIR):
    full_dir = os.path.join(RESULTS_DIR, dirname)

    if not os.path.isdir(full_dir):
        continue

    # -------------------------
    # DIRECTORY FILTERING
    # -------------------------
    if USE_PREFIX:
        # Expect <TEMP>_<polymer_id>
        if not dirname.startswith(f"{TEMP}_"):
            continue
        polymer_id = dirname.split("_", 1)[1]
    else:
        # Only accept dirs WITHOUT temp prefix
        if "_" in dirname:
            continue
        polymer_id = dirname

    # -------------------------
    # LOAD ANALYSIS
    # -------------------------
    analysis_file = os.path.join(full_dir, "analysis.json")

    if not os.path.exists(analysis_file):
        print(f"[SKIP] {dirname} (no analysis.json)")
        continue

    try:
        with open(analysis_file) as f:
            data = json.load(f)
    except json.JSONDecodeError:
        print(f"[WARN] {dirname} (invalid JSON)")
        continue

    analysis = data.get("analysis", {})

    row = {"polymer_id": polymer_id}

    # -------------------------
    # COLUMN NAMING
    # -------------------------
    for k, v in analysis.items():
        if USE_PREFIX:
            new_key = f"{k}_md_{TEMP}"
        else:
            new_key = f"{k}_md"

        row[new_key] = v
        all_keys.add(new_key)

    rows.append(row)

# =============================
# WRITE CSV
# =============================

fieldnames = ["polymer_id"] + sorted(all_keys)

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    for row in rows:
        writer.writerow(row)

# =============================
# DONE
# =============================

print(f"[OK] Saved CSV to {OUTPUT_CSV}")
print(f"[INFO] {len(rows)} polymers collected")
