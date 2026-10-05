import os
import json
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--smiles", required=True)
    args = parser.parse_args()

    REPO_PATH = "./BIO-SUSHY-workflows/FOOD_v2"

    input_dir = os.path.join(REPO_PATH, "inputs", "polymers")
    os.makedirs(input_dir, exist_ok=True)

 
    polymer_dir = os.path.join(REPO_PATH, "inputs", "polymers", args.name)   
    os.makedirs(polymer_dir, exist_ok=True)

    json_path = os.path.join(polymer_dir, "input.json")

    # write json if not exists
    if not os.path.exists(json_path):
        print(f"[INFO] Creating input JSON: {json_path}")

        data = {
            "name": args.name,
            "smiles_A": args.smiles,
            "polymer_type": "homopolymer"
        }

        with open(json_path, "w") as f:
            json.dump(data, f, indent=4)

    else:
        print(f"[INFO] Using existing JSON: {json_path}")


    # Run master runner
    subprocess.run(
        ["python", "master_runner.py", args.name],
        cwd=REPO_PATH,
        check=True
    )



if __name__ == "__main__":
    main()

