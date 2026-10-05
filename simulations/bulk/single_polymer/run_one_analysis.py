from pathlib import Path
import sys
import json
import subprocess

from rotacf_utils import run_rotacf, generate_backbone_ndx_from_folder


def run_cmd(cmd, cwd=None):
    print("Running:", " ".join(str(c) for c in cmd))

    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("Command failed")

    return result.stdout


def xvg_to_json(xvg_path, json_path):
    angles = []
    probs = []

    with open(xvg_path) as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            a, p = map(float, line.split())
            angles.append(a)
            probs.append(p)

    with open(json_path, "w") as f:
        json.dump({
            "angle": angles,
            "probability": probs
        }, f, indent=2)


if __name__ == "__main__":
    polymer_dir = Path(sys.argv[1]).resolve()

    if len(sys.argv) > 2:
        dummy_mdp = Path(sys.argv[2]).resolve()
    else:
        dummy_mdp = Path.cwd() / "dummy.mdp"

    print("Processing:", polymer_dir)

    md_dir = polymer_dir / "MD"

    # ======================================================
    # 1. BACKBONE + INDEX
    # ======================================================
    generate_backbone_ndx_from_folder(polymer_dir)

    # ======================================================
    # 2. ROTACF → TAU
    # ======================================================
    tau = run_rotacf(polymer_dir, dummy_mdp)

    print("✅ Tau:", tau)

    rotacf_file = polymer_dir / "rotacf.json"

    with open(rotacf_file, "w") as f:
        json.dump({"tau": tau}, f, indent=2)

    # ======================================================
    # 3. DIHEDRAL DISTRIBUTION
    # ======================================================
    print("Running dihedral distribution...")

    run_cmd([
        "gmx", "angle",
        "-f", str(md_dir / "traj.xtc"),
        "-n", str(md_dir / "dihedrals.ndx"),
        "-type", "dihedral",
        "-od", str(md_dir / "dihedral_dist.xvg")
    ])

    # ======================================================
    # 4. SAVE JSON DIHEDRAL
    # ======================================================
    xvg_file = md_dir / "dihedral_dist.xvg"
    json_file = md_dir / "dihedral.json"

    xvg_to_json(xvg_file, json_file)

    print("✅ Saved dihedral JSON")

    print("✅ DONE:", polymer_dir.name)

