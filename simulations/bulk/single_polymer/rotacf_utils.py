import sys
import subprocess
import mdtraj as md
import numpy as np
from pathlib import Path
from bondsList import generate_backbone_ndx_from_folder
import json
import csv

# ==============================
# Config
# ==============================

BOX_SIZE = 5.0  # nm (cubic box assigned to all frames)

script_path = Path(__file__).resolve().parent
default_dummy_mdp = script_path / "dummy.mdp"


# ==============================
# Utilities
# ==============================

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


# ==============================
# Needed conversions
# ==============================

def prepare_box(ff_dir: Path):
    """
    Create conf_box.gro with a cubic box.
    This is needed because many OpenMM exports end up with no/invalid box for GROMACS analysis.
    """
    run_cmd([
        "gmx", "editconf",
        "-f", "conf.gro",
        "-o", "conf_box.gro",
        "-bt", "cubic",
        "-box", str(BOX_SIZE), str(BOX_SIZE), str(BOX_SIZE)
    ], cwd=ff_dir)


def generate_tpr(ff_dir: Path, dummy_mdp: Path):
    """
    Generate a minimal .tpr needed by gmx rotacf, starting from a .top
    No MD is run with GROMACS: nsteps=0, this is only for generating the .tpr
    """
    run_cmd([
        "gmx", "grompp",
        "-f", str(dummy_mdp),
        "-c", "conf_box.gro",
        "-p", "topol.top",
        "-o", "topol_box.tpr",
        "-maxwarn", "1"
    ], cwd=ff_dir)


def convert_dcd_to_xtc(md_dir: Path):
    """
    Convert OpenMM DCD -> XTC and manually assign a box to every frame.
    This prevents GROMACS from seeing a zero box (which causes it to fail).
    """
    dcd = md_dir / "polymer_vac_trajectory.dcd"
    pdb = md_dir / "polymer_vac_final.pdb"
    out = md_dir / "traj.xtc"

    if not dcd.exists():
        raise FileNotFoundError(f"DCD not found: {dcd}")
    if not pdb.exists():
        raise FileNotFoundError(f"PDB not found: {pdb}")

    print("Converting DCD -> XTC with box...")

    t = md.load_dcd(str(dcd), top=str(pdb))

    # assign cubic box to each frame (nm)
    t.unitcell_lengths = np.tile([BOX_SIZE, BOX_SIZE, BOX_SIZE], (t.n_frames, 1))
    t.unitcell_angles  = np.tile([90.0, 90.0, 90.0], (t.n_frames, 1))

    t.save_xtc(str(out))

    print("Saved:", out)


# ==============================
# Actual execution
# ==============================

def run_rotacf(polymer_dir: Path, dummy_mdp: Path):
    """
    Pipeline for one polymer folder:
    - prepare box in force_field/
    - build .tpr with grompp 
    - convert trajectory in MD/
    - run rotacf 
    """
    md_dir = polymer_dir / "MD"
    ff_dir = polymer_dir / "force_field"

    if not md_dir.exists():
        raise FileNotFoundError(f"Missing MD dir: {md_dir}")
    if not ff_dir.exists():
        raise FileNotFoundError(f"Missing force_field dir: {ff_dir}")

    bonds_ndx = md_dir / "bonds.ndx"
    if not bonds_ndx.exists():
        raise FileNotFoundError(f"Missing bonds.ndx: {bonds_ndx}")

    print("\nProcessing:", polymer_dir.name)

    prepare_box(ff_dir)
    generate_tpr(ff_dir, dummy_mdp)
    convert_dcd_to_xtc(md_dir)

    output = run_cmd([
        "gmx", "rotacf",
        "-d",
        "-s", str(ff_dir / "topol_box.tpr"),
        "-f", str(md_dir / "traj.xtc"),
        "-n", str(bonds_ndx),
        "-o", str(md_dir / "rotacf.xvg"),
        "-fitfn", "exp",
        "-beginfit", "0",
        "-endfit", "100",
        "-P", "2"
    ])
    tau = None

    """ for line in output.splitlines():
        print(line)
        if "tau" in line:
            tau = float(line.split("=")[1].strip().split()[0])
            break """
    line=output.splitlines()[-1]
    if not line:
        raise Exception("Invalid line")
    
    count = 0
    tau = None
    for l in line.split()[1:]:
        tau = float(l)
        if tau is None:
            raise Exception("Error: tau is none. Skipping ", polymer_dir.name)
        count+=1
    if count != 5:
        raise Exception("Substring aren't 5")

    print("Finished:", polymer_dir.name)
    print("Tau:",tau)
    print("-" * 60)
    return tau
