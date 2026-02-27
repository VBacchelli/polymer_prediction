import sys
import subprocess
import mdtraj as md
import numpy as np
from pathlib import Path
from bondsList import generate_backbone_ndx_from_folder

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
    """
    Run a command and raise an error if it fails.
    """
    print("Running:", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


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

    run_cmd([
        "gmx", "rotacf",
        "-d",  # IMPORTANT: treat index as bond vectors (pairs), not atom triplets
        "-s", str(ff_dir / "topol_box.tpr"),
        "-f", str(md_dir / "traj.xtc"),
        "-n", str(bonds_ndx),
        "-o", str(md_dir / "rotacf.xvg")
    ])

    print("Finished:", polymer_dir.name)
    print("-" * 60)


# ==============================
# Batch logic
# ==============================

def run_batch(root_dir: Path, dummy_mdp: Path):
    """
    Scan only first-level subfolders of root_dir.
    Run if both MD/ and force_field/ are present.
    Runs the module that generates the bonds.ndx index file with the bonds list first.
    """
    successes = []
    failures = {}

    for polymer_dir in root_dir.iterdir():
        if not polymer_dir.is_dir():
            continue

        md_dir = polymer_dir / "MD"
        ff_dir = polymer_dir / "force_field"

        if md_dir.is_dir():
            try:
                generate_backbone_ndx_from_folder(md_dir)
            except Exception as e:
                print(f"Skipped {polymer_dir.name}: {e}")

        if md_dir.exists() and ff_dir.exists():
            try:
                run_rotacf(polymer_dir, dummy_mdp)
                successes.append(polymer_dir.name)
            except Exception as e:
                failures[polymer_dir.name] = str(e)
                print("Error in:", polymer_dir.name)
                print(e)
                print("=" * 60)

    # ==============================
    # Final execution summary
    # ==============================

    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)

    print(f"\nTotal processed: {len(successes) + len(failures)}")
    print(f"Success: {len(successes)}")
    print(f"Failed: {len(failures)}")


if __name__ == "__main__":

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage:")
        print("  python run_rotacf_pipeline.py /path/to/root [optional/path/to/dummy.mdp]")
        sys.exit(1)

    ROOT_DIR = Path(sys.argv[1]).resolve()

    if len(sys.argv) == 3:
        dummy_mdp = Path(sys.argv[2]).resolve()
    else:
        dummy_mdp = default_dummy_mdp

    if not ROOT_DIR.exists():
        print("Invalid root directory:", ROOT_DIR)
        sys.exit(1)

    if not dummy_mdp.exists():
        print("dummy.mdp not found:", dummy_mdp)
        sys.exit(1)

    print("Starting batch from root:", ROOT_DIR)
    print("Using dummy.mdp:", dummy_mdp)

    run_batch(ROOT_DIR, dummy_mdp)