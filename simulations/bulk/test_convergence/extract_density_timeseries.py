import os
import subprocess
import sys

def extract_density(edr_file, out_file, t_start=None):

    cmd = ["gmx", "energy", "-f", edr_file, "-o", out_file]

    if t_start is not None:
        cmd.extend(["-b", str(t_start)])

    print(f"[INFO] Running: {' '.join(cmd)}")

    res = subprocess.run(
        cmd,
        input="22\n",  # Density
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if res.returncode != 0:
        print("[ERROR]")
        print(res.stderr)
    else:
        print("[OK] Done:", out_file)


if __name__ == "__main__":
    pid = sys.argv[1]

    base_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(base_dir)

    BASE = os.path.join(
        parent_dir,
        "BIO-SUSHY-workflows",
        "FOOD_v2"
    )

    edr = os.path.join(
        BASE,
        "simulations",
        pid,
        "bulk",
        "equilibration",
        "equilibration.edr"
    )

    out = os.path.join(
        BASE,
        "simulations",
        pid,
        "density.xvg"
    )

    if not os.path.exists(edr):
        print("[ERROR] Missing:", edr)
        sys.exit(1)

    extract_density(edr, out, t_start=0)
