import os
import subprocess
import sys

def extract_rg(tpr_file, xtc_file, out_file, t_start=None):
    """
    Run gmx gyrate to extract Rg vs time
    """

    cmd = ["gmx", "gyrate", "-f", xtc_file, "-s", tpr_file, "-o", out_file]

    if t_start is not None:
        cmd.extend(["-b", str(t_start)])

    print(f"[INFO] Running: {' '.join(cmd)}")

    res = subprocess.run(
        cmd,
        input="1\n",   
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


    xtc = os.path.join(
        BASE,
        "simulations",
        pid,
        "bulk",
        "equilibration",
        "equilibration.xtc"
    )

    tpr = os.path.join(
        BASE,
        "simulations",
        pid,
        "bulk",
        "equilibration",
        "equilibration.tpr"
    )

    out = os.path.join(
        BASE,
        "simulations",
        pid,
        "rg.xvg"
    )

    if not os.path.exists(xtc):
        print("[ERROR] Missing xtc:", xtc)
        sys.exit(1)

    if not os.path.exists(tpr):
        print("[ERROR] Missing tpr:", tpr)
        sys.exit(1)

    extract_rg(tpr, xtc, out, t_start=1000)

