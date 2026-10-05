from pathlib import Path
import sys
import json

from rotacf_utils import run_rotacf, generate_backbone_ndx_from_folder

if __name__ == "__main__":
    polymer_dir = Path(sys.argv[1]).resolve()
    if len(sys.argv) > 2:
        dummy_mdp = Path(sys.argv[2]).resolve()
    else:
        dummy_mdp = Path.cwd() / "dummy.mdp"

    print("Processing:", polymer_dir)

    generate_backbone_ndx_from_folder(polymer_dir)

    tau = run_rotacf(polymer_dir, dummy_mdp)

    print("✅ Tau:", tau)

    
    out_file = polymer_dir / "rotacf.json"

    data = {"tau": tau}

    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)
        f.flush() 

