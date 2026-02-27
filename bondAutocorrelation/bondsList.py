import os
from pathlib import Path
from bson import ObjectId
from rdkit import Chem
import MDAnalysis as mda
import networkx as nx
from dotenv import load_dotenv
from pymongo import MongoClient

# ------------------ Mongo ------------------

load_dotenv()
_client = MongoClient(os.getenv("MONGO_URI"))
_collection = _client["PolymerPrediction"]["biceranoPolymers"]

def get_smiles(polymer_id):
    doc = (_collection.find_one({"_id": polymer_id})
           or _collection.find_one({"_id": ObjectId(polymer_id)}))
    if not doc:
        raise ValueError(f"No document found for ID {polymer_id}")
    return doc["smiles"]

# ------------------ SMILES ------------------

def get_attachments(smiles):
    mol = Chem.MolFromSmiles(smiles)
    wc = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(wc) != 2:
        raise ValueError("Expected exactly two [*] in SMILES")

    def heavy_deg(a):
        return sum(n.GetAtomicNum() > 1 for n in a.GetNeighbors())

    return [(n.GetSymbol(), heavy_deg(n)) for n in (w.GetNeighbors()[0] for w in wc)]

# ------------------ PDB loader ------------------

def load_universe(path):
    path = Path(path)

    # if folder → search pdb
    if path.is_dir():
        pdbs = sorted(path.glob("*.pdb"))
        if not pdbs:
            raise ValueError(f"No .pdb found in {path}")
        path = pdbs[0]

    u = mda.Universe(str(path))

    # if no bonds, guess
    if not u.bonds:
        print("⚠ No bonds found — guessing bonds.")
        u = mda.Universe(str(path), guess_bonds=True)

    return u, path

# ------------------ Backbone ------------------

def generate_backbone_ndx_from_folder(folder_path, verbose=True):

    folder = Path(folder_path)

    u, pdb_path = load_universe(folder)

    smiles = get_smiles(folder.name)
    attachments = get_attachments(smiles)

    if verbose:
        print("Attachment info:", attachments)
        print("Using PDB:", pdb_path.name)

    # heavy atom graph
    heavy = {a.index for a in u.atoms if a.element != "H"}
    G = nx.Graph((b.atoms[0].index, b.atoms[1].index)
                 for b in u.bonds
                 if b.atoms[0].index in heavy and b.atoms[1].index in heavy)

    if not G.nodes:
        raise ValueError("No heavy-atom bonds found.")

    deg1 = [n for n in G.nodes if G.degree(n) == 1]

    def candidates(elem, deg):
        return [n for n in deg1
                if u.atoms[n].element == elem and G.degree(n) == deg]

    start_set = candidates(*attachments[0])
    end_set   = candidates(*attachments[1])

    if not start_set or not end_set:
        raise ValueError("Could not match attachment atoms in PDB")

    best = max(
        ((i, j, nx.shortest_path_length(G, i, j))
         for i in start_set for j in end_set
         if nx.has_path(G, i, j)),
        key=lambda x: x[2],
        default=None
    )

    if not best:
        raise ValueError("No valid backbone endpoints found.")

    backbone = nx.shortest_path(G, best[0], best[1])
    if len(backbone) < 3:
        raise ValueError("Backbone too short.")

    cc_bonds = [(u.atoms[i].index + 1, u.atoms[j].index + 1)
                for i, j in zip(backbone[:-1], backbone[1:])
                if u.atoms[i].element == "C" and u.atoms[j].element == "C"]

    ndx_path = folder / "bonds.ndx"
    with open(ndx_path, "w") as f:
        f.write("[ bonds ]\n")
        for i, j in cc_bonds:
            f.write(f"{i:6d} {j:6d}\n")

    if verbose:
        print(f"Written {len(cc_bonds)} C–C bonds → {ndx_path.name}")

    return ndx_path

# ------------------ CLI ------------------

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <polymer_folder_or_MD_folder>")
        sys.exit(1)

    generate_backbone_ndx_from_folder(sys.argv[1])