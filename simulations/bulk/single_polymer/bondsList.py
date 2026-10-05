import os
from pathlib import Path
from bson import ObjectId
from rdkit import Chem
import MDAnalysis as mda
import networkx as nx
from dotenv import load_dotenv
from pymongo import MongoClient

# ==========================================================
# Mongo
# ==========================================================

load_dotenv("mongo.env")
_client = MongoClient(os.getenv("MONGO_URI"))
_collection = _client["PolymerPrediction"]["biceranoPolymers"]

def get_smiles(polymer_id):
    doc = (_collection.find_one({"_id": polymer_id})
           or _collection.find_one({"_id": ObjectId(polymer_id)}))
    if not doc:
        raise ValueError(f"No document found for ID {polymer_id}")
    return doc["smiles"]

# ==========================================================
# SMILES → monomer backbone sequence
# ==========================================================

def get_monomer_backbone_sequence(smiles):
    mol = Chem.MolFromSmiles(smiles)
    wc = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]

    if len(wc) != 2:
        raise ValueError("Expected exactly two [*] in SMILES")

    a1 = wc[0].GetNeighbors()[0]
    a2 = wc[1].GetNeighbors()[0]

    path = Chem.rdmolops.GetShortestPath(mol, a1.GetIdx(), a2.GetIdx())
    sequence = [mol.GetAtomWithIdx(i).GetSymbol() for i in path]

    return sequence

# ==========================================================
# PDB loader
# ==========================================================

def load_universe(path):
    path = Path(path)

    if path.is_dir():
        pdbs = sorted(path.glob("*.pdb"))
        if not pdbs:
            raise ValueError(f"No .pdb found in {path}")
        path = pdbs[0]

    u = mda.Universe(str(path))

    if not u.bonds:
        print("⚠ No bonds found — guessing bonds.")
        u = mda.Universe(str(path), guess_bonds=True)

    return u, path

# ==========================================================
# Heavy graph
# ==========================================================

def build_heavy_graph(u):
    heavy = {a.index for a in u.atoms if a.element != "H"}

    G = nx.Graph(
        (b.atoms[0].index, b.atoms[1].index)
        for b in u.bonds
        if b.atoms[0].index in heavy and b.atoms[1].index in heavy
    )

    if not G.nodes:
        raise ValueError("No heavy-atom bonds found.")

    return G

# ==========================================================
# Sequence check
# ==========================================================

def sequence_matches(polymer_sequence, monomer_sequence):
    L = len(monomer_sequence)

    if len(polymer_sequence) % L != 0:
        return False

    for i in range(0, len(polymer_sequence), L):
        if polymer_sequence[i:i+L] != monomer_sequence:
            return False

    return True

# ==========================================================
# Backbone extraction
# ==========================================================

def extract_backbone(G, u, monomer_sequence):

    elem_A = monomer_sequence[0]
    elem_C = monomer_sequence[-1]

    A_nodes = [n for n in G.nodes if u.atoms[n].element == elem_A]
    C_nodes = [n for n in G.nodes if u.atoms[n].element == elem_C]

    best_path = None
    max_len = -1

    for i in A_nodes:
        for j in C_nodes:
            if nx.has_path(G, i, j):
                path = nx.shortest_path(G, i, j)

                seq = [u.atoms[n].element for n in path]

                if sequence_matches(seq, monomer_sequence):
                    if len(path) > max_len:
                        max_len = len(path)
                        best_path = path

    if best_path is None:
        raise ValueError("No backbone satisfying full monomer periodicity found.")

    return best_path

# ==========================================================
# Main
# ==========================================================

def generate_backbone_ndx_from_folder(folder_path, verbose=True):

    folder = Path(folder_path)
    md_folder = folder / "MD"

    if not md_folder.exists():
        raise FileNotFoundError(f"Missing MD folder: {md_folder}")

    u, pdb_path = load_universe(md_folder)

    smiles = get_smiles(folder.name)
    monomer_sequence = get_monomer_backbone_sequence(smiles)

    if verbose:
        print("Using PDB:", pdb_path.name)
        print("Monomer backbone sequence:", monomer_sequence)

    G = build_heavy_graph(u)
    backbone = extract_backbone(G, u, monomer_sequence)

    if len(backbone) < len(monomer_sequence):
        raise ValueError("Backbone too short.")

    # ======================================================
    # BONDS
    # ======================================================

    cc_bonds = [
        (u.atoms[i].index + 1, u.atoms[j].index + 1)
        for i, j in zip(backbone[:-1], backbone[1:])
        if u.atoms[i].element == "C" and u.atoms[j].element == "C"
    ]

    ndx_path = md_folder / "bonds.ndx"

    with open(ndx_path, "w") as f:
        f.write("[ bonds ]\n")
        for i, j in cc_bonds:
            f.write(f"{i:6d} {j:6d}\n")

    # ======================================================
    # DIHEDRALS 
    # ======================================================

    dihedrals = [
        (
            u.atoms[a].index + 1,
            u.atoms[b].index + 1,
            u.atoms[c].index + 1,
            u.atoms[d].index + 1,
        )
        for a, b, c, d in zip(
            backbone[:-3],
            backbone[1:-2],
            backbone[2:-1],
            backbone[3:]
        )
    ]

    dih_path = md_folder / "dihedrals.ndx"

    with open(dih_path, "w") as f:
        f.write("[ dihedrals ]\n")
        for i, j, k, l in dihedrals:
            f.write(f"{i:6d} {j:6d} {k:6d} {l:6d}\n")

    # ======================================================

    if verbose:
        print(f"Backbone atoms: {len(backbone)}")
        print(f"C–C bonds written: {len(cc_bonds)}")
        print(f"Dihedrals written: {len(dihedrals)}")   # ✅ nuovo
        print(f"Saved → {ndx_path}")
        print(f"Saved → {dih_path}")  # ✅ nuovo

    return ndx_path, dih_path


# ==========================================================
# CLI
# ==========================================================

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <polymer_folder>")
        sys.exit(1)

    generate_backbone_ndx_from_folder(sys.argv[1])

