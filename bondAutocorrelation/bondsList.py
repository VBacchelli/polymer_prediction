import MDAnalysis as mda
import networkx as nx
from pathlib import Path

def generate_backbone_cc_ndx(topology_file,
                             output_ndx="bonds.ndx"):
    """
    Identify backbone on heavy-atom network (no H),
    then write only C–C bonds belonging to that backbone in GROMACS .ndx format.

    topology_file: path to a .pdb OR a directory containing a .pdb
    """

    topology_file = Path(topology_file)

    # allow passing a directory like .../<polymer_id>/MD/
    if topology_file.is_dir():
        pdbs = sorted(topology_file.glob("*.pdb"))
        if not pdbs:
            raise ValueError(f"No .pdb found in {topology_file}")

        topology_file = pdbs[0]   # deterministic (usually polymer_vac_final.pdb)
        output_ndx = topology_file.parent / output_ndx

    # first load normally
    u = mda.Universe(str(topology_file))

    # if no bonds present → reload with guess_bonds=True
    if len(u.bonds) == 0:
        print("⚠ No bonds found — guessing bonds.")
        u = mda.Universe(str(topology_file), guess_bonds=True)

    # build graph with heavy atoms (no H) 
    heavy = {a.index for a in u.atoms if a.element != "H"}

    G = nx.Graph()
    for bond in u.bonds:
        a1, a2 = bond.atoms
        if a1.index in heavy and a2.index in heavy:
            G.add_edge(a1.index, a2.index)

    if len(G.nodes) == 0:
        raise ValueError("No heavy-atom bonds found (check bonds/guess_bonds).")

    # find largest connected component
    largest = max(nx.connected_components(G), key=len)
    subgraph = G.subgraph(largest).copy()

    # try to find real endpoints (degree 1 nodes)
    endpoints = [n for n in subgraph.nodes if subgraph.degree[n] == 1]

    if len(endpoints) == 2:
        backbone = nx.shortest_path(subgraph, endpoints[0], endpoints[1])
    else:
        # fallback to double BFS
        start = next(iter(subgraph.nodes))
        dist1 = nx.single_source_shortest_path_length(subgraph, start)
        far = max(dist1, key=dist1.get)

        dist2 = nx.single_source_shortest_path_length(subgraph, far)
        opp = max(dist2, key=dist2.get)

        backbone = nx.shortest_path(subgraph, far, opp)

    backbone_set = set(backbone)

    # extract C–C bonds that lie on the backbone (for autocorrelation of C-C bonds only)
    cc_pairs = []
    for bond in u.bonds:
        a1, a2 = bond.atoms
        if a1.index in backbone_set and a2.index in backbone_set:
            if a1.element == "C" and a2.element == "C":
                cc_pairs.append((a1.index + 1, a2.index + 1))  # GROMACS is 1-based

    if not cc_pairs:
        raise ValueError("No C–C backbone bonds found.")

    # write ndx
    with open(output_ndx, "w") as f:
        f.write("[ bonds ]\n")
        for i, j in cc_pairs:
            f.write(f"{i} {j}\n")

    print(f"Backbone atoms: {len(backbone)}")
    print(f"C–C bonds written: {len(cc_pairs)}")
    print(f"File written: {output_ndx}")

    return str(output_ndx)