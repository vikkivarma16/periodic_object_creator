"""
topology_builder.py

Functions:
- build_topology(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2))
    Detects bonded pairs, neighbor lists, bond angles (triplets), dihedrals (quartets) and improper dihedrals.
    Exports results to a text file and returns a dictionary with arrays.

Usage:
- input_object: list of atom records, each record is a list-like:
      [x, y, z, attr1, attr2, ..., id]
  If id_index is None, the function expects the last element to be the id (or you can provide id_index).
- bond_length: expected bond length (float)
- tolerance: distance tolerance (float). Bond if abs(distance - bond_length) <= tolerance
- id_index: index of the id field in each atom record. If None, function will try to use the last element.
- coord_indices: tuple with indices of (x,y,z) in each record (default (0,1,2))

Output (dictionary):
{
  "bonds": list of [id1, id2],
  "neighbors": {id: [neigh_id,...], ...},
  "angles": list of [idA, idB, idC],
  "dihedrals": list of [idA, idB, idC, idD],
  "impropers": list of [central_id, n1, n2, n3],
  "export_file": "<filename>.txt"
}
"""

from math import sqrt
from itertools import combinations, product
from collections import defaultdict

def _distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)

def _canonical_bond(a, b):
    # canonical representation of bond: sorted tuple
    return tuple(sorted((a, b)))

def _canonical_angle(a, b, c):
    # angle a-b-c and c-b-a are the same; store canonical with lexicographic min between (a,b,c) and (c,b,a)
    t1 = (a, b, c)
    t2 = (c, b, a)
    return min(t1, t2)

def _canonical_dihedral(a, b, c, d):
    # dihedral a-b-c-d equivalent to reverse d-c-b-a; choose lexicographically smaller of the two
    t1 = (a, b, c, d)
    t2 = (d, c, b, a)
    return min(t1, t2)

def _canonical_improper(center, n1, n2, n3):
    # central atom first, neighbors as a sorted tuple to canonicalize order
    neigh_sorted = tuple(sorted((n1, n2, n3)))
    return (center,) + neigh_sorted

def find_bonds(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2)):
    """
    find_bonds(...)
    Finds bonded neighbors according to bond_length +/- tolerance.

    Returns:
      bonds_set: set of canonical bond tuples (id1, id2)
      neighbors: dict mapping id -> list of neighbor ids (unsorted)
      id_to_pos: dict mapping id -> (x,y,z)
    """
    # determine id_index if None -> assume last element
    if id_index is None:
        id_index = len(input_object[0]) - 1

    # build id -> position mapping and list to iterate
    id_list = []
    id_to_pos = {}
    for atom in input_object:
        atom_list = list(atom)
        aid = atom_list[id_index]
        x = atom_list[coord_indices[0]]
        y = atom_list[coord_indices[1]]
        z = atom_list[coord_indices[2]]
        id_list.append(aid)
        id_to_pos[aid] = (float(x), float(y), float(z))

    bonds_set = set()
    neighbors = defaultdict(list)

    ids = id_list
    n = len(ids)
    for i in range(n):
        id_i = ids[i]
        pos_i = id_to_pos[id_i]
        for j in range(i+1, n):
            id_j = ids[j]
            pos_j = id_to_pos[id_j]
            d = _distance(pos_i, pos_j)
            if abs(d - bond_length) <= tolerance:
                b = _canonical_bond(id_i, id_j)
                if b not in bonds_set:
                    bonds_set.add(b)
                    # add directed neighbors for subsequent angle/dihedral generation
                    neighbors[b[0]].append(b[1])
                    neighbors[b[1]].append(b[0])

    return bonds_set, neighbors, id_to_pos

def build_angles(bonds_set, neighbors):
    """
    build_angles(...)
    From neighbor lists, build unique angles A-B-C (center = B).
    Returns list of canonical triplets [A,B,C].
    """
    angle_set = set()
    for center, neighs in neighbors.items():
        if len(neighs) < 2:
            continue
        # all unordered pairs of neighbors produce an angle at 'center'
        for a, c in combinations(neighs, 2):
            canon = _canonical_angle(a, center, c)
            angle_set.add(canon)
    # return as list of lists
    return [list(x) for x in sorted(angle_set)]

def build_dihedrals(bonds_set, neighbors):
    """
    build_dihedrals(...)
    For each bond B-C, consider neighbors A of B (A != C) and D of C (D != B).
    All combinations A-B-C-D form dihedrals (canonicalized).
    Returns list of canonical dihedral quartets.
    """
    dihedral_set = set()
    for b, c in bonds_set:
        # note bonds_set uses sorted order; need to consider both directions when picking neighbors
        # We'll treat bond endpoints as (B,C) in both orders to find all combos
        for (B, C) in ((b, c), (c, b)):
            neigh_B = [x for x in neighbors.get(B, []) if x != C]
            neigh_C = [x for x in neighbors.get(C, []) if x != B]
            if not neigh_B or not neigh_C:
                continue
            for A in neigh_B:
                for D in neigh_C:
                    canon = _canonical_dihedral(A, B, C, D)
                    dihedral_set.add(canon)
    return [list(x) for x in sorted(dihedral_set)]

def build_impropers(neighbors):
    """
    build_impropers(...)
    For each central atom with 3+ neighbors, build improper tuples:
    (central, n1, n2, n3) for each combination of 3 neighbors. Canonicalized so neighbors sorted.
    """
    improper_set = set()
    for central, neighs in neighbors.items():
        if len(neighs) < 3:
            continue
        for trio in combinations(neighs, 3):
            canon = _canonical_improper(central, *trio)
            improper_set.add(canon)
    return [list(x) for x in sorted(improper_set)]

def export_topology_text(filename_base, bonds, angles, dihedrals, impropers):
    """
    export_topology_text(filename_base, bonds, angles, dihedrals, impropers)
    Writes topology lists to <filename_base>.txt in a readable format.

    Each section header is printed, then each interaction on its own line as space-separated ids.
    """
    fname = f"{filename_base}.txt"
    with open(fname, "w") as fh:
        fh.write("BONDS (id1 id2)\n")
        for b in bonds:
            fh.write(f"{b[0]} {b[1]}\n")

        fh.write("\nANGLES (idA idB idC)\n")
        for a in angles:
            fh.write(f"{a[0]} {a[1]} {a[2]}\n")

        fh.write("\nDIHEDRALS (idA idB idC idD)\n")
        for d in dihedrals:
            fh.write(f"{d[0]} {d[1]} {d[2]} {d[3]}\n")

        fh.write("\nIMPROPERS (center n1 n2 n3)\n")
        for imp in impropers:
            fh.write(f"{imp[0]} {imp[1]} {imp[2]} {imp[3]}\n")

    return fname

def build_topology(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="topology"):
    """
    build_topology(...)
    High-level function that finds bonds, neighbors, angles, dihedrals and impropers.
    Exports results to a text file (export_base + '.txt') and returns a dict with arrays.

    Returns:
      {
        "bonds": [[id1,id2], ...],
        "neighbors": {id: [neighs...]},
        "angles": [[a,b,c], ...],
        "dihedrals": [[a,b,c,d], ...],
        "impropers": [[center,n1,n2,n3], ...],
        "export_file": "<export_base>.txt"
      }
    """
    bonds_set, neighbors, id_to_pos = find_bonds(
        input_object, bond_length, tolerance, id_index=id_index, coord_indices=coord_indices
    )

    # convert bond set to sorted list of pairs for consistent output
    bonds = [list(b) for b in sorted(bonds_set)]

    angles = build_angles(bonds_set, neighbors)
    dihedrals = build_dihedrals(bonds_set, neighbors)
    impropers = build_impropers(neighbors)

    export_file = export_topology_text(export_base, bonds, angles, dihedrals, impropers)

    return {
        "bonds": bonds,
        "neighbors": dict((k, sorted(v)) for k, v in neighbors.items()),
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers,
        "export_file": export_file
    }

# Example quick test (uncomment to run as script)
if __name__ == "__main__":
    # tiny example: two water molecules with ids at last index
    atoms = [
        [0.0, 0.0, 0.0, "O", 1],
        [0.96, 0.0, 0.0, "H", 2],   # bonded to O (~0.96A)
        [-0.24, 0.93, 0.0, "H", 3],
        [3.0, 0.0, 0.0, "O", 4],
        [3.96, 0.0, 0.0, "H", 5],
        [2.76, 0.93, 0.0, "H", 6],
    ]
    # bond length approx 0.96, tolerance 0.2
    topo = build_topology(atoms, bond_length=0.96, tolerance=0.2, id_index=4, export_base="example_topo")
    print("Bonds:", topo["bonds"])
    print("Angles:", topo["angles"])
    print("Dihedrals:", topo["dihedrals"])
    print("Impropers:", topo["impropers"])
    print("Exported to:", topo["export_file"])

