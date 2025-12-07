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

"""
Enhanced topology_builder.py
Now also exports:
- number of bonds, angles, dihedrals, impropers
- coordinates of each atom involved in every topology interaction
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
    return tuple(sorted((a, b)))


def _canonical_angle(a, b, c):
    t1 = (a, b, c)
    t2 = (c, b, a)
    return min(t1, t2)


def _canonical_dihedral(a, b, c, d):
    t1 = (a, b, c, d)
    t2 = (d, c, b, a)
    return min(t1, t2)


def _canonical_improper(center, n1, n2, n3):
    neigh_sorted = tuple(sorted((n1, n2, n3)))
    return (center,) + neigh_sorted


def find_bonds(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2)):
    if id_index is None:
        id_index = len(input_object[0]) - 1

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

    n = len(id_list)
    for i in range(n):
        id_i = id_list[i]
        pos_i = id_to_pos[id_i]
        for j in range(i+1, n):
            id_j = id_list[j]
            pos_j = id_to_pos[id_j]
            d = _distance(pos_i, pos_j)
            if abs(d - bond_length) <= tolerance:
                b = _canonical_bond(id_i, id_j)
                if b not in bonds_set:
                    bonds_set.add(b)
                    neighbors[b[0]].append(b[1])
                    neighbors[b[1]].append(b[0])

    return bonds_set, neighbors, id_to_pos


def build_angles(bonds_set, neighbors):
    angle_set = set()
    for center, neighs in neighbors.items():
        if len(neighs) < 2:
            continue
        for a, c in combinations(neighs, 2):
            angle_set.add(_canonical_angle(a, center, c))
    return [list(x) for x in sorted(angle_set)]


def build_dihedrals(bonds_set, neighbors):
    dihedral_set = set()
    for b, c in bonds_set:
        for (B, C) in ((b, c), (c, b)):
            neigh_B = [x for x in neighbors.get(B, []) if x != C]
            neigh_C = [x for x in neighbors.get(C, []) if x != B]
            if not neigh_B or not neigh_C:
                continue
            for A in neigh_B:
                for D in neigh_C:
                    dihedral_set.add(_canonical_dihedral(A, B, C, D))
    return [list(x) for x in sorted(dihedral_set)]


def build_impropers(neighbors):
    improper_set = set()
    for central, neighs in neighbors.items():
        if len(neighs) < 3:
            continue
        for trio in combinations(neighs, 3):
            improper_set.add(_canonical_improper(central, *trio))
    return [list(x) for x in sorted(improper_set)]


def export_topology_text(filename_base, bonds, angles, dihedrals, impropers, id_to_pos):
    fname = f"{filename_base}.txt"
    with open(fname, "w") as fh:

        # ===== SUMMARY =====
        fh.write("SUMMARY COUNTS\n")
        fh.write(f"TOTAL_BONDS: {len(bonds)}\n")
        fh.write(f"TOTAL_ANGLES: {len(angles)}\n")
        fh.write(f"TOTAL_DIHEDRALS: {len(dihedrals)}\n")
        fh.write(f"TOTAL_IMPROPERS: {len(impropers)}\n\n")

        # ===== BONDS =====
        fh.write("BONDS (id1 id2 | coords1 | coords2)\n")
        for b in bonds:
            a1, a2 = b
            x1, y1, z1 = id_to_pos[a1]
            x2, y2, z2 = id_to_pos[a2]
            fh.write(f"{a1} {a2} | {x1:.6f} {y1:.6f} {z1:.6f} | {x2:.6f} {y2:.6f} {z2:.6f}\n")

        # ===== ANGLES =====
        fh.write("\nANGLES (A B C | coordsA | coordsB | coordsC)\n")
        for A, B, C in angles:
            fh.write(
                f"{A} {B} {C} | "
                f"{id_to_pos[A][0]:.6f} {id_to_pos[A][1]:.6f} {id_to_pos[A][2]:.6f} | "
                f"{id_to_pos[B][0]:.6f} {id_to_pos[B][1]:.6f} {id_to_pos[B][2]:.6f} | "
                f"{id_to_pos[C][0]:.6f} {id_to_pos[C][1]:.6f} {id_to_pos[C][2]:.6f}\n"
            )

        # ===== DIHEDRALS =====
        fh.write("\nDIHEDRALS (A B C D | coordsA | coordsB | coordsC | coordsD)\n")
        for A, B, C, D in dihedrals:
            fh.write(
                f"{A} {B} {C} {D} | "
                f"{id_to_pos[A][0]:.6f} {id_to_pos[A][1]:.6f} {id_to_pos[A][2]:.6f} | "
                f"{id_to_pos[B][0]:.6f} {id_to_pos[B][1]:.6f} {id_to_pos[B][2]:.6f} | "
                f"{id_to_pos[C][0]:.6f} {id_to_pos[C][1]:.6f} {id_to_pos[C][2]:.6f} | "
                f"{id_to_pos[D][0]:.6f} {id_to_pos[D][1]:.6f} {id_to_pos[D][2]:.6f}\n"
            )

        # ===== IMPROPERS =====
        fh.write("\nIMPROPERS (center n1 n2 n3 | coords ...)\n")
        for center, n1, n2, n3 in impropers:
            fh.write(
                f"{center} {n1} {n2} {n3} | "
                f"{id_to_pos[center][0]:.6f} {id_to_pos[center][1]:.6f} {id_to_pos[center][2]:.6f} | "
                f"{id_to_pos[n1][0]:.6f} {id_to_pos[n1][1]:.6f} {id_to_pos[n1][2]:.6f} | "
                f"{id_to_pos[n2][0]:.6f} {id_to_pos[n2][1]:.6f} {id_to_pos[n2][2]:.6f} | "
                f"{id_to_pos[n3][0]:.6f} {id_to_pos[n3][1]:.6f} {id_to_pos[n3][2]:.6f}\n"
            )

    return fname


def build_topology(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="topology"):

    bonds_set, neighbors, id_to_pos = find_bonds(
        input_object, bond_length, tolerance, id_index=id_index, coord_indices=coord_indices
    )

    bonds = [list(b) for b in sorted(bonds_set)]
    angles = build_angles(bonds_set, neighbors)
    dihedrals = build_dihedrals(bonds_set, neighbors)
    impropers = build_impropers(neighbors)

    export_file = export_topology_text(
        export_base, bonds, angles, dihedrals, impropers, id_to_pos
    )

    return {
        "bonds": bonds,
        "neighbors": {k: sorted(v) for k, v in neighbors.items()},
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers,
        "export_file": export_file
    }

