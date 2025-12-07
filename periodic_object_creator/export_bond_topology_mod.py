# Updated topology_builder.py with body index export and element counts
from math import sqrt
from itertools import combinations
from collections import defaultdict


def _distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)


def _canonical_bond(a, b):
    return tuple(sorted((a, b)))


def _canonical_angle(a, b, c):
    return min((a, b, c), (c, b, a))


def _canonical_dihedral(a, b, c, d):
    return min((a, b, c, d), (d, c, b, a))


def _canonical_improper(center, n1, n2, n3):
    return (center,) + tuple(sorted((n1, n2, n3)))


def find_bonds(input_object, bond_length, tolerance, many_body, id_index, id_body_index, coord_indices):
    if id_index is None:
        id_index = len(input_object[0]) - 1

    positions = {}
    body = {}
    ids = []

    for atom in input_object:
        rec = list(atom)
        aid = rec[id_index]
        x, y, z = rec[coord_indices[0]], rec[coord_indices[1]], rec[coord_indices[2]]

        if many_body:
            if id_body_index is None:
                raise ValueError("many_body=True requires id_body_index to be specified.")
            bid = rec[id_body_index]
        else:
            bid = 1

        positions[aid] = (float(x), float(y), float(z))
        body[aid] = bid
        ids.append(aid)

    bonds = set()
    neighbors = defaultdict(list)

    for i in range(len(ids)):
        ai = ids[i]
        for j in range(i + 1, len(ids)):
            aj = ids[j]
            if body[ai] != body[aj]:
                continue
            d = _distance(positions[ai], positions[aj])
            if abs(d - bond_length) <= tolerance:
                b = _canonical_bond(ai, aj)
                if b not in bonds:
                    bonds.add(b)
                    neighbors[b[0]].append(b[1])
                    neighbors[b[1]].append(b[0])

    return bonds, neighbors, positions, body


def build_angles(bonds, neighbors):
    angle_set = set()
    for center, neighs in neighbors.items():
        if len(neighs) < 2:
            continue
        for a, c in combinations(neighs, 2):
            angle_set.add(_canonical_angle(a, center, c))
    return [list(x) for x in sorted(angle_set)]


def build_dihedrals(bonds, neighbors):
    dihedral_set = set()
    for b, c in bonds:
        for B, C in ((b, c), (c, b)):
            neigh_B = [x for x in neighbors.get(B, []) if x != C]
            neigh_C = [x for x in neighbors.get(C, []) if x != B]
            for A in neigh_B:
                for D in neigh_C:
                    dihedral_set.add(_canonical_dihedral(A, B, C, D))
    return [list(x) for x in sorted(dihedral_set)]


def build_impropers(neighbors):
    improper_set = set()
    for center, neighs in neighbors.items():
        if len(neighs) < 3:
            continue
        for trio in combinations(neighs, 3):
            improper_set.add(_canonical_improper(center, *trio))
    return [list(x) for x in sorted(improper_set)]


def export_topology_text(filename_base, bonds, angles, dihedrals, impropers, positions, body):
    fname = f"{filename_base}.txt"
    total_bodies = len(set(body.values()))
    with open(fname, "w") as fh:
        fh.write("SUMMARY COUNTS")
        fh.write(f"TOTAL_ELEMENTS: {len(positions)}  \n ")
        fh.write(f"TOTAL_BODIES: {total_bodies}\n")
        fh.write(f"TOTAL_BONDS: {len(bonds)}  \n")
        fh.write(f"TOTAL_ANGLES: {len(angles)}  \n ")
        fh.write(f"TOTAL_DIHEDRALS: {len(dihedrals)} \n")
        fh.write(f"TOTAL_IMPROPERS: {len(impropers)} \n")

        fh.write("\nBONDS (# id1 id2 body | coords1 | coords2)\n\n")
        i = 1
        for a1, a2 in bonds:
            x1, y1, z1 = positions[a1]
            x2, y2, z2 = positions[a2]
            fh.write(f"{i}   {a1} {a2} {body[a1]}    "
                     f"{x1:.6f} {y1:.6f} {z1:.6f}    "
                     f"{x2:.6f} {y2:.6f} {z2:.6f}\n")
            i = i+1 

        fh.write("\nANGLES (# A B C body | coordsA | coordsB | coordsC)\n\n")
        i = 1
        for A, B, C in angles:
            xA, yA, zA = positions[A]
            xB, yB, zB = positions[B]
            xC, yC, zC = positions[C]
            fh.write(f"{i}   {A} {B} {C} {body[B]}    "
                     f"{xA:.6f} {yA:.6f} {zA:.6f}    "
                     f"{xB:.6f} {yB:.6f} {zB:.6f}    "
                     f"{xC:.6f} {yC:.6f} {zC:.6f}\n")
            i = i+1

        fh.write("\nDIHEDRALS (# A B C D body | coords...)\n\n")
        i=1
        for A, B, C, D in dihedrals:
            xA, yA, zA = positions[A]
            xB, yB, zB = positions[B]
            xC, yC, zC = positions[C]
            xD, yD, zD = positions[D]
            fh.write(f"{i}   {A} {B} {C} {D} {body[B]}    "
                     f"{xA:.6f} {yA:.6f} {zA:.6f}    "
                     f"{xB:.6f} {yB:.6f} {zB:.6f}    "
                     f"{xC:.6f} {yC:.6f} {zC:.6f}    "
                     f"{xD:.6f} {yD:.6f} {zD:.6f}\n")
            i = i+1

        fh.write("\nIMPROPERS (# n1 center n2 n3 body | coords...)\n\n")
        i=1
        for center, n1, n2, n3 in impropers:
            xC, yC, zC = positions[center]
            x1, y1, z1 = positions[n1]
            x2, y2, z2 = positions[n2]
            x3, y3, z3 = positions[n3]
            fh.write(f"{i}   {n1} {center} {n2} {n3} {body[center]}    "
                     f"{x1:.6f} {y1:.6f} {z1:.6f}    "
                     f"{xC:.6f} {yC:.6f} {zC:.6f}    "
                     f"{x2:.6f} {y2:.6f} {z2:.6f}    "
                     f"{x3:.6f} {y3:.6f} {z3:.6f}\n")
            i = i+1


    return fname


def build_topology(input_object, bond_length, tolerance, id_index=None,
                   coord_indices=(0, 1, 2), export_base="topology",
                   many_body=False, id_body_index=None):

    bonds, neighbors, positions, body = find_bonds(
        input_object,
        bond_length,
        tolerance,
        many_body,
        id_index,
        id_body_index,
        coord_indices,
    )

    bonds_list = [list(b) for b in sorted(bonds)]
    angles = build_angles(bonds, neighbors)
    dihedrals = build_dihedrals(bonds, neighbors)
    impropers = build_impropers(neighbors)

    export_file = export_topology_text(
        export_base, bonds_list, angles, dihedrals, impropers, positions, body
    )

    return {
        "bonds": bonds_list,
        "neighbors": {k: sorted(v) for k, v in neighbors.items()},
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers,
        "body": body,
        "export_file": export_file,
    }

