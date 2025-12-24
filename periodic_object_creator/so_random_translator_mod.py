import random
import math
from collections import defaultdict
from copy import deepcopy

def randomize_positions(
    input_object,
    idx,
    box,
    xyz_indices=(0, 1, 2),
    seed=None,
    translation=True,
    rotation=False
):
    """
    Rigid-body randomization of objects sharing the same ID with PBC.

    Parameters
    ----------
    input_object : list
        List of element records (lists or tuples).
    idx : int
        Index in each record that defines the rigid object ID (e.g. mol_id).
    box : tuple/list of 3 floats
        Simulation box size (Lx, Ly, Lz).
    xyz_indices : tuple of int, optional
        Indices of x, y, z coordinates inside each record.
    seed : int, optional
        Random seed.
    translation : bool
        Whether to apply random translation.
    rotation : bool
        Whether to apply random rotation after translation.

    Returns
    -------
    output_object : list
        Same object structure, coordinates updated with PBC.
    """

    if seed is not None:
        random.seed(seed)

    Lx, Ly, Lz = box
    xi, yi, zi = xyz_indices

    # Work on a deep copy to preserve original object
    output_object = deepcopy(input_object)

    # -------------------------------------------------
    # Group element indices by rigid object ID
    # -------------------------------------------------
    groups = defaultdict(list)
    for i, rec in enumerate(output_object):
        groups[rec[idx]].append(i)

    # -------------------------------------------------
    # Helper: random unit vector
    # -------------------------------------------------
    def random_unit_vector():
        z = 2.0 * random.random() - 1.0
        t = 2.0 * math.pi * random.random()
        r = math.sqrt(1.0 - z * z)
        return [r * math.cos(t), r * math.sin(t), z]

    # -------------------------------------------------
    # Helper: rotate vector using Rodrigues' formula
    # -------------------------------------------------
    def rotate(v, axis, angle):
        c = math.cos(angle)
        s = math.sin(angle)
        dot = sum(v[i] * axis[i] for i in range(3))
        cross = [
            axis[1] * v[2] - axis[2] * v[1],
            axis[2] * v[0] - axis[0] * v[2],
            axis[0] * v[1] - axis[1] * v[0],
        ]
        return [
            v[i] * c + cross[i] * s + axis[i] * dot * (1.0 - c)
            for i in range(3)
        ]

    # -------------------------------------------------
    # Apply PBC
    # -------------------------------------------------
    def apply_pbc(coord, box_length):
        return coord % box_length

    # -------------------------------------------------
    # Process each rigid object
    # -------------------------------------------------
    for _, indices in groups.items():

        # ---------- compute COM ----------
        com = [0.0, 0.0, 0.0]
        for i in indices:
            rec = output_object[i]
            com[0] += rec[xi]
            com[1] += rec[yi]
            com[2] += rec[zi]
        n = len(indices)
        com = [c / n for c in com]

        # ---------- random translation ----------
        if translation:
            new_com = [
                random.random() * Lx,
                random.random() * Ly,
                random.random() * Lz,
            ]
            disp = [
                new_com[0] - com[0],
                new_com[1] - com[1],
                new_com[2] - com[2],
            ]
        else:
            disp = [0.0, 0.0, 0.0]
            new_com = com

        # ---------- apply translation with PBC ----------
        for i in indices:
            rec = output_object[i]
            rec[xi] = apply_pbc(rec[xi] + disp[0], Lx)
            rec[yi] = apply_pbc(rec[yi] + disp[1], Ly)
            rec[zi] = apply_pbc(rec[zi] + disp[2], Lz)

        # ---------- random rotation ----------
        if rotation:
            axis = random_unit_vector()
            angle = (2.0 * random.random() - 1.0) * math.pi

            for i in indices:
                rec = output_object[i]
                v = [
                    rec[xi] - new_com[0],
                    rec[yi] - new_com[1],
                    rec[zi] - new_com[2],
                ]
                vr = rotate(v, axis, angle)
                rec[xi] = apply_pbc(new_com[0] + vr[0], Lx)
                rec[yi] = apply_pbc(new_com[1] + vr[1], Ly)
                rec[zi] = apply_pbc(new_com[2] + vr[2], Lz)

    return output_object

