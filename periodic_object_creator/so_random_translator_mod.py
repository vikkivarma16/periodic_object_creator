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
    Rigid-body randomization of molecules under PBC.

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
    output_object = deepcopy(input_object)

    # Group atoms by molecule
    groups = defaultdict(list)
    for i, rec in enumerate(output_object):
        groups[rec[idx]].append(i)

    # Helper: random unit vector
    def random_unit_vector():
        z = 2.0 * random.random() - 1.0
        t = 2.0 * math.pi * random.random()
        r = math.sqrt(1.0 - z * z)
        return [r * math.cos(t), r * math.sin(t), z]

    # Helper: Rodrigues' rotation formula
    def rotate(v, axis, angle):
        c = math.cos(angle)
        s = math.sin(angle)
        dot = sum(v[i] * axis[i] for i in range(3))
        cross = [
            axis[1] * v[2] - axis[2] * v[1],
            axis[2] * v[0] - axis[0] * v[2],
            axis[0] * v[1] - axis[1] * v[0],
        ]
        return [v[i] * c + cross[i] * s + axis[i] * dot * (1.0 - c) for i in range(3)]

    # Helper: PBC wrapping
    def wrap(coord, box_length):
        return coord % box_length

    # Helper: minimum image
    def min_image(dx, box_length):
        if dx > 0.5 * box_length:
            dx -= box_length
        elif dx < -0.5 * box_length:
            dx += box_length
        return dx

    # Process each molecule
    for indices in groups.values():
        # ---------- unwrap molecule relative to first atom ----------
        ref = output_object[indices[0]]
        ref_pos = [ref[xi], ref[yi], ref[zi]]
        unwrapped = []
        for i in indices:
            rec = output_object[i]
            x = ref_pos[0] + min_image(rec[xi] - ref_pos[0], Lx)
            y = ref_pos[1] + min_image(rec[yi] - ref_pos[1], Ly)
            z = ref_pos[2] + min_image(rec[zi] - ref_pos[2], Lz)
            unwrapped.append([x, y, z])

        # ---------- compute COM ----------
        n = len(unwrapped)
        com = [sum(u[i] for u in unwrapped) / n for i in range(3)]

        # ---------- random translation ----------
        if translation:
            new_com = [random.random() * Lx, random.random() * Ly, random.random() * Lz]
            disp = [new_com[i] - com[i] for i in range(3)]
        else:
            disp = [0.0, 0.0, 0.0]
            new_com = com

        # ---------- apply translation ----------
        for i, pos in zip(indices, unwrapped):
            rec = output_object[i]
            x = pos[0] + disp[0]
            y = pos[1] + disp[1]
            z = pos[2] + disp[2]
            rec[xi] = wrap(x, Lx)
            rec[yi] = wrap(y, Ly)
            rec[zi] = wrap(z, Lz)

        # ---------- random rotation ----------
        if rotation:
            axis = random_unit_vector()
            angle = (2.0 * random.random() - 1.0) * math.pi
            for i, pos in zip(indices, unwrapped):
                rec = output_object[i]
                # translate to COM
                v = [pos[j] + disp[j] - new_com[j] for j in range(3)]
                vr = rotate(v, axis, angle)
                # translate back and wrap
                rec[xi] = wrap(new_com[0] + vr[0], Lx)
                rec[yi] = wrap(new_com[1] + vr[1], Ly)
                rec[zi] = wrap(new_com[2] + vr[2], Lz)

    return output_object

