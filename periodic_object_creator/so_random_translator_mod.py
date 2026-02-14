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
    rotation=False,
    periodic=True,
    max_trials=10000,
    eps=1.0e-8
):
    """
    Rigid-body randomization of molecules.

    Parameters
    ----------
    periodic : bool
        True  -> periodic boundary conditions
        False -> closed box (no atom may leave)
    """

    if seed is not None:
        random.seed(seed)

    Lx, Ly, Lz = box
    xi, yi, zi = xyz_indices
    output_object = deepcopy(input_object)

    # -------------------------------
    # Group atoms by molecule
    # -------------------------------
    groups = defaultdict(list)
    for i, rec in enumerate(output_object):
        groups[rec[idx]].append(i)

    # -------------------------------
    # Helper functions
    # -------------------------------
    def random_unit_vector():
        z = 2.0 * random.random() - 1.0
        t = 2.0 * math.pi * random.random()
        r = math.sqrt(1.0 - z * z)
        return [r * math.cos(t), r * math.sin(t), z]

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

    def wrap(x, L):
        return x % L

    def min_image(dx, L):
        if dx > 0.5 * L:
            dx -= L
        elif dx < -0.5 * L:
            dx += L
        return dx

    def inside_box(pos):
        return (
            eps < pos[0] < Lx - eps and
            eps < pos[1] < Ly - eps and
            eps < pos[2] < Lz - eps
        )

    # -------------------------------
    # Process each molecule
    # -------------------------------
    for indices in groups.values():

        # --- build molecule geometry ---
        ref = output_object[indices[0]]
        ref_pos = [ref[xi], ref[yi], ref[zi]]
        unwrapped = []

        for i in indices:
            rec = output_object[i]
            if periodic:
                x = ref_pos[0] + min_image(rec[xi] - ref_pos[0], Lx)
                y = ref_pos[1] + min_image(rec[yi] - ref_pos[1], Ly)
                z = ref_pos[2] + min_image(rec[zi] - ref_pos[2], Lz)
            else:
                x, y, z = rec[xi], rec[yi], rec[zi]

            unwrapped.append([x, y, z])

        # --- center of mass ---
        n = len(unwrapped)
        com = [sum(p[i] for p in unwrapped) / n for i in range(3)]

        # --- axis-wise molecule extent ---
        dx_max = max(abs(p[0] - com[0]) for p in unwrapped)
        dy_max = max(abs(p[1] - com[1]) for p in unwrapped)
        dz_max = max(abs(p[2] - com[2]) for p in unwrapped)

        # -------------------------------
        # Random placement
        # -------------------------------
        placed = False
        trials = 0

        while not placed:
            trials += 1
            if not periodic and trials > max_trials:
                raise RuntimeError("Failed to place molecule inside closed box")

            # --- translation ---
            if translation:
                if periodic:
                    new_com = [
                        random.random() * Lx,
                        random.random() * Ly,
                        random.random() * Lz,
                    ]
                else:
                    new_com = [
                        random.uniform(dx_max, Lx - dx_max),
                        random.uniform(dy_max, Ly - dy_max),
                        random.uniform(dz_max, Lz - dz_max),
                    ]
            else:
                new_com = com[:]

            disp = [new_com[i] - com[i] for i in range(3)]

            # --- rotation ---
            if rotation:
                axis = random_unit_vector()
                angle = (2.0 * random.random() - 1.0) * math.pi
            else:
                axis = None
                angle = 0.0

            trial_positions = []

            for pos in unwrapped:
                v = [pos[i] + disp[i] - new_com[i] for i in range(3)]
                if rotation:
                    v = rotate(v, axis, angle)
                p = [new_com[i] + v[i] for i in range(3)]
                trial_positions.append(p)

            # --- validity check ---
            if periodic:
                placed = True
            else:
                placed = all(inside_box(p) for p in trial_positions)

        # -------------------------------
        # Commit positions
        # -------------------------------
        for idx_i, p in zip(indices, trial_positions):
            rec = output_object[idx_i]
            if periodic:
                rec[xi] = wrap(p[0], Lx)
                rec[yi] = wrap(p[1], Ly)
                rec[zi] = wrap(p[2], Lz)
            else:
                rec[xi], rec[yi], rec[zi] = p

    return output_object
