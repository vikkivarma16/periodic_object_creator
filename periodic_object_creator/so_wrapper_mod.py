import copy
import math


def periodic_wrapper(objects, box):
    """
    Wrap all particles back into the primary simulation box.

    Parameters
    ----------
    objects : list
        Each element is like:
        [x, y, z, ...]
        First three entries MUST be coordinates.

    box : list or tuple
        Periodic box dimensions [Lx, Ly, Lz]

    Returns
    -------
    wrapped : list
        Deep-copied list with wrapped coordinates
    """

    wrapped = copy.deepcopy(objects)
    Lx, Ly, Lz = box

    for i, obj in enumerate(objects):
        x, y, z = obj[:3]

        wrapped[i][0] = x - math.floor(x / Lx) * Lx
        wrapped[i][1] = y - math.floor(y / Ly) * Ly
        wrapped[i][2] = z - math.floor(z / Lz) * Lz

    return wrapped

