from collections import defaultdict
import copy

def min_image(dx, box):
    if dx >  0.5 * box:
        dx -= box
    elif dx < -0.5 * box:
        dx += box
    return dx


def periodic_unwrapper(objects, box, mol_id_index):
    """
    Unwrap molecules using minimum-image convention.

    Parameters
    ----------
    objects : list
        Each element is like:
        [x, y, z, ..., mol_id, ...]
        First three entries MUST be coordinates.

    box : list or tuple
        Periodic box dimensions [Lx, Ly, Lz]

    mol_id_index : int
        Index in each object that contains the molecule ID

    Returns
    -------
    unwrapped : list
        Deep-copied list with only coordinates updated
    """

    # group particle indices by molecule id
    mol_indices = defaultdict(list)
    for i, obj in enumerate(objects):
        mol_indices[obj[mol_id_index]].append(i)

    # deep copy so original data is preserved
    unwrapped = copy.deepcopy(objects)

    # unwrap each molecule independently
    for mol_id, indices in mol_indices.items():

        if len(indices) <= 1:
            continue  # nothing to unwrap

        # reference particle (first in molecule)
        ref_idx = indices[0]
        ref = objects[ref_idx][:3]

        for idx in indices[1:]:
            x, y, z = objects[idx][:3]

            dxx = min_image(x - ref[0], box[0])
            dyy = min_image(y - ref[1], box[1])
            dzz = min_image(z - ref[2], box[2])

            unwrapped[idx][0] = ref[0] + dxx
            unwrapped[idx][1] = ref[1] + dyy
            unwrapped[idx][2] = ref[2] + dzz

    return unwrapped

