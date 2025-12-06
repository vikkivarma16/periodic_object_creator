def overlap_eleminator(input_object_1, input_object_2, delete_from='input_object_1', tolerance=1e-6):
    """
    Deletes coordinates in one object that overlap with another object,
    preserving all extra attributes.

    Parameters:
    -----------
    input_object_1 : list of [x, y, z, ...] points with optional extra attributes
    input_object_2 : list of [x, y, z, ...] points with optional extra attributes
    delete_from : str, 'input_object_1' or 'input_object_2', which object to remove overlapping points from
    tolerance : float, distance tolerance to consider points as overlapping

    Returns:
    --------
    new_obj1, new_obj2 : lists of points after deletion, with attributes preserved
    """
    from math import sqrt

    def is_overlap(p1, p2, tol):
        return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2) < tol

    if delete_from == 'input_object_1':
        new_obj1 = [p1 for p1 in input_object_1 if all(not is_overlap(p1, p2, tolerance) for p2 in input_object_2)]
        new_obj2 = [p2[:] for p2 in input_object_2]  # preserve all attributes
    elif delete_from == 'input_object_2':
        new_obj2 = [p2 for p2 in input_object_2 if all(not is_overlap(p2, p1, tolerance) for p1 in input_object_1)]
        new_obj1 = [p1[:] for p1 in input_object_1]  # preserve all attributes
    else:
        raise ValueError("delete_from must be 'input_object_1' or 'input_object_2'")

    return new_obj1, new_obj2

