def cm_calculator(input_object):
    """
    Compute the center of mass (CM) of an object.
    Each element in input_object must be [x, y, z, type].
    If all atoms have identical mass (e.g., all 'C'), CM is the coordinate average.
    """
    n = len(input_object)
    if n == 0:
        return None

    sum_x = sum(p[0] for p in input_object)
    sum_y = sum(p[1] for p in input_object)
    sum_z = sum(p[2] for p in input_object)

    cm = [sum_x / n, sum_y / n, sum_z / n]
    return cm



