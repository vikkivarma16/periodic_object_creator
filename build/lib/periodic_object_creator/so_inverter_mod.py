# periodic_object_creator/so_inverter

def inverter(input_object, inversion_point):
    """
    Performs inversion of a 3D object around a given inversion point,
    preserving all additional attributes beyond [x, y, z].
    
    Parameters:
    -----------
    input_object : list of points [x, y, z, ...]
    inversion_point : [px, py, pz] the point around which inversion occurs
    
    Returns:
    --------
    object_return : list of inverted points with all original attributes preserved
    """
    from math import sin, cos, acos, sqrt
    from . import so_translator_mod as tl

    piee = acos(-1.0)
    tpiee = 2.0 * piee

    # Step 1: translate inversion_point to origin
    tvector = [-inversion_point[0], -inversion_point[1], -inversion_point[2]]
    obj_translated = tl.translator(input_object, tvector)

    # Step 2: invert coordinates
    obj_inverted = []
    for p in obj_translated:
        x, y, z = p[0], p[1], p[2]
        r = sqrt(x*x + y*y + z*z)
        if r < 1e-12:  # avoid division by zero
            obj_inverted.append(p[:])  # keep attributes as is
            continue
        theta = acos(z / r)
        rxy = sqrt(x*x + y*y)
        if rxy < 1e-12:
            phi = 0.0
        else:
            phi = acos(x / rxy)
            if y < 0.0:
                phi = tpiee - phi
        r_new = -r
        x_new = r_new * sin(theta) * cos(phi)
        y_new = r_new * sin(theta) * sin(phi)
        z_new = r_new * cos(theta)
        # preserve extra attributes
        new_point = [x_new, y_new, z_new] + p[3:]
        obj_inverted.append(new_point)

    # Step 3: translate back
    tvector = [inversion_point[0], inversion_point[1], inversion_point[2]]
    object_return = tl.translator(obj_inverted, tvector)

    # preserve attributes already done
    return object_return

