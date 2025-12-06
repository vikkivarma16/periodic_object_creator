def reflector(input_object, plane_normal, plane_location):
    """
    Reflects a 3D object (list of points) across a plane, preserving all extra attributes.

    Parameters:
    -----------
    input_object : list of [x, y, z, ...] points with optional extra attributes
    plane_normal : [nx, ny, nz] normal vector of the plane
    plane_location : [px, py, pz] a point on the plane

    Returns:
    --------
    object_return : list of reflected points with preserved attributes
    """
    from math import sqrt

    # Normalize the plane normal
    mag = sqrt(plane_normal[0]**2 + plane_normal[1]**2 + plane_normal[2]**2)
    n = [plane_normal[0]/mag, plane_normal[1]/mag, plane_normal[2]/mag]

    object_return = []

    for pt in input_object:
        # vector from plane point to current point
        v = [pt[0] - plane_location[0],
             pt[1] - plane_location[1],
             pt[2] - plane_location[2]]

        # dot product
        dot = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]

        # reflection formula: p' = p - 2 * (v . n) * n
        reflected_coords = [
            pt[0] - 2*dot*n[0],
            pt[1] - 2*dot*n[1],
            pt[2] - 2*dot*n[2]
        ]

        # preserve all other attributes beyond x, y, z
        object_return.append(reflected_coords + pt[3:])

    return object_return

