def scissor(input_object, plane_origin, plane_normal, keep_side="negative"):
    """
    Cuts a 3D object using a plane defined by: a point on the plane and a normal vector.
    
    input_object: list of input_object like [x, y, z, type, id, ...]
    plane_origin: [x0, y0, z0] a point on the cutting plane
    plane_normal: [nx, ny, nz] normal vector of the plane
    keep_side: "negative" (default) keeps input_object behind plane
               "positive" keeps input_object in front of plane
               
    returns: list of kept full input_object (not only coordinates)
    """
    import numpy as np

    n = np.array(plane_normal, dtype=float)
    n_norm = np.linalg.norm(n)

    if n_norm < 1e-12:
        raise ValueError("Plane normal cannot be zero.")
    n = n / n_norm   # normalize normal vector

    p0 = np.array(plane_origin, dtype=float)

    kept_objects = []

    for obj in input_object:
        # extract only coordinates (first three entries)
        p = np.array(obj[:3], dtype=float)

        # signed distance from point to plane
        dist = np.dot((p - p0), n)

        # filtering
        if keep_side == "negative":
            if dist <= 0:
                kept_objects.append(obj)

        elif keep_side == "positive":
            if dist >= 0:
                kept_objects.append(obj)

        else:
            raise ValueError("keep_side must be 'negative' or 'positive'.")

    return kept_objects

