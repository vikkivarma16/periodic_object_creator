def wrapper_cylindrical(input_object, cylinder_radius, object_size=None):
    """
    Wraps a sheet-like object around a cylinder whose axis is the Z-axis.
    The sheet's X-direction becomes the angular direction around the cylinder.

    input_object : list of elements [x, y, z, ...]
    cylinder_radius : radius of the cylinder
    object_size : optional, used only to warn if sheet exceeds circumference
    """

    from math import sin, cos, pi

    # Cylinder circumference
    max_size = 2.0 * pi * cylinder_radius

    # Optional size check
    if object_size is not None and object_size > max_size:
        print("Warning: object_size exceeds cylinder circumference...")

    wrapped = []

    for item in input_object:
        x, y, z = item[0], item[1], item[2]

        # Convert x → angular coordinate (around z-axis)
        theta = x / cylinder_radius   # radians

        # Mapping to cylindrical shell around Z-axis
        new_x = cylinder_radius * cos(theta)
        new_y = cylinder_radius * sin(theta)
        new_z = z  # IMPORTANT: original y becomes z-direction sheet height

        # Preserve all attributes after index 3
        wrapped.append([new_x, new_y, new_z] + item[3:])

    return wrapped

