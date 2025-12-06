# periodic_object_creator/so_rotator

def rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree):

    from math import sin, cos, fabs, acos, atan, tan, sqrt
    from . import so_translator_mod as st

    # ------------------------------------------------------------------
    # 1) Translate rotation axis position to origin
    # ------------------------------------------------------------------
    tvector = [-x for x in ro_axis_posi]
    obj1 = st.translator(input_object, tvector)

    # Make a deep copy
    objct = [row[:] for row in obj1]

    # ------------------------------------------------------------------
    # 2) Build quaternion rotation matrix
    # ------------------------------------------------------------------
    pi = acos(-1.0)
    theta = ro_degree * pi / 180.0

    dx, dy, dz = ro_axis_orien
    q1 = cos(theta / 2.0)
    q2 = dx * sin(theta / 2.0)
    q3 = dy * sin(theta / 2.0)
    q4 = dz * sin(theta / 2.0)

    q11 = 1.0 - 2.0 * (q3**2 + q4**2)
    q21 = 2.0 * (q2*q3 + q4*q1)
    q31 = 2.0 * (q2*q4 - q3*q1)

    q12 = 2.0 * (q2*q3 - q4*q1)
    q22 = 1.0 - 2.0 * (q2**2 + q4**2)
    q32 = 2.0 * (q3*q4 + q2*q1)

    q13 = 2.0 * (q2*q4 + q3*q1)
    q23 = 2.0 * (q3*q4 - q2*q1)
    q33 = 1.0 - 2.0 * (q2**2 + q3**2)

    # ------------------------------------------------------------------
    # 3) Rotate ALL objects but preserve extra attributes
    # ------------------------------------------------------------------
    rotated = []
    for row in objct:
        x, y, z = row[0], row[1], row[2]

        # rotate coords
        dxr = x*q11 + y*q12 + z*q13
        dyr = x*q21 + y*q22 + z*q23
        dzr = x*q31 + y*q32 + z*q33

        # preserve attributes
        rotated.append([dxr, dyr, dzr] + row[3:])

    # ------------------------------------------------------------------
    # 4) Translate back to original axis position
    # ------------------------------------------------------------------
    obj2 = st.translator(rotated, ro_axis_posi)

    # deep copy final output
    final_obj = [row[:] for row in obj2]

    return final_obj

