# periodic_object_creator/so_spherical_wrapping.py

def wrapper_spherical(input_object, sphere_radius, object_size):
    
    from math import sin, cos, atan2, acos, sqrt, pi
    from . import so_translator_mod as tl

    # ------------------------------------------------------------
    # 1) Translate object so COM sits at sphere "north pole"
    # ------------------------------------------------------------
    
    # Compute center of mass (only coords)
    cm = [0.0, 0.0, 0.0]
    for p in input_object:
        cm[0] += p[0]
        cm[1] += p[1]
        cm[2] += p[2]
    cm = [c/len(input_object) for c in cm]

    # Step A: move COM to origin
    obj1 = tl.translator(input_object, [-cm[0], -cm[1], -cm[2]])

    # Step B: move origin → north pole
    obj2 = tl.translator(obj1, [0.0, 0.0, sphere_radius])

    # ------------------------------------------------------------
    # 2) Convert each point to spherical coordinates
    # ------------------------------------------------------------

    def cartesian_to_spherical(x, y, z):
        r  = sqrt(x*x + y*y + z*z)
        theta = acos(z / r)
        phi   = atan2(y, x)
        return r, theta, phi

    def spherical_to_cartesian(r, theta, phi):
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        return [x, y, z]

    # ------------------------------------------------------------
    # 3) Angular scaling
    # ------------------------------------------------------------

    delta_theta = object_size / sphere_radius

    # Compute θ values
    thetas = []
    for p in obj2:
        _, th, _ = cartesian_to_spherical(p[0], p[1], p[2])
        thetas.append(th)

    theta_center = sum(thetas) / len(thetas)

    min_theta = min(thetas)
    max_theta = max(thetas)
    original_span = max_theta - min_theta

    if original_span < 1e-10:
        scale = 1.0
    else:
        scale = delta_theta / original_span

    # ------------------------------------------------------------
    # 4) Apply wrapping and preserve attributes
    # ------------------------------------------------------------

    wrapped = []
    for p, p_original in zip(obj2, input_object):
        r, th, ph = cartesian_to_spherical(p[0], p[1], p[2])

        th_new = theta_center + scale * (th - theta_center)
        ph_new = ph
        r_new = sphere_radius

        new_coords = spherical_to_cartesian(r_new, th_new, ph_new)

        # keep attributes from original object
        new_entry = new_coords + p_original[3:]

        wrapped.append(new_entry)

    return wrapped

