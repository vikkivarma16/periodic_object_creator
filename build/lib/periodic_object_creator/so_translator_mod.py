# periodic_object_creator/so_translator_mod.py

def translator(input_object, tvector):
    """
    Translates only the coordinate part (x, y, z) of each entry in input_object.
    All additional attributes are preserved exactly.

    Parameters
    ----------
    input_object : list of lists
        Example element: [x, y, z, type, id1, id2, ...]
    tvector : list or tuple of length 3
        Translation vector [tx, ty, tz]

    Returns
    -------
    object_return : list of lists
        Same as input but with updated coordinates.
    """
    
    object_return = []

    tx, ty, tz = tvector

    for item in input_object:

        # extract coordinates
        x, y, z = float(item[0]), float(item[1]), float(item[2])

        # translate them
        new_x = x + tx
        new_y = y + ty
        new_z = z + tz

        # keep all extra attributes untouched
        new_item = [new_x, new_y, new_z] + item[3:]

        object_return.append(new_item)

    return object_return

