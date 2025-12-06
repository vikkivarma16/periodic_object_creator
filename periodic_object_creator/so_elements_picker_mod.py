#periodic_object_creator/elements_picker


# periodic_object_creator/elements_picker

def elements_picker(input_object, indices):
    """
    Selects elements from input_object based on provided indices.
    Preserves all attributes of each element.

    Parameters:
    -----------
    input_object : list of points [x, y, z, ...]
    indices : list of integer indices to select

    Returns:
    --------
    return_object : list of selected points with all attributes
    """
    return [input_object[i] for i in indices]

