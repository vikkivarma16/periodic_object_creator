def get_object_size(input_object, coord_indices=(0,1,2)):
    """
    Calculate the size and extension of a 3D object along x, y, z axes.

    Parameters:
    -----------
    input_object : list
        List of object elements: [[x, y, z, ...], ...]
    coord_indices : tuple of ints
        Indices of coordinates in the element list (default: (0,1,2))

    Returns:
    --------
    dict
        {
            'x': {'min': float, 'max': float, 'size': float},
            'y': {'min': float, 'max': float, 'size': float},
            'z': {'min': float, 'max': float, 'size': float}
        }
    """

    if not input_object:
        return {'x': {}, 'y': {}, 'z': {}}

    xs = [atom[coord_indices[0]] for atom in input_object]
    ys = [atom[coord_indices[1]] for atom in input_object]
    zs = [atom[coord_indices[2]] for atom in input_object]

    size_info = {
        'x': {'min': min(xs), 'max': max(xs), 'size': max(xs)-min(xs)},
        'y': {'min': min(ys), 'max': max(ys), 'size': max(ys)-min(ys)},
        'z': {'min': min(zs), 'max': max(zs), 'size': max(zs)-min(zs)},
    }

    return size_info

