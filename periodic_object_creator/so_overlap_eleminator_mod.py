import ctypes
import numpy as np
from ctypes import c_double, c_int, POINTER
import os

# Load shared library
here = os.path.dirname(__file__)
lib_path = os.path.join(here, "support_engine_overlap_finder.so")
lib = ctypes.CDLL(lib_path)

# Define C function signature
lib.overlap_eliminator.argtypes = [
    POINTER(c_double), c_int,          # object A
    POINTER(c_double), c_int,          # object B
    c_int,                             # delete_from
    c_double,                          # tolerance
    POINTER(c_int)                     # delete mask
]

def overlap_eliminator(input_object_1, input_object_2,
                       delete_from='input_object_1',
                       tolerance=1e-6,
                       xyz_indices=None):
    """
    Deletes overlapping points between input_object_1 and input_object_2
    using the C shared library. Preserves extra attributes.

    Parameters
    ----------
    input_object_1, input_object_2 : list of points
        Each point can have extra attributes (e.g., [x,y,z,'C'])
    delete_from : str
        'input_object_1' or 'input_object_2' – which object to delete overlaps from
    tolerance : float
        Distance tolerance for overlap
    xyz_indices : list or tuple of 3 ints, optional
        Indices in the points corresponding to x, y, z. Defaults to [0,1,2]

    Returns
    -------
    new_obj1, new_obj2 : lists with all attributes preserved
    """

    if xyz_indices is None:
        xyz_indices = [0, 1, 2]

    # Extract coordinates and convert to float
    def extract_coords(obj):
        arr = np.array(obj, dtype=object)
        if arr.shape[1] < max(xyz_indices)+1:
            raise ValueError("Points do not contain enough dimensions for the given xyz_indices")
        coords = np.array([ [float(pt[i]) for i in xyz_indices] for pt in obj], dtype=np.float64)
        return coords

    A_coords = extract_coords(input_object_1)
    B_coords = extract_coords(input_object_2)

    delete_flag = 0 if delete_from == 'input_object_1' else 1
    NA, NB = len(A_coords), len(B_coords)

    # Prepare mask
    delete_mask = np.zeros(NA if delete_flag == 0 else NB, dtype=np.int32)

    # Call C function
    lib.overlap_eliminator(
        A_coords.ctypes.data_as(POINTER(c_double)), NA,
        B_coords.ctypes.data_as(POINTER(c_double)), NB,
        delete_flag,
        c_double(tolerance),
        delete_mask.ctypes.data_as(POINTER(c_int))
    )

    # Construct new objects preserving all attributes
    if delete_flag == 0:
        new_A = [pt for pt, delm in zip(input_object_1, delete_mask) if delm == 0]
        new_B = input_object_2
    else:
        new_A = input_object_1
        new_B = [pt for pt, delm in zip(input_object_2, delete_mask) if delm == 0]

    return new_A, new_B

