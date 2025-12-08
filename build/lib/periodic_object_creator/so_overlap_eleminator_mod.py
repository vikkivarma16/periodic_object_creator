import ctypes
import numpy as np
from ctypes import c_double, c_int, POINTER

# Load shared library

import os
here = os.path.dirname(__file__)
lib_path = os.path.join(here, "support_engine_overlap_finder.so")
lib = ctypes.CDLL(lib_path)



# Define C function signature
lib.overlap_eliminator.argtypes = [
    POINTER(c_double), c_int,
    POINTER(c_double), c_int,
    c_int,      # delete_from
    c_double,   # tolerance
    POINTER(c_int)
]

def overlap_eliminator(input_object_1, input_object_2,
                       delete_from='input_object_1',
                       tolerance=1e-6):

    A = np.array(input_object_1, dtype=float)
    B = np.array(input_object_2, dtype=float)

    if A.shape[1] < 3 or B.shape[1] < 3:
        raise ValueError("Each point must have at least x,y,z")

    delete_flag = 0 if delete_from == 'input_object_1' else 1

    NA = len(A)
    NB = len(B)

    delete_mask = np.zeros(NA if delete_flag==0 else NB, dtype=np.int32)

    lib.overlap_eliminator(
        A[:,:3].astype(np.float64).ctypes.data_as(POINTER(c_double)), NA,
        B[:,:3].astype(np.float64).ctypes.data_as(POINTER(c_double)), NB,
        delete_flag,
        c_double(tolerance),
        delete_mask.ctypes.data_as(POINTER(c_int))
    )

    if delete_flag == 0:
        new_A = [pt for pt,delm in zip(input_object_1, delete_mask) if delm == 0]
        new_B = input_object_2
    else:
        new_A = input_object_1
        new_B = [pt for pt,delm in zip(input_object_2, delete_mask) if delm == 0]

    return new_A, new_B

