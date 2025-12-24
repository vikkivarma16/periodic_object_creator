def overlap_remover(input_object,
                    molid_idx,
                    particle_idx,
                    mol_type_idx,
                    particle_type_idx,
                    sigma_matrix,
                    moving_mol_id,
                    box,
                    cell_size,
                    iter_max,
                    translation_step,
                    rotation_step,
                    max_particles_per_cell=64,
                    grid_shifting_rate=100000):

    import numpy as np
    import ctypes
    from ctypes import c_int, c_double, POINTER
    import os
    import sys

    # -------------------------------
    # Extract data
    # -------------------------------
    N = len(input_object)

    coords = np.zeros((N, 3), dtype=np.float64)
    mol_id = np.zeros(N, dtype=np.int32)
    part_type = np.zeros(N, dtype=np.int32)

    for i, rec in enumerate(input_object):
        coords[i] = rec[0:3]
        mol_id[i] = rec[molid_idx]
        part_type[i] = rec[particle_type_idx]

    # -------------------------------
    # REORDER particles so molecules are contiguous
    # -------------------------------
    order = np.argsort(mol_id)
    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(N)

    coords = coords[order]
    mol_id = mol_id[order]
    part_type = part_type[order]

    # -------------------------------
    # Build molecule indexing (NOW SAFE)
    # -------------------------------
    unique_mol, mol_count = np.unique(mol_id, return_counts=True)
    n_mol = len(unique_mol)

    mol_start = np.zeros(n_mol, dtype=np.int32)
    mol_count = mol_count.astype(np.int32)

    for i in range(1, n_mol):
        mol_start[i] = mol_start[i-1] + mol_count[i-1]

    # -------------------------------
    # Movable molecules
    # -------------------------------
    mol_movable = np.zeros(n_mol, dtype=np.int32)
    moving_set = set(moving_mol_id)

    for i, m in enumerate(unique_mol):
        if m in moving_set:
            mol_movable[i] = 1

    # -------------------------------
    # FIXED sigma_matrix handling
    # -------------------------------
    sigma_matrix = np.asarray(sigma_matrix, dtype=np.float64)
    if sigma_matrix.ndim != 2 or sigma_matrix.shape[0] != sigma_matrix.shape[1]:
        raise ValueError("sigma_matrix must be a square 2D array")

    n_ptype = sigma_matrix.shape[0]
    sigma = sigma_matrix.ravel()

    # -------------------------------
    # Prepare C arguments
    # -------------------------------
    coords_flat = coords.reshape(-1)
    box_arr = np.asarray(box, dtype=np.float64)

    # -------------------------------
    # Load shared library (portable)
    # -------------------------------
    here = os.path.dirname(__file__)

    if sys.platform == "darwin":
        libname = "support_engine_overlap_remover.dylib"
    elif sys.platform == "win32":
        libname = "support_engine_overlap_remover.dll"
    else:
        libname = "support_engine_overlap_remover.so"

    lib_path = os.path.join(here, libname)
    lib = ctypes.CDLL(lib_path)

    # -------------------------------
    # EXACT C SIGNATURE (17 args)
    # -------------------------------
    lib.relax_spherical_particles.argtypes = [
        POINTER(c_double),  # coords
        POINTER(c_int),     # mol_id
        POINTER(c_int),     # part_type
        c_int,              # N
        c_int,              # n_mol
        POINTER(c_int),     # mol_start
        POINTER(c_int),     # mol_count
        POINTER(c_int),     # mol_movable
        POINTER(c_double),  # sigma
        c_int,              # n_ptype
        POINTER(c_double),  # box
        c_double,           # cell_size
        c_int,              # iter_max
        c_double,           # step_trans
        c_double,           # step_rot
        c_int,              # max_particles_per_cell
        c_int               # grid_shifting_rate
    ]

    # -------------------------------
    # Call C function
    # -------------------------------
    lib.relax_spherical_particles(
        coords_flat.ctypes.data_as(POINTER(c_double)),
        mol_id.ctypes.data_as(POINTER(c_int)),
        part_type.ctypes.data_as(POINTER(c_int)),
        c_int(N),
        c_int(n_mol),
        mol_start.ctypes.data_as(POINTER(c_int)),
        mol_count.ctypes.data_as(POINTER(c_int)),
        mol_movable.ctypes.data_as(POINTER(c_int)),
        sigma.ctypes.data_as(POINTER(c_double)),
        c_int(n_ptype),
        box_arr.ctypes.data_as(POINTER(c_double)),
        c_double(cell_size),
        c_int(iter_max),
        c_double(translation_step),
        c_double(rotation_step),
        c_int(max_particles_per_cell),
        c_int(grid_shifting_rate)
    )

    # -------------------------------
    # Restore original particle order
    # -------------------------------
    coords = coords_flat.reshape((N, 3))
    coords = coords[inv_order]

    for i in range(N):
        input_object[i][0] = coords[i, 0]
        input_object[i][1] = coords[i, 1]
        input_object[i][2] = coords[i, 2]

    return input_object

