def assign_group_ids(obj, group_size=3, start_id=1, id_index_in_element=None):
    """
    assign_group_ids(input_object, group_size=3, start_id=1, id_index=None)
    Assigns or updates group IDs for each consecutive block of atoms.

    input_object – list of atom records in the form:
                   [[x, y, z, attr1, attr2, ...], ...]
    group_size – number of atoms per molecule/group (e.g., 3 for H2O)
    start_id – ID assigned to the first group
    id_index – position of an existing ID inside the atom record:
                 * If None → no existing ID; append new ID as last element.
                 * If an integer → update the value at that index.

    Behavior:
    - Atoms are processed in order.
    - Every 'group_size' atoms share the same group ID.
    - If id_index is None:
          A new ID is appended to each atom record.
      If id_index is an integer:
          The ID at that index is replaced with the new group ID.

    Output:
    Returns a new list of atom records with group IDs assigned or updated.
    """

    new_obj = []
    current_id = start_id

    for i, atom in enumerate(obj):
        # increment id after each completed block
        if i % group_size == 0 and i > 0:
            current_id += 1

        atom_list = list(atom)

        if id_index_in_element is None:
            # append new ID
            atom_list.append(current_id)
        else:
            # ensure index exists
            while len(atom_list) <= id_index_in_element:
                atom_list.append(None)
            # update ID at index
            atom_list[id_index_in_element] = current_id

        new_obj.append(atom_list)

    return new_obj

