def assign_group_ids(obj, group_size=3, start_id=1):
    new_obj = []
    current_id = start_id

    for i, atom in enumerate(obj):
        if i % group_size == 0 and i > 0:
            current_id += 1

        x, y, z, ele, _ = atom
        new_obj.append([x, y, z, ele, current_id])

    return new_obj

