def filter_broken_group(obj, group_size=3):
    # Count atoms per id
    count = {}
    for atom in obj:
        idnum = atom[4]
        count[idnum] = count.get(idnum, 0) + 1

    # IDs that are complete molecules
    good_ids = {gid for gid, c in count.items() if c == group_size}

    # Keep only atoms whose ID is good
    filtered = [atom for atom in obj if atom[4] in good_ids]

    return filtered

