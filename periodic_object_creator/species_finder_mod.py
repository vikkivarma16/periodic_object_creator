#periodic_object_forger/species_finder

def species_finder(basis):
    species = []
    for i in range(len(basis)):
        flag = 1
        for j in range(len(species)):
            if species[j] == basis[i]:
                flag = 0
                break
        if flag == 1:
            species.append(basis[i])
    return species
