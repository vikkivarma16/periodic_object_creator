#periodic_object_forger/bond_combinator

#It finds all the bond configuration and set the bond length one by default

def bond_combinator(species):
    
    bond_length=[]
    for i in range(len(species)):
        for j in range(len(species)):
            bond_length.append([species[i], species[j], 1.0])
    
    return bond_length
