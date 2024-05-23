# periodic_object_creator/jason_file_creator

def json_writer(input_data_js, filename):
    
    basis=[]
    
    
    import json

    li = []
    particle = []

    for i in range(len(input_data_js)):
        for j in range(len(input_data_js[i])):
            ind_particle = { 
                "id": input_data_js[i][j][0],
                "species": input_data_js[i][j][1],
                "coordinate": str( [input_data_js[i][j][2][0], input_data_js[i][j][2][1], input_data_js[i][j][2][2]]),
                "bonded_neighbours": str(input_data_js[i][j][3])
            }
            particle.append(ind_particle)

    
    
    
    # calculate the molecular properties bonds_angle
    total_triangle=[]
    total_particle=[]
    total_bonds=[]
    for i in range(len(input_data_js)):
        permutation=[]
        particle_set=[]
        bond_set=[]
        for j in range(len(input_data_js[i])):
            particle_set.append(input_data_js[i][j][0])
            
            for k in range(len(input_data_js[i][j][3])):
                if (input_data_js[i][j][0]<input_data_js[i][j][3][k]):
                    bond_set.append([input_data_js[i][j][0], input_data_js[i][j][3][k]])
                for l in range( k+1, len(input_data_js[i][j][3]) ):
                    permutation.append([input_data_js[i][j][0], input_data_js[i][j][3][k], input_data_js[i][j][3][l]])
                    
        total_triangle.append(permutation)
        total_particle.append(particle_set)
        total_bonds.append(bond_set)
    
    molecule=[]
    for i in range(len(input_data_js)):
        ind_molecule = { 
            "beads": str(total_particle[i]),
            "bonds": str(total_bonds[i]),
            "bond_angle": str(total_triangle[i])
        }
        molecule.append(ind_molecule)

    
    
    systems = {
        "system": {
            "particle": particle,
            "molecule": molecule
        }
    }

    # Use the filename variable correctly and open the file
    with open(filename, "w") as fu:
        json.dump(systems, fu, indent=4)

    print("Data has been written to", filename)

    
    
    return 0
