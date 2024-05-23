#periodic_object_creator/ example_object_simple_ring

def example_object_simple_ring(n_beads, r_curvature):
    
    
    if (n_beads>=6):
        from math import sin, cos, fabs, acos, atan, tan, sqrt


        piee = acos(-1.0)

        
        import species_finder as sf
        # Decide the basis.. For example ring has curvy geometry... with angular length 360 degree...
        basis=["A"]
        ju=[]
        for i in range(len(basis)):
            ju.append(basis[i])
        i=0    
        while (len(basis)<n_beads):
            basis.append(ju[i])
            i = i + 1
            if i >= len(ju):
                i = 0
        species=sf.species_finder(basis)
        
        import bond_combinator as bc
        
        bond_length=bc.bond_combinator(species)
        
        
        
    
        
        # sets default bond length to 1
        # for the non default values please specify the bond lengh in the formate 'bond_length=[[A, A, 1.1], [A, B , 1.3], [B, B, 1.2]....]' with all the possible combination
        
        angular_length_str="periodic"
        
        
        angle=piee/float(n_beads) # half of the angle....
        
        bond_l=2.0*r_curvature*sin(angle)
        
        for i in range(len(bond_length)):
            bond_length[i][2]=bond_l
            
            
        print(bond_length)

        import basis_curvy as clb
        
        basis_angular_posi, basis_coordinate, bond_length, basis_bonds=clb.curvy_linear_basis(basis, species, bond_length, angular_length_str)
        
        

    
    else:
        basis=[]
        species=[]
        bond_length=[]
        basis_bonds=[]
        basis_corodinate=[]
        
        
        print ("...Ring has not generated aborting the process... ...with all the zero zero values")
        
    
    return basis, species, bond_length, basis_bonds, basis_coordinate
