#periodic_object_creator/curvy_linear_basis




def curvy_linear_basis(basis, species, bond_length, angular_length_str):

    from math import sin, cos, fabs, acos, atan, tan, sqrt, asin
    piee = acos(-1.0)
    
    EPSILON=0.001
    

    
    if (angular_length_str=="periodic" or angular_length_str=="p"):
        basis_bonds=[]
        angular_length=2.0*piee
        for i in range(len(basis)):
            j=i-1
            k=i+1
            if (j<0):
                j=j+len(basis)
            if (k>=len(basis)):
                k=0
            basis_bonds.append([j, k])
            
            
        bond_indi=[]
        

        for i in range(len(basis_bonds)):
            for j in range(len(basis_bonds[i])):
                if(basis_bonds[i][j]>i):
                    bond_indi.append([i, basis_bonds[i][j]])
        
        bond_normalization=0.0
        for i in range(len(bond_indi)):
            for j in range(len(bond_length)):
                if (basis[bond_indi[i][0]] in bond_length[j] and basis[bond_indi[i][1]] in bond_length[j]):
                    
                    current_bl=bond_length[j][2]     
            bond_normalization=bond_normalization+current_bl
        bond_factor=bond_length[0][2]/bond_normalization
        
        
        angular_projection=angular_length*bond_factor
        angle=0.5*angular_projection
        half_bond=bond_length[0][2]*0.5
        curvature_radius=half_bond/sin(angle)
        
        
        basis_coordinate=[]
        
        linear_length=curvature_radius*angular_length
        
        # compatibility check....
        total_basis_length=0.0
        total_angular_length=0.0
        for i in range(len(basis)):
            for j in range(len(bond_length)):
                if (basis[basis_bonds[i][0]] in bond_length[j] and basis[basis_bonds[i][1]] in bond_length[j]):
                    current_bl=bond_length[j][2]
          
            angular_bl=2.0*asin(current_bl*0.5/curvature_radius)
            total_angular_length=total_angular_length+angular_bl
        
        parameter=fabs(total_angular_length-angular_length)
        
        
        
        if  (parameter>EPSILON):
            universal_bond_length=2.0*curvature_radius*sin((angular_length/len(basis))*0.5)
            for i in range(len(bond_length)):
                bond_length[i][2]=universal_bond_length
            print("...Universal bond length used as the given bond length was not compatible...")
            
        current_angular_posi=0.0
        
        basis_angular_posi=[]
        z=0  # linear curve is generated always in the XY-plain 
        for i in range(len(basis)):
            x=curvature_radius*cos(current_angular_posi)
            y=curvature_radius*sin(current_angular_posi)
            
            basis_angular_posi.append(current_angular_posi)
            if (i<len(basis)-1):
                for j in range(len(bond_length)):
                    if ( basis[i] in bond_length[j] and basis[basis_bonds[i][1]]):
                        current_bl=bond_length[j][2]
                angular_bond_length=2.0*asin(current_bl*0.5/curvature_radius)   
                current_angular_posi=current_angular_posi+angular_bond_length
                
            basis_coordinate.append([x,y,z])

            
    else:
        basis_bonds=[]
        angular_length=float(angular_length_str)
        for i in range(len(basis)):
            j=i-1
            k=i+1
            if (i==0):
                basis_bonds.append([k])
            elif (i==len(basis)-1):
                basis_bonds.append([j])
            else: 
                basis_bonds.append([j, k])
                
                
        bond_indi=[]
    
        for i in range(len(basis_bonds)):
            for j in range(len(basis_bonds[i])):
                if(basis_bonds[i][j]>i):
                    bond_indi.append([i, basis_bonds[i][j]])
        
        bond_normalization=0.0
        for i in range(len(bond_indi)):
            for j in range(len(bond_length)):
                if (bond_indi[i][0] in bond_length[j] and bond_indi[i][1] in bond_length[j]):
                    
                    current_bl=bond_length[j][2]        
            bond_normalization=bond_normalization+current_bl
        bond_factor=bond_length[0][2]/bond_normalization
        
        angular_projection=angular_length*bond_factor
        angle=0.5*angular_projection
        half_bond=bond_length[0][2]*0.5
        curvature_radius=half_bond/sin(angle)
        
        basis_coordinate=[]
        
        linear_length=curvature_radius*angular_length
        
        # compatibility check....
        total_basis_length=0.0
        total_angular_length=0.0
        for i in range(len(basis)):
            for j in range(len(bond_length)):
                if (basis[basis_bonds[i][0]] in bond_length[j] and basis[basis_bonds[i][1]] in bond_length[j]):
                    current_bl=bond_lendth[j][2]
                    
            angular_bl=2.0*asin(current_bl*0.5/curvature_radius)
            total_angular_length=total_angular_length+angular_bl
            
        if  (fabs(total_angular_length-angular_length)>EPSILON):
            universal_bond_length=2.0*curvature_radius*sin((angular_length/len(basis))*0.5)
            for i in range(len(bond_length)):
                bond_length[i][2]=universal_bond_length
            print("...Universal bond length used as the given bond length was not compatible...")
        
        
        basis_angular_posi=[]
        current_angular_posi=0.0
        z=0  # linear curve is generated always in the XY-plain 
        for i in range(len(basis)):
            x=curvature_radius*cos(current_angular_posi)
            y=curvature_radius*sin(current_angular_posi)
            
            basis_angular_posi.append(current_angular_posi)
            if (i<len(basis)-1):
                for j in range(len(bond_length)):
                    if ( basis[i] in bond_length[j] and basis[basis_bonds[i][1]]):
                        current_bl=bond_length[j][2]
                angular_bond_length=2.0*asin(current_bl*0.5/curvature_radius)   
                current_angular_posi=current_angular_posi+angular_bond_length
                
            basis_coordinate.append([x,y,z])
            
    return basis_angular_posi, basis_coordinate, bond_length, basis_bonds

    #...... particles coordinate lies over the curvature.....
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
