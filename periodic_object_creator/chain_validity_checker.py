#periodic_object_creator/chain_configuration_validity_checker


def chain_validator(rings_size, contact_bond_length, isometric_dif_factor):
    
    
    from math import sin, cos, acos, asin, atan, tan, sqrt
    piee=acos(-1.0)
    tpiee=2.0*piee
    # even rings are horizontal
    # odd rings are vertical
    # minimum chain radius is decided by the horizontal rings
    size_normalizer=0.0
    for i in range(0,len(rings_size), 2):
        size_normalizer=size_normalizer+rings_size[i]
    
    factor=[]
    
    for i in range(len(rings_size)):
        factor.append(rings_size[i]/size_normalizer)
        
    angular_width=[]
    
    for i in range(len(rings_size)):
        angular_width.append(tpiee*factor[i])
    
    #print("minimum factor", angular_width, factor)
        
    r_minimum=0.5*(0.5*(rings_size[0]+rings_size[2]))/sin(0.5*((angular_width[0]+angular_width[2])*0.5))
    
   
    
    size_normalizer=0.0
    for i in range(len(rings_size)):
        size_normalizer=size_normalizer+rings_size[i]-contact_bond_length[i][0]-contact_bond_length[i][1]
    factor=[]
    
    for i in range(len(rings_size)):
        factor.append((rings_size[i]-contact_bond_length[i][0]-contact_bond_length[i][1])/size_normalizer)
        
    angular_width=[]
    
    for i in range(len(rings_size)):
        angular_width.append(tpiee*factor[i])
    
    if (isometric_dif_factor==1):
        red=contact_bond_length[1][0]+contact_bond_length[1][1]+contact_bond_length[2][0]+contact_bond_length[2][1]+contact_bond_length[3][0]+contact_bond_length[3][1]
        max_curvature=rings_size[1]+rings_size[2]+rings_size[3]-red
        angular_length=angular_width[1]+angular_width[2]+angular_width[3]
        r_maximum=max_curvature/angular_length
    
    else:
        max_curvature=rings_size[2]-contact_bond_length[2][0]-contact_bond_length[2][1]
        angular_length=angular_width[2]
        r_maximum=0.5*max_curvature/sin(0.5*angular_length)
        
        print("this is being fired")
        
    print("\n\n...the returned R_minimum is given as...", r_minimum, "\n\n")
    
    print("\n\n...the returned R_maximum is given as...", r_maximum, "\n\n")
    
        
        
    return r_minimum, r_maximum
    
