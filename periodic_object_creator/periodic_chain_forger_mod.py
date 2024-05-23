# periodic_object_creator/periodic_ring_chain_forge



def periodic_ring_chain_forge(chain_radius, rings_coordinate, rings_size, chain_configuration_mode, specific_factor):
    
    
    from . import so_translator_mod as tl
    
    from . import so_rotator_mod as rt
    from math import sin, cos, fabs, acos, atan, tan, sqrt
    piee=acos(-1.0)
    tpiee=piee*2.0
    if (chain_configuration_mode==0):
        
        size_normalizer=0.0
        for i in range(len(rings_size)):
            size_normalizer=size_normalizer+rings_size[i]
        factor=[]
        
        for i in range(len(rings_size)):
            factor.append(rings_size[i]/size_normalizer)
        angular_position=[]
      
        ap=0.0
        for i in range(len(rings_size)):
            ap=ap+tpiee*factor[i]
            angular_position.append(ap)
        
        
        chain_coordinate=[]
        
        
        
        # specific factor can include any specification to a particular ring which is part of the chain. e.g, in this case it is about displacing the ring to avoid any hard core repulsion
        displacement_factor=[]
        for i in range(len(rings_size)):
            if (specific_factor==0):
            
                displacement_factor.append(0.0)
            
            
            else:
                val=rings_size[i]*0.5
                displacement_factor.append(sqrt(chain_radius*chain_radius-val*val)-chain_radius)
                
                
        
        for i in range(len(rings_coordinate)):
            input_object=[]
            for j in range(len(rings_coordinate[i])):
                input_object.append(rings_coordinate[i][j])
            tvector=[chain_radius+displacement_factor[i],0.0, 0.0]
            
            object_return=tl.translator(input_object, tvector)
            
           
            
            input_object=[]
            for j in range(len(object_return)):
                input_object.append(object_return[j])
            
           
            
            ro_axis_orien = [0.0, 0.0, 1.0]
            ro_axis_posi = [0.0, 0.0, 0.0]
            ro_degree = angular_position[i]*180/piee  # in degree
        
            
            object_return = rt.rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)
            
            input_object=[]
            for j in range(len(object_return)):
                input_object.append(object_return[j])
            
            
            chain_coordinate.append(input_object)
            
            
    else:
        chain_coordinate=[]
            
            
            
    return chain_coordinate       
