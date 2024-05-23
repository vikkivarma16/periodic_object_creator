#periodic_object_creator/so_rotator



def rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree):

    from math import sin, cos, fabs, acos, atan, tan, sqrt
    
    import so_translator as st
    
    
   
    tvector = []
    piee = acos(-1.0)

    for i in range(len(ro_axis_posi)):
        tvector.append(-ro_axis_posi[i])
    object_return = st.translator(input_object, tvector)
    objct = []
    for i in range(len(object_return)):
        objct.append(object_return[i])
    theta = ro_degree*piee/180.0
    dx = ro_axis_orien[0]
    dy = ro_axis_orien[1]
    dz = ro_axis_orien[2]
    q1 = cos(theta/2.0)
    q2 = dx*sin(theta/2.0)
    q3 = dy*sin(theta/2.0)
    q4 = dz*sin(theta/2.0)
    q11 = 1.0-2.0*((q3**2.0)+(q4**2.0))
    q21 = 2.0*((q2*q3)+(q4*q1))
    q31 = 2.0*((q2*q4)-(q3*q1))
    q12 = 2.0*((q2*q3)-(q4*q1))
    q22 = 1.0-2.0*((q2**2.0)+(q4**2.0))
    q32 = 2.0*((q3*q4)+(q2*q1))
    q13 = 2.0*((q2*q4)+(q3*q1))
    q23 = 2.0*((q3*q4)-(q2*q1))
    q33 = 1.0-2.0*((q2**2.0)+(q3**2.0))
    objct_return = []
    for i in range(len(objct)):
        dx = objct[i][0]*q11+objct[i][1]*q12+objct[i][2]*q13
        dy = objct[i][0]*q21+objct[i][1]*q22+objct[i][2]*q23
        dz = objct[i][0]*q31+objct[i][1]*q32+objct[i][2]*q33
        objct_return.append([dx, dy, dz])
    input_object = []
    for i in range(len(objct_return)):
        input_object.append(objct_return[i])
    tvector = []
    for i in range(len(ro_axis_posi)):
        tvector.append(ro_axis_posi[i])
    object_return = st.translator(input_object, tvector)
    
    input_object=[]
    for i in range(len(object_return)):
        input_object.append(object_return[i])
    object_return=[]
    for i in range(len(input_object)):
        object_return.append(input_object[i])
        
  
    
    return object_return

