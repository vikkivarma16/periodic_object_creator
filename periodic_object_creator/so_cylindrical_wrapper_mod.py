# periodic_object_creator/so_cylindrical_wrapping


#... it wraps a sheet like object with normal along the x axis around a cylinder along z axis ...

def wrap(input_object, cylinder_radius, object_size):
    
    from math import sin, cos, tan, atan, acos, asin, sqrt, fabs
    
    import math
    
    
    piee=acos(-1.0)
    
    maximum_os=2.0*piee*cylinder_radius
    if (maximum_os < object_size):
        print("...Warning:  Sorry! Empty object as been returned as the object size is greater than the size of the surface of the cylinders...")
        object_return=[]
        
    else:
        from . import so_translator_mod as tl
    
        tvector=[cylinder_radius, 0.0, 0.0]
        objct_return =tl.translator(input_object, tvector)
    

        objct_unit_2 = []

        for i in range(len(objct_return)):
            objct_unit_2.append(objct_return[i])

        cm_x = 0.0
        cm_y = 0.0
        cm_z = 0.0

        for i in range(len(objct_unit_2)):
            cm_x = cm_x+objct_unit_2[i][0]
            cm_y = cm_y+objct_unit_2[i][1]
            cm_z = cm_z+objct_unit_2[i][2]

        cm_x = cm_x/float(len(objct_unit_2))
        cm_y = cm_y/float(len(objct_unit_2))
        cm_z = cm_z/float(len(objct_unit_2))

        modxy = sqrt(cm_x*cm_x+cm_y*cm_y)
        theta_center = acos(cm_x/modxy)

        if (cm_y < 0.0):
            theta_center = 2.0*piee-theta_center

        if (fabs(theta_center-2.0*piee) < 0.000001):
            theta_center = 0.0
        # this work only when you have the center of mass of the object placed at x axis and its surface normal is along the x axis
    
        theta = 2.0*atan(object_size*0.5/cylinder_radius) # projection angle is being calculated....

       # print("theta_center", theta_center)

        length_curvature = theta*cylinder_radius

        factor = object_size/length_curvature

    

        for i in range(len(objct_unit_2)):
            x = objct_unit_2[i][0]
            y = objct_unit_2[i][1]
            z = objct_unit_2[i][2]
            modxy = sqrt(x*x+y*y)
            theta = acos(x/modxy)

            if (y < 0.0):
                theta = 2.0*piee-theta

            if (theta < piee):
                theta_diff = theta-theta_center
                theta_diff = factor*theta_diff
                theta_new = theta_diff+theta_center

            else:
                theta_diff = theta-theta_center
                theta_diff = 2.0*piee-theta_diff
                theta_diff = factor*theta_diff
                theta_new_diff = 2.0*piee-theta_diff
                theta_new = theta_new_diff+theta_center

            nx = cylinder_radius*cos(theta_new)
            ny = cylinder_radius*sin(theta_new)
            nz = z

            objct_unit_2[i][0] = nx
            objct_unit_2[i][1] = ny
            objct_unit_2[i][2] = nz

        input_object = []
        for i in range(len(objct_unit_2)):
            input_object.append(objct_unit_2[i])

        tvector = [-cylinder_radius, 0.0, 0.0]
        objct_return =tl.translator(input_object, tvector)
        #print("\n\n\n...Object have been wrapped successfully...\n\n\n")
        
        
        input_object = []
        for i in range(len(objct_return)):
            input_object.append(objct_return[i])
        
        object_return=[]
        for i in range(len(input_object)):
            object_return.append(input_object[i])

        
    

    return object_return
