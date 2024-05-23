# periodic_object_creator/so_reflector


def reflection(input_object, plain_normal, plain_location):

    from math import sin, cos, fabs, acos, atan, tan, sqrt
    piee = acos(-1.0)

    for i in range(len(input_object)):
        dx=input_object[i][0]-plain_location[0]
        dy=input_object[i][1]-plain_location[1]
        dz=input_object[i][2]-plain_location[2]
        magnitude=sqrt(dx*dx+dy*dy+dz*dz)
        dirx=dx/magnitude
        diry=dy/magintude
        dirz=dz/magnitude
        theta=acos(dirx*plain_normal[0] + dirx*plain_normal[1] + dirx*plain_normal[2])
        #incomplete....
    
    
    return object_return
