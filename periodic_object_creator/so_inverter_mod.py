# periodic_object_creator/so_inverter


def inversion(input_object, inversion_point):
    # Decide the basis.. For example ring has curvy geometry...
    
    from math import sin, cos, fabs, acos, atan, tan, sqrt

    from . import so_translator_mod as tl
    
    piee = acos(-1.0)
    tpiee=2.0*piee
    
    tvector=[-inversion_point[0], -inversion_point[1], -inversion_point[2]]
    object_return=tl.translator(input_object, tvector)
    
    input_object=[]
    
    for i in range(len(object_return)):
        input_object.append(object_return[i])
    
    for i in range(len(input_object)):
        x=input_object[i][0]
        y=input_object[i][1]
        z=input_object[i][2]
        
        r=sqrt(x*x+y*y+z*z)
        theta=acos(z/r)
        rxy=sqrt(x*x+y*y)
        phi=acos(x/rxy)
        if (y<0.0):
            phi=tpiee-phi
            
        r=-r
        input_object[i][0]=r*sin(theta)*cos(phi)
        input_object[i][1]=r*sin(theta)*sin(phi)
        input_object[i][2]=r*cos(theta)
    
    tvector=[inversion_point[0], inversion_point[1], inversion_point[2]]
    
    object_return=tl.translator(input_object, tvector)
    
    
    input_object=[]
    
    for i in range(object_return):
        input_object.append(object_return[i])
        
    
    object_return=[]
    
    for i in range(input_object):
        object_return.append(input_object[i])
        
        
    return object_return
            
      
