

import example_object_simple_ring as esr
import vtk_contour as vc
import vtk_particle as vp

import so_rotator as rt
import so_translator as tl 
import so_cylindrical_wrapper as wrapper

n_beads=20
r_curvature=15
basis, species, bond_length, basis_bonds, basis_coordinate= esr.example_object_simple_ring(n_beads, r_curvature)



isometric_dif_factor=1

chain_radius=20.0

flag=0
input_object=[]
for j in range(len(basis_coordinate)):
    input_object.append(basis_coordinate[j])

print(basis_coordinate)

ro_axis_orien = [0.0, 1.0, 0.0]
ro_axis_posi = [0.0, 0.0, 0.0]
ro_degree = 90  # in degree
object_return = rt.rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)
if (isometric_dif_factor==0):

    for j in range(len(object_return)):
        basis_coordinate[j][0]=object_return[j][0]
        basis_coordinate[j][1]=object_return[j][1]
        basis_coordinate[j][2]=object_return[j][2]
    print("\n\n\n\n",object_return)

else:
    input_object=[]
    for j in range(len(object_return)):
        input_object.append(object_return[j])
    cylinder_radius=chain_radius 
    
    object_size=2.0*r_curvature+bond_length[0][2]
    object_return=wrapper.wrap(input_object, cylinder_radius, object_size)
    for j in range(len(object_return)):
        basis_coordinate[j][0]=object_return[j][0]
        basis_coordinate[j][1]=object_return[j][1]
        basis_coordinate[j][2]=object_return[j][2]
    





input_data_cs=[]
end_connection_flag=[]
input_data_cs.append(basis_coordinate)
    
end_connection_flag.append(1)
    
    
input_data_ps=[]
for i in range(len(basis_coordinate)):
    input_data_ps.append([basis_coordinate[i][0], basis_coordinate[i][1], basis_coordinate[i][2], basis[i]])
    
    
    
    
    
    
filename="Countour_single_ring_"        # without extension
value=vc.contour_vis(input_data_cs, end_connection_flag, filename)


filename="Single_ring_"
value=vp.particle_vis(input_data_ps, filename)


