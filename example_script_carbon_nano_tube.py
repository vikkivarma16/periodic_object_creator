# periodic_object_creator/an_example_assignment
import math
from math import sin, cos, tan, asin, acos, atan


 
from periodic_object_creator.assign_mol_id_mod import assign_group_ids
from periodic_object_creator.so_cm_calculator_mod import cm_calculator
from periodic_object_creator.export_coordinate_particle_mod import export_xyz
from periodic_object_creator.export_bond_topology_mod import  build_topology
from periodic_object_creator.export_object_size_mod import get_object_size
from periodic_object_creator.so_elements_picker_mod import elements_picker
from periodic_object_creator.so_inverter_mod import inverter
from periodic_object_creator.so_overlap_eliminator_mod import overlap_eliminator
from periodic_object_creator.so_reflector_mod import reflector
from periodic_object_creator.so_replicator_mod import replicator
from periodic_object_creator.so_rotator_mod import rotator
from periodic_object_creator.so_scissor_mod import scissor
from periodic_object_creator.so_translator_mod import translator
from periodic_object_creator.so_wrapper_cylindrical_mod import wrapper_cylindrical
from periodic_object_creator.so_wrapper_spherical_mod import wrapper_spherical
from periodic_object_creator.vtk_particle_mod import particle_vis

#def assign_group_ids(obj, group_size=3, start_id=1, id_index_in_element=None):
#def export_xyz (cnt, "cord_cnt"):
#def get_object_size(input_object, coord_indices=(0,1,2)):
#def inverter(input_object, inversion_point):
#def build_topology(input_object, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="topology_data", many_body=True, id_body_index = 6):
#def overlap_eliminator(input_object_1, input_object_2, delete_from='obj1', tolerance=1e-6):
#def particle_vis(input_data_ps, filename):
#def picker(input_object, indices):
#def reflector(input_object, plane_normal, plane_location):
#def replicator(input_object):
#def rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree):
#def scissor(input_object, plane_origin, plane_normal, keep_side="negative"):
#def translator(input_object, tvector):
#def wrapper_cylindrical(input_object, cylinder_radius, object_size):
#def wrapper_spherical(input_object, sphere_radius, object_size):





bond_length =  1.42  # in Angstrom unit

lattice_constant =  2.46
unit_cells_x  =  10  # repetition along the x axis which is the first unit vector
unit_cells_y  =  10  # repetition along the other unit vector




basis_object =  [[0.0, 0.0, 0.0, "C"]  ] # this example contains just one element where the element contain coordinate and the particle type id minimum it can have more attributes
input_object  =   basis_object
tvector  =  [bond_length, 0.0, 0.0]
returned_object =  translator (input_object, tvector)


# A Hexagon is formed by rotating a point 6 times through 60 degree around the origin
basis_object_1 =  []
ro_axis_posi = [0.0, 0.0, 0.0]
ro_axis_orien = [0.0, 0.0, 1.0]
ro_degree = 60 
input_object = returned_object
for i in range (6):
    new_object  = rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)
    ro_degree =  ro_degree+60
    #for j in range (len(new_object)):
    basis_object_1.extend(new_object)
basis_object_1 =  rotator(basis_object_1, ro_axis_orien, ro_axis_posi, 30)







#Now the same hexagon is translated across x axis through the lattic constant to get a linear array of hexagones.

tvector  =  [lattice_constant, 0.0, 0.0]

basis_object_2 = []

old_object = basis_object_1
basis_object_2.extend(old_object)

for i in range(unit_cells_x):

    new_object = translator(old_object, tvector)

    remaining_object = overlap_eliminator(old_object, new_object, delete_from='input_object_2',   tolerance=0.1)
    old_object = new_object
    basis_object_2.extend(remaining_object[1])
    
    
    
# creating sheet of hexagons with coordinates like cnt


tvector  =  [lattice_constant, 0.0, 0.0]
pi  =  acos(-1.0000)
dx  =  cos(pi/3)
dy  =  sin (pi/3)
dz  = 0

tvector_1  =  [dx*lattice_constant, dy*lattice_constant, dz]

dx  =  cos(pi*2/3)
dy  =  sin(pi*2/3)
dz = 0
tvector_2  =  [dx*lattice_constant, dy*lattice_constant, dz]


basis_object_3 = []
old_object = basis_object_2
basis_object_3.extend(old_object)
flag = 0
for i in range(unit_cells_y):
    
    if (flag  == 0):
        tvector =  tvector_1
        flag  = 1
    else:
        tvector  = tvector_2
        flag =0
    
    new_object = translator(old_object, tvector)

    remaining_object = overlap_eliminator(
        old_object, new_object,
        delete_from='input_object_2',
        tolerance=0.1
    )

  
    old_object = new_object

    # FIXED
    basis_object_3.extend(remaining_object[1])
#particle_vis(basis_object_3, "basis_object_3")




# wrapping procedure along the cylinder 
object_size  =  (unit_cells_x+2)*lattice_constant
cylinder_radius  =  object_size/(2*pi)



# making normal along the y axis
ro_axis_posi = cm_calculator(basis_object_3)
ro_axis_orien = [1.0, 0.0, 0.0]
ro_degree = 90

basis_object_4  = rotator(basis_object_3, ro_axis_orien, ro_axis_posi, ro_degree)
cm = cm_calculator(basis_object_4)




# nutralizing the position of the center of mass of the sheet
tvec =  [-cm[0], -cm[1], -cm[2]]
basis_object_5  = translator(basis_object_4, tvec)
#particle_vis(basis_object_5, "basis_object_5")






cnt  =  wrapper_cylindrical(basis_object_5, cylinder_radius, object_size)

cm = cm_calculator(cnt)




# nutralizing the position of the center of mass of the sheet
tvec =  [-cm[0], -cm[1], -cm[2]]
cnt_new  = translator(cnt, tvec)

cnt  =  cnt_new

particle_vis(cnt, "cnt")
export_xyz (cnt, "cord_cnt")

current_id = 1
group_size = len(basis_object)
id_index = None
new_obj = assign_group_ids(cnt, group_size, current_id, id_index)

cnt = new_obj
bond_length =  1.42 
tolerance  = 0.2 



build_topology(cnt, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="cnt_full_topology", many_body=False)

size  = get_object_size(cnt, coord_indices=(0,1,2))
print (size)


current_id = 1
group_size = len(basis_object)
remove_existing_trailing_id=True
id_index = None
new_obj = assign_group_ids(basis_object_1, group_size, current_id, id_index)

hexagon = new_obj
bond_length =  1.42 
tolerance  = 0.2 

build_topology(hexagon, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="hexagon_full_topology", many_body=False)


