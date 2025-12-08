# periodic_object_creator/an_example_assignment
import math
from math import sin, cos, tan, asin, acos, atan



from periodic_object_creator.assign_mol_id_mod import assign_group_ids
from periodic_object_creator.export_bond_topology_mod import  build_topology
from periodic_object_creator.export_coordinate_particle_mod import export_xyz
from periodic_object_creator.filter_broken_mol_mod import filter_broken_group
from periodic_object_creator.so_cm_calculator_mod import cm_calculator
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
#def filter_broken_group(remaining_water_lattic, group_size, group_id_indices ):
#def inverter(input_object, inversion_point):
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
# particle_vis(cnt, "cnt")
# export_xyz (cnt, "cord_cnt")


current_id = 1
group_size = len(basis_object)
remove_existing_trailing_id=True
id_index = None
new_obj = assign_group_ids(cnt, group_size, current_id, id_index)

cnt = new_obj
bond_length =  1.42 
tolerance  = 0.2 
build_topology(cnt, bond_length, tolerance, id_index=None, coord_indices=(0,1,2), export_base="cnt_full_topology", many_body=False)








# Filling up water molecules.

basis = [[0.0000000, 0.000000,  0.00000, "O", 1], [0.8164904, 0.5773590, 0.00000, "H", 1], [-0.8164904, 0.5773590, 0.00000, "H", 1]]
tvector  = [0.8164904, 0.0, 0.0]
basis = translator (basis, tvector)
box_size =  [30.0, 30.0, 30.0]

rho  =  0.05

N = rho * box_size[0] * box_size[1] * box_size[2]

N_x  = int ( (N**(1.0/3.0)) * box_size[0]/ (box_size[0] * box_size[1] * box_size[2])**(1.0/3.0))
N_y  = int ( (N**(1.0/3.0)) * box_size[1]/ (box_size[0] * box_size[1] * box_size[2])**(1.0/3.0))
N_z  = int ( (N**(1.0/3.0)) * box_size[2]/ (box_size[0] * box_size[1] * box_size[2])**(1.0/3.0))


lattic_constant_x = box_size[0]/ N_x
lattic_constant_y = box_size[1]/ N_y
lattic_constant_z = box_size[2]/ N_z

basis_object_1 = []
tvector =  [lattic_constant_x, 0.0, 0.0]
old_object  =  basis
basis_object_1.extend(old_object)
for i in range(N_x-1):
    new_object  =  translator(old_object, tvector)
    old_object  =  new_object
    basis_object_1.extend(old_object)
    
basis_object_2 = []
tvector =  [0.0, lattic_constant_y, 0.0]
old_object  =  basis_object_1
basis_object_2.extend(old_object)
for i in range(N_x-1):
    new_object  =  translator(old_object, tvector)
    old_object  =  new_object
    basis_object_2.extend(old_object)
    
basis_object_3 = []
tvector =  [0.0, 0.0,  lattic_constant_z]
old_object  =  basis_object_2
basis_object_3.extend(old_object)
for i in range(N_x-1):
    new_object  =  translator(old_object, tvector)
    old_object  =  new_object
    basis_object_3.extend(old_object)




water_lattic  =  basis_object_3    


cm = cm_calculator(cnt)
tvector = [-cm[0], -cm[1], -cm[2]]
cnt_new  =  translator(cnt, tvector)

new_tvector  =  [box_size[0]/2.0, box_size[1]/2.0, box_size[2]/2.0]

dispersed_cnt =  translator(cnt_new, new_tvector)

# giving unique id to each basis in the object
current_id = 1
group_size = len(basis)
remove_existing_trailing_id=True
id_index = None
new_obj = assign_group_ids(water_lattic, group_size, current_id, id_index)   
water_lattic =  new_obj



# merging cnt with water and deleting water molecules which are overlapping with the cnt

remaining_object = overlap_eliminator(
        dispersed_cnt, water_lattic,
        delete_from='input_object_2',
        tolerance=1.0
    )

remaining_water_lattic  =  remaining_object[1]


group_id_indices = 5

filtered = filter_broken_group(remaining_water_lattic, group_size, group_id_indices )
remaining_water_lattic= filtered







particle_vis(remaining_water_lattic, "water")
export_xyz (remaining_water_lattic, "cord_water")
particle_vis(dispersed_cnt, "dispersed_cnt")
export_xyz (dispersed_cnt, "dispersed_cnt")
    


