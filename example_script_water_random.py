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
from periodic_object_creator.so_overlap_remover_mod import  overlap_remover
from periodic_object_creator.so_overlap_eliminator_mod import overlap_eliminator
from periodic_object_creator.so_random_translator_mod import randomize_positions
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
#def overlap_remover(input_object, molid_idx, particle_idx, mol_type_idx, particle_type_idx, sigma_matrix, moving_mol_id, box, cell_size, iter_max, translation_step, rotation_step, max_particles_per_cell=64, grid_shifting_rate=100000):
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




# Filling up water molecules.

basis = [[0.0000000, 0.000000,  0.00000, "O", 1, 1], [0.8164904, 0.5773590, 0.00000, "H", 2, 1], [-0.8164904, 0.5773590, 0.00000, "H", 2, 1]]

#basis = [[0.0000000, 0.000000,  0.00000, "O", 1, 1]] #, [0.8164904, 0.5773590, 0.00000, "H", 2, 1], [-0.8164904, 0.5773590, 0.00000, "H", 2, 1]]

tvector  = [0.8164904, 0.0, 0.0]
basis = translator (basis, tvector)
box_size =  [100.0, 100.0, 100.0]

particle_size = 3.15

particle_volume  =  particle_size**3 * acos(-1)/6


packing_fraction  = 0.4
rho  =  packing_fraction/particle_volume

print("\n\nComputed density for the water packing fraction,", packing_fraction, "is given as:  ",   rho, "\n\n")

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


current_id = 1
group_size = 1
remove_existing_trailing_id=True
id_index = None
new_obj = assign_group_ids(water_lattic, group_size, current_id, id_index)   
water_lattic =  new_obj


current_id = 1
group_size = len (basis)
remove_existing_trailing_id=True
id_index = None
new_obj = assign_group_ids(water_lattic, group_size, current_id, id_index)   
water_lattic =  new_obj


random_object =  randomize_positions (water_lattic, idx  =  7,  box  = box_size,  seed=3840983, rotation =  False)


particle_vis(random_object, "water_random")
export_xyz (random_object, "water_random")



water = random_object


input_object  =  water


particle_idx=  6
mol_idx  =  7

mol_type_idx = 5
particle_type_idx = 5

sigma_matrix =  [[1.0, 1.0], [1.0, 1.0]]

moving_mol_id  = []
for i in range(len(input_object)):
    moving_mol_id.append (input_object [i][mol_idx])
    
box  =  box_size

cell_size =  1.0

iter_max  =  1000000
translation_step  =  0.05
rotation_step  =  0.01


water_final = overlap_remover(input_object, mol_idx, particle_idx, mol_type_idx, particle_type_idx, sigma_matrix, moving_mol_id, box, cell_size, iter_max, translation_step, rotation_step)



particle_vis(water_final, "water_final")
export_xyz (water_final, "water_final")

