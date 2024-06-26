# periodic_object_creator/an_example_assignment
import math
from math import sin, cos, tan, asin, acos, atan


 

#from periodic_object_creator import example_object_sample_ring as sr

import periodic_object_creator.example_object_simple_ring_mod as sr
import periodic_object_creator.chain_validity_checker_mod as cvc
import periodic_object_creator.so_cylindrical_wrapper_mod as wrapper
import periodic_object_creator.so_rotator_mod as rt
import periodic_object_creator.so_translator_mod as tl
import periodic_object_creator.periodic_chain_forger_mod as pcf
import periodic_object_creator.vtk_particle_mod as vp 
import periodic_object_creator.vtk_contour_mod as vc
import periodic_object_creator.jason_file_creator_mod as js




# chain properties

chain_configuration_mode=0              # simplest chain made with horizontal and vertical rings
N_rings=20                              # keep it always even half for horizontal and half for vertical links
isometric_dif_factor=0                 # for the isometric deformed rings... wrapped around cylinder
chain_radius=4.0


# rings properties
n_beads=20  # minimum six beads are allowed
r_curvature=1.0


print("\n\n\n...given number of beads is...", n_beads, "\n\n\n")
print("\n\n\n... given ring radius is...", r_curvature, "\n\n\n")







if (N_rings % 2!=0):
    N_rings=N_rings+1




specific_factor=[]  # some special correction operation would be performed defined in the callable function
rings_coordinate=[]
contact_bond_length=[]
rings_size=[]
rings_particle_size=[]
flag=0
rings_bonds=[]
rings_ids=[]


# generates identical rings by calling example object ring function
id_counter=0
for i in range(N_rings):
     basis, species, bond_length, basis_bonds, basis_coordinate=sr.example_object_simple_ring(n_beads, r_curvature)
     rings_coordinate.append(basis_coordinate)
     rings_ids.append(basis)
     
     for j in range(len(basis_bonds)):
        for k in range(len(basis_bonds[j])):
            basis_bonds[j][k]=basis_bonds[j][k]+id_counter
     rings_bonds.append(basis_bonds)    
     rings_size.append(2.0*r_curvature+bond_length[0][2])
     id_counter=id_counter+n_beads
     
     # each rings have two contact points here the contact properties are set as default which also decides the size of the rings 
     contact_bond_length.append([bond_length[0][2], bond_length[0][2]])
     
     # Check if the number of beads and the radius of the ring is compatible... otherwise high bond lengths between the beads will not let the rings fit into each other to link..
     if (2.0*r_curvature<3.5*bond_length[0][2]):
        print("\n\n\n...bond length is given as", bond_length[0][2], "it represents a thick toroid which would not fit in each other to form the chain...\n\n")
        print("\n\n...beads are placed to too far increasing the bond length and making the torridal ring's inner radius smaller to fit the two rings in each link...\n\n")
        print("\n\n...please add the more number of beads...\n\n")
        print("\n\n...the process is aborted...\n\n")
        exit(0)
        
        
        
        

# particles_size is equal to the bond length or less than the bond length... if different than the bond length then ask for few more changes to implement validity checker properly
print("\n\n\nChosen particle's size is equal to contact bond length:", contact_bond_length[0][0])  # each ring have two contact points






# validity of the chain is being checked....
r_minimum, r_maximum = cvc.chain_validator(rings_size, contact_bond_length, isometric_dif_factor)






print("\n\n\nInput chain size is given as: ", chain_radius, "\n\n\n")




if (r_minimum> r_maximum):
    print("Warning:   ...incompatible configurations are given...Result: Broken or overlapping rings would be generated...")
    
if (chain_radius <r_minimum or chain_radius> r_maximum):
    print("Warning:   ...incompatible configurations are given...Result: Broken or overlapping rings would be generated...")
    chain_radius=0.5*(r_minimum+r_maximum)
    
    
    

print("\n\n\nAccepted chain size is given as: ", chain_radius, "\n\n\n")



# all the odd number rings are rotated vertically....


if (chain_configuration_mode==0):

    flag=0
    for i in range(len(rings_coordinate)):
        
        if (flag==0):
            flag=1
            
           
        else :
            flag=0
            input_object=[]
            for j in range(len(rings_coordinate[i])):
                input_object.append(rings_coordinate[i][j])
            
            ro_axis_orien = [0.0, 1.0, 0.0]
            ro_axis_posi = [0.0, 0.0, 0.0]
            ro_degree = 90  # in degree
            object_return = rt.rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)
            if (isometric_dif_factor==0):
            
                for j in range(len(object_return)):
                    rings_coordinate[i][j][0]=object_return[j][0]
                    rings_coordinate[i][j][1]=object_return[j][1]
                    rings_coordinate[i][j][2]=object_return[j][2]
        
            
            else:
                input_object=[]
                for j in range(len(object_return)):
                    input_object.append(object_return[j])
                cylinder_radius=chain_radius 
                
                object_size=rings_size[i]
                object_return=wrapper.wrap(input_object, cylinder_radius, object_size)
                for j in range(len(object_return)):
                    rings_coordinate[i][j][0]=object_return[j][0]
                    rings_coordinate[i][j][1]=object_return[j][1]
                    rings_coordinate[i][j][2]=object_return[j][2]
else:

    print("\n\n\n....Choose default configuration mode which is, 0 .... \n\n\n")
    print("\n\n\n...process aborted...\n\n\n")
        



# specific factor is used to displace the rings radially in case of overlap... you can do it explicitly and manually... 
# specific factor can also be used for other kind of operations...

specific_factor=[]
for i in range(len(rings_size)):
    specific_factor.append(0)
    
    
    
    
    
# here the final chain is being generated...    

chain_coordinate=pcf.periodic_ring_chain_forge(chain_radius, rings_coordinate, rings_size, chain_configuration_mode, specific_factor)


input_data_js=[]


id_counter=0
for i in range(len(chain_coordinate)):
    molecule=[]
    for j in range(len(chain_coordinate[i])):
        dum=[]
        dum.append(id_counter)
        dum.append(rings_ids[i][j])
        dum.append(chain_coordinate[i][j])
        dum.append(rings_bonds[i][j] )

        id_counter=id_counter+1
        molecule.append(dum)
    input_data_js.append(molecule)


filename="System.json"

value=js.json_writer(input_data_js, filename)





# it prints the contour lines for paraview... each contour line points are placed in a single list then the all the list are merged....
input_data_cs=[]
end_connection_flag=[]
for i in range(len(chain_coordinate)):
    input_data_cs.append(chain_coordinate[i])
    end_connection_flag.append(1)
    
filename="Countour"        # without extension
value=vc.contour_vis(input_data_cs, end_connection_flag, filename)


    
# it prints individual files with different species for paraview...
input_data_ps=[]
for i in range(len(chain_coordinate)):
    for j in range(len(chain_coordinate[i])):
        input_data_ps.append([chain_coordinate[i][j][0], chain_coordinate[i][j][1], chain_coordinate[i][j][2], rings_ids[i][j]])

filename="Chain_with_"
value=vp.particle_vis(input_data_ps, filename)





