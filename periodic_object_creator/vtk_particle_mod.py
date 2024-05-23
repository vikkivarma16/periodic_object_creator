# periodic_object_creator/vtk_particle

def particle_vis(input_data_ps, filename):
    
    basis=[]
    
    for i in range(len(input_data_ps)):
        if (input_data_ps[i][3] in basis):
            flag=1
        else:
            basis.append(input_data_ps[i][3])
    
    
    for i in range(len(basis)):
        count=0
        for j in range(len(input_data_ps)):
            if (input_data_ps[j][3]==basis[i]):
                count=count+1
        filename=filename+"particle_type_"+basis[i]+".vtk"
        fu=open(filename, "w")
        fu.write("# vtk DataFile Version 3.0\n")
        fu.write("Random data to test tensors\n")
        fu.write("ASCII\n")
        fu.write("DATASET POLYDATA\n")
        fu.write("POINTS  ")
        fu.write(str(count))
        fu.write("  float\n")
        for j in range(len(input_data_ps)):
	        if (input_data_ps[j][3]==basis[i]):
	            fu.write(str(input_data_ps[j][0]))
	            fu.write("  ")
	            fu.write(str(input_data_ps[j][1]))
	            fu.write("  ")
	            fu.write(str(input_data_ps[j][2]))
	            fu.write("\n")
	            
			
			
			
        
