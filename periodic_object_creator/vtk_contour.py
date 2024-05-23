# periodic_object_creator/vtk_contour

def contour_vis(input_data_cs, end_connection_flag, filename):
    filename=filename+".vtk"
    fu = open(filename, "w")
    
    
    import math
    from math import sin, cos, acos, atan, fabs
    count=0
    for i in range(len(input_data_cs)):
        count=count+len(input_data_cs[i])
    
    fu.write("# vtk DataFile Version 4.2\n")
    fu.write("Beads Data\n")
    fu.write("ASCII\n\n\n")
    fu.write("DATASET POLYDATA\n")
    fu.write("POINTS  ")
    fu.write(str(count))          
    fu.write("  float\n\n\n")

    i = 0

    for jt in range(len(input_data_cs)):

        for it in range(len(input_data_cs[jt])):
            i = i+1
            #print("here is the problem",input_data_cs[jt][it][0])
            if (fabs(input_data_cs[jt][it][0]) < 0.00000001):
                input_data_cs[jt][it][0] = 0.0

            if (fabs(input_data_cs[jt][it][1]) < 0.00000001):
                input_data_cs[jt][it][1] = 0.0

            if (fabs(input_data_cs[jt][it][2]) < 0.00000001):
                input_data_cs[jt][it][2] = 0.0

            fu.write(str(input_data_cs[jt][it][0]))
            fu.write("   ")
            fu.write(str(input_data_cs[jt][it][1]))
            fu.write("   ")
            fu.write(str(input_data_cs[jt][it][2]))
            fu.write("   ")
            fu.write("\n")
    fu.write("\n\n")   
    fu.write("LINES ")
    fu.write(str(len(input_data_cs)))
    fu.write("  ")
    
    total_contour_seg=0
    ind_cont=[]
    for i in range(len(end_connection_flag)):
        if (end_connection_flag[i]==0):
            total_contour_seg=total_contour_seg+(len(input_data_cs[i])-1)*2
            ind_cont.append((len(input_data_cs[i])-1)*2)
        else:
            total_contour_seg=total_contour_seg+(len(input_data_cs[i]))*2
            ind_cont.append((len(input_data_cs[i]))*2)
        total_contour_seg=total_contour_seg+1
    
    fu.write(str(total_contour_seg))
    fu.write("\n\n")
    

    i=1
    j=0
    for jt in range(len(input_data_cs)):
        fu.write(str(ind_cont[jt]))
        fu.write("  ")
        if (end_connection_flag[jt]==1):
            for it in range(len(input_data_cs[jt])):
                if (it==len(input_data_cs[jt])-1):
                    fu.write(str(i-1))
                    fu.write("  ")
                    fu.write(str(j))
                    fu.write("  ")
                    
                else:
                    fu.write(str(i-1))
                    fu.write("  ")
                    fu.write(str(i))
                    fu.write("  ")
                i=i+1
        else:
            for it in range(len(input_data_cs[jt])):
                fu.write(str(i-1))
                fu.write("  ")
                fu.write(str(i))
                fu.write("  ")
                i=i+1
            
        fu.write("\n")
        
        j=j+len(input_data_cs[jt])

    fu.close()
    
    return 0
