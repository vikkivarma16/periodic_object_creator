#periodic_object_creator/elements_picker


def picker(input_object, indices):

    return_object=[]
    for i in range(len(input_object)):
        dum=[]
        for j in range(len(indices)):
            dum.append(object[i][indices[j]])
            
        return_object.append(dum)


    return return_object
    
    
    
    
