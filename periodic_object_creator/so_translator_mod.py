# periodic_object_creator/so_translator

def translator(input_object, tvector):
    object_return = []
    for i in range(len(input_object)):
        dx = float(input_object[i][0]) + float(tvector[0])
        dy = float(input_object[i][1]) + float(tvector[1])
        dz = float(input_object[i][2]) + float(tvector[2])
        object_return.append([dx, dy, dz])
    return object_return
