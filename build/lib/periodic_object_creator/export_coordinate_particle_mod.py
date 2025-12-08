# periodic_object_creator/so_exporter_mod.py

def export_xyz(input_object, filename):
    """
    Exports coordinates + element type into a standard .xyz file.
    The first line = number of atoms
    The second line = comment
    Then: type  x  y  z
    """

    with open(filename+".xyz", "w") as f:
      
        for item in input_object:
            element = item[3] if len(item) > 3 else "X"
            x, y, z = item[0], item[1], item[2]
            f.write(f" {x:.10f}  {y:.10f}  {z:.10f}\n")


    with open(filename+".txt", "w") as f:
        f.write("# Full attribute export\n")
        f.write("# Format: index  x  y  z  attr1  attr2 ...\n\n")

        for i, item in enumerate(input_object):
            line = f"{i+1}  "
            line += " ".join(str(v) for v in item)
            f.write(line + "\n")

