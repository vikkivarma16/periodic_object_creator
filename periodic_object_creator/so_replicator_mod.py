# periodic_object_creator/so_replicator

# periodic_object_creator/so_replicator_mod.py

def replicator(input_object):
    return [p[:] for p in input_object]  # deep copy preserving all attributes

