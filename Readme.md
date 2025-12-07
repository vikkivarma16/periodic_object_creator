# Periodic object creator

A Python package to create 3D and 2D objects with periodic properties....
Author: Vikki Anand Varma
Email: vikkivarma16@gmail.com
Bio: PhD in Physics from IIT Delhi. Vikki Anand Varma specializes in computational modeling, molecular simulations, and development of Python-based tools for 3D object manipulation and visualization.



## Installation

```bash
pip install git+https://github.com/vikkivarma16/periodic_object_creator.git


# Periodic Object Creator

**Description:**  
`periodic_object_creator` is a Python package for manipulating 3D objects (particles, sheets, or molecular structures) using geometric operations such as rotation, translation, inversion, reflection, wrapping around cylinders or spheres, and more. The package also includes visualization utilities.

---

# Input Structure

All functions accept an `input_object` which is a list of elements. Each element is a list containing coordinates and additional properties. For example:

```python
[[x, y, z, id, type other_properties], ...]
```

- `x, y, z` : coordinates of the particle/atom/point  
- `type` : particle type or element symbol (e.g., `"C"` for carbon)  
- `id` : unique identifier of the coordinate 
- `other_properties` : optional additional information  

while the minimum attributes are "x, y, z, type" for each elements in the objects.

All functions operate on the coordinates but keep the other attributes unchanged unless specified.

---

## Functions and Their Functionalities

assign_group_ids(input_object, group_size=3, start_id=1, id_index=None)

Assigns a group ID to every consecutive block of `group_size` atoms.

**Parameters**
- `input_object` – list of atom records: `[[x, y, z, attr1, attr2, ...], ...]`
- `group_size` – number of atoms per group (e.g., 3 for H₂O)
- `start_id` – ID assigned to the first group
- `id_index` – where to place the ID:
  - `None` → append ID as the last element  
  - integer (e.g., 4) → overwrite/update ID at that index

**Behavior**
Atoms are processed in order. Every `group_size` atoms receive the same ID.  
IDs automatically increment for each new group.  
If `id_index` is given, the function updates the ID at that index;  
if not, the ID is appended.

**Output**
Returns a new object with group IDs added or updated for each atom.



export_xyz(input_object, filename)
Writes the atomic coordinates to an .xyz file and all additional attributes to a .txt file using the given filename.
input_object – list of elements, each in the form -[[x, y, z, attr1, attr2, ...],...], which is your input object
filename – base filename (without extension)

Output:
filename.xyz containing only coordinates
filename.txt containing element attributes



filter_broken_group(input_object, group_size=3, group_id_index=3)

Removes all atoms belonging to incomplete groups (broken molecules).

**Parameters**
- `input_object` – list of atom records: `[[x, y, z, attr1, attr2, ...], ...]`
- `group_size` – expected number of atoms per complete group
- `group_id_index` – index where the group ID is stored in each atom record

**Behavior**
Counts how many atoms share each group ID.  
Any ID that does not appear exactly `group_size` times is considered broken.  
All atoms belonging to broken groups are removed.

**Output**
Returns a filtered object containing only atoms from complete, intact groups.



build_topology( input_object, bond_length=0.96, tolerance=0.2, id_index=4, export_base="example_topo" )


This Python module automatically generates molecular topology from atomic coordinates.  
It detects **bonds, angles, dihedrals, and improper dihedrals** purely from geometry, and exports a complete topology file with coordinates.

---

## Features

- Detects **bonds** using a bond length and tolerance.
- Builds **neighbor lists** for each atom.
- Generates **bond angles** (triplets) from neighbors.
- Generates **dihedrals** (quartets) from bonded chains.
- Generates **improper dihedrals** for atoms with ≥3 neighbors (planarity enforcement optional).
- Exports **full topology** to a text file including coordinates.
- Returns all data structures as Python arrays for further processing.

---

## Usage


from topology_builder import build_topology

# Example: simple water molecule
atoms = [
    [0.0, 0.0, 0.0, "O", 1],
    [0.96, 0.0, 0.0, "H", 2],
    [-0.24, 0.93, 0.0, "H", 3],
]





elements_picker(input_object, indices)
Returns a subset of the input object containing only the elements at the specified indices.
Input:
input_object: list of object elements.
indices: list of integer indices to pick.
Output:
A list of selected elements with all attributes preserved.

inverter(input_object, inversion_point)
Performs geometrical inversion of the object through a specified point. Only coordinates are changed; other attributes remain unchanged.
Input:
input_object: list of object elements.
inversion_point: 3D point [x, y, z] about which inversion is performed.

overlap_eleminator(input_object_1, input_object_2, delete_from='obj1', tolerance=1e-6)
Deletes coordinates in one object that overlap with another object within a specified distance tolerance.
Input:
input_object_1, input_object_2: lists of object elements.
delete_from: string, either 'obj1' or 'obj2', specifies which object to remove overlapping elements from.
tolerance: float, distance threshold to consider elements as overlapping. After deletion, returns both the objects

particle_vis
Provides basic visualization of particle objects in 3D. Accepts a list of elements and renders them using size and type attributes.

reflector(input_object, plane_normal, plane_location)
Reflects a 3D object across a plane defined by a point and a normal vector. Only coordinates are changed; other attributes are preserved.
Input:
input_object: list of object elements.
plane_normal: 3D vector [nx, ny, nz] representing the plane normal.
plane_location: 3D point [px, py, pz] lying on the plane.

replicator(input_object)
Returns a copy of the input object with all coordinates and attributes intact.

rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)
Rotates the object around an arbitrary axis passing through a given point. Other attributes remain unchanged.
Input:
ro_axis_orien: 3D vector defining the rotation axis.
ro_axis_posi: 3D point through which the rotation axis passes.
ro_degree: rotation angle in degrees.

scissor(input_object, plane_origin, plane_normal, keep_side="negative")
Cuts a 3D object using a plane. Retains only the points on the specified side of the plane.
Input:
plane_origin: a point [x, y, z] on the cutting plane.
plane_normal: normal vector [nx, ny, nz] of the plane.
keep_side: 'negative' (default) or 'positive'.

translator(input_object, tvector)
Translates an object by a specified vector. All other attributes remain intact.
Input:
tvector: 3D translation vector [dx, dy, dz].

wrapper_cylindrical(input_object, cylinder_radius, object_size)
Wraps a sheet-like object around a cylinder along the z-axis. Compression is applied along the cylinder’s circumference.
Input:
cylinder_radius: radius of the cylinder.
object_size: length of the object along the wrapping direction.
Caution: Sheet object's center of mass must be at origin and normal of the object must be along y axis and the object size must be measured as the extension of the same in x axis.

wrapper_spherical(input_object, sphere_radius, object_size)
Wraps a sheet-like object onto a spherical surface. Compression is applied along the polar angle (θ) while keeping radius fixed.
Input:
sphere_radius: radius of the sphere.
object_size: approximate length of the object along the polar angle.

---

## Calling Methods

You can import and call these functions as follows:

```python
from periodic_object_creator.assign_mol_id_mod import assign_group_ids
from periodic_object_creator.export_coordinate_particle_mod import export_xyz
from periodic_object_creator.export_bond_topology_mod import  build_topology
from periodic_object_creator.filter_broken_mol_mod import filter_broken_group
from periodic_object_creator.so_cm_calculator_mod import cm_calculator
from periodic_object_creator.so_elements_picker_mod import elements_picker
from periodic_object_creator.so_inverter_mod import inverter
from periodic_object_creator.so_overlap_eleminator_mod import overlap_eleminator
from periodic_object_creator.so_reflector_mod import reflector
from periodic_object_creator.so_replicator_mod import replicator
from periodic_object_creator.so_rotator_mod import rotator
from periodic_object_creator.so_scissor_mod import scissor
from periodic_object_creator.so_translator_mod import translator
from periodic_object_creator.so_wrapper_cylindrical_mod import wrapper_cylindrical
from periodic_object_creator.so_wrapper_spherical_mod import wrapper_spherical
from periodic_object_creator.vtk_particle_mod import particle_vis
```

Example:

```python
# Translate an object
translated_obj = translator(input_object, [1.0, 0.0, 0.0])

# Rotate an object
rotated_obj = rotator(input_object, [0,0,1], [0,0,0], 90)

# Wrap on a sphere
wrapped_obj = wrapper_spherical(input_object, sphere_radius=10, object_size=2)
```

---

## Notes

- All functions preserve additional attributes of the input object (type, id, properties). Only the coordinates are modified.  
- Wrapping functions assume small patches to prevent excessive compression.  
- Visualization is provided through `particle_vis`.



