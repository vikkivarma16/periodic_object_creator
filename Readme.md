
````markdown
# Periodic Object Creator

A Python package to create, manipulate, and visualize 3D and 2D objects with periodic properties.  

**Author:** Vikki Anand Varma  
**Email:** vikkivarma16@gmail.com  
**Bio:** PhD in Physics from IIT Delhi. Specializes in computational modeling, molecular simulations, and Python-based tools for 3D object manipulation, geometry operations, and visualization in scientific computing and materials modeling.

---

## Installation

Install directly from GitHub using pip:

```bash
pip install git+https://github.com/vikkivarma16/periodic_object_creator.git
````

Ensure that you have **Python 3.7+** and the necessary dependencies installed. Visualization requires **Paraview** or any **VTK-compatible viewer**.

---

## Description

`periodic_object_creator` is a comprehensive Python package for creating, manipulating, and visualizing 3D objects such as particles, sheets, and molecular structures. The package provides advanced geometric operations including:

* Rotation about arbitrary axes
* Translation along vectors
* Geometrical inversion about points
* Reflection across planes
* Wrapping sheets around cylinders and spheres
* Detection and construction of molecular topologies (bonds, angles, dihedrals, improper dihedrals)

It also provides utilities for visualization, overlap elimination, and element filtering.

---

## Input Structure

All functions operate on an `input_object`, which is a **list of elements**, each element being a list of coordinates and properties:

```python
[[x, y, z, id, type, other_properties], ...]
```

* `x, y, z` : spatial coordinates of a particle, atom, or point
* `type` : particle type or element symbol (e.g., `"C"` for carbon)
* `id` : unique identifier of the particle/atom
* `other_properties` : optional additional information (charge, mass, etc.)

**Minimum attributes required per element**: `x, y, z, type`.
Functions **operate on coordinates** while keeping other attributes unchanged unless explicitly modified.

---

## Functions and Their Functionalities

### `assign_group_ids(input_object, group_size=3, start_id=1, id_index=None)`

Assigns a group ID to consecutive blocks of `group_size` atoms.

**Parameters**:

* `input_object`: list of atom records `[[x, y, z, attr1, attr2, ...], ...]`
* `group_size`: number of atoms per group (e.g., 3 for H₂O)
* `start_id`: ID for the first group
* `id_index`: location to place ID:

  * `None` → append ID
  * integer → overwrite existing attribute

**Behavior**:

* Processes atoms sequentially
* Assigns the same group ID to every `group_size` consecutive atoms
* Automatically increments group IDs for each new group

**Output**: returns a new object with group IDs added or updated.

---

### `export_xyz(input_object, filename)`

Exports atomic coordinates and attributes:

* `.xyz` file containing coordinates
* `.txt` file containing all attributes

**Parameters**:

* `input_object`: list of elements
* `filename`: base name (without extension)

---

### `get_object_size(input_object, coord_indices=[0,1,2])`

Calculates the **size and extension** of a 3D object along x, y, and z axes.

**Input:**

* `input_object` : List of elements with coordinates, e.g., `[[x, y, z, ...], ...]`

**Output:**

Dictionary with `min`, `max`, and `size` for each axis:

```python
{
  'x': {'min':..., 'max':..., 'size':...},
  'y': {'min':..., 'max':..., 'size':...},
  'z': {'min':..., 'max':..., 'size':...}
}
```

**Example:**

```python
size_info = get_object_size(obj)
print(size_info)
# {'x': {'min':-1, 'max':1, 'size':2}, ...}
```

---

### `filter_broken_group(input_object, group_size=3, group_id_index=3)`

Removes incomplete or broken groups.

**Parameters**:

* `input_object`: list of atom records
* `group_size`: expected number of atoms per complete group
* `group_id_index`: index storing group ID

**Behavior**:

* Counts atoms per group ID
* Removes atoms from groups with incomplete counts

**Output**: filtered object with only complete groups.

---

### `build_topology(input_object, bond_length=0.96, tolerance=0.2, id_index=4, export_base="example_topo", many_body=False, id_body_index=None)`

Automatically generates molecular topology with optional multi-body support:

* Detects **bonds** within same body (bond length ± tolerance)
* For many-body cases, searches topology only within the body specified by `id_body_index`
* Exports combined topology data including **body IDs** and coordinates
* Builds **neighbor lists**
* Generates **bond angles** (triplets) within same body
* Generates **dihedrals** (quartets) within same body
* Generates **improper dihedrals** (≥3 neighbors) within same body
* Returns data structures as Python arrays for further processing

**Usage Example**:

```python
from topology_builder import build_topology

atoms = [
    [0.0, 0.0, 0.0, "O", 1, 1],
    [0.96, 0.0, 0.0, "H", 2, 1],
    [-0.24, 0.93, 0.0, "H", 3, 1],
]

topology = build_topology(
    atoms,
    bond_length=0.96,
    tolerance=0.2,
    id_index=4,
    export_base="example_topo",
    many_body=True,
    id_body_index=5
)

print(topology["bonds"])
print(topology["angles"])
print(topology["dihedrals"])
print(topology["impropers"])
```

---

#### `elements_picker(input_object, indices)`

Selects and returns a subset of elements from the input object based on provided indices.  
**Parameters**:  
- `input_object`: list of object elements `[[x, y, z, id, type, ...], ...]`  
- `indices`: list of integers specifying which elements to pick  

**Behavior**:  
- Preserves all attributes of the selected elements  
- Returns a new list containing only the chosen elements  

**Example**:

```python
subset = elements_picker(input_object, [0, 2, 5])
print(subset)
```

---

#### `inverter(input_object, inversion_point)`

Performs a geometrical inversion of the object with respect to a specified point.
**Parameters**:

* `input_object`: list of object elements
* `inversion_point`: 3D coordinates `[x, y, z]` around which inversion occurs

**Behavior**:

* Only coordinates are modified
* Other attributes remain unchanged

**Example**:

```python
inverted_obj = inverter(input_object, [0.0, 0.0, 0.0])
```

---

#### `overlap_eliminator(input_object_1, input_object_2, delete_from='obj1', tolerance=1e-6)`

Removes overlapping atoms between two objects within a distance tolerance.
**Parameters**:

* `input_object_1, input_object_2`: lists of object elements
* `delete_from`: `'obj1'` or `'obj2'` (specifies which object to remove overlaps from)
* `tolerance`: maximum distance to consider atoms overlapping

**Behavior**:

* Detects overlapping atoms using Euclidean distance
* Removes overlaps from the specified object
* Returns both objects with overlaps removed

**Example**:

```python
clean_obj1, clean_obj2 = overlap_eliminator(obj1, obj2, delete_from='obj1', tolerance=1e-6)
```

---

#### `particle_vis(input_object, filename)`

Generates a 3D visualization of particle objects and exports a VTK file.
**Parameters**:

* `input_object`: list of object elements with coordinates and attributes
* `filename`: name of the output `.vtk` file

**Behavior**:

* Uses size and type attributes for rendering
* File can be visualized in Paraview or other VTK-compatible software

**Example**:

```python
particle_vis(input_object, "particles_output.vtk")
```

---

#### `reflector(input_object, plane_normal, plane_location)`

Reflects an object across a plane defined by a point and a normal vector.
**Parameters**:

* `input_object`: list of elements
* `plane_normal`: 3D vector `[nx, ny, nz]` representing the plane normal
* `plane_location`: 3D point `[px, py, pz]` lying on the plane

**Behavior**:

* Coordinates are updated based on reflection formula
* Other attributes remain unchanged

**Example**:

```python
reflected_obj = reflector(input_object, [0, 0, 1], [0, 0, 0])
```

---

#### `replicator(input_object)`

Creates an exact copy of the input object.
**Parameters**:

* `input_object`: list of elements

**Behavior**:

* Preserves all attributes and coordinates
* Useful for creating multiple instances of the same structure

**Example**:

```python
copy_obj = replicator(input_object)
```

---

#### `rotator(input_object, ro_axis_orien, ro_axis_posi, ro_degree)`

Rotates the object around an arbitrary axis passing through a specified point.
**Parameters**:

* `ro_axis_orien`: 3D vector defining rotation axis `[x, y, z]`
* `ro_axis_posi`: 3D point through which axis passes `[x, y, z]`
* `ro_degree`: rotation angle in degrees

**Behavior**:

* Rotates all coordinates along the defined axis
* Preserves other attributes

**Example**:

```python
rotated_obj = rotator(input_object, [0, 0, 1], [0, 0, 0], 90)
```

---

#### `scissor(input_object, plane_origin, plane_normal, keep_side="negative")`

Cuts a 3D object along a plane, retaining only points on the specified side.
**Parameters**:

* `plane_origin`: point `[x, y, z]` on the cutting plane
* `plane_normal`: normal vector `[nx, ny, nz]` of the plane
* `keep_side`: `'negative'` or `'positive'`

**Behavior**:

* Uses plane equation to determine which side to keep
* Returns only the points on the chosen side

**Example**:

```python
scissored_obj = scissor(input_object, [0, 0, 0], [0, 0, 1], keep_side="positive")
```

---

#### `translator(input_object, tvector)`

Translates the object along a specified vector.
**Parameters**:

* `tvector`: 3D translation vector `[dx, dy, dz]`

**Behavior**:

* Adds vector to all coordinates
* Other attributes remain intact

**Example**:

```python
translated_obj = translator(input_object, [1.0, 0.0, 0.0])
```

---

#### `wrapper_cylindrical(input_object, cylinder_radius, object_size)`

Wraps a sheet-like object around a cylinder along the z-axis.
**Parameters**:

* `cylinder_radius`: radius of the cylinder
* `object_size`: length of the object along the wrapping direction

**Behavior**:

* Compresses along circumference to fit cylinder
* Assumes sheet’s center of mass is at origin
* Normal of object must align along y-axis
* Object size measured along x-axis

**Example**:

```python
wrapped_cyl = wrapper_cylindrical(input_object, cylinder_radius=5, object_size=2)
```

---

#### `wrapper_spherical(input_object, sphere_radius, object_size)`

Wraps a sheet-like object onto a spherical surface.
**Parameters**:

* `sphere_radius`: radius of the sphere
* `object_size`: approximate length along polar angle θ

**Behavior**:

* Compresses along polar angle while keeping radius fixed
* Assumes small patches to prevent distortion

**Example**:

```python
wrapped_sphere = wrapper_spherical(input_object, sphere_radius=10, object_size=2)
```

## Calling Methods

```python
from periodic_object_creator.assign_mol_id_mod import assign_group_ids
from periodic_object_creator.export_coordinate_particle_mod import export_xyz
from periodic_object_creator.export_bond_topology_mod import build_topology
from periodic_object_creator.export_object_size import get_object_size
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
```

**Usage Examples**:

```python
# Translate an object
translated_obj = translator(input_object, [1.0, 0.0, 0.0])

# Rotate an object
rotated_obj = rotator(input_object, [0,0,1], [0,0,0], 90)

# Wrap on a sphere
wrapped_obj = wrapper_spherical(input_object, sphere_radius=10, object_size=2)

# Invert around a point
inverted_obj = inverter(input_object, [0,0,0])

# Remove overlapping atoms
clean_obj1, clean_obj2 = overlap_eliminator(obj1, obj2, delete_from='obj1', tolerance=1e-6)

# Visualize particles
particle_vis(input_object, "particles_output.vtk")
```

---

## Notes

* Functions **preserve all attributes** (type, id, properties); only coordinates are modified.
* Wrapping functions assume **small patches** to prevent excessive compression or distortion.
* Visualization outputs can be viewed in **Paraview** or any VTK-compatible software.
* Always verify the **center of mass and orientation** before applying cylindrical or spherical wrapping.
* Group ID assignment, filtering, and topology generation ensure **consistency and completeness** of molecular structures.

---
