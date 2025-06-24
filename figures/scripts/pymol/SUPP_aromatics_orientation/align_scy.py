from pymol import cmd
import math

# ----- Edit these selections -----
# Name to assign to your molecule object
mol_name = 'scy'
# Define the two atoms for the oriented vector (using PyMOL selection syntax)
atom1 = 'name CG'
atom2 = 'name CZ'
# Define the first three atoms of the aromatic ring to compute its plane normal
ring_atoms = ['name CG', 'name CD1', 'name CD2', 'name CE1', 'name CE2', 'name CZ']
# ----------------------------------

# 1) Load and center the molecule
cmd.load('./pdb/scy.pdb', mol_name)
cmd.center(mol_name)

# Helper functions
def get_vector(sel1, sel2):
    c1 = cmd.get_atom_coords(sel1)
    c2 = cmd.get_atom_coords(sel2)
    return [c2[i] - c1[i] for i in range(3)]

def normalize(v):
    norm = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return [v[i]/norm for i in range(3)]

def cross(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1],
            v1[2]*v2[0] - v1[0]*v2[2],
            v1[0]*v2[1] - v1[1]*v2[0]]

def rotate_to_align(vec, target=(0, 0, 1), selection=mol_name):
    v = normalize(vec)
    t = normalize(list(target))
    # Compute rotation axis
    axis = cross(v, t)
    axis_len = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    if axis_len < 1e-6:
        return
    axis = [axis[i]/axis_len for i in range(3)]
    # Compute rotation angle
    angle = math.degrees(math.acos(sum(v[i]*t[i] for i in range(3))))
    # Apply rotation around axis at origin
    cmd.rotate(axis, angle, selection, origin=(0, 0, 0))

# 2) Align aromatic ring plane with the XY-plane
p1 = cmd.get_atom_coords(ring_atoms[0])
p2 = cmd.get_atom_coords(ring_atoms[1])
p3 = cmd.get_atom_coords(ring_atoms[2])
v1 = [p2[i] - p1[i] for i in range(3)]
v2 = [p3[i] - p1[i] for i in range(3)]
normal = cross(v1, v2)
rotate_to_align(normal)

# 3) Align the oriented vector between atom1 and atom2 with the Z-axis
vec12 = get_vector(atom1, atom2)
rotate_to_align(vec12)

# 4) Translate the molecule so its center lies at Z = 11.398 Å A CHANGER SELON LE RESIDU
cmd.translate([0, 0, 11.398], mol_name)

# 5) Create a grid (plane) of oxygen pseudoatoms at Z = 20 Å
# Adjust grid_spacing and grid_extent as needed for density and size
grid_spacing = 1  # Å between oxygens
grid_extent = 10   # maximum absolute X/Y coordinate
for x in range(-grid_extent, grid_extent+1, grid_spacing):
    for y in range(-grid_extent, grid_extent+1, grid_spacing):
        cmd.pseudoatom('O_plane', elem='O', resn='PLN', pos=(x, y, 20))
# Show and color the oxygen plane
cmd.show('spheres', 'O_plane')
cmd.color('red', 'O_plane')

cmd.extend('rotate_to_align', rotate_to_align)

# End of script




