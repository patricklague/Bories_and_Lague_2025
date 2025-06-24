from pymol import cmd
import numpy as np
import math

# ------------------ Simplified Script ------------------
# Usage in PyMOL:
# 1) run pymol_orient_molecule.py
# 2) select_mol scy.pdb  # charge and select molecule
# 3) select_ring "CG CD1 CD2 CE1 CE2 CZ"  # define aromatic ring selection
# 4) orient_ring_to_xy()  # align ring to XY plane
# 5) add_oxygen_plane 19.5  # add oxygen plane at z=19.5 Å
# 6) translate_mol 0 0 5  # translate molecule if needed
# 7) set_angles 60 30 "CG CD2"  # set theta1 and theta2 on current selections

# ------------------ Commands ------------------

def select_mol(pdb_path):
    """
    Load a PDB and select it as 'mol'.
    Usage: select_mol pdb_path
    """
    cmd.load(pdb_path, 'mol')
    cmd.select('mol_sel', 'mol')


def select_ring(atom_names='CG CD1 CD2 CE1 CE2 CZ'):
    """
    Define 'ring_sel' selection for aromatic ring atoms in 'mol_sel'.
    Usage: select_ring "CG CD1 CD2 CE1 CE2 CZ"
    """
    names = atom_names.split()
    expr = 'mol_sel and (' + ' or '.join(f'name {n}' for n in names) + ')'
    cmd.select('ring_sel', expr)


def orient_ring_to_xy():
    """
    Rotate 'mol_sel' so that 'ring_sel' lies in the XY plane.
    Usage: orient_ring_to_xy
    """
    model = cmd.get_model('ring_sel')
    if not model.atom:
        print('Error: ring_sel is empty')
        return
    pts = np.array([a.coord for a in model.atom])
    C = pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts - C)
    normal = vh[-1]
    target = np.array([0.0, 0.0, 1.0])
    axis = np.cross(normal, target)
    if np.linalg.norm(axis) < 1e-6:
        return
    axis = axis / np.linalg.norm(axis)
    angle = math.degrees(math.acos(np.clip(np.dot(normal, target), -1, 1)))
    cmd.rotate(axis.tolist(), angle, 'mol_sel', origin=C.tolist())


def add_oxygen_plane(z=19.5, spacing=2.0):
    """
    Create a grid of pseudo-oxygen atoms at z-plane.
    Usage: add_oxygen_plane z [spacing]
    """
    cmd.delete('O_plane*')
    # get bounding box of molecule
    min_bb, max_bb = cmd.get_extent('mol_sel')[:3], cmd.get_extent('mol_sel')[3:]
    x_vals = np.arange(min_bb[0]-spacing, max_bb[0]+spacing, spacing)
    y_vals = np.arange(min_bb[1]-spacing, max_bb[1]+spacing, spacing)
    for ix, x in enumerate(x_vals):
        for iy, y in enumerate(y_vals):
            cmd.pseudoatom(f'O_plane_{ix}_{iy}', pos=[x, y, z], elem='O')
    cmd.group('O_plane', 'O_plane_*')
    cmd.show('spheres', 'O_plane')
    cmd.set('sphere_scale', 0.5, 'O_plane')


def translate_mol(dx, dy, dz):
    """
    Translate the molecule selection by dx, dy, dz.
    Usage: translate_mol dx dy dz
    """
    cmd.translate([dx, dy, dz], 'mol_sel')


def set_angles(theta1, theta2, atom_pair='CG CD2'):
    """
    Set angles for current molecule:
      theta1 = angle(atom_pair, Z-axis)
      theta2 = angle(ring normal, Z-axis)
    Usage: set_angles theta1 theta2 "atom1 atom2"
    """
    # theta1
    a1, a2 = atom_pair.split()
    p1 = cmd.get_atom_coords(f'mol_sel and name {a1}')
    p2 = cmd.get_atom_coords(f'mol_sel and name {a2}')
    v = np.array(p2) - np.array(p1)
    v /= np.linalg.norm(v)
    axis1 = np.cross(v, [0,0,1])
    if np.linalg.norm(axis1) > 1e-6:
        curr1 = math.degrees(math.acos(np.clip(np.dot(v, [0,0,1]), -1, 1)))
        d1 = theta1 - curr1
        axis1 /= np.linalg.norm(axis1)
        cmd.rotate(axis1.tolist(), d1, 'mol_sel', origin='mol_sel')
    # theta2
    model = cmd.get_model('ring_sel')
    pts = np.array([a.coord for a in model.atom])
    C = pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts - C)
    normal = vh[-1]
    axis2 = np.cross(normal, [0,0,1])
    if np.linalg.norm(axis2) > 1e-6:
        curr2 = math.degrees(math.acos(np.clip(np.dot(normal, [0,0,1]), -1, 1)))
        d2 = theta2 - curr2
        axis2 /= np.linalg.norm(axis2)
        cmd.rotate(axis2.tolist(), d2, 'mol_sel', origin=C.tolist())

# register commands
cmd.extend('select_mol', select_mol)
cmd.extend('select_ring', select_ring)
cmd.extend('orient_ring_to_xy', orient_ring_to_xy)
cmd.extend('add_oxygen_plane', add_oxygen_plane)
cmd.extend('translate_mol', translate_mol)
cmd.extend('set_angles', set_angles)



