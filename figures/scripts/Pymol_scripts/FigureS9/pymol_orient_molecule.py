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
 
# --- DEBUT FONCTIONS ORIENT---      
def orient(obj='scy',
           mode='atoms',
           atom1=None,
           atom2=None,
           ring_sel=None,
           angle=0.0,
           target_axis=(0,0,1)):
    """
    Oriente l'objet pour que l'angle entre un vecteur défini
    (atoms ou ring) et target_axis soit égal à 'angle',
    en tournant AUTOUR DU CENTRE du vecteur.

    obj         : nom de l'objet dans PyMOL
    mode        : 'atoms' ou 'ring'
    atom1,atom2 : sélections PyMOL pour mode 'atoms'
    ring_sel    : sélection PyMOL des atomes de l'anneau pour mode 'ring'
    angle       : angle désiré (°)
    target_axis : vecteur cible (tuple de 3 floats)
    """
    # 1) recentrer + origin sur 0
    cmd.center(obj)
    cmd.origin(obj)

    # 2) calcul du vecteur v ET du pivot
    if mode == 'atoms':
        if not atom1 or not atom2:
            raise ValueError("Pour mode='atoms', atom1 et atom2 sont requis.")
        p1 = np.array(cmd.get_atom_coords(atom1))
        p2 = np.array(cmd.get_atom_coords(atom2))
        v = p2 - p1
        pivot = (p1 + p2) / 2.0      # centre du vecteur
    elif mode == 'ring':
        if not ring_sel:
            raise ValueError("Pour mode='ring', ring_sel est requis.")
        model = cmd.get_model(ring_sel)
        coords = np.array([a.coord for a in model.atom])
        centroid = coords.mean(axis=0)
        coords_centered = coords - centroid
        _, _, vh = np.linalg.svd(coords_centered)
        v = vh[2]                    # vecteur normal
        pivot = centroid            # centroïde de l'anneau
    else:
        raise ValueError("mode doit être 'atoms' ou 'ring'.")

    # 3) normalisation et angle courant
    v = v / np.linalg.norm(v)
    a = np.array(target_axis, dtype=float)
    a /= np.linalg.norm(a)
    cos_theta = np.clip(np.dot(v, a), -1.0, 1.0)
    theta_curr = np.degrees(np.arccos(cos_theta))

    # 4) axe de rotation
    rot_axis = np.cross(v, a)
    norm_ra = np.linalg.norm(rot_axis)
    if norm_ra < 1e-6:
    # vecteur parallèle à l’axe cible
    # si on veut tout de même un angle ≠ 0 ou 180, on prend X comme axe de rotation
    	if abs(theta_curr - angle) > 1e-3:
        	rot_axis = np.array([1.0, 0.0, 0.0])
        	delta    = theta_curr - angle
        	cmd.rotate(list(rot_axis), delta, obj, origin=list(pivot))
        	print(f"Rotation « dégénérée » de {delta:.2f}° autour de X pour créer l’angle voulu.")
    	else:
        	print("Angle déjà correct (0° ou 180°).")
    return
    rot_axis /= norm_ra

    # 5) delta pour atteindre l'angle désiré
    delta = theta_curr - angle

    # 6) rotation autour de 'pivot'
    cmd.rotate(list(rot_axis),
               delta,
               obj,
               origin=list(pivot))
    print(f"[orient] {obj} tourné de {delta:.2f}° autour de {pivot.tolist()}")
    
def orient2(obj='scy',
            atom1=None, atom2=None, angle_atoms=0.0,
            ring_sel=None,           angle_ring=0.0):
    """
    obj         : nom de l'objet PyMOL
    atom1,atom2: sélections pour le vecteur atomique
    angle_atoms: angle (°) voulu entre atom2–atom1 et l'axe Z
    ring_sel    : sélection des atomes de l'anneau aromatique
    angle_ring  : angle (°) voulu entre normale d'anneau et Z

    Rotation 1 (anneau) : autour de X (plan YZ)
    Rotation 2 (atoms): autour de Y (plan XZ)
    """
    # 1) centre + origine
    cmd.center(obj)
    cmd.origin(obj)

    # 2) normale et centroïde de l'anneau
    model = cmd.get_model(ring_sel)
    coords = np.array([a.coord for a in model.atom])
    cen_ring = coords.mean(axis=0)
    coords_c = coords - cen_ring
    _,_,vh = np.linalg.svd(coords_c)
    n_ring = vh[2] / np.linalg.norm(vh[2])

    # angle courant anneau
    cos_r = np.clip(np.dot(n_ring, [0,0,1]), -1.0, 1.0)
    th_ring = np.degrees(np.arccos(cos_r))
    delta_r = th_ring - angle_ring

    # rotation autour de l'axe X
    cmd.rotate([1,0,0], delta_r, obj, origin=list(cen_ring))
    print(f"[anneau] Δ={delta_r:.2f}° autour de X en {cen_ring}")

    # 3) vecteur atom1→atom2 et son milieu
    p1 = np.array(cmd.get_atom_coords(atom1))
    p2 = np.array(cmd.get_atom_coords(atom2))
    mid = (p1 + p2) / 2.0
    v_at = (p2 - p1)
    v_at /= np.linalg.norm(v_at)

    # angle courant atoms
    cos_a = np.clip(np.dot(v_at, [0,0,1]), -1.0, 1.0)
    th_at = np.degrees(np.arccos(cos_a))
    delta_a = th_at - angle_atoms

    # rotation autour de l'axe Y
    cmd.rotate([0,1,0], delta_a, obj, origin=list(mid))
    print(f"[atoms]  Δ={delta_a:.2f}° autour de Y en {mid}")

# --- FIN FONCTIONS ORIENT---

# register commands
cmd.extend('select_mol', select_mol)
cmd.extend('select_ring', select_ring)
cmd.extend('orient_ring_to_xy', orient_ring_to_xy)
cmd.extend('add_oxygen_plane', add_oxygen_plane)
cmd.extend('translate_mol', translate_mol)
cmd.extend('set_angles', set_angles)



