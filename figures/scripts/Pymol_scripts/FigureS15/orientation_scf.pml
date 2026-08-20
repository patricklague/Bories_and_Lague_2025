#cd commande below should indicate the current directory
cd ./
run align_scf.py
set sphere_scale, 0.5
set orthoscopic, on
orient
run add_O_plane_10.py
#rotate z, 4
run pymol_orient_molecule.py

#theta1=angle_ring et theta2=angle_atoms
orient2('scf', atom1='name CG', atom2='name CZ', angle_atoms=37.851, ring_sel='name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ', angle_ring=91.826)

rotate z, 0
rotate x, -90
set sphere_transparency, 0.1
rotate y, 180

ray 1000
png ../../../plot/FigureS15/plan_scf_v1.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/FigureS15/plan_scf_v1_wo.png

show spheres, O_plane or O_plane2
rotate y, -90
ray 1000
png ../../../plot/FigureS15/plan_scf_v2.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/FigureS15/plan_scf_v2_wo.png

