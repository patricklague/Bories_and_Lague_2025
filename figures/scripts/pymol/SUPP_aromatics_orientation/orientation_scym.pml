cd Dropbox/SC_article/Bories_and_Lague_2025/figures/scripts/pymol/SUPP_aromatics_orientation/
run align_scym.py
set sphere_scale, 0.5
set orthoscopic, on
orient
run add_O_plane_10.py
#rotate z, 4
run orient_mol.py      

#theta1=angle_ring et theta2=angle_atoms
orient2('scym', atom1='name CG', atom2='name CZ', angle_atoms=27.053, ring_sel='name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ', angle_ring=102.550)

rotate z, 4
rotate x, -90
set sphere_transparency, 0.1

ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scym_v1.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scym_v1_wo.png

show spheres, O_plane or O_plane2
rotate y, 90
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scym_v2.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scym_v2_wo.png

