cd Dropbox/SC_article/Bories_and_Lague_2025/figures/scripts/pymol/SUPP_aromatics_orientation/
run align_scw.py
set sphere_scale, 0.5
set orthoscopic, on
orient
run add_O_plane_10.py
#rotate z, 4
run orient_mol.py

#theta1=angle_ring et theta2=angle_atoms
orient2('scw', atom1='name CZ3', atom2='name CE2', angle_atoms=41.368, ring_sel='name CD2 or name CZ2 or name CZ3 or name CD1 or name CD2 or name CH2 or name CE3', angle_ring=88.241)

rotate z, 2
rotate x, -90
set sphere_transparency, 0.1

ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scw_v1.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scw_v1_wo.png

show spheres, O_plane or O_plane2
rotate y, 90
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scw_v2.png
hide everything, O_plane or O_plane2
ray 1000
png ../../../plot/SUPP_aromatics_orientation/plan_scw_v2_wo.png

