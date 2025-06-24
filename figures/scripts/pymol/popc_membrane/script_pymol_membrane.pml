############################################################
#
#  PyMOL script to generate a view of a POPC bilayer 
#  immersed in water, similar to the provided image.
#
#
############################################################

# -------------- 1) Loading the structure --------------
load popc_bilayer.pdb, membrane

# Set white background
bg_color white

# Hide all default representations
hide everything

# -------------- 2) Water molecule representation --------------
# Assumes water molecules are named HOH (or WAT)
# and contain “O” atoms for oxygen and “H” for hydrogen.
# You can adjust the selector if your file uses a different residue name.

# Select all water molecules
select waters, resn TIP3

# Show water as spheres
show spheres, waters

# Adjust the water sphere size (tweak if needed)
set sphere_scale, 0.5, waters

# Color water oxygens red
color red, (waters and name O)

# Color water hydrogens white (if present)
color white, (waters and name H*)

# -------------- 3) Hydrocarbon tail representation --------------
# Select tails: all non-hydrogen atoms in POPC lipids
select lipid_tails, resn POPC and not (name H*)

# Show tails as sticks
show sticks, lipid_tails
set sticks_radius, 0.3, lipid_tails

# Color the chains gold
# ("gold" is a standard PyMOL color; you could also use "yellow" or RGB)
color gray60, lipid_tails
# "gray60" is also a nice color

# -------------- 4) Phosphate / choline headgroup representation --------------
# Show P and N atoms of the headgroups to highlight the interface

# Select phosphate (P) atoms in POPC
select P_atoms, resn POPC and name P

# Show P as spheres and color orange
show spheres, P_atoms
set sphere_scale, 0.9, P_atoms
color orange, P_atoms

# Select nitrogen (N) atoms in POPC (choline group)
select N_atoms, resn POPC and name N

# Show N as spheres and color blue
show spheres, N_atoms
set sphere_scale, 0.8, N_atoms
color blue, N_atoms

# -------------- 5) General display options --------------
# Disable residue/atom labels
#hide labels

# Disable specular lighting for a more matte finish
#set specular, 0

# Adjust quality of spheres and sticks
#set sphere_quality, 1        # from 1 (low) to 6 (high)
#set line_smooth, on
#set stick_quality, 15

# -------------- 6) Camera positioning --------------
# We want a side view of the bilayer (edge-on):
#  - First, orient the whole object
orient membrane

#  - Then rotate 90° around X axis to view bilayer edge
rotate x, 90
rotate x, 40
# You can fine-tune the angle (e.g. rotate y, angle) if needed:
#rotate y, 60
#rotate z, 30

# Center the view and zoom in (adjust value if needed)
zoom membrane, -30

# -------------- 7) Ray-tracing settings --------------
# Advanced tracing mode to avoid “black edges” when rendering
set ray_trace_mode, 0

# Anti-aliasing (smooth edges)
#set antialias, 2

# Depth cueing (to enhance 3D effect) -- optional
#set depth_cue, 1
#set depth_cue_power, 1

# -------------- 8) Final render (ray) --------------
# Image dimensions (adjust as needed, e.g., 1920×1080)
ray 1000

# Save as high-resolution PNG (300 dpi)
png popc_render.png, dpi=1000

# End of script
