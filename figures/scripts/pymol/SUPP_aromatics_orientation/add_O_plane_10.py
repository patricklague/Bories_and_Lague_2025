from pymol import cmd
import math

# 5) Create a grid (plane) of oxygen pseudoatoms at Z = 10 Å
# Adjust grid_spacing and grid_extent as needed for density and size
grid_spacing = 1  # Å between oxygens
grid_extent = 10   # maximum absolute X/Y coordinate
for x in range(-grid_extent, grid_extent+1, grid_spacing):
    for y in range(-grid_extent, grid_extent+1, grid_spacing):
        cmd.pseudoatom('O_plane2', elem='O', resn='PLN', pos=(x, y, 10))
# Show and color the oxygen plane
cmd.show('spheres', 'O_plane2')
cmd.color('blue', 'O_plane2')

# End of script




