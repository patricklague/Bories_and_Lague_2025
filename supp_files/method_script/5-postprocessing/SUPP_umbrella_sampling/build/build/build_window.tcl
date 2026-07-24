# build_window.tcl
#
# Build a single umbrella-sampling window:
#   - Take the equilibrated POPC bilayer (membrane/charmm-gui/namd/step5_input.{psf,pdb})
#   - Insert two copies of one analog (SCC / SCD / SCR) at z = z1 and z = z1 + 37 A
#   - Remove waters / ions overlapping with the inserted solutes
#   - Tag atoms in the beta column so colvars can pick them up:
#        beta = 1.0  -> solute 1  (the one near z1)
#        beta = 2.0  -> solute 2  (the one near z1 + 37)
#        beta = 3.0  -> POPC P atoms (bilayer reference)
#        beta = 0.0  -> everything else
#   - Write   system.psf  system.pdb  ref.pdb  (ref.pdb keeps the original
#     coords, used as cons-ref / fix-ref during equilibration)
#
# Usage (called by setup_all.sh):
#   vmd -dispdev text -e build_window.tcl -args <ANALOG> <Z1> <OUTDIR> <MEMBDIR> <ANALOGDIR>
#
#       ANALOG     : SCC | SCD | SCR
#       Z1         : center of solute 1 along z (in Angstrom, relative to box center)
#       OUTDIR     : directory in which to write system.psf / system.pdb / ref.pdb
#       MEMBDIR    : path to membrane/charmm-gui/namd  (contains step5_input.*)
#       ANALOGDIR  : path to analog/                   (contains scc.pdb, scd.pdb, scr.pdb
#                                                       and toppar-namd/top_all36_sidechains.str)

package require psfgen

if { [llength $argv] != 5 } {
    puts "Usage: vmd -dispdev text -e build_window.tcl -args ANALOG Z1 OUTDIR MEMBDIR ANALOGDIR"
    exit 1
}

set ANALOG    [string toupper [lindex $argv 0]]
set Z1        [lindex $argv 1]
set OUTDIR    [lindex $argv 2]
set MEMBDIR   [lindex $argv 3]
set ANALOGDIR [lindex $argv 4]

# solute 2 is placed 37 A BELOW solute 1 along z (MacCallum-style):
#   window  0 : Z1 =  0   Z2 = -37
#   window 37 : Z1 = 37   Z2 =   0
set DZ        -37.0
set Z2        [expr {$Z1 + $DZ}]

file mkdir $OUTDIR

# ---------------------------------------------------------------------------
# 1. load the analog topology
# ---------------------------------------------------------------------------
# top_all36_sidechains.str references atom types (HS, NH1, CT2, ...) defined
# in the master CHARMM36 protein topology, so that file must be loaded first.
# top_all36_cgenff.rtf brings in the few cgenff types referenced by SCP-style
# residues. Use toppar/ (full CHARMM stream with RESI records), NOT
# toppar-namd/ (which is the parameter-only version used by NAMD at run time).
topology $MEMBDIR/../toppar/top_all36_prot.rtf
topology $MEMBDIR/../toppar/top_all36_cgenff.rtf
topology $ANALOGDIR/toppar/top_all36_sidechains.str

# ---------------------------------------------------------------------------
# 2. read the equilibrated membrane PSF + PDB
# ---------------------------------------------------------------------------
readpsf  $MEMBDIR/step5_input.psf pdb $MEMBDIR/step5_input.pdb

# ---------------------------------------------------------------------------
# 3. place the two analog copies, dump their coords to temporary PDBs
#    (translated to xy=(0,0) and to the requested z)
# ---------------------------------------------------------------------------
mol new $ANALOGDIR/[string tolower $ANALOG].pdb type pdb waitfor all

set sel [atomselect top "all"]

# center on origin, then translate to (0, 0, z)
set com [measure center $sel weight none]
$sel moveby [vecscale -1.0 $com]

# solute 1
$sel moveby [list 0.0 0.0 $Z1]
$sel set segname "AN1"
$sel set resid   1
$sel writepdb    $OUTDIR/_an1.pdb
$sel moveby [list 0.0 0.0 [expr {-$Z1}]]

# solute 2
$sel moveby [list 0.0 0.0 $Z2]
$sel set segname "AN2"
$sel set resid   1
$sel writepdb    $OUTDIR/_an2.pdb

mol delete top

# ---------------------------------------------------------------------------
# 4. add the two solute segments with psfgen and load their coords
# ---------------------------------------------------------------------------
segment AN1 {
    pdb $OUTDIR/_an1.pdb
}
coordpdb $OUTDIR/_an1.pdb AN1

segment AN2 {
    pdb $OUTDIR/_an2.pdb
}
coordpdb $OUTDIR/_an2.pdb AN2

# ---------------------------------------------------------------------------
# 5. write a TEMP system, then use VMD to detect waters/ions clashing with
#    the analogs, and delete them with psfgen
# ---------------------------------------------------------------------------
writepsf $OUTDIR/_tmp.psf
writepdb $OUTDIR/_tmp.pdb

mol new     $OUTDIR/_tmp.psf type psf waitfor all
mol addfile $OUTDIR/_tmp.pdb type pdb waitfor all

# waters / ions within 3.0 A of any analog heavy atom
set clash [atomselect top \
    "(resname TIP3 SOD CLA POT) and \
     (same residue as within 3.0 of (segname AN1 AN2))"]

set clashList {}
foreach seg [$clash get segname] res [$clash get resid] {
    lappend clashList [list $seg $res]
}
set clashList [lsort -unique $clashList]
puts "Removing [llength $clashList] overlapping water/ion residues."

mol delete top

foreach pair $clashList {
    delatom [lindex $pair 0] [lindex $pair 1]
}

# ---------------------------------------------------------------------------
# 6. write the final psf / pdb
# ---------------------------------------------------------------------------
writepsf $OUTDIR/system.psf
writepdb $OUTDIR/system.pdb

# ---------------------------------------------------------------------------
# 7. mark beta column for colvars (solute1=1, solute2=2, POPC P=3) and
#    write system.pdb (overwrite) + ref.pdb (cons-ref during equilibration)
# ---------------------------------------------------------------------------
mol new     $OUTDIR/system.psf type psf waitfor all
mol addfile $OUTDIR/system.pdb type pdb waitfor all

set all [atomselect top "all"]
$all set beta 0.0
$all set occupancy 0.0

set s1 [atomselect top "segname AN1"]
$s1 set beta 1.0

set s2 [atomselect top "segname AN2"]
$s2 set beta 2.0

# bilayer reference: phosphorus atoms of POPC
set memref [atomselect top "resname POPC and name P"]
$memref set beta 3.0

$all writepdb $OUTDIR/system.pdb

# reference file used as cons-ref / fix-ref during equilibration
# occupancy = 1.0 on heavy atoms of solutes and POPC P -> harmonically restrained
$all set occupancy 0.0
set rstr [atomselect top "(segname AN1 AN2 and not hydrogen) or (resname POPC and name P)"]
$rstr set occupancy 1.0
$all writepdb $OUTDIR/ref.pdb

# clean up temp files
file delete $OUTDIR/_an1.pdb
file delete $OUTDIR/_an2.pdb
file delete $OUTDIR/_tmp.psf
file delete $OUTDIR/_tmp.pdb

puts "DONE building window: $OUTDIR  (analog=$ANALOG  z1=$Z1  z2=$Z2)"
exit
