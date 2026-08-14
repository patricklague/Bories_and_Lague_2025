#!/usr/bin/env python3
"""
extract_metadata.py

Read NAMD umbrella-sampling trajectories for one analog and produce WHAM
metadata using the ALC-2007 / MacCallum et al. setup:

    38 windows spaced by 1 A
    solute 1 centers:  0, 1, ..., 37 A
    solute 2 centers: -37, -36, ..., 0 A
    umbrella K: 3000 kJ/mol/nm^2 = 7.17 kcal/mol/A^2

The generated files can be passed directly to get_pmf.py:

    <out-dir>/traj/win<NN>_s1.dat   time(ps)  z(A)        for solute 1
    <out-dir>/traj/win<NN>_s2.dat   time(ps)  z(A)        for solute 2
    <out-dir>/metadata_s1.dat       file  center  K       (solute 1 / +z leaflet)
    <out-dir>/metadata_s2.dat       file  center  K       (solute 2 / -z leaflet)

Data layout per window  (runs/<analog>/win<NN>/) :
    system.psf            topology
    system.pdb            occupancy/beta flags
                            beta=1.0 -> solute 1   atoms
                            beta=2.0 -> solute 2   atoms
                            beta=3.0 -> bilayer reference (POPC P atoms)
    out/section*.dcd      production trajectory (multiple consecutive 1 ns chunks)

z is computed per frame as
    z_solute - z_reference
with a minimum-image correction along Lz (taken from the DCD's box).

Force-constant convention
-------------------------
Both NAMD/colvars and recent Grossfield WHAM (default) use
    U = (1/2) K (x - x0)^2
so we write the colvars force constant (kcal/mol/A^2) directly into the
metadata file. If you rebuild WHAM with --no-convert (i.e. so that it
assumes  U = K (x - x0)^2 ), pass --kwham-half to halve the constant
before writing the metadata.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np

try:
    import MDAnalysis as mda
except ImportError:
    sys.exit("ERROR: MDAnalysis is required (pip install MDAnalysis)")


# ---------------------------------------------------------------------------
# defaults that match templates/colvars.in and templates/production.namd
# ---------------------------------------------------------------------------
DEFAULT_KCV         = 7.17     # 3000 kJ/mol/nm^2 = 7.17 kcal/mol/A^2
DEFAULT_DZ_STEP     = 1.0      # Z1 spacing between windows (A)
DEFAULT_DZ_PAIR     = -37.0    # Z2 = Z1 + DZ_PAIR
DEFAULT_Z_MIN       = 0.0      # window 00 -> Z1 = 0
DEFAULT_FRAME_DT_PS = 10.0     # dcdfreq=5000, timestep=2 fs -> 10 ps/frame
DEFAULT_MAX_STEP_A  = 2.0      # restrained coordinate continuity check


def numeric_key(path: str) -> int:
    """Sort section1.dcd, section2.dcd, ..., section10.dcd, ... numerically."""
    m = re.search(r'section(\d+)\.dcd$', path)
    return int(m.group(1)) if m else 0


def list_sections(win_out: Path,
                  min_section: int = 1,
                  max_section: int | None = None) -> list[str]:
    """Return non-empty section*.dcd files with `min_section <= N <= max_section`."""
    # The trajectory is split into section*.dcd chunks. Filtering here lets us
    # use, for example, only section1..section30 without moving files around.
    files = glob.glob(str(win_out / "section*.dcd"))
    files = [f for f in files if os.path.getsize(f) > 0]
    files = [f for f in files
             if numeric_key(f) >= min_section
             and (max_section is None or numeric_key(f) <= max_section)]
    return sorted(files, key=numeric_key)


def parse_pdb_beta(pdb: Path) -> np.ndarray:
    """Return per-atom beta (tempfactor) values, in PDB file order.

    PSF + PDB Universes in MDAnalysis don't always expose tempfactors, so we
    parse the PDB directly. ATOM/HETATM columns 61-66 (1-indexed) hold beta.
    The order matches the PSF as long as system.pdb was generated alongside
    system.psf (the usual NAMD/psfgen convention).
    """
    betas = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    betas.append(float(line[60:66]))
                except ValueError:
                    betas.append(0.0)
    return np.asarray(betas)


def select_indices(pdb: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Atom indices for solute 1 (beta=1), solute 2 (beta=2), reference (beta=3)."""
    b = parse_pdb_beta(pdb)
    s1 = np.where(np.isclose(b, 1.0))[0]
    s2 = np.where(np.isclose(b, 2.0))[0]
    ref = np.where(np.isclose(b, 3.0))[0]
    if not (s1.size and s2.size and ref.size):
        raise RuntimeError(
            f"Could not find beta-flagged atoms in {pdb}: "
            f"|s1|={s1.size} |s2|={s2.size} |ref|={ref.size}"
        )
    return s1, s2, ref


def periodic_com_z(group, box_z: float) -> float:
    """Mass-weighted z center after unwrapping atoms around one group atom."""
    z = group.positions[:, 2].astype(float)
    masses = group.masses.astype(float)
    if not np.isfinite(box_z) or box_z <= 0:
        return float(np.average(z, weights=masses))
    anchor = z[0]
    dz = z - anchor
    dz -= box_z * np.round(dz / box_z)
    return float(anchor + np.average(dz, weights=masses))


def extract_window(win_dir: Path, frame_dt_ps: float, skip_frames: int,
                   min_section: int = 1, max_section: int | None = None
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Return (t_ps, z1_A, z2_A) for one window, or None if nothing usable."""
    psf = win_dir / "system.psf"
    pdb = win_dir / "system.pdb"
    sections = list_sections(win_dir / "out", min_section, max_section)
    if not (psf.exists() and pdb.exists() and sections):
        return None

    # PDB beta flags define the two restrained solutes and the membrane
    # reference atoms. They are converted to integer atom indices once here.
    s1_idx, s2_idx, ref_idx = select_indices(pdb)
    u = mda.Universe(str(psf), sections)
    sol1 = u.atoms[s1_idx]
    sol2 = u.atoms[s2_idx]
    ref  = u.atoms[ref_idx]

    z1, z2 = [], []
    for ts in u.trajectory:
        Lz = ts.dimensions[2] if ts.dimensions is not None else 0.0
        # DCD coordinates were written with wrapAll on. A direct COM of the
        # wrapped POPC-P atoms can jump by several Angstrom when the bilayer
        # straddles a periodic boundary. Unwrap each group locally first.
        zr  = periodic_com_z(ref, Lz)
        z1f = periodic_com_z(sol1, Lz) - zr
        z2f = periodic_com_z(sol2, Lz) - zr
        # Keep the solute-reference distance in the nearest periodic image.
        # This avoids jumps when a solute crosses the periodic z boundary.
        if Lz > 0:                                      # minimum-image
            z1f -= Lz * round(z1f / Lz)
            z2f -= Lz * round(z2f / Lz)
        z1.append(z1f); z2.append(z2f)

    z1 = np.asarray(z1); z2 = np.asarray(z2)
    if z1.size <= skip_frames:
        return None
    z1 = z1[skip_frames:]; z2 = z2[skip_frames:]
    t  = np.arange(1, z1.size + 1, dtype=float) * frame_dt_ps
    return t, z1, z2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analog",   required=True, help="analog name, e.g. scc")
    ap.add_argument("--runs-dir", required=True, type=Path)
    ap.add_argument("--out-dir",  required=True, type=Path)
    ap.add_argument("--equil-ps", type=float, default=5000.0,
                    help="ps of production to discard per window (default 5000)")
    ap.add_argument("--frame-dt-ps", type=float, default=DEFAULT_FRAME_DT_PS,
                    help="time between consecutive dcd frames (default 10 ps)")
    ap.add_argument("--max-step-A", type=float, default=DEFAULT_MAX_STEP_A,
                    help="reject a window if either extracted z coordinate "
                         "jumps farther than this between frames (default 2 A)")
    ap.add_argument("--kcv", type=float, default=DEFAULT_KCV,
                    help="colvars force constant in kcal/mol/A^2 (default 7.17)")
    ap.add_argument("--kwham-half", action="store_true",
                    help="write K/2 instead of K (only for WHAM compiled with "
                         "--no-convert, i.e. U = K(x-x0)^2 convention)")
    ap.add_argument("--z-min",   type=float, default=DEFAULT_Z_MIN)
    ap.add_argument("--dz-step", type=float, default=DEFAULT_DZ_STEP)
    ap.add_argument("--dz-pair", type=float, default=DEFAULT_DZ_PAIR)
    ap.add_argument("--min-section", type=int, default=1,
                    help="lowest section index to include (default 1)")
    ap.add_argument("--max-section", type=int, default=None,
                    help="highest section index to include (default: no limit)")
    args = ap.parse_args()

    analog_dir = args.runs_dir / args.analog
    if not analog_dir.is_dir():
        sys.exit(f"ERROR: {analog_dir} not found")

    out_dir   = args.out_dir.resolve()
    traj_dir  = out_dir / "traj"
    traj_dir.mkdir(parents=True, exist_ok=True)
    print(f"  writing to {out_dir}")

    # Most WHAM builds and NAMD colvars use U = 1/2 K(x-x0)^2, so by default
    # the metadata K is the same value used in the colvars restraint.
    k_wham      = args.kcv / 2.0 if args.kwham_half else args.kcv
    skip_frames = int(round(args.equil_ps / args.frame_dt_ps))

    # discover win?? directories
    win_dirs = sorted(p for p in analog_dir.glob("win[0-9][0-9]") if p.is_dir())
    if not win_dirs:
        sys.exit(f"ERROR: no win?? directories under {analog_dir}")

    meta1_path = out_dir / "metadata_s1.dat"
    meta2_path = out_dir / "metadata_s2.dat"
    # truncate so a partial run leaves a consistent (possibly empty) file
    meta1_path.write_text("")
    meta2_path.write_text("")

    meta1_lines, meta2_lines = [], []
    kept = 0
    frame_counts = []
    for win in win_dirs:
        nn = int(win.name[3:])
        # Window centers reproduce the two-solute ALC layout: solute 1 scans
        # 0..+37 A while solute 2 scans -37..0 A in the opposite leaflet.
        z1c = args.z_min + nn * args.dz_step
        z2c = z1c + args.dz_pair

        try:
            res = extract_window(win, args.frame_dt_ps, skip_frames,
                                 args.min_section, args.max_section)
        except Exception as exc:                       # noqa: BLE001
            print(f"  ! {args.analog}/{win.name}: {type(exc).__name__}: {exc}")
            continue
        if res is None:
            print(f"  ! {args.analog}/{win.name}: no usable data, skipped")
            continue
        t, z1, z2 = res

        max_step = max(np.max(np.abs(np.diff(z1))),
                       np.max(np.abs(np.diff(z2)))) if z1.size > 1 else 0.0
        if max_step > args.max_step_A:
            sys.exit(
                f"ERROR: {args.analog}/{win.name} has a {max_step:.3f} A "
                "frame-to-frame jump after periodic centering; inspect the "
                "trajectory or increase --max-step-A only if justified"
            )

        f1 = traj_dir / f"{win.name}_s1.dat"
        f2 = traj_dir / f"{win.name}_s2.dat"
        np.savetxt(f1, np.column_stack([t, z1]), fmt="%10.2f %10.4f")
        np.savetxt(f2, np.column_stack([t, z2]), fmt="%10.2f %10.4f")

        rel1 = os.path.relpath(f1, out_dir)
        rel2 = os.path.relpath(f2, out_dir)
        # Metadata is deliberately simple: time-series file, restraint center,
        # force constant. Bin width is not stored here; get_pmf.py chooses it.
        meta1_lines.append(f"{rel1:<32s} {z1c:8.3f} {k_wham:8.3f}")
        meta2_lines.append(f"{rel2:<32s} {z2c:8.3f} {k_wham:8.3f}")
        # flush after every window so partial results survive a later crash
        meta1_path.write_text("\n".join(meta1_lines) + "\n")
        meta2_path.write_text("\n".join(meta2_lines) + "\n")
        kept += 1
        frame_counts.append(z1.size)
        print(f"  + {win.name}  z1={z1c:6.2f}  z2={z2c:7.2f}  "
              f"frames={z1.size}")

    if kept != 38:
        sys.exit(f"ERROR: expected 38 complete windows, extracted {kept}")
    if len(set(frame_counts)) != 1:
        sys.exit(f"ERROR: inconsistent frame counts across windows: "
                 f"{sorted(set(frame_counts))}")

    print(f"==> {args.analog}: {kept}/{len(win_dirs)} windows -> "
          f"{meta1_path}, {meta2_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
