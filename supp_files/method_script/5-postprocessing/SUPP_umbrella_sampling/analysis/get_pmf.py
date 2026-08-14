#!/usr/bin/env python3
"""
get_pmf.py

Take the WHAM metadata produced by extract_metadata.py for three
independent 50 ns blocks of trajectory data (section101-section150,
section151-section200, section201-section250) and write, per block,
    <in-dir>/<block>/pmf_s1.dat        per-solute WHAM PMF (z, PMF, ...)
    <in-dir>/<block>/pmf_s2.dat
    <in-dir>/<block>/pmf-us-<analog>-<block>.dat
                                       block PMF (leaflets folded together)
and, combined over the three blocks,
    <out-dir>/pmf-us-<analog>.dat
        |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks

Within each block, solute 1 (+z leaflet, centers 0..+37 A) and solute 2
(-z leaflet, centers -37..0 A) are run through WHAM separately, then the
solute-2 profile is mirrored to |z| and averaged with the solute-1 profile
(ALC-2007 / MacCallum et al. two-leaflet treatment). This gives one PMF
estimate per 50 ns block.

The final PMF average and standard error are computed *across the three
blocks* (std/sqrt(n), n up to 3), which treats each 50 ns block as an
independent replicate. See get_pmf2.py for the alternative approach (single
50 ns block, error estimated from the two leaflets) used for shorter
(100 ns) trajectories.

By default WHAM uses 200 bins across HIST_MIN..HIST_MAX. With the default
range -38..38 A, this gives a bin width of 0.38 A. You can override either
the number of bins with --num-bins or the bin width with --bin-A.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# WHAM grid. Metadata contains window centers and force constants only; the
# PMF histogram grid is chosen here when WHAM is run.
# ---------------------------------------------------------------------------
DEFAULT_HIST_MIN = -38.0
DEFAULT_HIST_MAX =  38.0
DEFAULT_NUM_BINS = 200
DEFAULT_TOL      = 1.0e-5
DEFAULT_TEMP     = 303.15
DEFAULT_NUM_PAD  = 0
DEFAULT_SEED     = 12345
DEFAULT_Z_BULK   = 30.0    # |z| >= Z_BULK averaged to zero the PMF
DEFAULT_MAX_ITER = 100000  # safety cap for the python WHAM iteration
KB_KCAL          = 0.0019872041   # Boltzmann constant in kcal/(mol K)


def _read_metadata(meta: Path) -> tuple[list[Path], np.ndarray, np.ndarray]:
    """Read Grossfield-style metadata: trajectory-file, window-center, K."""
    centers, kvals, files = [], [], []
    with open(meta) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            f, c, k = parts[0], float(parts[1]), float(parts[2])
            p = Path(f)
            if not p.is_absolute():
                p = (meta.parent / p).resolve()
            files.append(p)
            centers.append(c)
            kvals.append(k)
    return files, np.asarray(centers), np.asarray(kvals)


def wham_py(meta: Path, pmf: Path, log: Path,
            hist_min: float, hist_max: float, nbins: int,
            temp: float, tol: float, max_iter: int = DEFAULT_MAX_ITER) -> None:
    """Pure-numpy 1D WHAM. Assumes U_i(x) = (1/2) K_i (x - x0_i)^2
    (the colvars / NAMD convention; this is also the default in Grossfield
    WHAM >= 2.0.10).

    Writes a Grossfield-style PMF file:
        # coor  Free  +/-  Prob  +/-
        x_k     g_k   0    p_k   0
    and a short log file.
    """
    files, centers, kvals = _read_metadata(meta)
    nwin = len(files)
    kT = KB_KCAL * temp

    # The WHAM grid is independent of the umbrella window spacing. More bins
    # give a smoother/finer PMF output grid, but do not create new sampling.
    edges = np.linspace(hist_min, hist_max, nbins + 1)
    x = 0.5 * (edges[:-1] + edges[1:])                     # bin centres
    dx = (hist_max - hist_min) / nbins

    H = np.zeros((nwin, nbins), dtype=float)               # counts per window
    n = np.zeros(nwin, dtype=float)
    for i, f in enumerate(files):
        if not f.is_file() or f.stat().st_size == 0:
            continue
        data = np.loadtxt(f, comments=("#", "@"))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        vals = data[:, 1] if data.shape[1] >= 2 else data[:, 0]
        h, _ = np.histogram(vals, bins=edges)
        H[i] = h
        n[i] = h.sum()

    if nwin != 38:
        raise RuntimeError(f"expected 38 umbrella windows, found {nwin} in {meta}")
    if np.any(n == 0):
        empty = ", ".join(str(i) for i in np.where(n == 0)[0])
        raise RuntimeError(f"empty/out-of-range umbrella windows in {meta}: {empty}")
    if n.sum() == 0:
        raise RuntimeError(f"no samples landed in [{hist_min},{hist_max}]")

    normalized = H / n[:, None]
    adjacent_overlap = np.minimum(normalized[:-1], normalized[1:]).sum(axis=1)
    if np.any(adjacent_overlap == 0):
        pairs = ", ".join(
            f"{i}-{i + 1}" for i in np.where(adjacent_overlap == 0)[0]
        )
        raise RuntimeError(f"zero histogram overlap between windows {pairs} in {meta}")

    # bias matrix U_ik = (1/2) K_i (x_k - x0_i)^2 ;  W = exp(-U/kT)
    U = 0.5 * kvals[:, None] * (x[None, :] - centers[:, None]) ** 2
    W = np.exp(-U / kT)

    sum_H = H.sum(axis=0)                                  # total count per bin
    F = np.zeros(nwin)                                     # window free energies / kT
    rho = np.full(nbins, np.nan)

    occupied = sum_H > 0
    if not occupied.any():
        raise RuntimeError("all histogram bins are empty")

    with open(log, "w") as flog:
        flog.write(f"# python WHAM  T={temp}  bins={nbins}  "
                   f"range=[{hist_min},{hist_max}]  tol={tol}\n")
        flog.write(f"# {nwin} windows, total counts = {n.sum():.0f}\n")
        flog.write(f"# adjacent overlap: min={adjacent_overlap.min():.4f} "
               f"median={np.median(adjacent_overlap):.4f}\n")
        converged = False
        for it in range(1, max_iter + 1):
            # rho_k = sum_i H_ik  /  sum_i n_i exp(F_i) W_ik
            denom = (n * np.exp(F))[:, None] * W            # (nwin, nbins)
            denom_sum = denom.sum(axis=0)
            rho = np.where(occupied & (denom_sum > 0),
                           sum_H / np.where(denom_sum > 0, denom_sum, 1.0),
                           0.0)
            # F_i = -ln( sum_k W_ik rho_k dx )
            Z = (W * rho[None, :]).sum(axis=1) * dx
            with np.errstate(divide="ignore"):
                F_new = -np.log(np.where(Z > 0, Z, np.nan))
            F_new -= F_new[0]                              # gauge fix
            diff = np.nanmax(np.abs(F_new - F))
            F = F_new
            if it % 200 == 0 or diff < tol:
                flog.write(f"  iter {it:6d}  max|dF| = {diff:.3e}\n")
            if np.isfinite(diff) and diff < tol:
                flog.write(f"# converged in {it} iterations\n")
                converged = True
                break
        else:
            flog.write(f"# WARNING: not converged after {max_iter} iterations\n")

    if not converged:
        raise RuntimeError(f"WHAM did not converge for {meta} after {max_iter} iterations")

    g = np.where(rho > 0, -kT * np.log(np.where(rho > 0, rho, np.nan)), np.nan)
    if np.isfinite(g).any():
        g = g - np.nanmin(g)

    with open(pmf, "w") as fh:
        fh.write("#Coor\t\tFree\t+/-\t\tProb\t\t+/-\n")
        for xk, gk, pk in zip(x, g, rho):
            gs = f"{gk:.6f}" if np.isfinite(gk) else "nan"
            ps = f"{pk:.6e}" if np.isfinite(pk) and pk > 0 else "nan"
            fh.write(f"{xk:.6f}\t{gs}\t0\t{ps}\t0\n")


def run_wham(wham_bin: str, hist_min: float, hist_max: float, nbins: int,
             tol: float, temp: float, num_pad: int, meta: Path,
             pmf: Path, log: Path, nboot: int, seed: int) -> None:
    """Run an external Grossfield WHAM binary with the same grid settings."""
    if shutil.which(wham_bin) is None:
        sys.exit(f"ERROR: '{wham_bin}' not on PATH. "
                 "Get it from http://membrane.urmc.rochester.edu/?page_id=126")
    cmd = [wham_bin, f"{hist_min}", f"{hist_max}", f"{nbins}", f"{tol}",
           f"{temp}", f"{num_pad}", str(meta), str(pmf)]
    if nboot > 0:
        cmd += [f"{nboot}", f"{seed}"]
    with open(log, "w") as fh:
        subprocess.run(cmd, cwd=meta.parent, stdout=fh, stderr=subprocess.STDOUT,
                       check=True)


def load_pmf(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the coordinate and free-energy columns from a WHAM PMF file."""
    d = np.loadtxt(path, comments=("#", "@"))
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0], d[:, 1]


def shift_to_bulk(z: np.ndarray, g: np.ndarray, z_bulk: float) -> np.ndarray:
    """Set the zero of a PMF to the average value in the bulk-water region."""
    finite = np.isfinite(g)
    bulk = finite & (np.abs(z) >= z_bulk)
    if not bulk.any():
        bulk = finite
    return g - np.nanmean(g[bulk])


def fold_side(zside: np.ndarray, gside: np.ndarray,
              target: np.ndarray) -> np.ndarray:
    """Interpolate one leaflet estimate onto the positive |z| output grid."""
    order = np.argsort(zside)
    zside, gside = zside[order], gside[order]
    finite = np.isfinite(gside)
    if finite.sum() < 2:
        return np.full_like(target, np.nan, dtype=float)
    return np.interp(target, zside[finite], gside[finite],
                     left=np.nan, right=np.nan)


def combine_leaflets(p1: Path, p2: Path, out: Path, analog: str, block: str,
                     z_bulk: float, bin_width: float
                     ) -> tuple[np.ndarray, np.ndarray]:
    """Combine the two offset solutes with the ALC-2007 leaflet method.

    In our setup solute 1 is restrained at 0..+37 A and solute 2 at -37..0 A.
    The two WHAM profiles are therefore independent estimates of the same
    |z|-PMF from opposite leaflets for one 50 ns block. We mirror solute 2
    onto the positive grid and average the two estimates into a single
    per-block PMF (returned for later averaging across blocks). The
    leaflet-to-leaflet spread is written to `out` for diagnostic purposes
    only; the reported uncertainty of the final PMF comes from the
    block-to-block spread (see combine_blocks).
    """
    z1, g1 = load_pmf(p1)
    z2, g2 = load_pmf(p2)
    if not np.allclose(z1, z2):
        sys.exit(f"{p1} and {p2} are on different z grids")
    z = z1
    g1 = shift_to_bulk(z, g1, z_bulk)
    g2 = shift_to_bulk(z, g2, z_bulk)

    # Keep the positive side from solute 1 and mirror the negative side from
    # solute 2, so both leaflets describe the same |z| coordinate.
    pos = z >= 0
    neg = z < 0
    zp = z[pos]
    zn = -z[neg]

    g_plus_leaflet = g1[pos]
    g_minus_leaflet = fold_side(zn, g2[neg], zp)

    stack = np.vstack([g_plus_leaflet, g_minus_leaflet])     # 2 x nbin
    n     = np.sum(np.isfinite(stack), axis=0)
    with np.errstate(invalid="ignore"):
        mean = np.divide(np.nansum(stack, axis=0), n,
                         out=np.full_like(zp, np.nan, dtype=float),
                         where=n > 0)
    with np.errstate(invalid="ignore"):
        std = np.where(n >= 2,
                       np.nanstd(stack, axis=0, ddof=1),
                       np.nan)
    se = std / np.sqrt(np.maximum(n, 1))

    finite = np.isfinite(mean)
    bulk   = finite & (zp >= z_bulk)
    ref    = (np.nanmean(mean[bulk]) if bulk.any()
              else (np.nanmean(mean[finite]) if finite.any() else 0.0))
    mean   = mean - ref

    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}', block '{block}' "
             f"(ALC-2007 two-leaflet averaging, {bin_width:g} A bins)\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_leaflets\n")
        fh.write(f"# zeroed on |z| >= {z_bulk:g} A; SE = std/sqrt(n)\n")
        fh.write("# leaflet estimates: solute 1 (+z) and mirrored solute 2 (-z)\n")
        fh.write("# diagnostic only -- final uncertainty comes from block-to-block spread\n")
        for zi, gi, si, ni in zip(zp, mean, se, n):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {int(ni):4d}\n")
    return zp, mean


def combine_blocks(z: np.ndarray, block_means: list[np.ndarray],
                   block_names: list[str], out: Path, analog: str
                   ) -> int:
    """Average the per-block (leaflet-folded) PMFs and estimate the
    standard error from the spread between the 50 ns blocks.
    """
    stack = np.vstack(block_means)                          # nblocks x nbin
    n = np.sum(np.isfinite(stack), axis=0)
    with np.errstate(invalid="ignore"):
        mean = np.divide(np.nansum(stack, axis=0), n,
                         out=np.full(z.shape, np.nan, dtype=float),
                         where=n > 0)
    with np.errstate(invalid="ignore"):
        std = np.where(n >= 2, np.nanstd(stack, axis=0, ddof=1), np.nan)
    se = std / np.sqrt(np.maximum(n, 1))

    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}'  "
                 f"(averaged over {len(block_names)} 50 ns blocks: "
                 f"{', '.join(block_names)})\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks\n")
        fh.write("# each block PMF folds solute 1 (+z) and mirrored solute 2 (-z) leaflets;\n")
        fh.write("# SE here = std/sqrt(n) across blocks\n")
        for zi, gi, si, ni in zip(z, mean, se, n):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {int(ni):4d}\n")
    return int(np.isfinite(mean).sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analog",   required=True)
    ap.add_argument("--in-dir",   required=True, type=Path,
                    help="parent directory containing one subdirectory per "
                         "50 ns block (see --blocks), each with "
                         "metadata_s1.dat/metadata_s2.dat from extract_metadata.py")
    ap.add_argument("--out-dir",  type=Path, default=None,
                    help="where to write the combined pmf-us-<analog>.dat "
                         "(default: same as --in-dir)")
    ap.add_argument("--blocks", nargs="+",
                    default=["block1", "block2", "block3"],
                    help="block subdirectory names under --in-dir, in order "
                         "(default: block1 block2 block3, matching "
                         "section101-150 / section151-200 / section201-250)")
    ap.add_argument("--backend", choices=["python", "grossfield"],
                    default="python",
                    help="WHAM implementation (default: python, no external dep)")
    ap.add_argument("--wham-bin", default="wham",
                    help="path/name of Grossfield wham binary (--backend grossfield)")
    ap.add_argument("--hist-min", type=float, default=DEFAULT_HIST_MIN)
    ap.add_argument("--hist-max", type=float, default=DEFAULT_HIST_MAX)
    ap.add_argument("--num-bins", type=int, default=DEFAULT_NUM_BINS,
                    help="number of WHAM histogram bins (default 200)")
    ap.add_argument("--bin-A",    type=float, default=None,
                    help="bin width in A; overrides --num-bins when set")
    ap.add_argument("--tol",      type=float, default=DEFAULT_TOL)
    ap.add_argument("--temp",     type=float, default=DEFAULT_TEMP)
    ap.add_argument("--num-pad",  type=int,   default=DEFAULT_NUM_PAD)
    ap.add_argument("--z-bulk",   type=float, default=DEFAULT_Z_BULK)
    ap.add_argument("--nboot",    type=int,   default=0)
    ap.add_argument("--seed",     type=int,   default=DEFAULT_SEED)
    args = ap.parse_args()

    out_dir = args.out_dir or args.in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    hist_range = args.hist_max - args.hist_min
    if hist_range <= 0:
        sys.exit("ERROR: --hist-max must be greater than --hist-min")
    if args.bin_A is not None:
        if args.bin_A <= 0:
            sys.exit("ERROR: --bin-A must be positive")
        nbins = int(round(hist_range / args.bin_A))
        bin_width = hist_range / nbins
    else:
        if args.num_bins <= 0:
            sys.exit("ERROR: --num-bins must be positive")
        nbins = args.num_bins
        bin_width = hist_range / nbins

    block_z: np.ndarray | None = None
    block_means: list[np.ndarray] = []
    used_blocks: list[str] = []
    for block in args.blocks:
        block_dir = args.in_dir / block
        meta1 = block_dir / "metadata_s1.dat"
        meta2 = block_dir / "metadata_s2.dat"
        missing = [m for m in (meta1, meta2) if not m.is_file() or m.stat().st_size == 0]
        if missing:
            print(f"  ! block '{block}': missing/empty metadata "
                  f"({', '.join(str(m) for m in missing)}), skipped")
            continue

        pmf1 = block_dir / "pmf_s1.dat"
        pmf2 = block_dir / "pmf_s2.dat"
        log1 = block_dir / "wham_s1.log"
        log2 = block_dir / "wham_s2.log"

        print(f"  [{block}] WHAM s1: bins={nbins} ({bin_width:g} A) "
              f"range=[{args.hist_min},{args.hist_max}]  T={args.temp}  "
              f"backend={args.backend}")
        if args.backend == "python":
            wham_py(meta1, pmf1, log1, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            wham_py(meta2, pmf2, log2, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
        else:
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta1, pmf1, log1, args.nboot, args.seed)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta2, pmf2, log2, args.nboot, args.seed)

        block_pmf = block_dir / f"pmf-us-{args.analog}-{block}.dat"
        zp, mean = combine_leaflets(pmf1, pmf2, block_pmf, args.analog, block,
                                    args.z_bulk, bin_width)
        if block_z is None:
            block_z = zp
        elif not np.allclose(block_z, zp):
            sys.exit(f"ERROR: block '{block}' PMF grid differs from earlier blocks")
        block_means.append(mean)
        used_blocks.append(block)
        print(f"  [{block}] -> {block_pmf}")

    if not block_means:
        sys.exit("ERROR: no usable blocks found under " + str(args.in_dir))

    combined = out_dir / f"pmf-us-{args.analog}.dat"
    nbin_ok = combine_blocks(block_z, block_means, used_blocks, combined, args.analog)
    print(f"==> {combined}   ({nbin_ok} usable bins, {len(used_blocks)} blocks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

=======
#!/usr/bin/env python3
"""
get_pmf.py

Take the WHAM metadata produced by extract_metadata.py for three
independent 50 ns blocks of trajectory data (section101-section150,
section151-section200, section201-section250) and write, per block,
    <in-dir>/<block>/pmf_s1.dat        per-solute WHAM PMF (z, PMF, ...)
    <in-dir>/<block>/pmf_s2.dat
    <in-dir>/<block>/pmf-us-<analog>-<block>.dat
                                       block PMF (leaflets folded together)
and, combined over the three blocks,
    <out-dir>/pmf-us-<analog>.dat
        |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks

Within each block, solute 1 (+z leaflet, centers 0..+37 A) and solute 2
(-z leaflet, centers -37..0 A) are run through WHAM separately, then the
solute-2 profile is mirrored to |z| and averaged with the solute-1 profile
(ALC-2007 / MacCallum et al. two-leaflet treatment). This gives one PMF
estimate per 50 ns block.

The final PMF average and standard error are computed *across the three
blocks* (std/sqrt(n), n up to 3), which treats each 50 ns block as an
independent replicate. See get_pmf2.py for the alternative approach (single
50 ns block, error estimated from the two leaflets) used for shorter
(100 ns) trajectories.

By default WHAM uses 200 bins across HIST_MIN..HIST_MAX. With the default
range -38..38 A, this gives a bin width of 0.38 A. You can override either
the number of bins with --num-bins or the bin width with --bin-A.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# WHAM grid. Metadata contains window centers and force constants only; the
# PMF histogram grid is chosen here when WHAM is run.
# ---------------------------------------------------------------------------
DEFAULT_HIST_MIN = -38.0
DEFAULT_HIST_MAX =  38.0
DEFAULT_NUM_BINS = 200
DEFAULT_TOL      = 1.0e-5
DEFAULT_TEMP     = 303.15
DEFAULT_NUM_PAD  = 0
DEFAULT_SEED     = 12345
DEFAULT_Z_BULK   = 35.0    # |z| >= Z_BULK averaged to zero the PMF
DEFAULT_MAX_ITER = 100000  # safety cap for the python WHAM iteration
KB_KCAL          = 0.0019872041   # Boltzmann constant in kcal/(mol K)


def _read_metadata(meta: Path) -> tuple[list[Path], np.ndarray, np.ndarray]:
    """Read Grossfield-style metadata: trajectory-file, window-center, K."""
    centers, kvals, files = [], [], []
    with open(meta) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            f, c, k = parts[0], float(parts[1]), float(parts[2])
            p = Path(f)
            if not p.is_absolute():
                p = (meta.parent / p).resolve()
            files.append(p)
            centers.append(c)
            kvals.append(k)
    return files, np.asarray(centers), np.asarray(kvals)


def wham_py(meta: Path, pmf: Path, log: Path,
            hist_min: float, hist_max: float, nbins: int,
            temp: float, tol: float, max_iter: int = DEFAULT_MAX_ITER) -> None:
    """Pure-numpy 1D WHAM. Assumes U_i(x) = (1/2) K_i (x - x0_i)^2
    (the colvars / NAMD convention; this is also the default in Grossfield
    WHAM >= 2.0.10).

    Writes a Grossfield-style PMF file:
        # coor  Free  +/-  Prob  +/-
        x_k     g_k   0    p_k   0
    and a short log file.
    """
    files, centers, kvals = _read_metadata(meta)
    nwin = len(files)
    kT = KB_KCAL * temp

    # The WHAM grid is independent of the umbrella window spacing. More bins
    # give a smoother/finer PMF output grid, but do not create new sampling.
    edges = np.linspace(hist_min, hist_max, nbins + 1)
    x = 0.5 * (edges[:-1] + edges[1:])                     # bin centres
    dx = (hist_max - hist_min) / nbins

    H = np.zeros((nwin, nbins), dtype=float)               # counts per window
    n = np.zeros(nwin, dtype=float)
    for i, f in enumerate(files):
        if not f.is_file() or f.stat().st_size == 0:
            continue
        data = np.loadtxt(f, comments=("#", "@"))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        vals = data[:, 1] if data.shape[1] >= 2 else data[:, 0]
        h, _ = np.histogram(vals, bins=edges)
        H[i] = h
        n[i] = h.sum()

    if nwin != 38:
        raise RuntimeError(f"expected 38 umbrella windows, found {nwin} in {meta}")
    if np.any(n == 0):
        empty = ", ".join(str(i) for i in np.where(n == 0)[0])
        raise RuntimeError(f"empty/out-of-range umbrella windows in {meta}: {empty}")
    if n.sum() == 0:
        raise RuntimeError(f"no samples landed in [{hist_min},{hist_max}]")

    normalized = H / n[:, None]
    adjacent_overlap = np.minimum(normalized[:-1], normalized[1:]).sum(axis=1)
    if np.any(adjacent_overlap == 0):
        pairs = ", ".join(
            f"{i}-{i + 1}" for i in np.where(adjacent_overlap == 0)[0]
        )
        raise RuntimeError(f"zero histogram overlap between windows {pairs} in {meta}")

    # bias matrix U_ik = (1/2) K_i (x_k - x0_i)^2 ;  W = exp(-U/kT)
    U = 0.5 * kvals[:, None] * (x[None, :] - centers[:, None]) ** 2
    W = np.exp(-U / kT)

    sum_H = H.sum(axis=0)                                  # total count per bin
    F = np.zeros(nwin)                                     # window free energies / kT
    rho = np.full(nbins, np.nan)

    occupied = sum_H > 0
    if not occupied.any():
        raise RuntimeError("all histogram bins are empty")

    with open(log, "w") as flog:
        flog.write(f"# python WHAM  T={temp}  bins={nbins}  "
                   f"range=[{hist_min},{hist_max}]  tol={tol}\n")
        flog.write(f"# {nwin} windows, total counts = {n.sum():.0f}\n")
        flog.write(f"# adjacent overlap: min={adjacent_overlap.min():.4f} "
               f"median={np.median(adjacent_overlap):.4f}\n")
        converged = False
        for it in range(1, max_iter + 1):
            # rho_k = sum_i H_ik  /  sum_i n_i exp(F_i) W_ik
            denom = (n * np.exp(F))[:, None] * W            # (nwin, nbins)
            denom_sum = denom.sum(axis=0)
            rho = np.where(occupied & (denom_sum > 0),
                           sum_H / np.where(denom_sum > 0, denom_sum, 1.0),
                           0.0)
            # F_i = -ln( sum_k W_ik rho_k dx )
            Z = (W * rho[None, :]).sum(axis=1) * dx
            with np.errstate(divide="ignore"):
                F_new = -np.log(np.where(Z > 0, Z, np.nan))
            F_new -= F_new[0]                              # gauge fix
            diff = np.nanmax(np.abs(F_new - F))
            F = F_new
            if it % 200 == 0 or diff < tol:
                flog.write(f"  iter {it:6d}  max|dF| = {diff:.3e}\n")
            if np.isfinite(diff) and diff < tol:
                flog.write(f"# converged in {it} iterations\n")
                converged = True
                break
        else:
            flog.write(f"# WARNING: not converged after {max_iter} iterations\n")

    if not converged:
        raise RuntimeError(f"WHAM did not converge for {meta} after {max_iter} iterations")

    g = np.where(rho > 0, -kT * np.log(np.where(rho > 0, rho, np.nan)), np.nan)
    if np.isfinite(g).any():
        g = g - np.nanmin(g)

    with open(pmf, "w") as fh:
        fh.write("#Coor\t\tFree\t+/-\t\tProb\t\t+/-\n")
        for xk, gk, pk in zip(x, g, rho):
            gs = f"{gk:.6f}" if np.isfinite(gk) else "nan"
            ps = f"{pk:.6e}" if np.isfinite(pk) and pk > 0 else "nan"
            fh.write(f"{xk:.6f}\t{gs}\t0\t{ps}\t0\n")


def run_wham(wham_bin: str, hist_min: float, hist_max: float, nbins: int,
             tol: float, temp: float, num_pad: int, meta: Path,
             pmf: Path, log: Path, nboot: int, seed: int) -> None:
    """Run an external Grossfield WHAM binary with the same grid settings."""
    if shutil.which(wham_bin) is None:
        sys.exit(f"ERROR: '{wham_bin}' not on PATH. "
                 "Get it from http://membrane.urmc.rochester.edu/?page_id=126")
    cmd = [wham_bin, f"{hist_min}", f"{hist_max}", f"{nbins}", f"{tol}",
           f"{temp}", f"{num_pad}", str(meta), str(pmf)]
    if nboot > 0:
        cmd += [f"{nboot}", f"{seed}"]
    with open(log, "w") as fh:
        subprocess.run(cmd, cwd=meta.parent, stdout=fh, stderr=subprocess.STDOUT,
                       check=True)


def load_pmf(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the coordinate and free-energy columns from a WHAM PMF file."""
    d = np.loadtxt(path, comments=("#", "@"))
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0], d[:, 1]


def shift_to_bulk(z: np.ndarray, g: np.ndarray, z_bulk: float) -> np.ndarray:
    """Set the zero of a PMF to the average value in the bulk-water region."""
    finite = np.isfinite(g)
    bulk = finite & (np.abs(z) >= z_bulk)
    if not bulk.any():
        bulk = finite
    return g - np.nanmean(g[bulk])


def fold_side(zside: np.ndarray, gside: np.ndarray,
              target: np.ndarray) -> np.ndarray:
    """Interpolate one leaflet estimate onto the positive |z| output grid."""
    order = np.argsort(zside)
    zside, gside = zside[order], gside[order]
    finite = np.isfinite(gside)
    if finite.sum() < 2:
        return np.full_like(target, np.nan, dtype=float)
    return np.interp(target, zside[finite], gside[finite],
                     left=np.nan, right=np.nan)


def combine_leaflets(p1: Path, p2: Path, out: Path, analog: str, block: str,
                     z_bulk: float, bin_width: float
                     ) -> tuple[np.ndarray, np.ndarray]:
    """Combine the two offset solutes with the ALC-2007 leaflet method.

    In our setup solute 1 is restrained at 0..+37 A and solute 2 at -37..0 A.
    The two WHAM profiles are therefore independent estimates of the same
    |z|-PMF from opposite leaflets for one 50 ns block. We mirror solute 2
    onto the positive grid and average the two estimates into a single
    per-block PMF (returned for later averaging across blocks). The
    leaflet-to-leaflet spread is written to `out` for diagnostic purposes
    only; the reported uncertainty of the final PMF comes from the
    block-to-block spread (see combine_blocks).
    """
    z1, g1 = load_pmf(p1)
    z2, g2 = load_pmf(p2)
    if not np.allclose(z1, z2):
        sys.exit(f"{p1} and {p2} are on different z grids")
    z = z1
    g1 = shift_to_bulk(z, g1, z_bulk)
    g2 = shift_to_bulk(z, g2, z_bulk)

    # Keep the positive side from solute 1 and mirror the negative side from
    # solute 2, so both leaflets describe the same |z| coordinate.
    pos = z >= 0
    neg = z < 0
    zp = z[pos]
    zn = -z[neg]

    g_plus_leaflet = g1[pos]
    g_minus_leaflet = fold_side(zn, g2[neg], zp)

    stack = np.vstack([g_plus_leaflet, g_minus_leaflet])     # 2 x nbin
    n     = np.sum(np.isfinite(stack), axis=0)
    with np.errstate(invalid="ignore"):
        mean = np.divide(np.nansum(stack, axis=0), n,
                         out=np.full_like(zp, np.nan, dtype=float),
                         where=n > 0)
    with np.errstate(invalid="ignore"):
        std = np.where(n >= 2,
                       np.nanstd(stack, axis=0, ddof=1),
                       np.nan)
    se = std / np.sqrt(np.maximum(n, 1))

    finite = np.isfinite(mean)
    bulk   = finite & (zp >= z_bulk)
    ref    = (np.nanmean(mean[bulk]) if bulk.any()
              else (np.nanmean(mean[finite]) if finite.any() else 0.0))
    mean   = mean - ref

    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}', block '{block}' "
             f"(ALC-2007 two-leaflet averaging, {bin_width:g} A bins)\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_leaflets\n")
        fh.write(f"# zeroed on |z| >= {z_bulk:g} A; SE = std/sqrt(n)\n")
        fh.write("# leaflet estimates: solute 1 (+z) and mirrored solute 2 (-z)\n")
        fh.write("# diagnostic only -- final uncertainty comes from block-to-block spread\n")
        for zi, gi, si, ni in zip(zp, mean, se, n):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {int(ni):4d}\n")
    return zp, mean


def combine_blocks(z: np.ndarray, block_means: list[np.ndarray],
                   block_names: list[str], out: Path, analog: str
                   ) -> int:
    """Average the per-block (leaflet-folded) PMFs and estimate the
    standard error from the spread between the 50 ns blocks.
    """
    stack = np.vstack(block_means)                          # nblocks x nbin
    n = np.sum(np.isfinite(stack), axis=0)
    with np.errstate(invalid="ignore"):
        mean = np.divide(np.nansum(stack, axis=0), n,
                         out=np.full(z.shape, np.nan, dtype=float),
                         where=n > 0)
    with np.errstate(invalid="ignore"):
        std = np.where(n >= 2, np.nanstd(stack, axis=0, ddof=1), np.nan)
    se = std / np.sqrt(np.maximum(n, 1))

    with open(out, "w") as fh:
        fh.write(f"# PMF for analog '{analog}'  "
                 f"(averaged over {len(block_names)} 50 ns blocks: "
                 f"{', '.join(block_names)})\n")
        fh.write("# columns: |z|(A)   PMF(kcal/mol)   SE(kcal/mol)   n_blocks\n")
        fh.write("# each block PMF folds solute 1 (+z) and mirrored solute 2 (-z) leaflets;\n")
        fh.write("# SE here = std/sqrt(n) across blocks\n")
        for zi, gi, si, ni in zip(z, mean, se, n):
            fh.write(f"{zi:8.3f} {gi:12.5f} {si:12.5f} {int(ni):4d}\n")
    return int(np.isfinite(mean).sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analog",   required=True)
    ap.add_argument("--in-dir",   required=True, type=Path,
                    help="parent directory containing one subdirectory per "
                         "50 ns block (see --blocks), each with "
                         "metadata_s1.dat/metadata_s2.dat from extract_metadata.py")
    ap.add_argument("--out-dir",  type=Path, default=None,
                    help="where to write the combined pmf-us-<analog>.dat "
                         "(default: same as --in-dir)")
    ap.add_argument("--blocks", nargs="+",
                    default=["block1", "block2", "block3"],
                    help="block subdirectory names under --in-dir, in order "
                         "(default: block1 block2 block3, matching "
                         "section101-150 / section151-200 / section201-250)")
    ap.add_argument("--backend", choices=["python", "grossfield"],
                    default="python",
                    help="WHAM implementation (default: python, no external dep)")
    ap.add_argument("--wham-bin", default="wham",
                    help="path/name of Grossfield wham binary (--backend grossfield)")
    ap.add_argument("--hist-min", type=float, default=DEFAULT_HIST_MIN)
    ap.add_argument("--hist-max", type=float, default=DEFAULT_HIST_MAX)
    ap.add_argument("--num-bins", type=int, default=DEFAULT_NUM_BINS,
                    help="number of WHAM histogram bins (default 200)")
    ap.add_argument("--bin-A",    type=float, default=None,
                    help="bin width in A; overrides --num-bins when set")
    ap.add_argument("--tol",      type=float, default=DEFAULT_TOL)
    ap.add_argument("--temp",     type=float, default=DEFAULT_TEMP)
    ap.add_argument("--num-pad",  type=int,   default=DEFAULT_NUM_PAD)
    ap.add_argument("--z-bulk",   type=float, default=DEFAULT_Z_BULK)
    ap.add_argument("--nboot",    type=int,   default=0)
    ap.add_argument("--seed",     type=int,   default=DEFAULT_SEED)
    args = ap.parse_args()

    out_dir = args.out_dir or args.in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    hist_range = args.hist_max - args.hist_min
    if hist_range <= 0:
        sys.exit("ERROR: --hist-max must be greater than --hist-min")
    if args.bin_A is not None:
        if args.bin_A <= 0:
            sys.exit("ERROR: --bin-A must be positive")
        nbins = int(round(hist_range / args.bin_A))
        bin_width = hist_range / nbins
    else:
        if args.num_bins <= 0:
            sys.exit("ERROR: --num-bins must be positive")
        nbins = args.num_bins
        bin_width = hist_range / nbins

    block_z: np.ndarray | None = None
    block_means: list[np.ndarray] = []
    used_blocks: list[str] = []
    for block in args.blocks:
        block_dir = args.in_dir / block
        meta1 = block_dir / "metadata_s1.dat"
        meta2 = block_dir / "metadata_s2.dat"
        missing = [m for m in (meta1, meta2) if not m.is_file() or m.stat().st_size == 0]
        if missing:
            print(f"  ! block '{block}': missing/empty metadata "
                  f"({', '.join(str(m) for m in missing)}), skipped")
            continue

        pmf1 = block_dir / "pmf_s1.dat"
        pmf2 = block_dir / "pmf_s2.dat"
        log1 = block_dir / "wham_s1.log"
        log2 = block_dir / "wham_s2.log"

        print(f"  [{block}] WHAM s1: bins={nbins} ({bin_width:g} A) "
              f"range=[{args.hist_min},{args.hist_max}]  T={args.temp}  "
              f"backend={args.backend}")
        if args.backend == "python":
            wham_py(meta1, pmf1, log1, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            wham_py(meta2, pmf2, log2, args.hist_min, args.hist_max, nbins,
                    args.temp, args.tol)
        else:
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta1, pmf1, log1, args.nboot, args.seed)
            print(f"  [{block}] WHAM s2: bins={nbins}")
            run_wham(args.wham_bin, args.hist_min, args.hist_max, nbins, args.tol,
                     args.temp, args.num_pad, meta2, pmf2, log2, args.nboot, args.seed)

        block_pmf = block_dir / f"pmf-us-{args.analog}-{block}.dat"
        zp, mean = combine_leaflets(pmf1, pmf2, block_pmf, args.analog, block,
                                    args.z_bulk, bin_width)
        if block_z is None:
            block_z = zp
        elif not np.allclose(block_z, zp):
            sys.exit(f"ERROR: block '{block}' PMF grid differs from earlier blocks")
        block_means.append(mean)
        used_blocks.append(block)
        print(f"  [{block}] -> {block_pmf}")

    if not block_means:
        sys.exit("ERROR: no usable blocks found under " + str(args.in_dir))

    combined = out_dir / f"pmf-us-{args.analog}.dat"
    nbin_ok = combine_blocks(block_z, block_means, used_blocks, combined, args.analog)
    print(f"==> {combined}   ({nbin_ok} usable bins, {len(used_blocks)} blocks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

