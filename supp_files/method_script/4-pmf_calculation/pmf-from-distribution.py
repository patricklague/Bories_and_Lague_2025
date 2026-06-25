import numpy as np
import os

# ============================================================
# USER PARAMETERS
# ============================================================

analogs = [
    "sca", "scv", "scl", "sci", "scc", "scm", "scs", "sct", "scq",
    "scn", "scf", "scy", "scw", "scp", "glyd", "schd", "sche", "scdn",
    "scen", "sckn", "scrn", "schp", "sccm", "scym", "scd", "sce", "sck", "scr",
    "scrn-1", "scw-1"
]
n_trajectories = 3
n_batches = 3  # columns in each trajectory file (400-600, 600-800, 800-1000 ns)

# Distribution mode directory (e.g. "monomer_4.5A", "total", ...)
mode = "total"  # used to find distribution files and name output dir

# PMF constants
k_B = 0.008314  # kJ/(mol*K)
T = 303.15
epsilon = 1e-5

# PMF recalibration window (angstroms)
recal_lo = 40.0
recal_hi = 50.0

# Z-bin width (for extending the grid)
z_bin = 1.0

# Paths
distribution_dir = os.path.join("..", "distribution_data", mode)
batch_labels = ["400-600ns", "600-800ns", "800-1000ns"]

# ============================================================
# Helper: density -> symmetrised PMF on |z| grid
# ============================================================

def density_to_pmf(z, density, z_pos):
    """
    Symmetrise density(z) around z=0 onto |z| grid, compute PMF,
    recalibrate so <PMF> = 0 in the [recal_lo, recal_hi] window.
    Returns pmf array aligned with z_pos.
    """
    # Map density onto |z| grid: average neg and pos sides
    sym_density = np.zeros_like(z_pos, dtype=float)
    counts = np.zeros_like(z_pos, dtype=float)

    for zi, di in zip(z, density):
        az = abs(zi)
        idx = np.argmin(np.abs(z_pos - az))
        sym_density[idx] += di
        counts[idx] += 1.0

    # Average where we have contributions from both sides
    nonzero = counts > 0
    sym_density[nonzero] /= counts[nonzero]

    # Re-normalise so sum = 1 on the half-grid
    total = sym_density.sum()
    if total > 0:
        sym_density /= total

    # PMF = -kT ln(p + epsilon)
    pmf_raw = -k_B * T * np.log(sym_density + epsilon)

    # Recalibrate: shift so mean PMF in [recal_lo, recal_hi] = 0
    mask_win = (z_pos >= recal_lo) & (z_pos <= recal_hi)
    if mask_win.any():
        pmf_raw -= pmf_raw[mask_win].mean()

    return pmf_raw

# ============================================================
# Processing
# ============================================================

for analog in analogs:
    dist_analog_dir = os.path.join(distribution_dir, analog)
    out_analog_dir = os.path.join(mode, analog)
    os.makedirs(out_analog_dir, exist_ok=True)

    # --- Pass 1: read all trajectory files, build common |z| grid ---
    traj_data = {}  # traj -> (z_array, density_matrix with n_batches cols)
    all_abs_z = set()

    for traj in range(1, n_trajectories + 1):
        fpath = os.path.join(dist_analog_dir, f"trajectory{traj}.dat")
        if not os.path.isfile(fpath):
            print(f"WARNING: {fpath} not found, skipping.")
            continue
        arr = np.loadtxt(fpath, skiprows=1)  # z, batch1, batch2, batch3
        z = arr[:, 0]
        densities = arr[:, 1:]  # shape (n_z, 3)
        traj_data[traj] = (z, densities)
        all_abs_z.update(np.abs(z).tolist())

    if not traj_data:
        print(f"WARNING: No trajectory data for {analog}, skipping.")
        continue

    # Common |z| grid (positive only), extended to 50 A
    z_pos = np.array(sorted(all_abs_z))
    z_pos = z_pos[z_pos >= 0]
    if z_pos[-1] < 50:
        extra = np.arange(z_pos[-1] + z_bin, 50 + z_bin, z_bin)
        z_pos = np.concatenate([z_pos, extra])
    n_z = len(z_pos)

    # --- Pass 2: compute PMF per batch, write per-trajectory files ---
    all_pmfs = []  # collect all 9 PMF arrays for the summary

    for traj in range(1, n_trajectories + 1):
        if traj not in traj_data:
            continue
        z, densities = traj_data[traj]
        traj_pmfs = []
        for b in range(n_batches):
            pmf = density_to_pmf(z, densities[:, b], z_pos)
            traj_pmfs.append(pmf)
            all_pmfs.append(pmf)

        # Write trajectory file: z | PMF_batch1 | PMF_batch2 | PMF_batch3
        output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
        header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
        np.savetxt(
            output_file,
            np.column_stack([z_pos] + traj_pmfs),
            fmt="%12.6f",
            header=header,
            comments="",
        )
        print(f"Written {output_file}")

    # --- Summary: mean and SE of PMF across all 9 batches ---
    all_pmfs = np.array(all_pmfs)  # shape (up to 9, n_z)
    mean_pmf = np.mean(all_pmfs, axis=0)
    se_pmf = np.std(all_pmfs, axis=0, ddof=1) / np.sqrt(all_pmfs.shape[0])

    summary_file = os.path.join(out_analog_dir, f"pmf_{analog}.dat")
    header = f"{'z':>12s} {'mean':>12s} {'se':>12s}"
    np.savetxt(
        summary_file,
        np.column_stack([z_pos, mean_pmf, se_pmf]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {summary_file}")

print("Done.")

