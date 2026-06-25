import numpy as np
import os

# ============================================================
# USER PARAMETERS — modify these as needed
# ============================================================

# List of analogs to process
analogs = [
    "sca", "scv", "scl", "sci", "scc", "scm", "scs", "sct", "scq",
    "scn", "scf", "scy", "scw", "scp", "glyd", "schd", "sche", "scdn",
    "scen", "sckn", "scrn", "schp", "sccm", "scym", "scd", "sce", "sck", "scr", "scrn-1", "scw-1"
    ]
# Number of trajectories per analog
n_trajectories = 3

# Counting mode: 'monomer', 'multimer', or 'total'
mode = "total"

# Cutoff column to use (used for 'monomer' and 'multimer' modes)
# Options: "4.5A_cutoff"
cutoff = "4.5A_cutoff"

# Z-bin width in angstroms
z_bin = 1.0

# Batch definitions: (label, start_frame, end_frame exclusive)
# Each 200 ns = 20,000 frames; raw data contains frames 0–59999 (400–1000 ns)
frames_per_batch = 20000
batches = [
    ("400-600ns",  0,     20000),
    ("600-800ns",  20000, 40000),
    ("800-1000ns", 40000, 60000),
]

# Input / output directories
raw_data_dir = "raw_data"
output_dir = f"{mode}_{cutoff.split('_')[0]}" if mode != "total" else mode

# ============================================================
# Column mapping
# ============================================================
cutoff_col_map = {
    "4.5A_cutoff": 5,
}

frame_col = 0  # column index of frame number
z_col = 4      # column index of z in the data file

# ============================================================
# Processing
# ============================================================

cutoff_idx = cutoff_col_map[cutoff]

for analog in analogs:
    out_analog_dir = os.path.join(output_dir, analog)
    os.makedirs(out_analog_dir, exist_ok=True)

    # Collect z-values per (trajectory, batch) for two-pass processing
    all_z_data = {}  # (traj, batch_idx) -> z_values array

    # --- Pass 1: read all trajectories & store filtered z-values per batch ---
    for traj in range(1, n_trajectories + 1):
        input_file = os.path.join(raw_data_dir, analog, f"{analog}_contacts_{traj}.dat")
        if not os.path.isfile(input_file):
            print(f"WARNING: {input_file} not found, skipping.")
            continue

        data = np.loadtxt(input_file, skiprows=1)
        frames = data[:, frame_col]
        z_all = data[:, z_col]

        # Mode filter
        if mode == "monomer":
            mode_mask = data[:, cutoff_idx] == 0
        elif mode == "multimer":
            mode_mask = data[:, cutoff_idx] != 0
        else:
            mode_mask = np.ones(len(data), dtype=bool)

        for b_idx, (_, b_start, b_end) in enumerate(batches):
            batch_mask = (frames >= b_start) & (frames < b_end)
            all_z_data[(traj, b_idx)] = z_all[mode_mask & batch_mask]

    if not all_z_data:
        print(f"WARNING: No data for {analog}, skipping.")
        continue

    # --- Common z-grid across all trajectories and batches ---
    non_empty = [z for z in all_z_data.values() if len(z) > 0]
    if not non_empty:
        print(f"WARNING: All batches empty for {analog}, skipping.")
        continue
    global_z_min = np.floor(min(z.min() for z in non_empty))
    global_z_max = np.ceil(max(z.max() for z in non_empty))
    bins = np.arange(global_z_min, global_z_max + z_bin, z_bin)
    bin_centers = (bins[:-1] + bins[1:]) / 2.0
    n_bins = len(bin_centers)

    # --- Pass 2: histogram each batch & write per-trajectory files ---
    all_densities = []  # collect all 9 density arrays for the summary

    for traj in range(1, n_trajectories + 1):
        traj_densities = []
        for b_idx in range(len(batches)):
            key = (traj, b_idx)
            if key in all_z_data and len(all_z_data[key]) > 0:
                counts, _ = np.histogram(all_z_data[key], bins=bins)
                total = counts.sum()
                density = counts / total if total > 0 else np.zeros(n_bins)
            else:
                density = np.zeros(n_bins)
            traj_densities.append(density)
            all_densities.append(density)

        # trajectory file: z | batch1 | batch2 | batch3
        output_file = os.path.join(out_analog_dir, f"trajectory{traj}.dat")
        header = f"{'z':>12s} {'400-600ns':>12s} {'600-800ns':>12s} {'800-1000ns':>12s}"
        np.savetxt(
            output_file,
            np.column_stack([bin_centers] + traj_densities),
            fmt="%12.6f",
            header=header,
            comments="",
        )
        print(f"Written {output_file}")

    # --- Summary: mean and SE across all 9 batches ---
    all_densities = np.array(all_densities)  # shape (9, n_bins)
    mean_density = np.mean(all_densities, axis=0)
    se_density = np.std(all_densities, axis=0, ddof=1) / np.sqrt(all_densities.shape[0])

    summary_file = os.path.join(out_analog_dir, f"summary_{analog}.dat")
    header = f"{'z':>12s} {'mean':>12s} {'se':>12s}"
    np.savetxt(
        summary_file,
        np.column_stack([bin_centers, mean_density, se_density]),
        fmt="%12.6f",
        header=header,
        comments="",
    )
    print(f"Written {summary_file}")

print("Done.")
