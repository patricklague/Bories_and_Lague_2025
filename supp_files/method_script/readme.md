# Method Scripts

This directory contains all the scripts used to generate, analyze, and process the simulation data.
The workflow is divided into five sequential steps, each in its own subdirectory.

The three production trajectories per system are stored on Borealis (DOI unavailable) as
`POPC-<AA>/trajectory{1,2,3}/section{401-500,501-600,601-700,701-800,801-900,901-1000}.dcd`
plus the matching `popc-<aa>.psf`. See [Running on Borealis-formatted data](#running-on-borealis-formatted-data) at the end of this file for the exact command substitutions required to drive steps 2–4 from the public archive.

The full list of dependencies (NAMD, VMD, catdcd, Packmol, Python packages, …) for both the pipeline below and the figure notebooks is given in the top-level [README.md](../../README.md#dependencies).

---

## 1-system_generation/

**Purpose:** Build the bilayer–solute systems (PSF/PDB files) for NAMD simulations.

**Entry point:** `build/run.sh`

`run.sh` loops over all amino acid analogs and calls `main.sh` for each one. `main.sh` orchestrates the full build in three steps:
1. Distribute solutes in the simulation box using `packmol-POPC.sh` (with solute PDB files from `solutes/`)
2. Combine solutes with the bilayer system using `psfgen.inp` (VMD/psfgen)
3. Generate the extra-bonds file using `extraBondsFile.py`

**Contents:**
- `build/` — Build scripts (`run.sh`, `main.sh`, `packmol-POPC.sh`, `psfgen.inp`, `extraBondsFile.py`, `update_step7.sh`)
- `solutes/` — PDB files for all amino acid side-chain analogs and associated topology/parameter files
- `charmm-gui.tgz` — Base POPC bilayer system from CHARMM-GUI

---

## 2-system_analysis/

**Purpose:** Analyze the MD trajectories: cell dimensions, membrane thickness, density profiles, lipid order parameters, and aromatic-ring orientations.

This step is split into two sibling pipelines:

### `membrane_parm_analysis/`

**Entry point:** `analysis-run.sh`

`analysis-run.sh` loops over all analogs and trajectories, substituting placeholders in `analysis_per_system.sh` and running it. `analysis_per_system.sh` performs:
1. Trajectory centering on the bilayer mid-plane (uses `indexNoWat.vmd`, `center.vmd`)
2. Cell-dimension extraction from per-section `.xsc` files (not shared)
3. Membrane thickness (uses `thickness.vmd`)
4. POPC component density profiles (uses `densityProfiles-popc.vmd` or `densityProfiles-water-total.vmd`)
5. Lipid acyl-chain deuterium order parameters (uses `orders-popc.vmd`)

When batch analysis is enabled (`batches=1`), the production window 401–1000 ns is split into the three 200 ns batches `401-600 / 601-800 / 801-1000`. Each batch produces its own centered DCD (`batch401-600.dcd`, `batch601-800.dcd`, `batch801-1000.dcd`) and the density / order-parameter outputs are tagged with the same suffix.

The density-profile VMD scripts depend on `density_profile.tcl` (VMD Density Profile Tool plugin, sourced at runtime).

### `aromatic_analysis/`

**Entry point:** `run.sh`

Computes the per-frame aromatic-ring atom coordinates (SCF, SCY, SCW) needed for the ring-orientation analysis. `run.sh` drives `get_trajectory.sh` (water/membrane stripping + bilayer centering, same recipe as step 3) and the post-processing scripts:

- `ring_orientation_analysis.py` — produces one `orientation_<aa>_traj{1,2,3}.csv` per trajectory (ring-normal and pair-vector angles vs Z, plus centroid depth) under `figures/data/aromatics_orientation/raw_data/total/<AA>/vector_orientations/`.
- `freq_angle_analysis.py` — averages a 2D depth/angle histogram over the 9 batches (3 trajectories × 3 frame windows) and writes `figures/data/aromatics_orientation/freq_angle_<aa>.dat`.
- `extract_top_angles.sh`, `Ring_xyz_preparation.ipynb` — auxiliary helpers used to prepare the input atom-coordinate files.

---

## 3-distribution_extraction/

**Purpose:** Extract per-frame solute positions from the centered trajectories and bin them into z-density distributions, identifying monomers from multimers via a contact cutoff.

**Entry point:** `counting_script.ipynb`

`counting_script.ipynb` is the master notebook that drives the whole step. It launches `get_trajectory.sh` for each analog/trajectory, which in turn:
1. Removes water/ions (`indexNoWat.vmd`)
2. Centers the bilayer (`center.vmd`)
3. Removes the membrane (`indexNoMemb.vmd`)
4. Writes the per-frame solute positions and per-frame contact counts at the 4.5 Å cutoff.

Once the raw contact files are produced, `raw_to_distribution.py` bins them into the per-trajectory z-density distributions consumed by step 4. Outputs land under `figures/data/distribution_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/trajectory{1,2,3}.dat`.

---

## 4-pmf_calculation/

**Purpose:** Compute potentials of mean force (PMFs) from the density distributions extracted in step 3.

**Entry point:** `pmf-from-distribution.py`

`pmf-from-distribution.py` loops over all analogs, copies the relevant distribution files from step 3, and stores the output PMF data. The script handles the `total`, `monomer_4.5A`, and `multimer_4.5A` distribution modes through its parameters. Outputs land under `figures/data/pmf_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/`.

---

## 5-postprocessing/

**Purpose:** Take the raw outputs of steps 1–4 and turn them into the artifacts that are actually shipped (Borealis trajectory bundle + the curated `figures/data/` tree). This step is the only one allowed to write outside `supp_files/method_script/`.

**Contents:**
- `SUPP_membrane_parm/` — Per-system extraction and aggregation scripts that turn the raw `analysis_per_system.sh` outputs into the canonical `figures/data/SUPP_membrane_parm/` tables.
  - `thickness/`, `area_per_lipid/`, `order_parameter/`, `densityProfiles/` each contain a `run.sh` driver and a `get_*.py` extractor that merge the three per-trajectory raw files into a single per-analog `<aa>-<param>.dat`. Each `run.sh` keeps the original local source path as the active `DIR=...` and exposes a commented reviewer-facing alternative pointing at `../../../../results/POPC-aa/POPC-$aa/`. The merged file is moved to `figures/data/SUPP_membrane_parm/<subdir>/`.
  - `compute_thickness.py`, `compute_apl.py`, `compute_acm.py`, `compute_density_deviation.py`, `compute_order_deviation.py` — top-level aggregators that read the per-analog tables under `figures/data/SUPP_membrane_parm/<subdir>/` and produce the `computed_*.csv` summary files used by the figure notebooks.

All scripts in this step apply the canonical naming conventions defined in `correctives/list_of_script.md` and emit headers compatible with `correctives/fix_header.py` so that no post-hoc header rewriting is needed.

