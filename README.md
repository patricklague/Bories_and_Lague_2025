# Bories and Lagüe, 2026

Supporting data and scripts for the article.

---
## 1. Data Source and Citation

- Repository: Github
- DCD trajectory location : <Borealis_link> (DOI : unavalaible)
- Companion article:
  - Authors: `S. Bories and P. Lague`
  - Title: `Amino Acid Insertion Energetics in a POPC Bilayer from Unbiased Molecular Dynamics`
  - Journal / preprint: `...`
  - Year: `2026`

If you use these data, please cite both the Borealis dataset and the companion article above.

## figures/

Data, scripts, and outputs for all figures in the article.

- **data/** — Raw and processed data used for plotting
  - `pmf_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/trajectory{1,2,3}.dat` — per-trajectory PMFs
  - `distribution_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/trajectory{1,2,3}.dat` — per-trajectory z-density distributions, plus `raw_data/` (per-frame contact counts)
  - `densityProfile-popc/` — POPC component density profiles
  - `aromatics_orientation/` — ring-orientation raw data, per-trajectory `vector_orientations/`, and the aggregated `freq_angle_<aa>.dat` 2D histograms
  - `SUPP_membrane_parm/` — per-analog `thickness/`, `area_per_lipid/`, `order_parameter/`, `densityProfiles/` tables and the top-level `computed_*.csv` summaries
  - `SUPP_monomer/` — `monomer_rates_45A_9batches.dat` (mean ± std monomer/multimer rates over 9 batches)
  - `SUPP_hydrophobicity/`, `extra_analysis/` — auxiliary tables for the supplementary figures
- **scripts/** — Jupyter notebooks that generate the figures (`PMF_plots.ipynb`, `density_popc_plot.ipynb`, `pKa_plot.ipynb`, `aromatics_orientation_plots.ipynb`, `SUPP_Membrane_parm.ipynb`, `SUPP_monomer_plots.ipynb`, `SUPP_hydrophobicity_plots.ipynb`, `SUPP_aromatics_orientation_plots.ipynb`) plus helper scripts under `aromatics_orientation/` and PyMOL scripts under `pymol/`
- **plot/** — Output figure files (PNG)
- **draw/** — Hand-drawn figures (Pages files for article and abstract figures, ring orientation diagrams)

---

## supp_files/

Supplementary files: simulation inputs, outputs, parameters, and method scripts.

- **input_systems/** — Input PSF/PDB files for each amino acid analog system (one subdirectory per analog: sca, scc, ..., scym) and the pure POPC bilayer
- **output_systems/** — Output trajectory frames for each system (one subdirectory per analog).
- **parameter_files/** — CHARMM and NAMD topology/parameter files for both standard lipids (`toppar-charmm/`, `toppar-namd/`) and side-chain analogs (`sidechain-toppar-charmm/`, `sidechain-toppar-namd/`)
- **popc/** — Pure POPC reference data: density distributions and PMFs
- **method_script/** — All analysis scripts organized in five sequential steps (see [method_script/readme.md](supp_files/method_script/readme.md)):
  1. `1-system_generation/` — Build bilayer–solute systems
  2. `2-system_analysis/` — Trajectory analysis split into `membrane_parm_analysis/` (centering, cell dimensions, thickness, density profiles, deuterium order parameters; supports the 401-600/601-800/801-1000 ns batch split) and `aromatic_analysis/` (per-trajectory ring-orientation CSVs and aggregated 2D depth/angle histograms for SCF/SCY/SCW)
  3. `3-distribution_extraction/` — Extract solute z-density distributions and identify monomers via the 4.5 Å contact cutoff
  4. `4-pmf_calculation/` — Compute potentials of mean force from the distributions
  5. `5-postprocessing/` — `format_borealis_trajectories.sh` (build the public Borealis trajectory bundle from the local NAMD output), helper scripts to extract output frames (`get_last_frames.sh`, `get_frame.tcl`), `SUPP_membrane_parm/` per-parameter `run.sh` extractors that write into `figures/data/SUPP_membrane_parm/`, and the top-level `compute_*.py` aggregators that produce the `computed_*.csv` summaries; also contains `SUPP_monomer/monomer_rate.py` (writes `monomer_rates_45A_9batches.dat`)

The trajectories themselves are distributed on Borealis as `POPC-<AA>/trajectory{1,2,3}/section{401-500,501-600,601-700,701-800,801-900,901-1000}.dcd` plus the matching `popc-<aa>.psf`. See [supp_files/method_script/readme.md](supp_files/method_script/readme.md#running-on-borealis-formatted-data) for the `catdcd` substitutions required to drive steps 2–4 directly from the Borealis archive.
