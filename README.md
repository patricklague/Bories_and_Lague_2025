# Bories and Lagüe, 2026

Supporting data and scripts for the article.

---
## 1. Data Source and Citation

- Repository: Github
- Companion article:
  - Authors: `S. Bories and P. Lague`
  - Title: `Amino Acid Insertion Energetics in a POPC Bilayer from Unbiased Molecular Dynamics`
  - Journal / preprint: `...`
  - Year: `2026`

If you use these data, please cite the companion article above.

## figures/

Data, scripts, and outputs for all figures in the article.

- **data/** — Raw and processed data used for plotting
  - `pmf_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/trajectory{1,2,3}.dat` — per-trajectory PMFs
  - `distribution_data/{total,monomer_4.5A,multimer_4.5A}/<aa>/trajectory{1,2,3}.dat` — per-trajectory z-density distributions, plus `raw_data/` (per-frame contact counts)
  - `densityProfile-popc/` — POPC component density profiles
  - `aromatics_orientation/` — ring-orientation raw data, per-trajectory `vector_orientations/`, and the aggregated `freq_angle_<aa>.dat` 2D histograms
  - `SUPP_membrane_parm/` — per-analog `thickness/`, `area_per_lipid/`, `order_parameter/`, `densityProfiles/` tables and the top-level `computed_*.csv` summaries
  - `SUPP_monomer/` — `monomer_rates_45A_9batches.dat` (mean ± std monomer/multimer rates)
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
  5. `5-postprocessing/` — helper scripts to extract output frames (`get_last_frames.sh`, `get_frame.tcl`), `SUPP_membrane_parm/` per-parameter `run.sh` extractors that write into `figures/data/SUPP_membrane_parm/`, and the top-level `compute_*.py` aggregators that produce the `computed_*.csv` summaries; also contains `SUPP_monomer/monomer_rate.py` (writes `monomer_rates_45A_9batches.dat`)

The dependencies required to re-run the analysis pipeline (NAMD, VMD, Packmol, …) are listed in the [Dependencies](#dependencies) section below.

---

## Dependencies

The versions listed below are the ones actually used to produce the published data (Python packages taken from the repository's `.venv`). Any recent compatible release should work.

### External tools

| Tool | Version | Used by | Purpose |
| --- | --- | --- | --- |
| **NAMD** | 2.14 (Alliance Canada) / NAMD3 (local) | `supp_files/method_script/1-system_generation/submit-alliancecan.sh`, `submit-local.sh` | Run the MD simulations |
| **VMD** | ≥ 1.9.4 | steps 1–3 (invoked as `vmd -dispdev text -e …`) | Trajectory manipulation (centering, atom selection), `psfgen`, and the Density Profile Tool plugin |
| **psfgen** | bundled with VMD ≥ 1.9.4 | `1-system_generation/build/psfgen.inp` | Combine solutes with the pre-built POPC bilayer into the final PSF/PDB |
| **VMD Density Profile Tool** (`density_profile.tcl`) | bundled with VMD | `2-system_analysis/membrane_parm_analysis/densityProfiles-*.vmd` | Compute POPC component and water density profiles |
| **catdcd** | 4.0 (bundled with VMD) | `analysis_per_system.sh`, `get_trajectory.sh` | Concatenate/stride the per-section production DCDs |
| **Packmol** | ≥ 20 | `1-system_generation/build/packmol-POPC.sh` | Place solutes in the simulation box before assembling the system |
| **CHARMM-GUI *Membrane Builder*** | output shipped as `1-system_generation/charmm-gui.tgz` | step 1 | Provide the base POPC bilayer system |
| **PyMOL** | ≥ 2.5 | `figures/scripts/Pymol_scripts/` (`.pml` scenes and their `.py` helpers) | Render the molecular scenes |

### Python

Tested with **Python 3.13.9** (Jupyter notebooks and all `.py` scripts). Python ≥ 3.10 is expected to work.

| Package | Version | Used by | Purpose |
| --- | --- | --- | --- |
| `numpy` | 2.3.5 | all figure notebooks, all analysis `.py` scripts, `counting_script.ipynb`, `pymol_orient_molecule.py` | Numerical arrays and math |
| `pandas` | 2.3.3 | all figure notebooks and every `5-postprocessing/**` script | Tabular I/O and aggregation |
| `matplotlib` | 3.10.6 | all `figures/scripts/*.ipynb` (incl. `matplotlib.ticker`, `matplotlib.lines`, `matplotlib.patches`) | Plot rendering |
| `scipy` | 1.16.3 | `figures/scripts/FigureS10_plot.ipynb` (`scipy.stats.linregress`); `5-postprocessing/SUPP_membrane_parm/compute_{density,order}_deviation.py` (`scipy.interpolate.interp1d`) | Regression and 1D interpolation |
| `scikit-learn` | 1.7.2 | `figures/scripts/FigureS10_plot.ipynb` (`LinearRegression`, `mean_squared_error`) | Linear regression diagnostics |
| `MDAnalysis` | any recent release | `supp_files/method_script/3-distribution_extraction/counting_script.ipynb` | Trajectory parsing and per-frame contact counts (`distance_array`) |
| `polars` | any recent release | `supp_files/method_script/3-distribution_extraction/counting_script.ipynb` | Fast DataFrame processing of the per-frame contact tables |
| `jupyter` / `jupyterlab` / `notebook` | 1.1.1 / 4.4.7 / 7.4.5 | all `*.ipynb` files under `figures/scripts/` and `supp_files/method_script/3-distribution_extraction/` | Run the notebooks |

A minimal Python install covering everything is, for example:

```bash
pip install numpy pandas matplotlib scipy scikit-learn MDAnalysis polars jupyter
```

Any HPC-specific pieces (SLURM directives in `submit-alliancecan.sh`, the `module load namd-multicore/2.14` line) can be adapted or ignored on a local workstation.
