# Bulk simulations

This directory contains the bulk-polymer workflows that were executed on a
cluster using SLURM.

## Job submission and simulation

- [`create_bulk_jobs.py`](./create_bulk_jobs.py): reads polymers from MongoDB
  and submits one bulk simulation job per polymer.
- [`bulk_job_template.sbatch`](./bulk_job_template.sbatch): SLURM template
  executed for each bulk simulation job.
- [`run_bulk_from_smiles.py`](./run_bulk_from_smiles.py): runs a bulk workflow
  starting from a SMILES input.
- [`run_bulk_analysis.py`](./run_bulk_analysis.py): runs the bulk analysis
  workflow on completed simulations.
- [`bulk_analyzer.py`](./bulk_analyzer.py): analysis implementation for bulk
  simulation outputs.
- [`create_analysis_jobs.py`](./create_analysis_jobs.py): finds completed
  bulk results and submits their analysis jobs.
- [`analysis_job_template.sbatch`](./analysis_job_template.sbatch): SLURM
  template for bulk analysis jobs.

## Data preparation and exploration

- [`polymer_properties.py`](./polymer_properties.py): utilities for extracting
  and computing polymer properties. Originally from BIO-SUSHY but edited for the specific workflow.
- [`make_csv.py`](./make_csv.py): creates CSV datasets from the available
  bulk results.
- [`data_exploration.py`](./data_exploration.py): exploratory inspection of
  bulk data.
- [`baseEnvironment.yml`](./baseEnvironment.yml): Conda environment used by
  the bulk workflow.
- [`mongo.env`](./mongo.env): local MongoDB configuration. It contains
  credentials/configuration and must remain untracked.


## Subdirectories

### `single_polymer/`

Scripts for running or analysing a single polymer within the bulk-related
workflow:

- [`simulate_polymer.py`](./single_polymer/simulate_polymer.py): simulation
  setup for one polymer.
- [`run_one_analysis.py`](./single_polymer/run_one_analysis.py): runs the
  analysis for one completed polymer simulation.
- [`run_one_rotacf.py`](./single_polymer/run_one_rotacf.py): runs rotational
  autocorrelation processing for one polymer.
- [`rotacf_utils.py`](./single_polymer/rotacf_utils.py): shared utilities for
  rotational autocorrelation processing.
- [`bondsList.py`](./single_polymer/bondsList.py): generates the bond index
  file used by the single-polymer analysis.
- [`make_dataset.py`](./single_polymer/make_dataset.py): assembles the
  resulting single-polymer measurements into a dataset.
- [`create_Rg_jobs.py`](./single_polymer/create_Rg_jobs.py): submits Rg jobs.
- [`create_rotacf_jobs.py`](./single_polymer/create_rotacf_jobs.py): submits
  rotational autocorrelation jobs.
- [`Rg_job_template.sbatch`](./single_polymer/Rg_job_template.sbatch): SLURM
  template for Rg jobs.
- [`rotacf_job_template.sbatch`](./single_polymer/rotacf_job_template.sbatch):
  SLURM template for rotational autocorrelation jobs.
- [`dummy.mdp`](./single_polymer/dummy.mdp): minimal GROMACS input used by the
  analysis workflow.

### `test_convergence/`

Scripts for checking the convergence of bulk features over time:

- [`extract_rg_timeseries.py`](./test_convergence/extract_rg_timeseries.py):
  extracts radius-of-gyration time series with GROMACS.
- [`extract_density_timeseries.py`](./test_convergence/extract_density_timeseries.py):
  extracts density time series with GROMACS.
- [`rg_ts.sbatch`](./test_convergence/rg_ts.sbatch): SLURM job template for
  radius-of-gyration time series.
- [`density_ts.sbatch`](./test_convergence/density_ts.sbatch): SLURM job
  template for density time series.
