# Polymer prediction

Repository for the prediction of polymer properties with various ML techniques. 
The project combines dataset preparation, molecular-dynamics simulations, feature generation, exploratory analysis, and machine-learning experiments.

## Project structure

```text
.
├── data/              # Versioned input datasets and generated feature tables
├── datasetToJson/     # Dataset ETL and MongoDB ingestion
├── notebooks/         # Notebooks containing model experiments and exploratory analysis
├── simulations/       # Single-chain and bulk molecular-dynamics workflows
├── results/           # Local generated results, ignored by Git
├── runs/              # Local simulation runs, ignored by Git
└── runlog.py          # Utility for recording or inspecting project runs
```

More detailed documentation is available in:

- [`datasetToJson/README.md`](./datasetToJson/README.md): ETL pipeline and
  MongoDB storage.
- [`simulations/README.md`](./simulations/README.md): simulation workflows.
- [`simulations/single_chain/README.md`](./simulations/single_chain/README.md):
  single-chain files and commands.
- [`simulations/bulk/README.md`](./simulations/bulk/README.md): bulk workflows
  and cluster/SLURM files.

## Main workflows

### 1. Dataset preparation

The [`datasetToJson/`](./datasetToJson/) pipeline converts heterogeneous
datasets into a common JSON schema and stores the resulting documents in
MongoDB. Dataset-specific mappings are implemented in
[`datasetToJson/wrappers/`](./datasetToJson/wrappers/).

The database connection is configured through environment variables. Do not
commit credentials or local `.env` files.

### 2. Molecular-dynamics features

The [`simulations/`](./simulations/) directory contains two distinct
workflows:

- `single_chain/`: local or WSL2 workflow for individual polymer chains and
  chain-level properties such as radius of gyration.
- `bulk/`: bulk workflow originally executed on the CNR Ariadne cluster through SLURM. 

The generated feature tables are stored under
[`data/md_features/`](./data/md_features/).

### 3. Analysis and prediction

The [`notebooks/`](./notebooks/) directory contains experiments for:

- exploratory property distributions and statistical analysis;
- gradient-boosting baselines;
- Chemprop and GNN models;
- ChemBERTa embeddings;
- comparison of molecular-dynamics features with measured properties.

The notebooks are exploratory artifacts rather than a single automated
training pipeline. Check the notebook cells for the input dataset and output
locations used by each experiment.

## Environment setup

There is no single root-level environment file. Use the environment associated
with each specific workflow.

