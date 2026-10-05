# Single-chain simulations

This directory contains the workflow for building and simulating individual
polymer chains. The scripts are based on
[BIO-SUSHY](https://github.com/daimoners/BIO-SUSHY-tutorials).

## Files

### Simulation workflow

- [`simulate.py`](./simulate.py): shared implementation for loading
  BIO-SUSHY and running one single-chain simulation. It builds the polymer,
  generates the force field, runs the vacuum MD simulation, and stores the
  workflow state and results.
- [`main.py`](./main.py): example entrypoint for one polyethylene simulation.
  It uses a fixed SMILES, chain length, temperature, and number of steps, then
  prints the mean radius of gyration.
- [`simBatch.py`](./simBatch.py): batch entrypoint that reads polymers from the
  MongoDB `biceranoPolymers` collection and runs the same workflow for each
  polymer.
- [`environment.yml`](./environment.yml): Conda environment used for the
  single-chain workflow. It is intended for Linux or WSL2.
- `BIO-SUSHY-tutorials/`: local checkout of the BIO-SUSHY dependency. It is
  downloaded or updated by `simulate.py` and is ignored by Git.

The example entrypoints write generated outputs to `simulation_results/` and
`simulation_resultsHf/`. These directories are local run outputs and should
not be committed.

### Radius-of-gyration utilities

- [`bondAutocorrelation/bondsList.py`](./bondAutocorrelation/bondsList.py):
  identifies the polymer backbone and writes the bond index file needed by the
  single-chain feature workflow.
- [`bondAutocorrelation/run_rotacf_pipeline.py`](./bondAutocorrelation/run_rotacf_pipeline.py):
  processes completed single-chain results and runs the required GROMACS
  conversion and rotational-correlation steps.
- [`bondAutocorrelation/dummy.mdp`](./bondAutocorrelation/dummy.mdp):
  minimal GROMACS parameter file used to generate the topology needed by the
  post-processing command.

The `bondAutocorrelation/` name is retained for compatibility with the
existing scripts. Its utilities are part of the single-chain workflow and are
intended to support the future radius-of-gyration (Rg) feature.

## Running

Create the environment:

```bash
conda env create -f environment.yml
conda activate simulation
```

Run the example from this directory:

```bash
python main.py
```

Run the MongoDB batch workflow from this directory:

```bash
python simBatch.py
```

