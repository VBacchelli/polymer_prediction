This folder contains workflows to build and run simple molecular dynamics simulations of polymer systems.

The workflows are grouped by system type:

- `single_chain/`: scripts, configuration, and utilities for single-chain simulations.
- `bulk/`: workflows for bulk simulations executed on the CNR's Ariadne cluster.

More detailed file-by-file descriptions are available in the workflow
documentation:

- [`single_chain/README.md`](./single_chain/README.md)
- [`bulk/README.md`](./bulk/README.md)

The bulk workflow is intentionally kept close to its original cluster layout.
The scripts in `bulk/` use that directory as their working directory and
refer to local files such as `mongo.env`, the `.sbatch` templates, `results/`,
and `logs/`. They should therefore be launched from `simulations/bulk/`
without moving the root-level scripts or job templates.

The workflow follows the approach described in [BIO-SUSHY](https://github.com/daimoners/BIO-SUSHY-tutorials) (CNR Daimon team) for polymer construction,
force-field generation, and simulation setup.

These scripts are used to obtain chain-level properties (e.g. radius of gyration, end-to-end distance), which are then collected and added to the dataset used to train the predictive model for polymer properties.

## Environment

The scripts are intended to be run on Linux or WSL2.
This avoids file-locking issues that can occur on native Windows systems
when handling temporary files during the simulation workflow.
