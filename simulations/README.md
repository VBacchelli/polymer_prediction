This folder contains scripts to build and run simple molecular dynamics simulations of polymer systems.

The workflow follows the approach described in [BIO-SUSHY](https://github.com/daimoners/BIO-SUSHY-tutorials) (CNR Daimon team) for polymer construction,
force-field generation, and simulation setup.

These scripts are used to obtain chain-level conformational properties (e.g. radius of gyration, end-to-end distance), which are then collected and added to the dataset used to train the predictive model for polymer properties.

## Environment

The scripts are intended to be run on Linux or WSL2.
This avoids file-locking issues that can occur on native Windows systems
when handling temporary files during the simulation workflow.

All required dependencies are specified in the provided [`environment.yml`](./simulations/environment.yml)
and should be installed using Conda with:

```bash
conda env create -f environment.yml
conda activate simulation
```
