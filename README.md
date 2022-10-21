# Cancer rREs paper analyses

This repository contains reproducible analyses for the paper titled "a genome-wide atlas of recurrent repeat expansions in human cancer." [Read the preprint on bioArxiv](https://www.biorxiv.org/content/10.1101/2022.08.24.505159v1).

## Usage
The `ehdn-lrdn-tool/` directory links to the ExpansionHunter denovo local read depth normalization tool. Instructions for using the tool are within the submodule.

All the analyses can be found in the `analysis/` directory. To run the analyses, first install the required packages using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
$ conda env create -f environment.yml
```

Then, activate the conda environment and add it as a Jupyter kernel:
```
$ conda activate cancer-rre
$ python -m ipykernel install --user --name=cancer-rre
```

Finally, start a Jupyter notebooks instance to run the analyses interactively:
```
$ jupyter notebooks
```

Make sure to select `cancer-rre` as the kernel.
