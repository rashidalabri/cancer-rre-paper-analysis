# LRDN Benchmark

This directory contains a Snakemake workflow for benchmarking the local read depth normalization method. First, run the workflow using the `snakemake` command. Then, use the `Benchmark Plots.ipynb` notebook to reproduce the plots found in Extended Data Fig. 4 .

**Note**
Before running the workflow, update the `SAMTOOLS_PATH` variable in `scripts/extract.py`. If the executable is found in your default environment, then you can change it to `samtools` instead of the path.
