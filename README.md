# Description
Understanding the spatial distribution of gene expression in the pancreas is essential for advancing knowledge of its function and islet biology. 
The 10x Genomics Visium platform offers a robust method for quantifying gene expression within a spatial context. 
Gene expression data from Visium consists of gene-specific unique molecular identifier (UMI) counts, measured across thousands of 55m-diameter “spots” with spot-specific barcodes on a 6.5 x 6.5 mm Visium v1 capture slide. 
However, transcript bleed from nearby spots can contaminate these data.
This pipeline uses 10x Genomics Space Ranger v2.1.1 for preprocessing of Visium spatial gene expression data, and SpotClean v1.4.1 to correct for transcript bleed.

## Quickstart
Quickstart for deploying this pipeline locally and on a high performance compute cluster.

### 1. Environment
Create a conda environment from the `env/environment.yml` file using:
```
conda env create -f env/environment.yml
```
then activate the environment with
```
conda activate environment
```

### 2. Input files/config
Prepare input files with one sample per line, named \<MODEL\>_samples.txt. The model name should match in `config.yaml`.
### 3. Run
Use `bash dryrun.sh` to conduct a dry-run (lists rules Snakemake will run) and `bash run.sh` to run the full pipeline.
