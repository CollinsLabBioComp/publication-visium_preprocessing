"""
Snakemake Visium preprocessing pipeline
----------------------

Snakemake for preprocessing and bleed correction of Visium v1 Spatial Gene Expression data
"""

import yaml
import pandas as pd
from snakemake.io import expand
import numpy as np

# Load config file
config_file = "config.yaml"
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

# Load from Visium input_file file
input_file = pd.read_csv(config['paths']['input_file'])


# Load variables from config files
FASTQ_DIR = list(config['fastq_dirs'])
EXPERIMENT_SAMPLE = list(input_file['experiment_sample'])
SPOTCLEAN_RADII = list(config['spotclean_radii'])
SPOTCLEAN_KERNELS = list(config['spotclean_kernels'])
GENES = config['plot_gene']
SLIDE = dict(zip(EXPERIMENT_SAMPLE, input_file['slide']))
AREA = dict(zip(EXPERIMENT_SAMPLE, input_file['area']))

rule all:
    input: 
      expand('results/plots/{experiment_sample}/kernels_{gene}_expression.png', experiment_sample = EXPERIMENT_SAMPLE, gene = GENES),
      expand('results/plots/{experiment_sample}/normalized_{gene}_expression.png', experiment_sample = EXPERIMENT_SAMPLE, gene = GENES)


################################################################################
# Step 1: Run Space Ranger
################################################################################
rule make_fastqs:
    input:
      run = '{fastq_dir}/Permeabilization_Transcripts',
      csv = '{fastq_dir}/samples_{fastq_dir}.csv'
    output: directory('results/fastqs/{fastq_dir}_fastq')
    envmodules:
        "spaceranger/2.1.1"
    shell:
      '''
      mkdir -p results/fastqs

      spaceranger mkfastq \
        --id={wildcards.fastq_dir}_fastq \
        --run={input.run} \
        --csv={input.csv} \
        --delete-undetermined \
        --output-dir=results/fastqs/{wildcards.fastq_dir}_fastq
      cd ..

      '''

rule spaceranger_count:
    input:
      transcriptome = config['paths']['transcriptome']
    output: directory('results/counts/{experiment_sample}/outs')
    params:
      slide=lambda wildcards: SLIDE[wildcards.experiment_sample],
      area=lambda wildcards: AREA[wildcards.experiment_sample]
    envmodules:
        "spaceranger"
    shell:
      '''
      cwd=$(pwd)
      mkdir -p results/counts
      cd results/counts
      
      spaceranger count --id={wildcards.experiment_sample} \
        --transcriptome={input.transcriptome} \
        --fastqs=../fastqs/HTAN_fastq \
        --sample={wildcards.experiment_sample} \
        --image=${{cwd}}/data/images/{wildcards.experiment_sample}.tif \
        --reorient-images true \
        --create-bam true \
        --slide {params.slide} \
        --area {params.area}
      
      cd ..
      '''


################################################################################
# Step 2: Run bleed correction using SpotClean
################################################################################
rule correct_bleed:
    input: 'results/counts/{experiment_sample}/outs'
    output: 'results/cleaned/{experiment_sample}/outs/kernel={kernel}/radius={radius}/cleaned_feature_bc_matrix.h5ad'
    params:
      kernel=lambda wildcards: SPOTCLEAN_KERNELS,
      radius=lambda wildcards: SPOTCLEAN_RADII
    shell:
      '''
      mkdir -p results/cleaned

      Rscript scripts/010-correct_bleed.R \
        --input {input} \
        --output results/cleaned/{wildcards.experiment_sample}/outs \
        --kernel {wildcards.kernel} \
        --radius {wildcards.radius}
      '''


################################################################################
# Step 3: Convert SpotClean anndata format to 10x H5 format
################################################################################
rule convert_anndata:
    input: 'results/cleaned/{experiment_sample}/outs/kernel={kernel}/radius={radius}/cleaned_feature_bc_matrix.h5ad'
    output: 'results/cleaned/{experiment_sample}/outs/kernel={kernel}/radius={radius}/cleaned_feature_bc_matrix.h5'
    shell:
      '''
      # Copy spatial and analysis directories from SpaceRanger outs
      cp -rn results/counts/{wildcards.experiment_sample}/outs/spatial results/cleaned/{wildcards.experiment_sample}/outs/kernel={wildcards.kernel}/radius={wildcards.radius}
      cp -rn results/counts/{wildcards.experiment_sample}/outs/analysis results/cleaned/{wildcards.experiment_sample}/outs/kernel={wildcards.kernel}/radius={wildcards.radius}

      python scripts/022-convert_anndata.py \
        --dataDir results/cleaned/{wildcards.experiment_sample}/outs/kernel={wildcards.kernel}/radius={wildcards.radius}/
      '''

rule expression_csv:
    input: expand('results/cleaned/{experiment_sample}/outs/kernel={kernel}/radius={radius}/cleaned_feature_bc_matrix.h5', experiment_sample = EXPERIMENT_SAMPLE, radius = SPOTCLEAN_RADII, gene = GENES, kernel = SPOTCLEAN_KERNELS)
    output: 'results/plots/{experiment_sample}/kernel={kernel}/spot_{gene}_expression.csv'
    params:
      radii = SPOTCLEAN_RADII
    shell:
      '''
      python scripts/030-expression_csv.py \
        --spacerangerDir results/counts/{wildcards.experiment_sample}/outs/ \
        --spotcleanDir results/cleaned/{wildcards.experiment_sample}/outs/kernel={wildcards.kernel}/ \
        --outDir results/plots/{wildcards.experiment_sample}/kernel={wildcards.kernel}/ \
        -r {params.radii} \
        -g {wildcards.gene}
      '''


################################################################################
# Step 4: Plot results
################################################################################
rule plot_expression:
    input: expand('results/plots/{experiment_sample}/kernel={kernel}/spot_{gene}_expression.csv', experiment_sample = EXPERIMENT_SAMPLE, kernel = SPOTCLEAN_KERNELS, gene = GENES)
    output: 'results/plots/{experiment_sample}/normalized_{gene}_expression.png'
    params:
      gene = GENES
    shell:
      '''
      Rscript scripts/040-plot_expression.R \
        -p 'results/plots/{wildcards.experiment_sample}/' \
        -g {wildcards.gene}
      '''


rule plot_expression_grid:
    input: expand('results/plots/{experiment_sample}/normalized_{gene}_expression.png', experiment_sample = EXPERIMENT_SAMPLE, gene = GENES)
    output: 'results/plots/{experiment_sample}/kernel={kernel}/combined_{gene}_expression.png'
    params:
      radii = SPOTCLEAN_RADII,
      gene = GENES
    shell:
      '''
      Rscript scripts/041-plot_expression_grid.R \
        -p 'results/plots/{wildcards.experiment_sample}/kernel={wildcards.kernel}/' \
        -r "{params.radii}" \
        -k {wildcards.kernel} \
        -g {wildcards.gene}
      '''
      
rule compare_kernels:
    input: expand('results/plots/{experiment_sample}/kernel={kernel}/combined_{gene}_expression.png', experiment_sample = EXPERIMENT_SAMPLE, kernel = SPOTCLEAN_KERNELS, gene = GENES)
    output: 'results/plots/{experiment_sample}/kernels_{gene}_expression.png'
    params:
      radii = SPOTCLEAN_RADII,
      gene = GENES,
      kernels = SPOTCLEAN_KERNELS
    shell:
      '''
      Rscript scripts/042-compare_kernels.R \
        -p 'results/plots/{wildcards.experiment_sample}/' \
        -g {wildcards.gene} \
        -k '{params.kernels}'
      '''