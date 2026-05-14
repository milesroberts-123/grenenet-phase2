# grenenet-phase2

This is the repository for the GrENE-net phase 2 project, including snakemake workflow and R notebooks to recreate analyses and simulations.

## Repository structure

```
config/             Configuration files (YAML/TSV) for simulation parameters, sample metadata, and experimental design
data/               Empirical data files (allele frequencies, genome pairings, climate site data, metadata CSVs)
resources/          Pre-processing scripts, R source code, simulation scripts, and helper functions
results/            Analysis outputs organized by analysis type and date (by-analysis/, by-date/)
workflow/           Snakemake pipeline including Snakefile, rules, environments, profiles, and notebooks
  envs/             Conda environment specifications (bcftools, plink2, msprime, fastp, bwa, kmc, etc.)
  notebooks/        R analysis notebooks (e.g., 00_cam5_deletion_allele.Rmd)
  profiles/         Cluster/profile configurations (default, msu, local)
  rules/             Snakemake rule modules (freqk, fst_sim, msprime, fastp, legacy tools)
  scripts/           Standalone helper scripts (Python, SLiM simulations, pyslim utilities)
```


