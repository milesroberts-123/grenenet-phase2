# GrENE-net phase 2 — Agent notes

Research repo for population-genetics simulations and pool-seq analyses. Primary languages: Snakemake, R, Python, SLiM.

## Entry points

- **Snakemake pipeline**: `workflow/Snakefile` is the main entry point.
  - `all_gt_sims` — SLiM genetic-drift simulations (`slim_results/{ID}_neutfreqs.tsv`).
  - `all_fst_sims` — msprime/GCTA Fst workflow (`msprime_results/{ID}.txt`, `gcta_results/{ID}.fst`).
  - `all_poolseq` — fastp + freqk empirical allele-frequency calls (`freqk_call_results/{sample}.txt`).
- **R analysis**: `resources/03_linked_selection.Rmd` is the current primary notebook.
- **Parameter generation**: `resources/06_sim_parameters.R` regenerates `config/gt_params.tsv` and `config/fst_params.tsv`.

## Running the pipeline

Run from inside `workflow/`, e.g.:

```bash
cd workflow
snakemake --profile profiles/default all_gt_sims --dry-run
snakemake --profile profiles/local all_gt_sims
```

Profiles:

- `profiles/default` — Slurm on Savio (account `co_moilab`, partition `savio3_htc`, QOS `savio_lowprio`).
- `profiles/local` — local execution, 20 jobs.
- `profiles/msu` — MSU Slurm profile (account `josephsnodes`).

The Snakefile pins `min_version("9.5.1")` and loads `config/config.yaml` from `../config/config.yaml`.

## R workflow conventions

- R source modules live in `resources/R/` (`gt.R`, `misc.R`, `plot.R`).
- Modules are loaded with `box::use()`, e.g. `box::use(./R/gt[...])`.
- Tests are in `resources/R/tests.R` and run via `testthat::test_file("./R/tests.R")` from within the Rmd.
- The main notebook writes outputs to `results/by-date/<YYYY-MM-DD>/` using `Sys.Date()` in the path.

## Important gotchas

- **Absolute paths in `config/config.yaml`**: they currently point to Savio-specific paths (`/global/scratch/users/...` and `/global/home/users/...`) and must be updated for a new machine.
- **`config/reads.tsv` also contains absolute FASTQ paths** and is filtered in the Snakefile to samples starting with `S`:

  ```python
  reads = reads[reads['sample'].str.contains('^S')]
  ```

  Only those samples drive `all_poolseq`.
- **`reads.tsv` vs `samples_phase1_poolseq.tsv`**: the latter is a simple sample list; `reads.tsv` actually maps samples to R1/R2 paths.
- **`.gitignore` is aggressive**: it ignores `*.tsv`, `*.csv`, `*.txt`, `*.html`, `*.pdf`, `*.fa`, `*.fasta`, `*.vcf*`, all of `data/`, all of `results/`, and `.snakemake/`. Existing tracked files stay tracked, but new outputs and intermediate files will not show up as untracked.
- **No requirements/lock files**: dependencies come from the Snakemake conda envs in `workflow/envs/` and whatever R packages the notebooks `library()` load.
- **`config/sampling_times.tsv`**: loaded in `workflow/Snakefile` and used by included rules.
- **Simulation parameter tables are generated**: edit `resources/06_sim_parameters.R` if you need to change the design, then rerun it to rebuild `config/gt_params.tsv` and `config/fst_params.tsv`.
- **Do not set memory resources in `profiles/default`**: Savio (UC Berkeley) restricts memory requests on `savio3_htc`/`savio_lowprio`. Memory resources belong in `profiles/msu` (e.g. `mem_mb_per_cpu`) or `profiles/local` only.
- **Legacy rules** are commented out in the Snakefile under `rules/legacy/`; they are not part of the current pipeline.
- **Three simulation types** in `slim.smk`: `unstruct`, `struct`, and `bank`. Each uses a different `.slim` script (`gt_expectations.slim`, `gt_expectations_structured.slim`, `gt_expectations_seed_bank.slim`).
