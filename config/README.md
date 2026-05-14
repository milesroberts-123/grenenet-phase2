# config/ — Configuration Reference

All configuration files for the GrENE-net Phase 2 pipeline live here. They are
grouped into four categories: pipeline settings, simulation parameters, sample
metadata, and time-point data.

---

## 1. Pipeline settings

### `config.yaml`

Snakemake's central config file, loaded at the top of `workflow/Snakefile`. It
defines paths and parameter defaults that every rule can access via `config["key"]`.

| Key | Type | Description |
|-----|------|-------------|
| `sim_fasta` | path | Pseudogenome reference FASTA used for SIMULATION |
| `sim_vcf` | path | Reference VCF used for SIMULATION |
| `reference` | path | Empirical read-mapping reference genome |
| `freqk_index` | path | Path to the freqk reference index |
| `phase1_poolseq_root` | path | Base directory for Phase 1 pool-seq FASTQ files |
| `phase1_poolseq_samples` | path | TSV listing Phase 1 pool-seq sample IDs |
| `phase1_poolseq_suffix` | list | Suffix appended to Phase 1 sample IDs (R1 / R2) |
| `adaptest_root` | path | Base directory for AdapTest deep-seq data |
| `adapteset_samples` | path | TSV listing AdapTest sample IDs |
| `seed_root` | path | Base directory for seed-mix FASTQ files |
| `seed_samples` | path | TSV listing seed-mix sample IDs |
| `seed_suffix` | list | Suffix appended to seed-mix sample IDs (.1_P.fq.gz / .2_P.fq.gz) |
| `pool_root` | path | Base directory for pooled samples |
| `pool_samples` | path | TSV listing pooled sample IDs |
| `phase2_poolseq_root` | path | Base directory for Phase 2 pool-seq FASTQ files |
| `phase2_poolseq_samples` | path | TSV listing Phase 2 pool-seq sample IDs |
| `phase2_poolseq_suffix` | list | Suffix appended to Phase 2 sample IDs (_R1_001.fastq.gz / _R2_001.fastq.gz) |
| `k` | int | K-mer size for KMC (default 31) |
| `contam_mincount` | int | Minimum k-mer count to flag as contaminant |
| `sample_mincount` | int | Minimum k-mer count to keep as sample |
| `unqualLimit` | int | fastp: maximum number of unqualified bases allowed |
| `qualThresh` | int | fastp: quality threshold for filtering |
| `windowLength` | int | fastp: sliding window length for QC |
| `nBaseLimit` | int | fastp: maximum N bases per read |
| `custom_contams` | list | Custom contaminant sources (currently `grenenet_pseudogenomes`) |

---

## 2. Simulation parameter tables

Both files share the **identical schema** (TSV, tab-delimited). They differ only
in the `adjust` column and the breadth of `N` values.

### Shared columns

| Column | Values | Description |
|--------|--------|-------------|
| `nmu` | 1e-08 | Per-base mutation rate |
| `tmu` | 1e-09 | Per-base recombination rate |
| `R` | 1e-08 | Migration rate between demes |
| `N` | 1000, 1025, 1333, 1904, 1980 | Effective population size (five levels: full and subsamples) |
| `L` | 5e+06 | Genome length in bp |
| `sigma` | 0, 0.05, 0.5, 0.95, 0.99 | Divergence/fission parameter across demes |
| `alpha` | 0, 0.01 | Selection coefficient on beneficial alleles |
| `gamma` | 5 | Selection coefficient on deleterious alleles / dominance |
| `tau` | 105 | Number of generations |
| `rep` | 1–32 | Replicate number |
| `adjust` | TRUE / FALSE | Whether to adjust effective diploid pool size for pool-seq simulation |
| `ID` | 1–1418 | Unique simulation run identifier |

- **`gt_params.tsv`** — Simulation parameters for **genetic drift (GT) curves**.  `adjust = TRUE` for IDs 1–600 (N varying across five levels), `adjust = FALSE` for IDs 601–1418 (N = 1000 only). 50 sigma × 2 alpha × 10 tau replications = 1418 rows total.
- **`parameters.tsv`** — Identical schema and data to `gt_params.tsv`.  Appears to be a duplicate or predecessor; used in earlier iterations.

### `fst_params.tsv`

Simulation parameters for **Fst scans** (msprime forward sim). Contains 50
parameter combinations (5 sigma × 2 alpha × 5 tau reps), all with `N = 1000`
and `adjust = TRUE`.

### `samples_adaptest.tsv`

List of AdapTest experimental sample IDs (MEAJM prefix, 769 rows, column
`sample`). Used to map AdapTest deep-seq libraries to the pipeline.

### `samples_phase1_poolseq.tsv`

List of Phase 1 pool-seq sample IDs (MLFH prefix, 2415 rows, column
`sample`). Cross-referenced by `config.yaml` (`phase1_poolseq_samples`) and the
`reads.tsv` file.

### `samples_phase2_poolseq.tsv`

List of Phase 2 pool-seq sample IDs (MEAJM prefix, 471 rows + "Undetermined",
column `sample`). Mapped from root `adaptest_root` with `phase2_poolseq_suffix`.

### `samples_pool.tsv`

List of pooled sample IDs (MEAJM prefix, 483 rows + trailing blank, column
`sample`). Corresponds to the `pool_root` path in `config.yaml`.

### `samples_seed_mix.tsv`

List of seed-mix individual samples (9 rows: S1–S8, column `sample`). Mapped
from root `seed_root` with `seed_suffix`.

---

## 3. Time points

### `sampling_times.csv`

Empirical sampling times for historical herbarium material and modern accessions.
460 rows.

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `year` | Collection year (range ~1817–2007) |
| `type` | `historical` (herbarium) or `modern` (recently collected) |
| `time` | Relative time metric for simulation (in generations?) |

### `sampling_times.tsv`

Simulation-ready sampling time points. 140 rows (124 `historical` + 16
`modern`).

| Column | Description |
|--------|-------------|
| `type` | `historical` or `modern` |
| `time` | Relative generation time |
| `n` | Number of samples at this time point |
| `ID` | Simulation run ID |

---

## 4. Read / file paths

### `reads.tsv`

Mapping of sample IDs to read1 and read2 FASTQ paths.  224 rows.

| Column | Description |
|--------|-------------|
| `sample` | MLFH sample ID |
| `read1` | Absolute path to R1 FASTQ |
| `read2` | Absolute path to R2 FASTQ |

---

## 5. Notes

- All `samples_*.tsv` files contain a single `sample` column with one ID per
  line. They are cross-referenced by `config.yaml` via `~//config/` paths.
- Absolute paths in `config.yaml` and `reads.tsv` are specific to the compute
  environment. Before copying this repo to another machine, update those paths.
- `fst_params.tsv` and `parameters.tsv` currently contain identical data.  If
  one is meant to diverge, clarify its intended role.
- `sampling_times.csv` and `sampling_times.tsv` serve different purposes: the
  CSV holds empirical collection metadata; the TSV holds simulation-ready time
  points. They should not be conflated.
