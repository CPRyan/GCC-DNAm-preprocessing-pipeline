---
output:
  html_document: default
  pdf_document: default
header-includes:
  - |
    ```{=latex}
    \usepackage{fvextra}
    \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
    \renewcommand{\_}{\textunderscore\allowbreak{}}
    ```
pagetitle: "GCC DNAm Preprocessing Pipeline"
---
# GCC DNAm Preprocessing Pipeline

This repository, written in R, provides a clean, chunk-aware pipeline for Illumina DNA methylation array preprocessing from raw IDAT files through sample QC, Noob normalization, probe filtering, and BMIQ normalization, with downstream epigenetic clocks, cell-type estimates, summary figures and a PDF report.

The pipeline is designed for large datasets that may not fit in memory as a single object. Chunks are supported throughout the workflow: samples are split into chunks, each step processes one chunk at a time, and the chunk-level outputs are merged at the end. All steps can be run with one command or one at a time.

The pipeline assumes users begin with a sample sheet and a directory of paired IDAT files. AWS/S3 downloading and study-specific data access are intentionally excluded from the public pipeline.

**Tested on:** EPICv2 arrays (dried blood spots, 48 samples). The pipeline detects the array type (450k, EPICv1 or EPICv2) from the data and has code for all three, but only EPICv2 has been tested with this version. On 450k and EPICv1 data, check the results with extra care: the SNP outlier and sex checks, the probe filtering and the downstream modules have not been tested there, and the 450k probe lists are not included (see [Probe filtering](#probe-filtering)).

**Using it takes four steps:** [1. Install the R packages](#step-1-install-the-r-packages), [2. Prepare your data](#step-2-prepare-your-data), [3. Clone the pipeline and fill in the config](#step-3-clone-the-pipeline-and-fill-in-the-config), [4. Run the pipeline](#step-4-run-the-pipeline). Then open the report (see [What you get](#what-you-get)).

## Step 1: Install the R packages

R 4.0 or newer is needed. The example config switches on every result (cell counts, PC clocks, SystemsAge, DunedinPACE and FraminghamPACE), so a default run needs all of the packages below. Install them in R:

```r
# from CRAN
install.packages(c("yaml", "ggplot2", "RPMM", "remotes", "BiocManager"))

# from Bioconductor
BiocManager::install(c("minfi", "Biobase", "wateRmelon", "FlowSorted.Blood.EPIC"))

# from GitHub
remotes::install_github("hhhh5/ewastools")
remotes::install_github("MorganLevineLab/methylCIPHER")
remotes::install_github("danbelsky/DunedinPACE")
remotes::install_github("will-marella/FraminghamPACE")
```

**Array annotation.** Install the packages for your array:

```r
# EPICv2
BiocManager::install(c("IlluminaHumanMethylationEPICv2manifest",
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))

# EPICv1
BiocManager::install(c("IlluminaHumanMethylationEPICmanifest",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))

# 450k
BiocManager::install(c("IlluminaHumanMethylation450kmanifest",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
```

What the packages are used for:

| Package | Used for |
|---|---|
| `minfi`, `Biobase`, `wateRmelon`, `RPMM` | loading the IDATs, QC, Noob and BMIQ normalization |
| `ewastools` | the EWAS tools control metrics, the SNP outlier check and the sex check |
| `FlowSorted.Blood.EPIC` | cell counts |
| `methylCIPHER` | PC clocks and SystemsAge |
| `DunedinPACE` | DunedinPACE |
| `FraminghamPACE` | FraminghamPACE |
| `yaml`, `ggplot2` | reading the config; the figures and the report |

Notes:

- **Skipping a result.** If you switch a result off in the config (`noob_downstream` and `bmiq_downstream`), you do not need its package. If a package is missing for a result that is switched on, that module fails with a warning and the other modules still run.
- **Cell-count reference.** The default is `FlowSorted.Blood.EPIC`. The extended reference (12 cell types) is restricted: it needs a licence agreement signed by an authorized representative of your institution (see the LICENSE on the package's GitHub page, github.com/immunomethylomics/FlowSorted.BloodExtended.EPIC). To use it, install it with `remotes::install_github("immunomethylomics/FlowSorted.BloodExtended.EPIC")` and set `noob_downstream: cell_reference: FlowSorted.BloodExtended.EPIC` in the config.

## Step 2: Prepare your data

Prepare the IDAT files and the sample sheet wherever you like, on your computer or on a server. In Step 3 you point the config to them.

**IDAT files.** One folder with a pair of files for every sample, `..._Red.idat` and `..._Grn.idat`. Subfolders are searched too.

```text
203546000006_R01C01_Red.idat
203546000006_R01C01_Grn.idat
```

**Sample sheet.** A CSV file with one row per sample and these columns:

| Column | Used for | What to put in it |
|---|---|---|
| `Sample_Name` | your name or ID for the sample (for example the participant ID), shown next to the `IDATid` in the logs, the report and the output tables | required; put it in the first column |
| `Sentrix_ID` | finding the IDAT files | the chip ID, 12 digits, e.g. `203546000006` |
| `Sentrix_Position` | finding the IDAT files | the position on the chip, e.g. `R01C01` |
| `IDATid` | the sample ID in all outputs | unique for every sample; usually `Sentrix_ID` + `_` + `Sentrix_Position` |
| `Age` | clocks and plots against age | age in years (a number) |
| `Female` | sex check, clocks, plots by sex | `1` = female, `0` = male |
| `SubjectID` | keeping samples of one person together | only needed when the samples are split into chunks (`n_chunks` above 1); samples of the same person share the same value |

`Age` and `Female` must have exactly these names, because the clocks look for them. The other column names can be changed in the config (section `sample_sheet`).

An example with three samples:

```text
Sample_Name,Sentrix_ID,Sentrix_Position,SubjectID,IDATid,Female,Age
sample_001,203546000006,R01C01,subj_001,203546000006_R01C01,1,34
sample_002,203546000006,R02C01,subj_002,203546000006_R02C01,0,41
sample_003,203546000006,R03C01,subj_003,203546000006_R03C01,1,29
```

Notes:

- Spreadsheet programs can turn the 12-digit `Sentrix_ID` into `2.03546E+11` when you save. Open the CSV in a text editor to check that the digits are intact.
- An Illumina sample sheet with a `[Data]` section also works.

## Step 3: Clone the pipeline and fill in the config

### Clone the repository

Clone the repository from a terminal and go into its folder:

```bash
git clone https://github.com/CPRyan/GCC-DNAm-preprocessing-pipeline.git
cd GCC-DNAm-preprocessing-pipeline
```

All the commands in this step and in Step 4 are run from this repository folder. In a new terminal, go back to it first (`cd /path/to/GCC-DNAm-preprocessing-pipeline`); `pwd` should print a path that ends in `GCC-DNAm-preprocessing-pipeline`.

### Fill in the config

Copy the example config and edit the copy:

```bash
cp config/example_config.yml config/my_study.yml
```

At the top of the file, set the run name and the three paths (full paths, no `~`):

```yaml
study_name: My_study

paths:
  idat_dir: /home/me/data/idats
  sample_sheet: /home/me/data/sample_sheet.csv
  output_dir: /home/me/output
```

The `study_name` is used to name all outputs. You set only the output folder; the subfolders `data`, `figures` and `logs` are created inside it.

That is all you need to start. Everything below the line `ADVANCED SETTINGS` in the file has working defaults and can be left as it is. The advanced settings are grouped as follows:

| Section | What it controls |
|---|---|
| `sample_sheet` | the names of the sample sheet columns, and whether IDATs are searched in subfolders |
| `chunking` | number of chunks, grouping columns, random seed |
| `noob_downstream` | which results are calculated on the Noob data (cell counts and their reference, PC clocks, SystemsAge, DunedinPACE) |
| `qc_thresholds` | sample QC thresholds |
| `ewas_tools_qc` | the EWAS tools checks (must stay on) |
| `sample_filtering` | which failed checks remove a sample |
| `probe_filtering` | probe filters, their thresholds and the probe list file names |
| `bmiq` | number of cores and log detail for BMIQ |
| `bmiq_downstream` | FraminghamPACE, the report, and which figures are made |

## Step 4: Run the pipeline

### Run everything with one command

From the repository folder (see Step 3):

```bash
Rscript scripts/run_all_steps.R --config=config/my_study.yml
```

This runs steps 00 to 09 in order, for all chunks. Each step runs in its own R process, exactly as if you had started it yourself, so it writes its own log (see [Logs](#logs)) and frees its memory before the next step starts. The run stops at the first step that fails. At the end, a summary table shows the status, run time, number of warnings and number of failed modules of every step.

The first run downloads two reference files into the `data/` folder of the repository (about 6 GB in total), so it needs an internet connection and disk space. BMIQ is the slowest step.

When it finishes, look at:

- `output/data/Derived_Variables/My_study_report.pdf`: the report, with the sample and probe numbers and which samples were removed and why, ready for a methods section;
- `output/data/My_study_methods_summary.txt`: the same numbers as text;
- `output/data/Derived_Variables/`: the files meant for further analysis (see [What you get](#what-you-get));
- `output/logs/`: one log per step, with everything the step printed.

### Options

```bash
# check your inputs only (fast): validates the sample sheet and the IDAT files
Rscript scripts/run_all_steps.R --config=config/my_study.yml --to=00

# run only some steps, for example stop before BMIQ (the slowest step)
Rscript scripts/run_all_steps.R --config=config/my_study.yml --to=07

# restart from a step, for example after fixing an error
# (earlier steps are not rerun: their saved outputs are used)
Rscript scripts/run_all_steps.R --config=config/my_study.yml --from=05

# run a range of steps
Rscript scripts/run_all_steps.R --config=config/my_study.yml --from=05 --to=07
```

BMIQ took about half a minute per sample on one core in our test. Set `bmiq: n_cores` in the config to use more cores.

### Long runs on a server

To keep a run going after you close the terminal, start it in the background and follow the progress in the log:

```bash
nohup Rscript scripts/run_all_steps.R --config=config/my_study.yml > /dev/null 2>&1 &
tail -F output/logs/run_all_steps.log
```

Press Ctrl+C to stop following the log; the run itself continues. Each step also writes its own log in `output/logs/`.

### Run step by step

```bash
Rscript scripts/00_validate_inputs.R --config=config/my_study.yml
Rscript scripts/01_create_targets_and_chunks.R --config=config/my_study.yml
Rscript scripts/02_load_rgset_chunks.R --config=config/my_study.yml --chunk=all
Rscript scripts/03_run_raw_qc.R --config=config/my_study.yml --chunk=all
Rscript scripts/04_filter_samples.R --config=config/my_study.yml --chunk=all
Rscript scripts/05_noob_normalize.R --config=config/my_study.yml --chunk=all
Rscript scripts/06_noob_downstream.R --config=config/my_study.yml --chunk=all
Rscript scripts/07_probe_filter.R --config=config/my_study.yml --chunk=all
Rscript scripts/08_bmiq.R --config=config/my_study.yml --chunk=all
Rscript scripts/09_bmiq_downstream.R --config=config/my_study.yml
```

### Large studies: one chunk at a time

If the data do not fit in memory, set `chunking.n_chunks` in the config to split the samples (`1` means no chunking). Every step then works on one chunk at a time. `--chunk=all` runs all chunks in sequence and is the default; `--chunk=N` runs a single chunk.

Steps 02 to 08 can be run for one chunk at a time, and each chunk is independent of the others (see [Probe filtering](#probe-filtering) for what this means for the probes). Step 09 works on one chunk or on all of them: with `--chunk=N`, the merged table, the figures and the report cover that chunk only, and their names end in `_chunkN`. Run it without `--chunk`, after all chunks have finished step 08, to merge them.

```bash
# for each chunk N = 1, 2, ...: steps 02 to 08
Rscript scripts/02_load_rgset_chunks.R --config=config/my_study.yml --chunk=N
Rscript scripts/03_run_raw_qc.R --config=config/my_study.yml --chunk=N
Rscript scripts/04_filter_samples.R --config=config/my_study.yml --chunk=N
Rscript scripts/05_noob_normalize.R --config=config/my_study.yml --chunk=N
Rscript scripts/06_noob_downstream.R --config=config/my_study.yml --chunk=N
Rscript scripts/07_probe_filter.R --config=config/my_study.yml --chunk=N
Rscript scripts/08_bmiq.R --config=config/my_study.yml --chunk=N

# once: M values, FraminghamPACE, merged results, figures and report for all chunks
Rscript scripts/09_bmiq_downstream.R --config=config/my_study.yml
```

If a run is interrupted, rerun the failed step for the relevant chunk. Chunk-level output files are written separately, so completed chunks do not need to be regenerated unless upstream settings changed.

## What you get

The pipeline writes outputs to `paths.output_dir`: data files under `data`, figures under `figures`, and logs under `logs`.

### Derived variables

The outputs meant for further analysis and sharing are collected in `data/Derived_Variables`, which is created at the end of the analysis (step 09):

| File | Content |
|---|---|
| `{study_name}_BMIQ_betas.rds` | the BMIQ-normalized beta values: a matrix with the probes in the rows and the samples in the columns (named by the `IDATid`). All chunks together, on the probes they share. Read it with `readRDS()` |
| `{study_name}_BMIQ_Mvalues.rds` | the M values calculated from these betas, in the same layout |
| `{study_name}_report.pdf` | the report (see below) |
| `{study_name}_cell_composition_and_clocks.csv` | one row per sample: cell composition, PC clocks, SystemsAge, DunedinPACE and FraminghamPACE |
| `{study_name}_sample_qc_filtering_decisions.csv` | for every sample: the checks it failed, which of them were used to remove it, and whether it was removed; all chunks together |
| `{study_name}_sample_qc_pass_fail.csv` | for every sample: the result of every QC check; all chunks together |
| `{study_name}_control_probe_pca_scores.csv` | the control-probe PCA scores of every sample (technical variation); calculated before samples are removed, so removed samples are included |

### Report

At the end of step 09, the pipeline writes a report and a text file with the numbers behind it (switch both off with `bmiq_downstream: report: false`):

- **`data/Derived_Variables/{study_name}_report.pdf`**: a report of about ten pages, in A4 format:
  1. Overview: the samples kept and removed, the probes kept and removed, the checks and the removed samples.
  2. Samples: the identity and control checks, all sample QC checks with their rules, the removed samples, and the samples that failed a check that is not used for filtering.
  3. Sample QC plots: one point per sample for signal intensity, probes not detected, bead count, SNP outlier score, the sex check (X against Y signal) and the control-probe PCA, coloured by passed, flagged or removed.
  4. Probes: every filter with its number, percentage and reason, and a chart (for each chunk, if the chunks end up with different probes).
  5. Normalization: the distribution of the BMIQ beta values and of the M values.
  6. Cell composition: the estimated cell-type proportions.
  7. Epigenetic clocks: PC clocks against age.
  8. SystemsAge: the scores of the systems against age.
  9. Pace of aging: the DunedinPACE and FraminghamPACE scores.
  10. Settings and software: the settings of the run and the versions of R and the main packages.

  The plots in pages 5 to 9 are the same as the figures in the `figures` folder. If you switch a plot off in `bmiq_downstream: plots`, its page is left out of the report. A page that fails is replaced by a note, and the other pages are still made.
- **`data/{study_name}_methods_summary.txt`**: the same numbers as text for a methods section, plus the settings and the versions of R and the main packages.

### Figures

All figures are PDFs in `figures`, named `{study_name}_...`:

| File | Content |
|---|---|
| `{study_name}_cell_counts_boxplots_{reference}.pdf` | cell-type proportions |
| `{study_name}_clocks_vs_age.pdf` | PC clocks against chronological age |
| `{study_name}_systemsage_vs_age.pdf` | SystemsAge scores against chronological age |
| `{study_name}_pace_scores_distribution.pdf` | histograms of the DunedinPACE and FraminghamPACE scores |
| `{study_name}_pace_scores_by_sex.pdf` | pace scores by sex |
| `{study_name}_bmiq_beta_distribution.pdf` | density of the BMIQ betas, one curve per sample |
| `{study_name}_bmiq_mvalue_distribution.pdf` | density of the M values, one curve per sample |

### Data files

- **Validated manifest**: `{study_name}_sample_manifest_validated.csv`
- **Chunk assignments**: `{study_name}_targets_chunks.csv`
- **Raw RGsets**: `{study_name}_RGset_raw_chunk{chunk}.rds`
- **QC metrics**: `{study_name}_sample_qc_metrics_chunk{chunk}.csv`
- **QC pass/fail flags**: `{study_name}_sample_qc_pass_fail_chunk{chunk}.csv`
- **Filtering decisions**: `{study_name}_sample_qc_filtering_decisions_chunk{chunk}.csv`
- **Control-probe PCA (all chunks)**: `{study_name}_control_probe_pca_scores_all_chunks.csv`, `{study_name}_control_probe_pca_variance_all_chunks.csv` and `{study_name}_control_probe_pca_all_chunks.rds`
- **Good-sample RGsets**: `{study_name}_RGset_good_samples_chunk{chunk}.rds`
- **Noob-normalized GRsets**: `{study_name}_Noob_GRset_chunk{chunk}.rds`
- **Downstream results per chunk**: one table per module, for example `{study_name}_PCClocks_noob_chunk{chunk}.csv`
- **Probe-filtered GRsets**: `{study_name}_probe_filtered_GRset_chunk{chunk}.rds`
- **Probe filtering table**: `{study_name}_probe_filter_attrition_chunk{chunk}.csv`
- **BMIQ beta matrices**: `{study_name}_BMIQ_betas_chunk{chunk}.rds`
- **BMIQ M values**: `{study_name}_BMIQ_Mvalues_chunk{chunk}.rds`
- **Methods summary**: `{study_name}_methods_summary.txt`

QC outputs are intentionally separated into measured values and decisions:

- **Metrics table**: observed values such as mean bead count, low-bead fraction, average detection P-value, SNP outlier score, predicted sex, and EWAS tools values
- **Pass/fail table**: all QC flags, regardless of whether they are used for filtering
- **Filtering decisions table**: which failed checks were actually used to remove each sample

### Log files

`logs/{script name}.log`, one for every step (see [Logs](#logs)).

## What each step does

### `00_validate_inputs.R`

Validates the sample sheet and confirms that every sample has paired `_Red.idat` and `_Grn.idat` files. IDATs may be located directly in the IDAT directory or in recursive subdirectories.

### `01_create_targets_and_chunks.R`

Creates the target manifest and assigns samples to chunks. If grouping columns are provided, related samples are kept in the same chunk where possible. With `n_chunks: 1` there is a single chunk.

### `02_load_rgset_chunks.R`

Loads raw IDAT files into `minfi` `RGChannelSetExtended` objects, one RDS file per chunk.

### `03_run_raw_qc.R`

Runs raw sample QC and the mandatory EWAS tools checks: the Illumina control metrics, the SNP outlier check and the sex check (see [Sample QC checks](#sample-qc-checks)). It produces separate metrics, pass/fail and control-probe PCA outputs.

### `04_filter_samples.R`

Uses the config file to decide which failed checks remove samples. Each check can be switched on or off. The log lists the number of samples removed by each check, and every removed sample with its ID, name and reason.

### `05_noob_normalize.R`

Noob-normalizes good-sample RGsets and saves chunk-level normalized `GenomicRatioSet` objects.

### `06_noob_downstream.R`

Runs the downstream modules that use Noob-normalized data before probe filtering, as enabled in the config: cell counts (FlowSorted.Blood.EPIC or FlowSorted.BloodExtended.EPIC), PC clocks and SystemsAge (methylCIPHER), and DunedinPACE. Each module writes its own table per chunk, and a module that fails gives a warning without stopping the others.

The reference files for the PC clocks (about 2.1 GB) and SystemsAge (about 3.6 GB) are downloaded automatically into the `data/` folder of the repository the first time those modules run, so the first run needs an internet connection and disk space.

### `07_probe_filter.R`

Filters probes after Noob normalization (see [Probe filtering](#probe-filtering)). Each chunk is filtered on its own, with its own samples.

### `08_bmiq.R`

Runs BMIQ normalization on the probe-filtered data. BMIQ works on each sample separately, so chunking does not change its results. It can be time-consuming: use `--to=07` with `run_all_steps.R` to stop before it and inspect the earlier outputs first. The `bmiq` section of the config sets the number of cores and whether BMIQ's own progress output is shown.

### `09_bmiq_downstream.R`

Works on the BMIQ-normalized betas:

1. calculates M values (saved per chunk);
2. runs FraminghamPACE (ElasticNet and Ridge), when enabled;
3. merges all downstream results (cell counts, PC clocks, SystemsAge, DunedinPACE, FraminghamPACE) from all chunks into one table, `Derived_Variables/{study_name}_cell_composition_and_clocks.csv`;
4. creates the summary figures;
5. writes the PDF report and a text file with the numbers for a methods section, and the sample QC tables of all chunks (see [What you get](#what-you-get)).

It runs on all chunks by default. With `--chunk=N`, the merged table, the figures and the report cover that chunk only, and their file names end in `_chunkN`.

## Sample QC checks

Every sample goes through the checks below. A sample is removed if it fails a check that is switched on in `sample_filtering`. The thresholds are in `qc_thresholds`; the values shown are the defaults.

| Check | A sample fails if |
|---|---|
| Signal intensity | the median log2 methylated or unmethylated intensity is below 10.5 |
| Bead count | more than 5% of its probes have fewer than 5 beads, or its mean bead count is 2 or less |
| Detection p-value | more than 1% of its probes have a detection p-value above 0.05 |
| EWAS tools control metrics (17) | it fails one of the Illumina control metrics, which include the bisulfite conversion controls (the thresholds are set by the `ewastools` package) |
| SNP outliers | its SNP outlier score (`ewastools`) is above -4. The SNP probes of the array should fall into three genotype clusters; a poor fit suggests a contaminated or failed sample |
| Sex mismatch | the sex predicted from the X and Y chromosome probes (`ewastools`) differs from the sex in the sample sheet. This suggests a swapped or mislabelled sample. By default this check only flags the sample: it is listed in the log and the report and kept in the analysis |

Notes:

- **EWAS tools metrics.** `ewastools::read_idats()` and `ewastools::control_metrics()` are called directly; `GCC.DNAmPipeline` is not required. Each of the 17 metrics is listed under `sample_filtering.ewas_tools_metrics` in the config and can be switched off there.
- **Not removing a sample.** Set any check to `false` in `sample_filtering` to record it but not use it for sample removal (in the example config, only the sex mismatch is set to `false`). Samples that fail such a check are listed as "flagged but kept" in the log and in the report, so you can look at them.
- **The sex check** only flags by default (`sex_mismatch: false` in `sample_filtering`). Samples whose predicted sex differs from the reported sex are listed in the log and in the report, with the reported and the predicted sex, and they are kept in the analysis, so that you can check the sex in your sample sheet and decide. To remove them instead, set `sex_mismatch: true` and rerun from step 04 (the later steps use the samples that step 04 keeps). The check needs the reported sex in the `Female` column (1 or 0) and both sexes in the data. A sample with a missing or inconclusive sex is never flagged. If the sex cannot be predicted for more than 20% of the samples, or more than half of the samples mismatch (which usually means the column is coded the wrong way round), the check is not run and a warning is shown.
- **Rerunning step 03.** The EWAS tools checks take a while, so they are calculated once and saved as `{study_name}_ewas_tools_qc.rds` in `data`. If you change `sex_col` or `snp_outlier_threshold` after running step 03, delete this file before rerunning, so that it is calculated again.

## Probe filtering

Step 07 removes probes in this order, and prints and saves how many each filter removed and why (`{study_name}_probe_filter_attrition_chunk{chunk}.csv`):

1. probes without a genomic position in the array annotation;
2. SNP probes (a SNP at the CpG or single-base-extension site);
3. cross-reactive and polymorphic probes (from the probe lists);
4. probes on the X and Y chromosomes (`remove_xy`);
5. probes with too few beads (`bead_min_count`, `bead_fail_fraction`);
6. probes that are not detected (`detection_p_threshold`, `detection_p_fail_fraction`).

Each filter can be switched on or off in the `probe_filtering` section of the config. To keep the X and Y chromosome probes, for example if all your samples have the same sex or you want to analyse the sex chromosomes, set `remove_xy: false`. If your samples are of both sexes, keep them removed unless you analyse them, because these probes differ strongly between the sexes.

The probes are filtered chunk by chunk: the bead count and detection p-value criteria ("in more than X% of the samples") are judged on the samples of each chunk. This has two consequences:

- With fewer than 100 samples in a chunk, one failing sample is already more than 1%, so a probe is removed if it fails in a single sample. Chunks of different sizes therefore use effectively different thresholds.
- Different chunks can end up with different probes. Combine chunks only on the probes they share.

The log of step 07 says so when the first applies, and the report shows the probe numbers of each chunk when the chunks differ.

The array type (EPICv1, EPICv2 or 450k) is detected from the data, so it does not need to be set. `remove_cross_reactive` and `remove_polymorphic` act as one filter, because the probe lists contain both kinds of probe.

For EPIC arrays, the pipeline includes Pidsley and McCartney probe lists under `data/probe_lists/`. These lists were published for EPICv1, so probes that exist only on EPICv2 are not covered by them. For 450k arrays, the config includes expected Chen and Benton filenames, but those source files still need to be added if running 450k data. The bespoke `ReliableProbeList.csv` is not included or used by this pipeline.

## Logs

Every script writes a log with everything it printed (output, messages, warnings and errors) to `logs/{script name}.log`, for example `03_run_raw_qc.log`. The file names have no date or time, and a new run overwrites the log of the previous run. With `--chunk=N`, the name ends in `_chunkN`, so runs for different chunks do not overwrite each other. `run_all_steps.R` writes its own `run_all_steps.log` with the start and end of each step and the summary.

A warning that is repeated many times is written to the log three times, and the total number of repeats is written at the end of the log. Warnings that are known to be harmless carry a short note saying so.

## Analysis methods

The pipeline currently implements:

- Raw IDAT loading with `minfi::read.metharray()`
- Control-metric QC, SNP outlier detection and sex checks with `ewastools`
- Detection P-value QC with `minfi::detectionP()`
- Beadcount QC from extended RGsets or `wateRmelon`
- Signal intensity QC using `minfi::preprocessRaw()` followed by methylated/unmethylated intensity summaries
- Control-probe PCA from red/green control intensities
- Noob normalization with `minfi::preprocessNoob()`
- Cell-type proportions, PC clocks, SystemsAge and DunedinPACE on Noob-normalized data
- Probe filtering with `minfi`, beadcount, detection P-values, and external probe lists, with the beadcount and detection criteria judged within each chunk
- BMIQ normalization with `wateRmelon::BMIQ()`
- M values (`log2(beta / (1 - beta))`) and FraminghamPACE on the BMIQ-normalized betas

Cell counts, PC clocks, SystemsAge and DunedinPACE run on Noob-normalized data before probe filtering. FraminghamPACE runs on BMIQ-normalized betas, as recommended by the author of the package.

## Additional documentation

- `BLUEPRINT.md`: design decisions and intended workflow
- `TODO.md`: current implementation status and next steps
- `IMPLEMENTATION_NOTES.md`: smoke-test status and important implementation notes
- `docs/input_specification.md`: input format details
- `docs/output_specification.md`: output file descriptions
- `docs/qc_decisions.md`: QC metrics versus filtering decisions
- `docs/chunking.md`: chunking behavior
- `docs/package_requirements.md`: package requirements

## Data sharing and security

Do not commit raw IDATs, RDS/RData outputs, local configs, AWS keys, or study-specific restricted data. The repository `.gitignore` excludes common raw and intermediate data products by default. The downloaded reference files (`*.qs2`) should be ignored too.

AWS/S3 download and upload scripts should live outside this public pipeline or in a private study-specific repository.

## Contact

For bugs or issues, please submit an Issue in this repository.