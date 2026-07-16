# DNAm Preprocessing Pipeline

This repository, written in R, provides a clean, chunk-aware pipeline for
Illumina DNA methylation array preprocessing from raw IDAT files through
sample QC, Noob normalization, probe filtering, and BMIQ normalization.

The pipeline is designed for large datasets that may not fit in memory as
a single object. Chunks are supported throughout the workflow: users can
run the entire preprocessing sequence one chunk at a time, then merge
chunk-level outputs at the end.

The pipeline assumes users begin with a sample sheet and a directory of
paired IDAT files. AWS/S3 downloading and study-specific data access are
intentionally excluded from the public pipeline.

## Installation

Clone the repository:

```bash
git clone https://github.com/CPRyan/dnam_pipeline_clean.git
cd dnam_pipeline_clean
```

Install required R packages. At minimum, the pipeline uses `minfi`,
`ewastools`, `wateRmelon`, `yaml`, `ggplot2`, `UpSetR`, and array-specific
Illumina manifest/annotation packages.

`ewastools` is installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("hhhh5/ewastools")
```

Bioconductor packages can be installed with:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "minfi",
  "Biobase",
  "wateRmelon",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
))
```

Additional optional packages are needed for cell counts, PC clocks, and
DunedinPACE. See `docs/package_requirements.md`.

## Quick Example

Copy the example config and edit paths/columns for your study:

```bash
cp config/example_config.yml config/my_study.yml
```

Then run the pipeline step by step:

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
Rscript scripts/09_merge_chunks.R --config=config/my_study.yml
```

On memory-limited machines, run one chunk at a time:

```bash
Rscript scripts/02_load_rgset_chunks.R --config=config/my_study.yml --chunk=1
Rscript scripts/03_run_raw_qc.R --config=config/my_study.yml --chunk=1
Rscript scripts/04_filter_samples.R --config=config/my_study.yml --chunk=1
Rscript scripts/05_noob_normalize.R --config=config/my_study.yml --chunk=1
Rscript scripts/06_noob_downstream.R --config=config/my_study.yml --chunk=1
Rscript scripts/07_probe_filter.R --config=config/my_study.yml --chunk=1
Rscript scripts/08_bmiq.R --config=config/my_study.yml --chunk=1
```

Repeat for each chunk, then run `scripts/09_merge_chunks.R`.

## Pipeline Steps

### `00_validate_inputs.R`

Validates the sample sheet and confirms that every sample has paired
`_Red.idat` and `_Grn.idat` files. IDATs may be located directly in the
IDAT directory or in recursive subdirectories.

### `01_create_targets_and_chunks.R`

Creates the target manifest and assigns samples to chunks. If grouping
columns are provided, related samples are kept in the same chunk where
possible.

### `02_load_rgset_chunks.R`

Loads raw IDAT files into `minfi` `RGChannelSetExtended` objects, one RDS
file per chunk.

### `03_run_raw_qc.R`

Runs raw sample QC and mandatory EWAS tools control-metric QC. This step
produces separate metrics, pass/fail, filtering-decision, control-probe
PCA, and QC plot outputs.

### `04_filter_samples.R`

Uses the config file to decide which failed QC metrics remove samples.
EWAS tools metrics can be filtered independently, and sex mismatch is
collected but disabled for filtering by default.

### `05_noob_normalize.R`

Noob-normalizes good-sample RGsets and saves chunk-level normalized
`GenomicRatioSet` objects.

### `06_noob_downstream.R`

Runs downstream modules that should use Noob-normalized data before probe
filtering, including cell counts, PC clocks, and DunedinPACE when enabled.

### `07_probe_filter.R`

Filters probes after Noob normalization. Supported filters include SNP
probes, cross-reactive/polymorphic probes, XY probes, beadcount-failed
probes, and detection-P-failed probes.

### `08_bmiq.R`

Runs BMIQ normalization as a standard but explicit step. It is not hidden
inside the automatic workflow because it can be time-consuming and users
may want to inspect prior outputs first.

### `09_merge_chunks.R`

Merges standard chunk-level CSV outputs after all chunks have completed.

## Configuration

The pipeline is controlled by one YAML config file. The key value is
`study_name`, which is used to name all outputs:

```yaml
study_name: HRS
```

Example output names:

```text
HRS_RGset_raw_chunk1.rds
HRS_sample_qc_metrics_chunk1.csv
HRS_sample_qc_pass_fail_chunk1.csv
HRS_RGset_good_samples_chunk1.rds
HRS_Noob_GRset_chunk1.rds
HRS_probe_filtered_GRset_chunk1.rds
HRS_BMIQ_betas_chunk1.rds
```

Important config sections:

- **`paths`**: sample sheet, IDAT directory, output directory, probe-list directory
- **`sample_sheet`**: sample ID column, IDAT basename column, Sentrix columns, recursive IDAT search
- **`chunking`**: number of chunks, grouping columns, random seed
- **`qc_thresholds`**: sample-level QC thresholds
- **`sample_filtering.standard_metrics`**: standard QC flags used for sample removal
- **`sample_filtering.ewas_tools_metrics`**: EWAS tools metrics used for sample removal
- **`noob_downstream`**: optional downstream modules run on Noob-normalized data
- **`probe_filtering`**: probe filtering choices and probe-list filenames
- **`bmiq`**: BMIQ settings

## EWAS Tools QC

EWAS tools control-metric QC is mandatory. The pipeline calls
`ewastools::read_idats()` and `ewastools::control_metrics()` directly;
it does not require `GCC.DNAmPipeline`.

Each EWAS tools metric can be used or ignored during sample filtering:

```yaml
sample_filtering:
  ewas_tools_metrics:
    Restoration: true
    Staining Green: true
    Staining Red: true
    Extension Green: true
    Extension Red: true
    Hybridization High/Medium: true
    Hybridization Medium/Low: true
    Target Removal 1: true
    Target Removal 2: true
    Bisulfite Conversion I Green: true
    Bisulfite Conversion I Red: true
    Bisulfite Conversion II: true
    Specificity I Green: true
    Specificity I Red: true
    Specificity II: true
    Non-polymorphic Green: true
    Non-polymorphic Red: true
```

Set any metric to `false` to collect it but not use it for sample removal.

## Return Values and Outputs

The pipeline writes outputs to `paths.output_dir`, with data files under
`Output/data` and figures under `Output/figures` by default.

Key outputs include:

- **Validated manifest**: `{study_name}_sample_manifest_validated.csv`
- **Chunk assignments**: `{study_name}_targets_chunks.csv`
- **Raw RGsets**: `{study_name}_RGset_raw_chunk{chunk}.rds`
- **QC metrics**: `{study_name}_sample_qc_metrics_chunk{chunk}.csv`
- **QC pass/fail flags**: `{study_name}_sample_qc_pass_fail_chunk{chunk}.csv`
- **Filtering decisions**: `{study_name}_sample_qc_filtering_decisions_chunk{chunk}.csv`
- **Good-sample RGsets**: `{study_name}_RGset_good_samples_chunk{chunk}.rds`
- **Noob-normalized GRsets**: `{study_name}_Noob_GRset_chunk{chunk}.rds`
- **Probe-filtered GRsets**: `{study_name}_probe_filtered_GRset_chunk{chunk}.rds`
- **BMIQ beta matrices**: `{study_name}_BMIQ_betas_chunk{chunk}.rds`

QC outputs are intentionally separated into measured values and decisions:

- **Metrics table**: observed values such as mean bead count, low-bead fraction, average detection P-value, predicted sex, and EWAS tools values
- **Pass/fail table**: all QC flags, regardless of whether they are used for filtering
- **Filtering decisions table**: which failed metrics were actually used to remove each sample

## Chunking

For large DNAm datasets, the full dataset may exceed available memory.
This pipeline supports chunked processing across all major steps, not just
IDAT loading.

Use:

```bash
--chunk=all
```

to run all chunks in sequence, or:

```bash
--chunk=1
```

to run a single chunk.

If a run is interrupted, rerun the failed step for the relevant chunk.
Chunk-level output files are written separately, so completed chunks do
not need to be regenerated unless upstream settings changed.

## Data Format Requirements

### IDAT Files

The IDAT directory must contain paired Red and Green files for each
sample:

```text
203546000006_R01C01_Red.idat
203546000006_R01C01_Grn.idat
```

Files may be located in recursive subdirectories if
`sample_sheet.recursive_idat_search` is `true`.

### Sample Sheet

The sample sheet must be a CSV or Illumina sample sheet with a `[Data]`
section. It must contain either:

- A basename column matching IDAT filenames without `_Red.idat` or `_Grn.idat`
- Sentrix ID and Sentrix position columns that can be combined into the IDAT basename

Required config values:

```yaml
sample_sheet:
  sample_id_col: IDATid
  basename_col: null
  sentrix_id_col: Sentrix_ID
  sentrix_position_col: Sentrix_Position
  recursive_idat_search: true
```

Example sample sheet structure:

```text
Sample_Name,Sentrix_ID,Sentrix_Position,SubjectID,IDATid
sample_001,203546000006,R01C01,subj_001,203546000006_R01C01
sample_002,203546000006,R02C01,subj_002,203546000006_R02C01
```

If your sample sheet already has a basename column, set:

```yaml
sample_sheet:
  sample_id_col: barcode
  basename_col: barcode
```

## Probe Filtering

Probe filtering is configured under `probe_filtering`:

```yaml
probe_filtering:
  array_type: EPICv1
  remove_snp_probes: true
  remove_cross_reactive: true
  remove_polymorphic: true
  remove_xy: true
  remove_failed_detection: true
  remove_failed_beadcount: true
  bead_min_count: 3
  bead_fail_fraction: 0.05
  detection_p_threshold: 0.05
  detection_p_fail_fraction: 0.01
```

For EPICv1 arrays, the pipeline can use Pidsley and McCartney probe
lists. For 450k arrays, it can use Chen and Benton lists. The config
stores the expected filenames under `probe_filtering.cross_reactive_files`.

## Analysis Methods

The pipeline currently implements:

- Raw IDAT loading with `minfi::read.metharray()`
- Mandatory control-metric QC with `ewastools`
- Detection P-value QC with `minfi::detectionP()`
- Beadcount QC from extended RGsets or `wateRmelon`
- Signal intensity QC using `minfi::preprocessRaw()` followed by methylated/unmethylated intensity summaries
- Control-probe PCA from red/green control intensities
- Noob normalization with `minfi::preprocessNoob()`
- Probe filtering with `minfi`, beadcount, detection P-values, and external probe lists
- BMIQ normalization with `wateRmelon::BMIQ()`

Cell counts, PC clocks, and DunedinPACE are designed to run on
Noob-normalized data before probe filtering. The cell-count wrapper is
included; PC clock and DunedinPACE wrappers require local configuration
of their external dependencies before use.

## Additional Documentation

- `BLUEPRINT.md`: design decisions and intended workflow
- `TODO.md`: current implementation status and next steps
- `IMPLEMENTATION_NOTES.md`: smoke-test status and important implementation notes
- `docs/input_specification.md`: input format details
- `docs/output_specification.md`: output file descriptions
- `docs/qc_decisions.md`: QC metrics versus filtering decisions
- `docs/chunking.md`: chunking behavior
- `docs/package_requirements.md`: package requirements

## Data Sharing and Security

Do not commit raw IDATs, RDS/RData outputs, local configs, AWS keys, or
study-specific restricted data. The repository `.gitignore` excludes
common raw and intermediate data products by default.

AWS/S3 download and upload scripts should live outside this public
pipeline or in a private study-specific repository.

## Contact

For bugs or issues, please submit an Issue in this repository.
