# DNAm Pipeline Clean Blueprint

This directory is the clean, shareable version of the DNAm preprocessing pipeline. The original `Code/` and `Data/` directories are intentionally left untouched and can be used as references during development.

## Core Assumptions

- Users provide a sample sheet and a directory containing IDAT files.
- IDAT files may be nested in subdirectories.
- AWS/S3 loading is not part of the public pipeline.
- Example IDAT data used during development must not be committed or distributed.
- The pipeline is study-name driven: set `study_name` once, and outputs use that prefix.
- Chunking is pipeline-wide. Large datasets may need every downstream step run per chunk before chunk-level outputs are merged.
- EWAS tools QC is mandatory to run.
- EWAS tools filtering is configurable per metric.
- Sex mismatch is collected but not necessarily used for filtering.
- Cell counts, PC clocks, and DunedinPACE are run on Noob-normalized data before probe filtering.
- Probe filtering removes SNP, cross-reactive/polymorphic, XY, beadcount-failed, and detection-P-failed probes according to config.
- BMIQ is a standard step but must be run explicitly by the user.

## Intended Workflow

```text
Validate inputs
Create targets and chunk assignments
Load raw RGset chunks
Run raw QC and mandatory EWAS tools QC
Create QC metrics and pass/fail tables
Filter failed samples
Noob normalize good samples
Run Noob downstream modules: cell counts, PC clocks, DunedinPACE
Probe filtering
Run BMIQ explicitly
Merge chunk outputs
```

## Directory Design

```text
dnam_pipeline_clean/
├── README.md
├── BLUEPRINT.md
├── TODO.md
├── IMPLEMENTATION_NOTES.md
├── config/
├── R/
├── scripts/
├── docs/
└── data-raw/
```

## Output Naming

All generated filenames derive from `study_name` and, where needed, `chunk`:

```text
{study_name}_targets.csv
{study_name}_targets_chunks.csv
{study_name}_RGset_raw_chunk{chunk}.rds
{study_name}_sample_qc_metrics_chunk{chunk}.csv
{study_name}_sample_qc_pass_fail_chunk{chunk}.csv
{study_name}_sample_qc_filtering_decisions_chunk{chunk}.csv
{study_name}_RGset_good_samples_chunk{chunk}.rds
{study_name}_Noob_GRset_chunk{chunk}.rds
{study_name}_Noob_cell_counts_chunk{chunk}.csv
{study_name}_Noob_PC_clocks_chunk{chunk}.csv
{study_name}_Noob_DunedinPACE_chunk{chunk}.csv
{study_name}_probe_filtered_GRset_chunk{chunk}.rds
{study_name}_BMIQ_betas_chunk{chunk}.rds
```

## EWAS Tools QC Metrics

EWAS tools QC is run for all samples. Filtering can be toggled independently for each metric:

```text
Restoration
Staining Green
Staining Red
Extension Green
Extension Red
Hybridization High/Medium
Hybridization Medium/Low
Target Removal 1
Target Removal 2
Bisulfite Conversion I Green
Bisulfite Conversion I Red
Bisulfite Conversion II
Specificity I Green
Specificity I Red
Specificity II
Non-polymorphic Green
Non-polymorphic Red
```

## QC Table Contract

Each chunk produces three sample-level QC tables:

- Metrics table: numeric/string values such as mean bead count, fraction low bead count, average detection P-value, predicted sex, donor ID, and EWAS tools values.
- Pass/fail table: pass/fail or TRUE/FALSE decisions for all metrics.
- Filtering decisions table: records which failed metrics were actually used to remove samples.

## Continuity Note

If development continues in a new conversation, start by reading `BLUEPRINT.md`, `TODO.md`, and `IMPLEMENTATION_NOTES.md`.
