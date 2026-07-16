# Implementation Notes

## Chunking

Chunking is not only an IDAT loading strategy. It is supported through the full pipeline. Scripts should accept either `--chunk=all` or `--chunk=N`.

## Naming

Use `study_name` from the config for all outputs. Avoid literal study names such as `HRS` inside reusable code.

## QC Semantics

Keep metrics and filtering decisions separate:

- Metrics tables contain measured values.
- Pass/fail tables contain whether a metric passed its threshold.
- Filtering decisions tables show whether each failed metric was used for sample removal.

## EWAS Tools

EWAS tools control-metric QC is implemented locally in `R/ewas_tools_qc.R`. It uses `ewastools::read_idats()` and `ewastools::control_metrics()` directly, based on the validated manifest basenames. The pipeline no longer depends on `GCC.DNAmPipeline`.

## Local Test Data

Local test IDATs can be referenced by `config/local_test_config.yml`, but that file and all IDAT data are ignored by git.

## Smoke Test Status

- `scripts/00_validate_inputs.R --config=config/local_test_config.yml.example` succeeded on the local 24-sample test dataset.
- `scripts/01_create_targets_and_chunks.R --config=config/local_test_config.yml.example` succeeded and created 2 chunks.
- `scripts/02_load_rgset_chunks.R --config=config/local_test_config.yml.example --chunk=1` succeeded and wrote a raw RGset chunk.
- Installed `ewastools` from `hhhh5/ewastools` GitHub into the local R library.
- `scripts/03_run_raw_qc.R --config=config/local_test_config.yml.example --chunk=1` now succeeds and writes chunk 1 QC tables.
- During raw-QC testing, `calculate_signal_metrics()` was corrected to call `minfi::preprocessRaw()` before `minfi::getMeth()`/`minfi::getUnmeth()`, because those methods do not work directly on `RGChannelSetExtended`.
