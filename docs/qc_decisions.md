# QC Decisions

The pipeline separates measured QC values from filtering decisions.

## Metrics Table

Contains numeric/string values such as detection P-values, mean bead count, low-bead fraction, predicted sex, and donor ID.

## Pass/Fail Table

Contains all QC flags, including EWAS tools control metrics.

## Filtering Decisions Table

Uses the config to determine which failed flags actually remove samples.

Sex mismatch should often be collected but not used automatically, because it can be wrong or reflect metadata issues.
