# Input Specification

## Required Inputs

- `sample_sheet`: CSV or Illumina sample sheet with a `[Data]` section.
- `idat_dir`: directory containing paired IDAT files.

## Required Sample Sheet Columns

At minimum, the sample sheet needs:

- A sample identifier column configured as `sample_sheet.sample_id_col`.
- Either a basename column configured as `sample_sheet.basename_col`, or both Sentrix ID and position columns.

The basename should match IDAT filenames without `_Red.idat` or `_Grn.idat`.

## Recursive IDAT Search

If `sample_sheet.recursive_idat_search` is `true`, IDATs may be nested in subdirectories under `paths.idat_dir`.
