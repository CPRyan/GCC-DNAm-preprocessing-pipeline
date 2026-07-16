# TODO

## Current Build Stage

- [x] Create clean pipeline directory.
- [x] Add persistent blueprint and implementation notes.
- [x] Add example config.
- [x] Add initial chunk-aware R modules and scripts.
- [ ] Port full legacy plotting functions with cleaned inputs/outputs.
- [x] Port initial array-aware probe filtering in detail.
- [ ] Validate probe filtering outputs against the legacy pipeline on a known dataset.
- [ ] Add package lockfile with `renv`.
- [ ] Run full test workflow on local non-shareable IDAT subset.
- [ ] Write final GitHub-ready README after test workflow succeeds.

## Important Decisions

- Keep original pipeline untouched.
- Do not commit IDATs, RDS outputs, AWS credentials, or local configs.
- EWAS tools QC is mandatory, but each EWAS tools metric is independently filterable.
- BMIQ is standard but explicitly run by the user.
- Cell counts, PC clocks, and DunedinPACE run on Noob-normalized data before probe filtering.
