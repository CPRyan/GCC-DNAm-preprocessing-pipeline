# Output Specification

Outputs are written to `paths.output_dir` with subdirectories:

- `data/`: RDS and CSV data outputs.
- `figures/`: QC plots.
- `logs/`: reserved for logs.

Key output categories:

- Validated manifest and chunk assignment files.
- Raw RGsets by chunk.
- Sample QC metrics, pass/fail, and filtering decisions by chunk.
- Good-sample RGsets by chunk.
- Noob-normalized GRsets by chunk.
- Noob downstream outputs by chunk.
- Probe-filtered GRsets by chunk.
- BMIQ outputs by chunk.
- Merged final tables.
