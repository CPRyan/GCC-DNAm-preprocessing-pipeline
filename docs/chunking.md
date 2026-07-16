# Chunking

Chunking is supported throughout the pipeline.

Use `--chunk=all` to process all chunks or `--chunk=N` to process one chunk. After all chunk-level steps finish, run `scripts/09_merge_chunks.R`.

Grouping columns can be set in config so related samples remain in the same chunk when possible.
