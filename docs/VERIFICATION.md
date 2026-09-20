# Verification status

## File-level provenance

The 32 copied scientific/workflow files were read back from `main` and compared against their expected upstream Git blob SHA-1 values.

**Result: 32 / 32 exact matches.**

This establishes byte-level identity of the imported frozen source files with the source commits recorded in `UPSTREAM_SOURCES.json`.

## Runtime status

This verification does **not** claim that the full 10,000-mock G0-G19 battery or 5,000-mock H0-H6 hardening battery has already executed successfully inside this new repository.

The repository contains the exact workflows and pinned numerical dependencies needed to run those batteries. Runtime execution and generated artifact hashes are a separate validation layer.
