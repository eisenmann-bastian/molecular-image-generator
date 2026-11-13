# Project instructions

- Keep input processing and rendering in `molecular_image_generator/`; preserve
  named and SMARTS selection, CSV row isolation, and no-overwrite behavior.
- Configure `molecular_utilities.logging` only in the CLI using the
  `molecular_image_generator` namespace; do not configure on import.
- Update tests, examples, and the [README](README.md) for behavior changes.
  Run `./examples/run-tests.sh` with an installed environment.
- Example name lookups use PubChem and need network access; do not overwrite
  existing output or environments.
