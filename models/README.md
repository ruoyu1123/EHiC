# Empirical Hi-C Models

`hicreate` can use empirical matrix models as reusable Hi-C contact templates.

Each model lives in its own directory:

```text
models/<model-name>/
  manifest.tsv
  matrix.tsv   # or matrix.hcm
  offset.tsv
```

The manifest is a simple key/value file:

```text
name    human_cell_40mb
matrix  matrix.tsv
offset  offset.tsv
format  dense
```

Supported matrix formats are:

- `sparse`: text rows `bin1 bin2 value`
- `dense`: headerless square numeric matrix
- `binary`: compact binary sparse `.hcm` format

When a model is selected with `--empirical-model`:

- without `--matrix`, the model is remapped to the target reference and used as the full contact matrix;
- with `--matrix`, the input matrix keeps its cis contacts, while model trans contacts are remapped, normalized, and used to replace the input trans contacts.

When `--empirical-model` is omitted, `--species-model` selects the default empirical model. The bundled preset maps `auto`, `human`, `homo_sapiens`, `mammal`, and `chm13` to `human_cell_40mb`.

This lets users keep sample-specific cis structure while swapping trans architecture from a curated model.

## Included Starter Model

`human_cell_40mb` is a small starter model at 40 Mb bins, derived from the existing human-cell test matrix in this workspace. Replace or extend this directory with real curated matrices for production analyses.
