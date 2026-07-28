# hicreate

`hicreate` is a small C++ command-line program for generating 150 bp paired-end Hi-C reads from a reference genome and an optional contact matrix.

The current workflow is:

1. Read a reference FASTA with one or more contigs.
2. Read a Hi-C contact matrix in sparse or dense form, or use an empirical model matrix if no matrix is supplied.
3. Read an optional offset file that maps matrix bins back to source contigs.
4. Perform an in silico restriction digest.
5. Sample fragment-fragment ligation events according to matrix contact frequencies.
6. Build the two ends of each virtual ligation molecule from restriction-fragment coordinates.
7. Write paired FASTQ files with built-in Illumina-like quality and substitution errors.

## Build

This project now uses lightweight compile scripts for day-to-day builds so it does not create extra `build*` directories during local development.

Windows:

```powershell
powershell -ExecutionPolicy Bypass -File build.ps1
```

Linux:

```bash
sh build.sh
```

Outputs:

- Windows: `hicreate.exe`
- Linux: `hicreate`

## Command Line

```text
hicreate ref.fa 1000 \
         [--matrix matrix.tsv] \
         [--matrix-format auto|sparse|dense|binary] \
         [--empirical-model human_cell_40mb] [--model-dir models] \
         [--coverage X | --pairs 100000] [--output-prefix sim] \
         [--offset offset.tsv] [--enzyme-site AAGCTT] [--seed 42] [--threads 4] \
         [--trans-ratio 0.12] [--force-contig-reuse] [--species-model auto]
```

Required arguments:

- `ref.fa`: input reference FASTA as the first positional argument
- `1000`: genomic bin size as the second positional argument

The older flag form is still supported for compatibility:

```text
hicreate --reference ref.fa --bin-size 1000 [options]
```

Optional arguments:

- `-m`, `--matrix`: input Hi-C matrix. If omitted, `hicreate` uses the empirical model linked to `--species-model`
- `--matrix-format`: matrix parser mode (`auto`, `sparse`, `dense`, `binary`), default `auto`
- `--empirical-model`, `--model`: empirical matrix model name or manifest path. Without `--matrix`, use the model as the full matrix. With `--matrix`, preserve input cis contacts and replace input trans contacts with the model's remapped trans contacts. If omitted, the input matrix is used as-is.
- `--model-dir`: empirical model directory, default `models`
- `-S`, `--species-model`: empirical preset selector. Default `auto`; currently `auto`, `human`, `homo_sapiens`, `mammal`, and `chm13` resolve to `human_cell_40mb` unless `--empirical-model` is provided
- `-f`, `--offset`: contig-to-global-bin mapping file
- `-c`, `--coverage`: target read depth over the reference genome; pairs are computed as `ceil(coverage * reference_bases / 300)`
- `-p`, `--pairs`: number of 150 bp paired-end read pairs to write, default `100000` when `--coverage` is omitted
- `-o`, `--output-prefix`: prefix for output files, default `sim`
- `-e`, `--enzyme-site`: restriction enzyme motif, default `AAGCTT`
- `-s`, `--seed`: random seed
- `-j`, `--threads`: worker threads for read generation, default `1`; use `0` to auto-detect hardware threads
- `-t`, `--trans-ratio`: target fraction of trans-chromosomal interaction mass. This rescales existing empirical trans contacts when an empirical model is in use
- `--force-contig-reuse`: allow an offset with fewer source contigs than the reference. Reused-target trans blocks are synthesized from real source trans blocks and the final trans mass is normalized

## Input Files

Reference FASTA:

- Standard FASTA with one or more contigs
- Sequence lines are concatenated per contig
- Header lines start with `>`

Sparse matrix:

```text
bin1    bin2    value
0       0       1000
0       20      47.619
40      60      47.619
```

- `bin1` and `bin2` are global bin indices
- `value` is the contact weight or contact frequency
- Rows with non-positive weights are ignored
- Header rows such as `bin1 bin2 value` are accepted

Dense matrix:

- Headerless square numeric matrix; use `--matrix-format dense` for 3x3 dense matrices to avoid sparse/dense ambiguity
- Mainly kept for compatibility with earlier versions
- If needed, it is resized to the expected global bin count

Offset file:

```text
contig  start_bin   end_bin
chr1    0           2000
chr2    2000        3500
```

- `start_bin` is inclusive
- `end_bin` is exclusive
- This lets one global sparse matrix describe multiple contigs or chromosomes
- hicmap-style offsets are also accepted directly:

```text
name    bin_offset  length      bin_num
chr1    0           40000000    2000
chr2    2000        30000000    1500
```

- For hicmap-style offsets, `start_bin = bin_offset` and `end_bin = bin_offset + bin_num`
- The `length` column is accepted for compatibility and is not used for bin interval calculation
	

Empirical model directory:

```text
models/human_cell_40mb/
  manifest.tsv
  matrix.tsv
  offset.tsv
```

`manifest.tsv` is a key/value file:

```text
name    human_cell_40mb
matrix  matrix.tsv
offset  offset.tsv
format  dense
```

- Supported formats are `sparse`, `dense`, and compact binary sparse `binary`/`.hcm`.
- If no `--matrix` and no `--empirical-model` are supplied, `--species-model` chooses the default empirical model.
- If `--empirical-model MODEL` is used without `--matrix`, the empirical model is remapped to the target reference and used as the full contact matrix.
- If `--empirical-model MODEL` is used with `--matrix`, input cis contacts are kept unchanged, while input trans contacts are replaced by the model's trans contacts after resizing/remapping and normalization.
- If `--trans-ratio` is also provided in trans-replacement mode, the replacement trans contacts are scaled to the requested final trans fraction. Otherwise they are scaled to the original input trans total.

## Example

Windows:

```powershell
.\hicreate.exe data\ref.fa 1000 -m data\matrix.tsv -f data\offset.tsv `
  -p 1000 -o data\sim
```

Linux:

```bash
./hicreate data/ref.fa 1000 -m data/matrix.tsv -f data/offset.tsv \
  -p 1000 -o data/sim
```

## Output

Main output:

- `<prefix>_R1.fastq`: first reads
- `<prefix>_R2.fastq`: second reads

`--pairs` is the exact number of records written to each FASTQ file.
If `--coverage` is provided, it replaces `--pairs` and computes the FASTQ record count from the reference size. `--coverage` and `--pairs` are mutually exclusive; using both is an error. There is no default 30x coverage unless you explicitly pass `--coverage 30`.

## Simulation Notes

- The simulation is driven by the input Hi-C matrix contact frequencies or by an empirical model selected with `--empirical-model`/`--species-model`.
- If `--matrix` is provided, input cis contacts are preserved and trans contacts are replaced by the selected empirical model.
- If both `--matrix` and an empirical model are provided, the input matrix provides cis contacts and the empirical model provides trans contacts.
- If `--matrix` is provided without `--empirical-model`, the input matrix is used as-is after remapping to the target reference.
- If `--matrix` and `--offset` exactly match the target reference contig names and bin ranges, the matrix is used directly and resize/remap is skipped.
- If `--matrix` and `--offset` do not perfectly match the target reference, same-name contigs are paired first; remaining contigs are paired one-to-one by offset order. Cis is scaled within each paired contig and trans is scaled between the corresponding paired contigs.
- By default, a matrix offset must contain at least as many contigs as the target reference. With `--force-contig-reuse`, source contigs may be reused; cis remains inside each target copy, while contacts between copies use a real trans block between that source and another source contig. The resulting trans mass is normalized to the pre-synthesis total.
- When a source contig and target contig have different bin spans, remapping is done per contig for cis and per contig-pair for trans, and each source-bin contact is distributed across neighboring target bins instead of collapsing the full weight into a single target bin.
- If `--matrix` is provided without `--offset`, the matrix is resized globally to the target genome bin count.
- If `--matrix` is omitted, the selected empirical model goes through the same offset check: exact offset matches are used directly, and only mismatched offsets are remapped to the target reference.
- Restriction digestion is modeled explicitly from the reference sequence.
- Common enzyme cut offsets are recognized for motifs such as HindIII and DpnII/MboI.
- Ligation products include an explicit fill-in ligation junction, but the program does not materialize full ligation molecules per read pair.
- For each sampled ligation event, reads are sliced from the ligated restriction-fragment ends outward. Short fragments are extended across the virtual ligation junction mathematically, without materializing long ligated sequences.
- Read start positions are sampled at random offsets within the virtual ligation templates, so repeated sampling from the same ligation event can yield different read windows.
- Each run reports the matrix cis/trans contact counts and weight fractions before FASTQ generation, which helps diagnose whether an unexpected map comes from the contact model or from downstream alignment/filtering.
- Read pairs are generated in bounded FASTQ blocks and streamed to disk, so memory does not scale with `--pairs` or `--coverage`.
- `--threads` parallelizes read-pair sampling, template construction, quality simulation, and FASTQ block formatting; a single writer preserves ordered output and avoids file-write lock contention.
- FASTQ qualities use an Illumina-like positional profile: high Q values near the start of each read, gradually lower Q values toward the 3' end, and Phred-derived substitution error probabilities.
- Only 150 bp paired-end output is currently supported.

## Included Example Data

The `data/` directory currently contains:

- `ref.fa`: example reference
- `matrix.tsv`: sparse example matrix
- `offset.tsv`: example offset file

These files are suitable for a quick paired FASTQ test run.

