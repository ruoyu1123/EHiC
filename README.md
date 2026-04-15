# hicreate

`hicreate` is a small C++ command-line program for generating a simple Hi-C ligation library from a reference genome and a contact matrix.

The current workflow is:

1. Read a reference FASTA with one or more contigs.
2. Read a Hi-C contact matrix in sparse or dense form.
3. Read an optional offset file that maps global bins back to contigs.
4. Perform an in silico restriction digest.
5. Sample fragment-fragment ligation events according to matrix contact frequencies.
6. Write ligation products to FASTA.
7. Optionally call `art_illumina` for read simulation.

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
hicreate --reference ref.fa [--matrix matrix.tsv] --bin-size 1000 \
         --read-length 150 --pairs 100000 --output-prefix sim \
         [--offset offset.tsv] [--enzyme-site AAGCTT] [--seed 42] \
         [--insert-mean 150] [--insert-std 25] [--skip-art] \
         [--trans-ratio 0.10] [--synthetic-contacts 200000] \
         [--cis-decay-alpha 1.0] [--max-cis-distance-bins 200] \
         [--species-model generic_plant] [--arrangement-model auto] \
         [--collision-randomness 0.35]
```

Required arguments:

- `--reference`: input reference FASTA
- `--bin-size`: genomic bin size used by the matrix
- `--read-length`: read length for ART
- `--pairs`: number of ligation products to sample
- `--output-prefix`: prefix for output files

Optional arguments:

- `--matrix`: input Hi-C matrix (if omitted, `hicreate` builds a synthetic matrix from the genome)
- `--offset`: contig-to-global-bin mapping file
- `--enzyme-site`: restriction enzyme motif, default `AAGCTT`
- `--seed`: random seed
- `--insert-mean`: ART insert mean
- `--insert-std`: ART insert standard deviation
- `--trans-ratio`: target fraction of trans-chromosomal interaction mass (default `0.10`)
- `--synthetic-contacts`: number of sparse contacts used to build synthetic matrix (auto if omitted)
- `--cis-decay-alpha`: cis distance-decay exponent (default `1.0`)
- `--max-cis-distance-bins`: max cis bin separation for synthetic matrix (default `200`)
- `--species-model`: synthetic matrix species preset (`generic_plant`, `arabidopsis`, `rice`, `maize`, `wheat`, `barley`)
- `--arrangement-model`: chromosome arrangement override (`auto`, `territory`, `rabl`, `rosette`, `nonrabl`)
- `--collision-randomness`: random collision baseline in trans sampling (0-1, default `0.35`)
- `--skip-art`: only write ligation FASTA, do not call ART

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

Dense matrix:

- Headerless square numeric matrix
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

## Example

Windows:

```powershell
.\hicreate.exe --reference data\ref.fa --matrix data\matrix.tsv --offset data\offset.tsv `
  --bin-size 1000 --read-length 150 --pairs 1000 --output-prefix data\sim --skip-art
```

Linux:

```bash
./hicreate --reference data/ref.fa --matrix data/matrix.tsv --offset data/offset.tsv \
  --bin-size 1000 --read-length 150 --pairs 1000 --output-prefix data/sim --skip-art
```

## Output

Main output:

- `<prefix>_fragments.fa`: sampled ligation product library

If ART is enabled and installed:

- paired FASTQ files produced by `art_illumina`

## Simulation Notes

- The simulation is driven by the input Hi-C matrix contact frequencies.
- If `--matrix` is provided, matrix weights are rebalanced to match the requested `--trans-ratio`.
- If `--matrix` is omitted, a sparse synthetic matrix is generated:
  - cis contacts follow a power-law distance decay (`1/(distance+1)^alpha`)
  - trans contacts are sampled across chromosomes according to `--trans-ratio`
  - species-aware defaults include chromosome arrangement effects (`Rabl`, `Rosette`, or territory-like)
  - random collision and chromosome arrangement are mixed via `--collision-randomness`
- Restriction digestion is modeled explicitly from the reference sequence.
- Common enzyme cut offsets are recognized for motifs such as HindIII and DpnII/MboI.
- Ligation products include an explicit fill-in ligation junction.
- Each sampled ligation event is written as an independent FASTA record instead of concatenating everything into one sequence.

## Included Example Data

The `data/` directory currently contains:

- `ref.fa`: example reference
- `matrix.tsv`: sparse example matrix
- `offset.tsv`: example offset file

These files are suitable for a quick `--skip-art` test run.
