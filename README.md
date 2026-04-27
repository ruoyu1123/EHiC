# hicreate

`hicreate` is a small C++ command-line program for generating 150 bp paired-end Hi-C reads from a reference genome and an optional contact matrix.

The current workflow is:

1. Read a reference FASTA with one or more contigs.
2. Read a Hi-C contact matrix in sparse or dense form, or build a synthetic one if no matrix is supplied.
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
         [--coverage X | --pairs 100000] [--output-prefix sim] \
         [--offset offset.tsv] [--enzyme-site AAGCTT] [--seed 42] [--threads 4] \
         [--trans-ratio 0.10] [--synthetic-contacts 200000] \
         [--cis-decay-alpha 1.0] [--max-cis-distance-bins 200] \
         [--species-model auto] [--arrangement-model auto] \
         [--trans-model auto] [--trans-hotspots 8] \
         [--collision-randomness 0.35]
```

Required arguments:

- `ref.fa`: input reference FASTA as the first positional argument
- `1000`: genomic bin size as the second positional argument

The older flag form is still supported for compatibility:

```text
hicreate --reference ref.fa --bin-size 1000 [options]
```

Optional arguments:

- `-m`, `--matrix`: input Hi-C matrix (if omitted, `hicreate` builds a synthetic matrix from the genome)
- `-f`, `--offset`: contig-to-global-bin mapping file
- `-c`, `--coverage`: target read depth over the reference genome; pairs are computed as `ceil(coverage * reference_bases / 300)`
- `-p`, `--pairs`: number of 150 bp paired-end read pairs to write, default `100000` when `--coverage` is omitted
- `-o`, `--output-prefix`: prefix for output files, default `sim`
- `-e`, `--enzyme-site`: restriction enzyme motif, default `AAGCTT`
- `-s`, `--seed`: random seed
- `-j`, `--threads`: worker threads for read generation, default `1`; use `0` to auto-detect hardware threads
- `-t`, `--trans-ratio`: target fraction of trans-chromosomal interaction mass (default `0.10`)
- `--synthetic-contacts`: number of sparse contacts used to build synthetic matrix (auto if omitted)
- `--cis-decay-alpha`: cis distance-decay exponent (default `1.0`)
- `--max-cis-distance-bins`: max cis bin separation for synthetic matrix (default `200`)
- `-S`, `--species-model`: synthetic matrix species preset (`auto`, `generic_plant`, `human`, `arabidopsis`, `rice`, `maize`, `wheat`, `barley`)
- `-A`, `--arrangement-model`: chromosome arrangement override (`auto`, `territory`, `rabl`, `rosette`, `nonrabl`)
- `-T`, `--trans-model`: trans interaction style (`auto`, `territory`, `random`, `telomere`, `centromere`, `compartment`, `hubs`)
- `--trans-hotspots`: number of trans hub bins for `--trans-model hubs`, default `8`
- `--collision-randomness`: random collision baseline in trans sampling (0-1, default `0.35`)

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

- The simulation is driven by the input Hi-C matrix contact frequencies.
- If `--matrix` is provided, matrix weights are rebalanced to match the requested `--trans-ratio`.
- If `--matrix` and `--offset` do not perfectly match the target reference, contacts are remapped onto the target contigs by name and relative position, then missing signal is filled with a synthetic Hi-C background.
- If `--matrix` is provided without `--offset`, the matrix is resized globally to the target genome bin count.
- If `--matrix` is omitted, a sparse synthetic matrix is generated:
  - cis contacts follow a power-law distance decay (`1/(distance+1)^alpha`)
  - trans contacts are sampled across chromosomes according to `--trans-ratio`
  - `territory` trans contacts emphasize chromosome-size and nuclear-distance effects
  - `random` trans contacts provide a uniform collision background
  - `telomere` trans contacts enrich different-chromosome telomere/subtelomere bins
  - `centromere` trans contacts enrich approximate pericentromeric/chromocenter-like bins
  - `compartment` trans contacts enrich same synthetic A/B-like compartment bins across chromosomes
  - `hubs` trans contacts create several point-like interchromosomal hot spots
  - species-aware defaults include chromosome arrangement effects (`Rabl`, `Rosette`, or territory-like)
  - random collision and chromosome arrangement are mixed via `--collision-randomness`
- Restriction digestion is modeled explicitly from the reference sequence.
- Common enzyme cut offsets are recognized for motifs such as HindIII and DpnII/MboI.
- Ligation products include an explicit fill-in ligation junction, but the program does not materialize full ligation molecules per read pair.
- For each sampled ligation event, reads are sliced from the ligated restriction-fragment ends outward. Short fragments are extended across the virtual ligation junction mathematically, without materializing long ligated sequences.
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
