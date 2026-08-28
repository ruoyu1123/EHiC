# hicreate

`hicreate` is a C++ command-line simulator for Hi-C, Pore-C, and CiFi chromosome-conformation reads from a reference genome and an optional contact matrix. Hi-C output is 150 bp paired-end Illumina-like FASTQ. Pore-C and CiFi output single-molecule, multi-fragment concatemer FASTQ plus segment-level truth tables.

The current workflow is:

1. Read a reference FASTA with one or more contigs.
2. Read a Hi-C contact matrix in sparse or dense form, or use an empirical model matrix if no matrix is supplied.
3. Read an optional offset file that maps matrix bins back to source contigs.
4. Perform an in silico restriction digest.
5. Sample fragment-fragment ligation events according to matrix contact frequencies.
6. For Hi-C, build the two ends of each virtual ligation molecule and write paired Illumina-like reads.
7. For Pore-C or CiFi, extend a multi-contact concatemer around a matrix-sampled anchor and write platform-specific long reads.

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
         [--assay hic|porec|cifi] \
         [--matrix matrix.tsv] \
         [--matrix-format auto|sparse|dense|binary] \
         [--empirical-model human_cell_40mb] [--model-dir models] \
         [--coverage X | --pairs 100000 | --reads 10000] [--output-prefix sim] \
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

- `--assay`: output assay, one of `hic`, `porec`, or `cifi`; default `hic`
- `-m`, `--matrix`: input Hi-C matrix. If omitted, `hicreate` uses the empirical model linked to `--species-model`
- `--matrix-format`: matrix parser mode (`auto`, `sparse`, `dense`, `binary`), default `auto`
- `--empirical-model`, `--model`: empirical matrix model name or manifest path. Without `--matrix`, use the model as the full matrix. With `--matrix`, preserve input cis contacts and replace input trans contacts with the model's remapped trans contacts. If omitted, the input matrix is used as-is.
- `--model-dir`: empirical model directory, default `models`
- `-S`, `--species-model`: empirical preset selector. Default `auto`; currently `auto`, `human`, `homo_sapiens`, `mammal`, and `chm13` resolve to `human_cell_40mb` unless `--empirical-model` is provided
- `-f`, `--offset`: contig-to-global-bin mapping file
- `-c`, `--coverage`: target depth; Hi-C pairs use `ceil(coverage * reference_bases / 300)`, while long-read counts use the configured mean read length
- `-p`, `--pairs`: number of 150 bp paired-end read pairs to write, default `100000` when `--coverage` is omitted
- `--reads`: number of single-molecule Pore-C or CiFi reads, default `10000` when `--coverage` is omitted
- `-o`, `--output-prefix`: prefix for output files, default `sim`
- `-e`, `--enzyme-site`: restriction enzyme motif; defaults are `AAGCTT` for Hi-C and `GATC` for Pore-C/CiFi
- `-s`, `--seed`: random seed
- `-j`, `--threads`: worker threads for read generation, default `1`; use `0` to auto-detect hardware threads
- `-t`, `--trans-ratio`: target fraction of trans-chromosomal interaction mass. This rescales existing empirical trans contacts when an empirical model is in use
- `--force-contig-reuse`: allow an offset with fewer source contigs than the reference. Reused-target trans blocks are synthesized from real source trans blocks and the final trans mass is normalized

Long-read-only arguments:

- `--long-read-mean`, `--long-read-min`, `--long-read-max`: override the assay read-length distribution
- `--long-read-sigma`: override the lognormal read-length sigma
- `--long-read-qv`: override mean long-read Phred QV
- `--max-segments`: maximum restriction-fragment segments in one concatemer, default `256`; the program warns if this limit truncates sampled molecules

Pore-C defaults are mean/min/max lengths of 10,000/3,000/100,000 bp, lognormal sigma 0.8, and mean QV 13 with Nanopore-like substitution, deletion, and insertion errors. CiFi defaults are 9,350/5,000/25,000 bp, sigma 0.65, mean QV 38, HiFi-like errors, and a 1.8% PCR-template duplicate rate.

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

Pore-C with DpnII digestion:

```bash
./hicreate ref.fa 1000 --assay porec -m matrix.tsv -f offset.tsv \
  --enzyme-site GATC --reads 10000 -j 8 -o porec_sim
```

CiFi with a target depth:

```bash
./hicreate ref.fa 1000 --assay cifi -m matrix.tsv -f offset.tsv \
  --enzyme-site GATC --coverage 5 -j 8 -o cifi_sim
```

## Output

Hi-C output:

- `<prefix>_R1.fastq`: first reads
- `<prefix>_R2.fastq`: second reads

Pore-C/CiFi output:

- `<prefix>_porec.fastq` or `<prefix>_cifi.fastq`: single-end concatemer reads
- `<prefix>_porec_truth.tsv` or `<prefix>_cifi_truth.tsv`: one row per source segment with contig coordinates, strand, matrix bin, and template coordinates

`--pairs` is the exact number of records written to each Hi-C FASTQ file.
If `--coverage` is provided, it replaces the assay's exact count option and computes the FASTQ record count from the reference size. `--coverage` cannot be combined with `--pairs` or `--reads`. There is no default 30x coverage unless you explicitly pass `--coverage 30`.

For Pore-C/CiFi, `--reads` is the exact FASTQ record count. Long-read coverage converts requested bases to a read count using the configured mean read length; because sampled lengths and sequencing indels vary, the final measured coverage can differ slightly and is reported at completion.

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
- Hi-C remains fixed at 150 bp paired-end output. Its fragment sampling and Illumina error implementation are separate from the new long-read modules.
- Pore-C and CiFi begin with one contact sampled exactly from the pairwise matrix. Additional segments are sampled from contact distributions conditioned on the first anchor bin. This anchor model avoids a trans contact turning into an artificial chromosome-scale random walk.
- A pairwise matrix does not uniquely determine a higher-order joint contact distribution. Consequently, all-vs-all pairs recovered from simulated concatemers should be close to, but are not mathematically guaranteed to reproduce, every input matrix weight exactly. Use the truth TSV to evaluate or recalibrate a downstream deconcatenation workflow.
- Segment orientation and restriction-fragment choice are random. Long-read sequences contain explicit ligation junctions and platform-specific substitutions and indels.
- CiFi template duplication models PCR amplification; independent HiFi sequencing errors are added to each duplicate copy.

## Protocol Basis

The Pore-C implementation follows the assay structure described by Deshpande et al., *Nanopore sequencing of DNA concatemers reveals higher-order features of chromatin structure* ([doi:10.1101/833590](https://doi.org/10.1101/833590)), and the size-selection/sequencing details reported for high-throughput Pore-C ([doi:10.1038/s41467-023-36899-x](https://doi.org/10.1038/s41467-023-36899-x)). These workflows crosslink chromatin, digest and proximity-ligate fragments, size-select multi-fragment DNA, and sequence individual concatemers with Oxford Nanopore.

The CiFi profile follows McGinty et al., *CiFi: accurate long-read chromosome conformation capture with low-input requirements* ([doi:10.1038/s41467-025-66918-y](https://doi.org/10.1038/s41467-025-66918-y)) and its library protocol ([doi:10.17504/protocols.io.4r3l21zxpg1y/v1](https://doi.org/10.17504/protocols.io.4r3l21zxpg1y/v1)). CiFi adds high-fidelity PCR enrichment, SMRTbell preparation, >5 kb size selection, and PacBio CCS/HiFi sequencing. The publication reports approximately 9.35 kb mean HiFi reads, median QV 38, and enzyme-dependent segment counts (DpnII median 17; HindIII median 2).

## Included Example Data

The `data/` directory currently contains:

- `ref.fa`: example reference
- `matrix.tsv`: sparse example matrix
- `offset.tsv`: example offset file

These files are suitable for a quick paired FASTQ test run.

