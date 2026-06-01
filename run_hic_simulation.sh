#!/bin/bash
# Hi-C Simulation and Visualization Pipeline
# This script demonstrates the complete workflow:
# 1. Prepare reference genome and contact matrix
# 2. Simulate Hi-C reads using hicreate
# 3. Align reads using minimap2
# 4. Build Hi-C contact matrix
# 5. Visualize the heatmap

set -e

echo "========================================="
echo "Hi-C Simulation and Visualization Pipeline"
echo "========================================="

# Configuration
WORKDIR="hic_simulation_demo"
REF_FA="$WORKDIR/reference.fa"
MATRIX_TSV="$WORKDIR/matrix.tsv"
OFFSET_TSV="$WORKDIR/offset.tsv"
OUTPUT_PREFIX="$WORKDIR/simulated"
R1_FASTQ="$OUTPUT_PREFIX"_R1.fastq
R2_FASTQ="$OUTPUT_PREFIX"_R2.fastq
ALIGNED_SAM="$WORKDIR/aligned.sam"
PAIRS_FILE="$WORKDIR/aligned.pairs"
COOLER_FILE="$WORKDIR/contacts.cool"
HEATMAP_PNG="$WORKDIR/heatmap.png"

BIN_SIZE=10000
NUM_PAIRS=50000
ENZYME_SITE="AAGCTT"  # HindIII
SEED=42

# Create working directory
mkdir -p "$WORKDIR"

echo ""
echo "Step 1: Creating synthetic reference genome..."
echo "----------------------------------------------"

# Create a small synthetic reference genome (3 chromosomes, ~1Mb each)
python3 << 'EOF'
import random
import sys

random.seed(42)

def generate_chromosome(name, length):
    """Generate a random chromosome sequence"""
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return name, sequence

def write_fasta(filename, chromosomes):
    """Write chromosomes to FASTA file"""
    with open(filename, 'w') as f:
        for name, seq in chromosomes:
            f.write(f">{name}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# Generate 3 chromosomes of ~1Mb each
chromosomes = [
    generate_chromosome("chr1", 1000000),
    generate_chromosome("chr2", 1200000),
    generate_chromosome("chr3", 800000),
]

write_fasta("hic_simulation_demo/reference.fa", chromosomes)
print(f"Created reference genome with {len(chromosomes)} chromosomes")
for name, seq in chromosomes:
    print(f"  {name}: {len(seq)} bp")
EOF

echo ""
echo "Step 2: Generating realistic Hi-C contact matrix..."
echo "---------------------------------------------------"

# Generate a contact matrix with realistic features:
# - Strong diagonal (cis contacts decrease with distance)
# - Some compartment structure
# - A few loops/TADs
python3 << 'EOF'
import numpy as np
import sys

np.random.seed(42)

# Chromosome sizes in bins (10kb bins)
chr_sizes = [100, 120, 80]  # chr1: 1Mb, chr2: 1.2Mb, chr3: 800kb
total_bins = sum(chr_sizes)

# Initialize contact matrix
matrix = np.zeros((total_bins, total_bins))

# Function to add cis contacts with distance decay
def add_cis_contacts(start_bin, size, strength=1000):
    """Add cis contacts with power-law distance decay"""
    for i in range(size):
        for j in range(i, size):
            distance = abs(j - i)
            if distance == 0:
                weight = strength * 10
            else:
                # Power law decay: P(s) ~ s^(-1)
                weight = strength / (distance ** 0.75)
                # Add some noise
                weight *= np.random.lognormal(0, 0.3)
            
            global_i = start_bin + i
            global_j = start_bin + j
            matrix[global_i, global_j] = weight
            matrix[global_j, global_i] = weight

# Function to add compartments (A/B compartments)
def add_compartments(start_bin, size, compartment_strength=2):
    """Add A/B compartment structure"""
    compartment = np.ones(size)
    # Create alternating A/B compartments (~200kb = 20 bins)
    compartment_size = 20
    for i in range(0, size, compartment_size * 2):
        end = min(i + compartment_size, size)
        compartment[i:end] = -1
    
    # Add compartment-compartment interactions
    for i in range(size):
        for j in range(i, size):
            if compartment[i] * compartment[j] > 0:  # Same compartment
                matrix[start_bin + i, start_bin + j] *= compartment_strength
                matrix[start_bin + j, start_bin + i] *= compartment_strength

# Function to add TAD-like structures
def add_tad(start_bin, tad_start, tad_end, loop_strength=50):
    """Add TAD boundary and internal loops"""
    tad_size = tad_end - tad_start
    
    # Strengthen internal contacts
    for i in range(tad_start, tad_end):
        for j in range(i, tad_end):
            matrix[start_bin + i, start_bin + j] *= 1.5
            matrix[start_bin + j, start_bin + i] *= 1.5
    
    # Add corner peak (loop)
    if tad_size > 10:
        loop_i = tad_start + tad_size // 3
        loop_j = tad_start + 2 * tad_size // 3
        matrix[start_bin + loop_i, start_bin + loop_j] += loop_strength
        matrix[start_bin + loop_j, start_bin + loop_i] += loop_strength

# Generate contacts for each chromosome
current_bin = 0
for chr_idx, size in enumerate(chr_sizes):
    print(f"Generating contacts for chr{chr_idx+1} ({size} bins)...")
    
    # Add basic cis contacts
    add_cis_contacts(current_bin, size, strength=500)
    
    # Add compartment structure
    add_compartments(current_bin, size, compartment_strength=1.8)
    
    # Add 2-3 TADs per chromosome
    num_tads = np.random.randint(2, 4)
    for _ in range(num_tads):
        tad_size = np.random.randint(15, 30)
        tad_start = np.random.randint(0, size - tad_size)
        add_tad(current_bin, tad_start, tad_start + tad_size, 
                loop_strength=np.random.uniform(30, 80))
    
    current_bin += size

# Add weak trans contacts (inter-chromosomal)
trans_strength = 5
for i in range(total_bins):
    for j in range(i+1, total_bins):
        # Check if on different chromosomes
        chr_i = 0 if i < chr_sizes[0] else (1 if i < chr_sizes[0] + chr_sizes[1] else 2)
        chr_j = 0 if j < chr_sizes[0] else (1 if j < chr_sizes[0] + chr_sizes[1] else 2)
        
        if chr_i != chr_j:
            weight = trans_strength * np.random.exponential(1)
            matrix[i, j] = weight
            matrix[j, i] = weight

# Convert to sparse format and write
print("Writing sparse matrix...")
with open("hic_simulation_demo/matrix.tsv", 'w') as f:
    f.write("bin1\tbin2\tvalue\n")
    count = 0
    for i in range(total_bins):
        for j in range(i, total_bins):
            if matrix[i, j] > 0.1:  # Filter very low values
                f.write(f"{i}\t{j}\t{matrix[i, j]:.3f}\n")
                count += 1

print(f"Wrote {count} contacts to matrix.tsv")
print(f"Matrix size: {total_bins} x {total_bins}")
print(f"Total contact weight: {np.sum(matrix):.2f}")

# Calculate cis/trans ratio
cis_weight = 0
trans_weight = 0
current_bin = 0
for chr_idx, size in enumerate(chr_sizes):
    chr_start = current_bin
    chr_end = current_bin + size
    cis_weight += np.sum(matrix[chr_start:chr_end, chr_start:chr_end])
    current_bin += size

trans_weight = np.sum(matrix) - cis_weight
print(f"Cis contacts: {cis_weight:.2f} ({100*cis_weight/(cis_weight+trans_weight):.1f}%)")
print(f"Trans contacts: {trans_weight:.2f} ({100*trans_weight/(cis_weight+trans_weight):.1f}%)")
EOF

echo ""
echo "Step 3: Creating offset file..."
echo "-------------------------------"

# Create offset file mapping bins to chromosomes
python3 << 'EOF'
chr_sizes = [100, 120, 80]  # in bins
current_bin = 0

with open("hic_simulation_demo/offset.tsv", 'w') as f:
    f.write("contig\tstart_bin\tend_bin\n")
    for i, size in enumerate(chr_sizes):
        f.write(f"chr{i+1}\t{current_bin}\t{current_bin + size}\n")
        current_bin += size

print("Created offset.tsv")
print("Contents:")
with open("hic_simulation_demo/offset.tsv", 'r') as f:
    print(f.read())
EOF

echo ""
echo "Step 4: Simulating Hi-C reads with hicreate..."
echo "----------------------------------------------"

# Check if hicreate is built
if [ ! -f "./hicreate" ]; then
    echo "Building hicreate..."
    bash build.sh
fi

# Run hicreate simulation
./hicreate "$REF_FA" $BIN_SIZE \
    --matrix "$MATRIX_TSV" \
    --offset "$OFFSET_TSV" \
    --pairs $NUM_PAIRS \
    --output-prefix "$OUTPUT_PREFIX" \
    --enzyme-site "$ENZYME_SITE" \
    --seed $SEED \
    --threads 4

echo ""
echo "Simulated FASTQ files:"
ls -lh "$R1_FASTQ" "$R2_FASTQ" 2>/dev/null || echo "Error: FASTQ files not created"

echo ""
echo "Step 5: Building reference index with minimap2..."
echo "-------------------------------------------------"

# Check if minimap2 is available
if ! command -v minimap2 &> /dev/null; then
    echo "Warning: minimap2 not found. Installing..."
    # Try to install minimap2
    if command -v conda &> /dev/null; then
        conda install -c bioconda minimap2 -y
    elif command -v apt-get &> /dev/null; then
        sudo apt-get update && sudo apt-get install -y minimap2
    else
        echo "Please install minimap2 manually: https://github.com/lh3/minimap2"
        exit 1
    fi
fi

# Build index
minimap2 -x map-pb -d "$WORKDIR/ref.mmi" "$REF_FA"
echo "Reference index created: $WORKDIR/ref.mmi"

echo ""
echo "Step 6: Aligning simulated reads with minimap2..."
echo "-------------------------------------------------"

# Align R1 and R2 separately, then merge
# For Hi-C data, we need to align both ends
minimap2 -ax map-pb --secondary=no "$WORKDIR/ref.mmi" "$R1_FASTQ" > "$WORKDIR/r1.sam" 2>"$WORKDIR/r1.log"
minimap2 -ax map-pb --secondary=no "$WORKDIR/ref.mmi" "$R2_FASTQ" > "$WORKDIR/r2.sam" 2>"$WORKDIR/r2.log"

echo "Alignment complete"
echo "R1 alignments: $(grep -c '^@' "$WORKDIR/r1.sam" || echo 0) headers, $(grep -cv '^@' "$WORKDIR/r1.sam" || echo 0) reads"
echo "R2 alignments: $(grep -c '^@' "$WORKDIR/r2.sam" || echo 0) headers, $(grep -cv '^@' "$WORKDIR/r2.sam" || echo 0) reads"

echo ""
echo "Step 7: Processing alignments into Hi-C pairs..."
echo "------------------------------------------------"

# Create a simple pairs file from the SAM alignments
# In production, you would use pairtools or similar
python3 << 'EOF'
import re
from collections import defaultdict

def parse_sam(sam_file):
    """Parse SAM file and extract alignment information"""
    alignments = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            
            # Only keep mapped reads with good quality
            if rname == '*' or mapq < 20:
                continue
            
            # Extract read name base (remove /1 or /2 suffix)
            base_name = re.sub(r'/[12]$', '', qname)
            read_num = 1 if qname.endswith('/1') else 2
            
            alignments[(base_name, read_num)] = {
                'chrom': rname,
                'pos': pos,
                'mapq': mapq,
                'cigar': cigar
            }
    
    return alignments

# Parse both SAM files
print("Parsing R1 alignments...")
r1_alignments = parse_sam("hic_simulation_demo/r1.sam")
print(f"  Found {len(r1_alignments)} R1 alignments")

print("Parsing R2 alignments...")
r2_alignments = parse_sam("hic_simulation_demo/r2.sam")
print(f"  Found {len(r2_alignments)} R2 alignments")

# Find paired alignments
paired_count = 0
valid_pairs = []

for base_name in r1_alignments:
    if base_name in r2_alignments:
        r1 = r1_alignments[base_name]
        r2 = r2_alignments[base_name]
        
        # Create a pair entry
        # Format: read_id \t chrom1 \t pos1 \t strand1 \t chrom2 \t pos2 \t strand2
        valid_pairs.append({
            'id': base_name,
            'chrom1': r1['chrom'],
            'pos1': r1['pos'],
            'chrom2': r2['chrom'],
            'pos2': r2['pos'],
        })
        paired_count += 1

print(f"\nFound {paired_count} valid read pairs")

# Write pairs file
with open("hic_simulation_demo/aligned.pairs", 'w') as f:
    f.write("#read_id\tchrom1\tpos1\tstrand1\tchrom2\tpos2\tstrand2\n")
    for pair in valid_pairs:
        # Simplified: assume both reads are forward strand
        f.write(f"{pair['id']}\t{pair['chrom1']}\t{pair['pos1']}\t+\t{pair['chrom2']}\t{pair['pos2']}\t+\n")

print(f"Wrote {len(valid_pairs)} pairs to aligned.pairs")

# Statistics
cis_count = sum(1 for p in valid_pairs if p['chrom1'] == p['chrom2'])
trans_count = len(valid_pairs) - cis_count
print(f"\nCis pairs: {cis_count} ({100*cis_count/len(valid_pairs):.1f}%)")
print(f"Trans pairs: {trans_count} ({100*trans_count/len(valid_pairs):.1f}%)")
EOF

echo ""
echo "Step 8: Building contact matrix from aligned pairs..."
echo "-----------------------------------------------------"

# Bin the pairs and create a contact matrix
python3 << 'EOF'
import numpy as np
from collections import defaultdict

# Parameters
bin_size = 10000  # 10kb bins

# Chromosome sizes (should match reference)
chrom_sizes = {
    'chr1': 1000000,
    'chr2': 1200000,
    'chr3': 800000,
}

# Calculate number of bins per chromosome
chrom_bins = {chrom: (size + bin_size - 1) // bin_size 
              for chrom, size in chrom_sizes.items()}

# Create bin offset mapping
chrom_offsets = {}
current_offset = 0
for chrom in ['chr1', 'chr2', 'chr3']:
    chrom_offsets[chrom] = current_offset
    current_offset += chrom_bins[chrom]

total_bins = sum(chrom_bins.values())
print(f"Total bins: {total_bins}")
print(f"Chromosome bins: {chrom_bins}")

# Initialize contact matrix
contact_matrix = np.zeros((total_bins, total_bins))

# Read pairs and bin them
print("Binning read pairs...")
pair_count = 0
cis_count = 0
trans_count = 0

with open("hic_simulation_demo/aligned.pairs", 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 7:
            continue
        
        chrom1 = fields[1]
        pos1 = int(fields[2])
        chrom2 = fields[4]
        pos2 = int(fields[5])
        
        # Convert to bin coordinates
        bin1 = chrom_offsets[chrom1] + pos1 // bin_size
        bin2 = chrom_offsets[chrom2] + pos2 // bin_size
        
        # Ensure bin indices are within bounds
        if bin1 >= total_bins or bin2 >= total_bins:
            continue
        
        # Add to matrix (symmetric)
        contact_matrix[bin1, bin2] += 1
        if bin1 != bin2:
            contact_matrix[bin2, bin1] += 1
        
        pair_count += 1
        if chrom1 == chrom2:
            cis_count += 1
        else:
            trans_count += 1

print(f"Processed {pair_count} pairs")
print(f"Cis: {cis_count}, Trans: {trans_count}")

# Save the contact matrix
np.save("hic_simulation_demo/contact_matrix.npy", contact_matrix)

# Also save as sparse format for compatibility
with open("hic_simulation_demo/reconstructed_matrix.tsv", 'w') as f:
    f.write("bin1\tbin2\tvalue\n")
    count = 0
    for i in range(total_bins):
        for j in range(i, total_bins):
            if contact_matrix[i, j] > 0:
                f.write(f"{i}\t{j}\t{contact_matrix[i, j]:.0f}\n")
                count += 1

print(f"Wrote {count} non-zero contacts to reconstructed_matrix.tsv")

# Print some statistics
print(f"\nContact matrix statistics:")
print(f"  Shape: {contact_matrix.shape}")
print(f"  Total contacts: {np.sum(contact_matrix):.0f}")
print(f"  Mean contacts per bin: {np.mean(contact_matrix):.2f}")
print(f"  Max contacts: {np.max(contact_matrix):.0f}")
EOF

echo ""
echo "Step 9: Visualizing Hi-C contact heatmap..."
echo "-------------------------------------------"

# Create visualization using Python/matplotlib
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-white')
sns.set_context("paper", font_scale=1.2)

# Load contact matrix
contact_matrix = np.load("hic_simulation_demo/contact_matrix.npy")

# Chromosome information
chrom_info = [
    ('chr1', 100),
    ('chr2', 120),
    ('chr3', 80),
]

# Apply log transformation for better visualization
# Add small pseudocount to avoid log(0)
log_matrix = np.log10(contact_matrix + 1)

# Create figure with multiple panels
fig = plt.figure(figsize=(16, 14))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# Panel 1: Full contact matrix heatmap
ax1 = fig.add_subplot(gs[0:2, 0:2])
im1 = ax1.imshow(log_matrix, cmap='YlOrRd', aspect='auto', 
                 norm=colors.Normalize(vmin=0, vmax=np.percentile(log_matrix, 99)))
ax1.set_title('Hi-C Contact Matrix (log10 scale)', fontsize=14, fontweight='bold')
ax1.set_xlabel('Genomic Position (bins)', fontsize=12)
ax1.set_ylabel('Genomic Position (bins)', fontsize=12)

# Add chromosome boundaries
current_pos = 0
for chrom_name, chrom_size in chrom_info:
    ax1.axhline(y=current_pos, color='blue', linewidth=1, alpha=0.5)
    ax1.axvline(x=current_pos, color='blue', linewidth=1, alpha=0.5)
    # Label chromosomes
    mid_pos = current_pos + chrom_size / 2
    ax1.text(mid_pos, -5, chrom_name, ha='center', va='top', 
             fontsize=10, rotation=45)
    current_pos += chrom_size

plt.colorbar(im1, ax=ax1, label='log10(contacts + 1)')

# Panel 2: Diagonal profile (average contacts vs genomic distance)
ax2 = fig.add_subplot(gs[0, 2])
distances = []
avg_contacts = []

max_distance = min(100, contact_matrix.shape[0] // 2)
for d in range(max_distance):
    diag_values = []
    for i in range(contact_matrix.shape[0] - d):
        if contact_matrix[i, i+d] > 0:
            diag_values.append(contact_matrix[i, i+d])
    if diag_values:
        distances.append(d)
        avg_contacts.append(np.mean(diag_values))

ax2.plot(distances, avg_contacts, 'b-', linewidth=2)
ax2.set_xlabel('Genomic Distance (bins)', fontsize=11)
ax2.set_ylabel('Average Contacts', fontsize=11)
ax2.set_title('Contact Decay Profile', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

# Panel 3: Cis vs Trans distribution
ax3 = fig.add_subplot(gs[1, 2])
cis_contacts = []
trans_contacts = []

current_bin = 0
for chrom_name, chrom_size in chrom_info:
    chrom_start = current_bin
    chrom_end = current_bin + chrom_size
    
    # Cis contacts (within chromosome)
    cis_sum = np.sum(contact_matrix[chrom_start:chrom_end, chrom_start:chrom_end])
    cis_contacts.append(cis_sum)
    
    # Trans contacts (with other chromosomes)
    trans_sum = 0
    for other_start, (other_name, other_size) in enumerate(chrom_info):
        if other_name != chrom_name:
            other_bin_start = sum(s for n, s in chrom_info[:other_start])
            other_bin_end = other_bin_start + other_size
            trans_sum += np.sum(contact_matrix[chrom_start:chrom_end, 
                                              other_bin_start:other_bin_end])
    trans_contacts.append(trans_sum)
    
    current_bin += chrom_size

x_pos = np.arange(len(chrom_info))
width = 0.35
bars1 = ax3.bar(x_pos - width/2, cis_contacts, width, label='Cis', color='steelblue')
bars2 = ax3.bar(x_pos + width/2, trans_contacts, width, label='Trans', color='coral')

ax3.set_xlabel('Chromosome', fontsize=11)
ax3.set_ylabel('Total Contacts', fontsize=11)
ax3.set_title('Cis vs Trans Contacts', fontsize=12, fontweight='bold')
ax3.set_xticks(x_pos)
ax3.set_xticklabels([c[0] for c in chrom_info])
ax3.legend()
ax3.set_yscale('log')

# Panel 4: Individual chromosome heatmaps
chrom_matrices = []
current_bin = 0
for chrom_name, chrom_size in chrom_info:
    chrom_mat = contact_matrix[current_bin:current_bin+chrom_size, 
                               current_bin:current_bin+chrom_size]
    chrom_matrices.append((chrom_name, chrom_mat))
    current_bin += chrom_size

for idx, (chrom_name, chrom_mat) in enumerate(chrom_matrices):
    ax = fig.add_subplot(gs[2, idx])
    log_chrom = np.log10(chrom_mat + 1)
    im = ax.imshow(log_chrom, cmap='YlOrRd', aspect='auto',
                   norm=colors.Normalize(vmin=0, vmax=np.percentile(log_chrom, 99)))
    ax.set_title(f'{chrom_name}', fontsize=11, fontweight='bold')
    ax.set_xlabel('Position (bins)', fontsize=9)
    ax.set_ylabel('Position (bins)', fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# Save figure
plt.savefig("hic_simulation_demo/heatmap.png", dpi=300, bbox_inches='tight')
print("Heatmap saved to hic_simulation_demo/heatmap.png")

# Also create a simplified version focusing on chr1
fig2, axes = plt.subplots(1, 2, figsize=(14, 6))

# Original matrix (chr1 only)
chr1_mat = chrom_matrices[0][1]
log_chr1 = np.log10(chr1_mat + 1)
im1 = axes[0].imshow(log_chr1, cmap='YlOrRd', aspect='auto')
axes[0].set_title('Original Contact Matrix (chr1)', fontsize=13, fontweight='bold')
axes[0].set_xlabel('Genomic Position (10kb bins)', fontsize=11)
axes[0].set_ylabel('Genomic Position (10kb bins)', fontsize=11)
plt.colorbar(im1, ax=axes[0], label='log10(contacts + 1)')

# Reconstructed matrix (simulate what we got from alignment)
# Load the reconstructed matrix
reconstructed = np.load("hic_simulation_demo/contact_matrix.npy")
chr1_recon = reconstructed[0:100, 0:100]
log_chr1_recon = np.log10(chr1_recon + 1)
im2 = axes[1].imshow(log_chr1_recon, cmap='YlOrRd', aspect='auto')
axes[1].set_title('Reconstructed from Aligned Reads (chr1)', fontsize=13, fontweight='bold')
axes[1].set_xlabel('Genomic Position (10kb bins)', fontsize=11)
axes[1].set_ylabel('Genomic Position (10kb bins)', fontsize=11)
plt.colorbar(im2, ax=axes[1], label='log10(contacts + 1)')

plt.tight_layout()
plt.savefig("hic_simulation_demo/chr1_comparison.png", dpi=300, bbox_inches='tight')
print("Chromosome 1 comparison saved to hic_simulation_demo/chr1_comparison.png")

# Print summary statistics
print("\n" + "="*50)
print("SIMULATION SUMMARY")
print("="*50)
print(f"Reference genome: 3 chromosomes, {sum(s for _, s in chrom_info)} bins")
print(f"Simulated read pairs: 50,000")
print(f"Bin size: 10 kb")
print(f"Enzyme: HindIII (AAGCTT)")
print(f"\nContact matrix properties:")
print(f"  Total contacts: {np.sum(contact_matrix):.0f}")
print(f"  Non-zero bins: {np.count_nonzero(contact_matrix)}")
print(f"  Sparsity: {100*(1 - np.count_nonzero(contact_matrix)/contact_matrix.size):.1f}%")
print(f"\nVisualization files created:")
print(f"  - hic_simulation_demo/heatmap.png")
print(f"  - hic_simulation_demo/chr1_comparison.png")
print("="*50)
EOF

echo ""
echo "========================================="
echo "Pipeline Complete!"
echo "========================================="
echo ""
echo "Output files:"
echo "  Reference:     $REF_FA"
echo "  Input matrix:  $MATRIX_TSV"
echo "  Offset file:   $OFFSET_TSV"
echo "  Simulated R1:  $R1_FASTQ"
echo "  Simulated R2:  $R2_FASTQ"
echo "  Alignments:    $WORKDIR/r1.sam, r2.sam"
echo "  Pairs file:    $PAIRS_FILE"
echo "  Contact matrix: $WORKDIR/contact_matrix.npy"
echo "  Heatmap:       $HEATMAP_PNG"
echo ""
echo "To view the heatmap:"
echo "  display $HEATMAP_PNG"
echo "  # or open with any image viewer"
echo ""
