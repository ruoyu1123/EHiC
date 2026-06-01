#!/usr/bin/env python3
"""
Hi-C Simulation and Visualization Demo
This script demonstrates:
1. Creating a synthetic reference genome
2. Generating a realistic Hi-C contact matrix
3. Simulating Hi-C reads using hicreate
4. Visualizing the original and simulated contact matrices
"""

import subprocess
import os
import sys
import numpy as np

def run_command(cmd, description):
    """Run a shell command and report progress"""
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"{'='*60}")
    print(f"Running: {cmd}\n")
    result = subprocess.run(cmd, shell=True, capture_output=False)
    if result.returncode != 0:
        print(f"Warning: Command failed with return code {result.returncode}")
        return False
    return True

def main():
    workdir = "hic_simulation_demo"
    os.makedirs(workdir, exist_ok=True)
    
    print("="*60)
    print("Hi-C Simulation and Visualization Pipeline")
    print("="*60)
    
    # Step 1: Create synthetic reference genome
    print("\n[Step 1/5] Creating synthetic reference genome...")
    ref_fa = f"{workdir}/reference.fa"
    
    create_ref_script = f'''
import random
random.seed(42)

def generate_chromosome(name, length):
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return name, sequence

chromosomes = [
    generate_chromosome("chr1", 1000000),
    generate_chromosome("chr2", 1200000),
    generate_chromosome("chr3", 800000),
]

with open("{ref_fa}", 'w') as f:
    for name, seq in chromosomes:
        f.write(f">{{name}}\\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\\n")

print(f"Created reference with {{len(chromosomes)}} chromosomes")
for name, seq in chromosomes:
    print(f"  {{name}}: {{len(seq)}} bp")
'''
    
    with open(f"{workdir}/create_ref.py", 'w') as f:
        f.write(create_ref_script)
    
    run_command(f"python3 {workdir}/create_ref.py", "Creating reference genome")
    
    # Step 2: Generate realistic contact matrix
    print("\n[Step 2/5] Generating Hi-C contact matrix with realistic features...")
    matrix_tsv = f"{workdir}/matrix.tsv"
    offset_tsv = f"{workdir}/offset.tsv"
    
    create_matrix_script = f'''
import numpy as np
np.random.seed(42)

# Chromosome sizes in bins (10kb bins)
chr_sizes = [100, 120, 80]  # chr1: 1Mb, chr2: 1.2Mb, chr3: 800kb
total_bins = sum(chr_sizes)
matrix = np.zeros((total_bins, total_bins))

def add_cis_contacts(start_bin, size, strength=1000):
    for i in range(size):
        for j in range(i, size):
            distance = abs(j - i)
            if distance == 0:
                weight = strength * 10
            else:
                weight = strength / (distance ** 0.75)
                weight *= np.random.lognormal(0, 0.3)
            
            global_i = start_bin + i
            global_j = start_bin + j
            matrix[global_i, global_j] = weight
            matrix[global_j, global_i] = weight

def add_compartments(start_bin, size, compartment_strength=2):
    compartment = np.ones(size)
    compartment_size = 20
    for i in range(0, size, compartment_size * 2):
        end = min(i + compartment_size, size)
        compartment[i:end] = -1
    
    for i in range(size):
        for j in range(i, size):
            if compartment[i] * compartment[j] > 0:
                matrix[start_bin + i, start_bin + j] *= compartment_strength
                matrix[start_bin + j, start_bin + i] *= compartment_strength

def add_tad(start_bin, tad_start, tad_end, loop_strength=50):
    tad_size = tad_end - tad_start
    for i in range(tad_start, tad_end):
        for j in range(i, tad_end):
            matrix[start_bin + i, start_bin + j] *= 1.5
            matrix[start_bin + j, start_bin + i] *= 1.5
    
    if tad_size > 10:
        loop_i = tad_start + tad_size // 3
        loop_j = tad_start + 2 * tad_size // 3
        matrix[start_bin + loop_i, start_bin + loop_j] += loop_strength
        matrix[start_bin + loop_j, start_bin + loop_i] += loop_strength

current_bin = 0
for chr_idx, size in enumerate(chr_sizes):
    print(f"Generating contacts for chr{{chr_idx+1}} ({{size}} bins)...")
    add_cis_contacts(current_bin, size, strength=500)
    add_compartments(current_bin, size, compartment_strength=1.8)
    
    num_tads = np.random.randint(2, 4)
    for _ in range(num_tads):
        tad_size = np.random.randint(15, 30)
        tad_start = np.random.randint(0, size - tad_size)
        add_tad(current_bin, tad_start, tad_start + tad_size, 
                loop_strength=np.random.uniform(30, 80))
    
    current_bin += size

# Add weak trans contacts
trans_strength = 5
for i in range(total_bins):
    for j in range(i+1, total_bins):
        chr_i = 0 if i < chr_sizes[0] else (1 if i < chr_sizes[0] + chr_sizes[1] else 2)
        chr_j = 0 if j < chr_sizes[0] else (1 if j < chr_sizes[0] + chr_sizes[1] else 2)
        
        if chr_i != chr_j:
            weight = trans_strength * np.random.exponential(1)
            matrix[i, j] = weight
            matrix[j, i] = weight

# Write sparse matrix
with open("{matrix_tsv}", 'w') as f:
    f.write("bin1\\tbin2\\tvalue\\n")
    count = 0
    for i in range(total_bins):
        for j in range(i, total_bins):
            if matrix[i, j] > 0.1:
                f.write(f"{{i}}\\t{{j}}\\t{{matrix[i, j]:.3f}}\\n")
                count += 1

print(f"Wrote {{count}} contacts to matrix.tsv")

# Calculate statistics
cis_weight = 0
trans_weight = 0
current_bin = 0
for chr_idx, size in enumerate(chr_sizes):
    chr_start = current_bin
    chr_end = current_bin + size
    cis_weight += np.sum(matrix[chr_start:chr_end, chr_start:chr_end])
    current_bin += size

trans_weight = np.sum(matrix) - cis_weight
print(f"Cis: {{cis_weight:.2f}} ({{100*cis_weight/(cis_weight+trans_weight):.1f}}%)")
print(f"Trans: {{trans_weight:.2f}} ({{100*trans_weight/(cis_weight+trans_weight):.1f}}%)")

# Save matrix for later use
np.save("{workdir}/original_matrix.npy", matrix)

# Write offset file
with open("{offset_tsv}", 'w') as f:
    f.write("contig\\tstart_bin\\tend_bin\\n")
    current_bin = 0
    for i, size in enumerate(chr_sizes):
        f.write(f"chr{{i+1}}\\t{{current_bin}}\\t{{current_bin + size}}\\n")
        current_bin += size

print("Created offset.tsv")
'''
    
    with open(f"{workdir}/create_matrix.py", 'w') as f:
        f.write(create_matrix_script)
    
    run_command(f"python3 {workdir}/create_matrix.py", "Generating contact matrix")
    
    # Step 3: Simulate Hi-C reads using hicreate
    print("\n[Step 3/5] Simulating Hi-C reads with hicreate...")
    
    output_prefix = f"{workdir}/simulated"
    cmd = (f"./hicreate {ref_fa} 10000 "
           f"--matrix {matrix_tsv} "
           f"--offset {offset_tsv} "
           f"--pairs 50000 "
           f"--output-prefix {output_prefix} "
           f"--enzyme-site AAGCTT "
           f"--seed 42 "
           f"--threads 4")
    
    success = run_command(cmd, "Running hicreate simulation")
    
    if not success:
        print("\nError: hicreate simulation failed!")
        print("Please check that hicreate is properly built.")
        return
    
    # Check output files
    r1_fastq = f"{output_prefix}_R1.fastq"
    r2_fastq = f"{output_prefix}_R2.fastq"
    
    if os.path.exists(r1_fastq) and os.path.exists(r2_fastq):
        print(f"\n✓ Successfully generated:")
        print(f"  {r1_fastq}: {os.path.getsize(r1_fastq)} bytes")
        print(f"  {r2_fastq}: {os.path.getsize(r2_fastq)} bytes")
    else:
        print("\n✗ Error: FASTQ files not created!")
        return
    
    # Step 4: Simulate alignment and create reconstructed matrix
    print("\n[Step 4/5] Simulating alignment and reconstructing contact matrix...")
    
    simulate_alignment_script = f'''
import numpy as np
np.random.seed(42)

# Load original matrix
original_matrix = np.load("{workdir}/original_matrix.npy")

# Simulate imperfect reconstruction from aligned reads
# Add sampling noise and reduce coverage
coverage_factor = 0.3  # Simulate ~30% of contacts captured
noise_level = 0.2

# Sample contacts based on original weights
reconstructed = np.zeros_like(original_matrix)
total_original = np.sum(original_matrix)

if total_original > 0:
    # Normalize to probabilities
    probs = original_matrix / total_original
    
    # Sample contacts
    num_samples = int(total_original * coverage_factor)
    indices = np.random.choice(len(probs.flat), size=num_samples, p=probs.flat)
    
    for idx in indices:
        i, j = np.unravel_index(idx, probs.shape)
        reconstructed[i, j] += 1
        if i != j:
            reconstructed[j, i] += 1
    
    # Add some noise
    noise = np.random.poisson(noise_level, reconstructed.shape)
    reconstructed = reconstructed + noise

# Save reconstructed matrix
np.save("{workdir}/reconstructed_matrix.npy", reconstructed)

# Statistics
print(f"Original matrix total contacts: {{np.sum(original_matrix):.0f}}")
print(f"Reconstructed matrix total contacts: {{np.sum(reconstructed):.0f}}")
print(f"Capture efficiency: {{100*np.sum(reconstructed)/np.sum(original_matrix):.1f}}%")
print(f"Non-zero bins: {{np.count_nonzero(reconstructed)}} / {{reconstructed.size}}")
'''
    
    with open(f"{workdir}/simulate_alignment.py", 'w') as f:
        f.write(simulate_alignment_script)
    
    run_command(f"python3 {workdir}/simulate_alignment.py", "Simulating alignment process")
    
    # Step 5: Visualize results
    print("\n[Step 5/5] Creating visualization...")
    
    visualize_script = f'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

plt.style.use('seaborn-v0_8-white')
sns.set_context("paper", font_scale=1.2)

# Load matrices
original = np.load("{workdir}/original_matrix.npy")
reconstructed = np.load("{workdir}/reconstructed_matrix.npy")

# Chromosome information
chrom_info = [
    ('chr1', 100),
    ('chr2', 120),
    ('chr3', 80),
]

# Apply log transformation
log_original = np.log10(original + 1)
log_reconstructed = np.log10(reconstructed + 1)

# Create comprehensive figure
fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.35)

# Panel 1: Original full matrix
ax1 = fig.add_subplot(gs[0, 0:2])
im1 = ax1.imshow(log_original, cmap='YlOrRd', aspect='auto',
                 norm=colors.Normalize(vmin=0, vmax=np.percentile(log_original, 99)))
ax1.set_title('Original Contact Matrix (Input)', fontsize=13, fontweight='bold')
ax1.set_xlabel('Genomic Position (bins)', fontsize=11)
ax1.set_ylabel('Genomic Position (bins)', fontsize=11)

# Add chromosome boundaries
current_pos = 0
for chrom_name, chrom_size in chrom_info:
    ax1.axhline(y=current_pos, color='blue', linewidth=1.5, alpha=0.6)
    ax1.axvline(x=current_pos, color='blue', linewidth=1.5, alpha=0.6)
    mid_pos = current_pos + chrom_size / 2
    ax1.text(mid_pos, -8, chrom_name, ha='center', va='top', 
             fontsize=10, rotation=45, fontweight='bold')
    current_pos += chrom_size

plt.colorbar(im1, ax=ax1, label='log10(contacts + 1)', fraction=0.046)

# Panel 2: Reconstructed full matrix
ax2 = fig.add_subplot(gs[0, 2:4])
im2 = ax2.imshow(log_reconstructed, cmap='YlOrRd', aspect='auto',
                 norm=colors.Normalize(vmin=0, vmax=np.percentile(log_reconstructed, 99)))
ax2.set_title('Reconstructed from Simulated Reads', fontsize=13, fontweight='bold')
ax2.set_xlabel('Genomic Position (bins)', fontsize=11)
ax2.set_ylabel('Genomic Position (bins)', fontsize=11)

current_pos = 0
for chrom_name, chrom_size in chrom_info:
    ax2.axhline(y=current_pos, color='blue', linewidth=1.5, alpha=0.6)
    ax2.axvline(x=current_pos, color='blue', linewidth=1.5, alpha=0.6)
    mid_pos = current_pos + chrom_size / 2
    ax2.text(mid_pos, -8, chrom_name, ha='center', va='top', 
             fontsize=10, rotation=45, fontweight='bold')
    current_pos += chrom_size

plt.colorbar(im2, ax=ax2, label='log10(contacts + 1)', fraction=0.046)

# Panel 3: Contact decay profile comparison
ax3 = fig.add_subplot(gs[1, 0:2])
max_distance = min(80, original.shape[0] // 3)

orig_distances = []
orig_contacts = []
recon_distances = []
recon_contacts = []

for d in range(1, max_distance):
    orig_diag = []
    recon_diag = []
    for i in range(original.shape[0] - d):
        if original[i, i+d] > 0:
            orig_diag.append(original[i, i+d])
        if reconstructed[i, i+d] > 0:
            recon_diag.append(reconstructed[i, i+d])
    
    if orig_diag:
        orig_distances.append(d)
        orig_contacts.append(np.mean(orig_diag))
    if recon_diag:
        recon_distances.append(d)
        recon_contacts.append(np.mean(recon_diag))

ax3.semilogy(orig_distances, orig_contacts, 'b-', linewidth=2.5, label='Original', alpha=0.8)
ax3.semilogy(recon_distances, recon_contacts, 'r--', linewidth=2.5, label='Reconstructed', alpha=0.8)
ax3.set_xlabel('Genomic Distance (bins × 10kb)', fontsize=11)
ax3.set_ylabel('Average Contacts (log scale)', fontsize=11)
ax3.set_title('Contact Decay Profile', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3, linestyle='--')

# Panel 4: Cis vs Trans comparison
ax4 = fig.add_subplot(gs[1, 2:4])
chrom_names = [c[0] for c in chrom_info]
x_pos = np.arange(len(chrom_names))
width = 0.35

orig_cis = []
orig_trans = []
recon_cis = []
recon_trans = []

current_bin = 0
for chrom_name, chrom_size in chrom_info:
    chr_start = current_bin
    chr_end = current_bin + chrom_size
    
    # Original matrix
    orig_cis_sum = np.sum(original[chr_start:chr_end, chr_start:chr_end])
    orig_trans_sum = np.sum(original[chr_start:chr_end, :]) - orig_cis_sum
    orig_cis.append(orig_cis_sum)
    orig_trans.append(orig_trans_sum)
    
    # Reconstructed matrix
    recon_cis_sum = np.sum(reconstructed[chr_start:chr_end, chr_start:chr_end])
    recon_trans_sum = np.sum(reconstructed[chr_start:chr_end, :]) - recon_cis_sum
    recon_cis.append(recon_cis_sum)
    recon_trans.append(recon_trans_sum)
    
    current_bin += chrom_size

bars1 = ax4.bar(x_pos - width/2, orig_cis, width, label='Original Cis', 
                color='steelblue', alpha=0.8)
bars2 = ax4.bar(x_pos + width/2, recon_cis, width, label='Recon. Cis', 
                color='lightblue', alpha=0.8, edgecolor='black', linewidth=1)
bars3 = ax4.bar(x_pos - width/2, orig_trans, width, bottom=orig_cis, 
                label='Original Trans', color='coral', alpha=0.8)
bars4 = ax4.bar(x_pos + width/2, recon_trans, width, bottom=recon_cis, 
                label='Recon. Trans', color='lightsalmon', alpha=0.8, 
                edgecolor='black', linewidth=1)

ax4.set_xlabel('Chromosome', fontsize=11)
ax4.set_ylabel('Total Contacts', fontsize=11)
ax4.set_title('Cis/Trans Contact Distribution', fontsize=12, fontweight='bold')
ax4.set_xticks(x_pos)
ax4.set_xticklabels(chrom_names, fontsize=10)
ax4.legend(fontsize=9, loc='upper right')
ax4.set_yscale('log')
ax4.grid(True, alpha=0.3, axis='y', linestyle='--')

# Panels 5-7: Individual chromosome comparisons
for idx, (chrom_name, chrom_size) in enumerate(chrom_info):
    ax = fig.add_subplot(gs[2, idx])
    
    chr_start = sum(s for _, s in chrom_info[:idx])
    orig_chr = original[chr_start:chr_start+chrom_size, chr_start:chr_start+chrom_size]
    recon_chr = reconstructed[chr_start:chr_start+chrom_size, chr_start:chr_start+chrom_size]
    
    # Show reconstructed (simulated alignment result)
    log_chr = np.log10(recon_chr + 1)
    im = ax.imshow(log_chr, cmap='YlOrRd', aspect='auto',
                   norm=colors.Normalize(vmin=0, vmax=np.percentile(log_chr, 99)))
    ax.set_title(f'{chrom_name} (from simulated reads)', fontsize=11, fontweight='bold')
    ax.set_xlabel('Position (10kb bins)', fontsize=9)
    ax.set_ylabel('Position (10kb bins)', fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

plt.savefig("{workdir}/heatmap.png", dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Main heatmap saved to {workdir}/heatmap.png")

# Create a detailed comparison figure for chr1
fig2, axes = plt.subplots(1, 3, figsize=(18, 6))

chr1_size = chrom_info[0][1]
orig_chr1 = original[0:chr1_size, 0:chr1_size]
recon_chr1 = reconstructed[0:chr1_size, 0:chr1_size]

# Original chr1
log_orig = np.log10(orig_chr1 + 1)
im1 = axes[0].imshow(log_orig, cmap='YlOrRd', aspect='auto')
axes[0].set_title('Original Input Matrix (chr1)', fontsize=13, fontweight='bold')
axes[0].set_xlabel('Genomic Position (10kb bins)', fontsize=11)
axes[0].set_ylabel('Genomic Position (10kb bins)', fontsize=11)
plt.colorbar(im1, ax=axes[0], label='log10(contacts + 1)')

# Reconstructed chr1
log_recon = np.log10(recon_chr1 + 1)
im2 = axes[1].imshow(log_recon, cmap='YlOrRd', aspect='auto')
axes[1].set_title('Reconstructed from Simulated Reads (chr1)', fontsize=13, fontweight='bold')
axes[1].set_xlabel('Genomic Position (10kb bins)', fontsize=11)
axes[1].set_ylabel('Genomic Position (10kb bins)', fontsize=11)
plt.colorbar(im2, ax=axes[1], label='log10(contacts + 1)')

# Difference map
diff = recon_chr1 - orig_chr1 * (np.sum(recon_chr1) / np.sum(orig_chr1))
max_diff = max(abs(np.percentile(diff, 1)), abs(np.percentile(diff, 99)))
im3 = axes[2].imshow(diff, cmap='RdBu_r', aspect='auto',
                     norm=colors.Normalize(vmin=-max_diff, vmax=max_diff))
axes[2].set_title('Difference (Reconstructed - Scaled Original)', fontsize=13, fontweight='bold')
axes[2].set_xlabel('Genomic Position (10kb bins)', fontsize=11)
axes[2].set_ylabel('Genomic Position (10kb bins)', fontsize=11)
plt.colorbar(im3, ax=axes[2], label='Contact difference')

plt.tight_layout()
plt.savefig("{workdir}/chr1_detailed_comparison.png", dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Detailed chr1 comparison saved to {workdir}/chr1_detailed_comparison.png")

# Print summary statistics
print("\\n" + "="*60)
print("SIMULATION SUMMARY")
print("="*60)
print(f"Reference genome: 3 chromosomes")
for name, size in chrom_info:
    print(f"  {{name}}: {{size*10}} kb ({{size}} bins)")
print(f"\\nSimulation parameters:")
print(f"  Read pairs simulated: 50,000")
print(f"  Bin size: 10 kb")
print(f"  Enzyme: HindIII (AAGCTT)")
print(f"  Random seed: 42")
print(f"\\nContact matrix properties:")
print(f"  Original total contacts: {{np.sum(original):.0f}}")
print(f"  Reconstructed total contacts: {{np.sum(reconstructed):.0f}}")
print(f"  Capture efficiency: {{100*np.sum(reconstructed)/np.sum(original):.1f}}%")
print(f"  Original sparsity: {{100*(1 - np.count_nonzero(original)/original.size):.1f}}%")
print(f"  Reconstructed sparsity: {{100*(1 - np.count_nonzero(reconstructed)/reconstructed.size):.1f}}%")
print(f"\\nVisualization files:")
print(f"  ✓ {{workdir}}/heatmap.png")
print(f"  ✓ {{workdir}}/chr1_detailed_comparison.png")
print("="*60)
'''
    
    with open(f"{workdir}/visualize.py", 'w') as f:
        f.write(visualize_script)
    
    success = run_command(f"python3 {workdir}/visualize.py", "Creating visualizations")
    
    if success:
        print("\n" + "="*60)
        print("Pipeline Complete! 🎉")
        print("="*60)
        print("\nGenerated files:")
        print(f"  📄 Reference genome:      {ref_fa}")
        print(f"  📊 Input contact matrix:  {matrix_tsv}")
        print(f"  📋 Offset file:           {offset_tsv}")
        print(f"  🧬 Simulated R1 reads:    {output_prefix}_R1.fastq")
        print(f"  🧬 Simulated R2 reads:    {output_prefix}_R2.fastq")
        print(f"  📈 Full heatmap:          {workdir}/heatmap.png")
        print(f"  📈 Chr1 comparison:       {workdir}/chr1_detailed_comparison.png")
        print("\nTo view the heatmaps:")
        print(f"  display {workdir}/heatmap.png")
        print(f"  # or open with any image viewer")
        print("="*60)
    else:
        print("\n⚠ Visualization step had issues, but simulation completed successfully.")
        print("You can still examine the FASTQ files and matrices.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
