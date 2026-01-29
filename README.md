# UMI-40mer-Processor
This script is a specialized bioinformatics tool designed for the high-throughput processing of Large-scale pooled lentiviral UMI libraries. It is specifically tailored to recover and validate 40-nt random sequences (UMIs/tags) from the specialized DNBSEQ 50-bp paired-end sequencing data.

Overview
The processing pipeline reconstructs a continuous 40-nt sequence from read pairs that contain specific 30-nt flanking motifs at positions 21–50. To ensure maximum sequence integrity, the script implements a strict quality filter, only collecting reconstructed 40-mers where all bases meet a minimum quality threshold (default Q30).

Scientific Logic
Library Structure: Two random 20-mers flanking a fixed 39-nt internal sequence.

Motif Validation:

R1 Motif: GGCATTCTGGATAGTGTCAAGGGAATAAAT (positions 21–50).

R2 Motif: GCAAATATCATTTATTCCCTTGACACTATC (positions 21–50).

Reconstruction: Concatenates the first 20-nt of Read 1 with the reverse-complement of the first 20-nt of Read 2.

Features
Memory Efficiency: Uses a spill-to-disk architecture and external merge-sort to maintain a bounded RAM footprint, allowing for the processing of exceptionally large datasets.

Multi-threaded Performance: Supports worker multiprocessing and utilizes pigz for high-speed decompression of .gz files.

Quality Splitting: Automatically separates counts into "Q-pass" (all bases ≥Q30) and "sub-Q" files.

Flexible Scanning: Allows for a user-defined shift (±N nt) in motif positions to account for sequencing indels or offsets.

Requirements
Python: 3.7+.

Dependencies: Standard Python libraries (argparse, concurrent.futures, gzip, etc.).

Optional: pigz for faster decompression performance.

Usage
Basic Command
Bash

python3 UMI_40mer_processing.py 'path/to/reads/BMR*_L*_1.fq.gz' \
    -o Output_Prefix \
    --fwd-motif GGCATTCTGGATAGTGTCAAGGGAATAAAT \
    --rev-motif GCAAATATCATTTATTCCCTTGACACTATC \
    --q-threshold 30 \
    --jobs 4
Key Arguments
r1_files: Positional argument for R1 files (supports globbing).

--fwd-motif / --rev-motif: The 30-nt sequences expected in R1 and R2.

--q-threshold: Phred score threshold for quality filtering (default: 30).

--flush-every-unique: Controls how often the script spills counts to disk (adjust based on available RAM).

Output Files
The script generates three primary outputs:

*_combined_counts.Q30.txt: Unique 40-nt sequences passing the quality threshold with their respective counts.

*_combined_counts.subQ30.txt: Unique 40-nt sequences that failed the quality threshold.

*_summary.txt: A detailed processing log including total read pairs, motif match rates, and unique sequence totals.
