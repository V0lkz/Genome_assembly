from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import subprocess

def generate_random_genome_sequence(length):
    bases = ['A', 'T', 'C', 'G']
    genome_sequence = ''.join(random.choice(bases) for _ in range(length))
    return Seq(genome_sequence)

def cut_genome_sequence(genome_sequence, num_reads, min_read_length, max_read_length, overlap):
    reads = []
    for _ in range(num_reads):
        read_length = random.randint(min_read_length, max_read_length)
        start_position = random.randint(0, len(genome_sequence) - read_length)
        end_position = start_position + read_length
        read = genome_sequence[start_position:end_position]
        reads.append(read)

    # Add overlap between consecutive reads
    for i in range(1, len(reads)):
        reads[i] = reads[i-1][-overlap:] + reads[i]

    return reads

def create_paired_reads(single_end_reads):
    paired_reads = []
    for i in range(0, len(single_end_reads), 2):
        if i + 1 < len(single_end_reads):
            read1 = single_end_reads[i]
            read2 = single_end_reads[i + 1]
            paired_reads.append((read1, read2))
    return paired_reads

# Specify parameters for cutting the genome sequence
num_reads = 400
min_read_length = 150
max_read_length = 200
overlap = 20

# Generate the random genome sequence
genome_length = 10000
random_genome_sequence = generate_random_genome_sequence(genome_length)

# Cut the genome sequence into reads
reads = cut_genome_sequence(random_genome_sequence, num_reads, min_read_length, max_read_length, overlap)

# Create SeqRecord objects for each read
read_records = [SeqRecord(Seq(read), id=f"Read_{i+1}", description="") for i, read in enumerate(reads)]

# Create a SeqRecord object for the original genome sequence
original_record = SeqRecord(random_genome_sequence, id="Original", description="")

# Write the reads to a FASTA file
fasta_file_path_reads = "reads.fasta"
SeqIO.write(read_records, fasta_file_path_reads, "fasta")
print(f"Reads written to {fasta_file_path_reads}")

# Write the original genome sequence to a FASTA file
fasta_file_path_original = "original.fasta"
SeqIO.write(original_record, fasta_file_path_original, "fasta")
print(f"Original sequence written to {fasta_file_path_original}")

# Generate single-end reads and paired-end reads
single_end_reads = cut_genome_sequence(random_genome_sequence, num_reads, min_read_length, max_read_length, overlap)
paired_end_reads = create_paired_reads(single_end_reads)

# Write paired-end reads to a FASTQ file
fastq_file_path_paired_end_reads = "paired_end_reads.fastq"
with open(fastq_file_path_paired_end_reads, "w") as fq_file:
    for i, (read1, read2) in enumerate(paired_end_reads):
        fq_file.write(f"@Pair_{i+1}_Read1\n{read1}\n+\n{'I' * len(read1)}\n")
        fq_file.write(f"@Pair_{i+1}_Read2\n{read2}\n+\n{'I' * len(read2)}\n")
print(f"Paired-end reads written to {fastq_file_path_paired_end_reads}")




