import subprocess
from Bio import SeqIO

def run_blast(query_file, subject_file, output_file):
    blast_cmd = ["blastn", "-query", query_file, "-subject", subject_file, "-out", output_file, "-outfmt", "6"]
    subprocess.run(blast_cmd)

def parse_blast_output(output_file):
    with open(output_file, "r") as file:
        lines = file.readlines()
        alignments = [line.split() for line in lines]
        return alignments

def main():
    # File paths
    fasta_file_path_reads = "original.fasta"
    contigs_file_path = "./output/contigs.fa"  # Replace with the actual path to your contigs file
    blast_output_file = "blast_output.txt"

    # Run BLAST
    run_blast(fasta_file_path_reads, contigs_file_path, blast_output_file)

    # Parse BLAST output
    alignments = parse_blast_output(blast_output_file)

    # Display results
    print("Top BLAST alignments:")
    for alignment in alignments:
        query_id, subject_id, percent_identity, alignment_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_value, bit_score = alignment
        print(f"Query: {query_id}, Subject: {subject_id}, Percent Identity: {percent_identity}%")

if __name__ == "__main__":
    main()


