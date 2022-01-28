import subprocess
import os
import sys
from Bio import SeqIO


def _create_bedtools_command(path_to_fasta, path_to_bed, out_path):
    """
    https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html

    Parameters
    ----------
    - path_to_fasta: path to fasta file with the genomic secuence.
    - path_to_bed: path to the bed file defining the "cuts".
    - out_path: path where the output FASTA file is written.

    Returns
    -------
    bed_tools_command: call for bedtools getfasta
    """
    return ['bedtools', 'getfasta', '-fi', path_to_fasta, '-bed', path_to_bed, '-name', '-fo', out_path]


def extract_genome(path_to_fasta, path_to_bed, out_path):
    """
    Extracts the nucleotide secuence according to a bed file using the getfasta function from  bedtools.
    It writes a new FASTA file with all the secuences, named according to the "Name" field in the .bed file.
    The stderr of bedtools is sent to stdout.

    Parameters
    ----------
    path_to_fasta: path to fasta file with the genomic secuence.
    path_to_bed: path to the bed file defining the "cuts".
    out_path: path where the output FASTA file is written.
    """
    if not os.path.exists(path_to_fasta):
        sys.stderr.write("Could not find fasta file: " + path_to_fasta)
        sys.exit(202)

    if not os.path.exists(path_to_fasta):
        sys.stderr.write("Could not find bed file: " + path_to_bed)
        sys.exit(203)

    bed_tools_command = _create_bedtools_command(path_to_fasta,
                                                 path_to_bed,
                                                 out_path)
    p = subprocess.Popen(bed_tools_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    my_stdout, my_stderr = p.communicate()
    print(my_stderr.decode("utf-8"))


def transform_the_wild_type_fasta(original_output_wild_fasta, input_vcf, output_wild_fasta):
    """
    Transform wild type fasta files into vcf files.

    Parameters
    ----------
    original_output_wild_fasta
    input_vcf
    output_wild_fasta
    """
    n_base_dictionary = {}
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            line = line.strip()
            if line[:1] != "#":
                line = line.strip().split('\t')
                alt_base = line[4]
                if alt_base.upper() == "N":
                    n_base_dictionary[line[2]] = line[13].split(":")[1]

    with open(output_wild_fasta, 'w') as out_fasta:
        for seq in SeqIO.parse(original_output_wild_fasta, 'fasta'):
            desc = seq.description
            snp_id = desc.split("|")[1].split(":")[1].strip()
            conn = seq.id
            gene_id = conn.split(":")[1]
            if snp_id not in n_base_dictionary:
                SeqIO.write(seq, out_fasta, "fasta-2line")
            else:
                if gene_id not in n_base_dictionary[snp_id]:
                    SeqIO.write(seq, out_fasta, "fasta-2line")
