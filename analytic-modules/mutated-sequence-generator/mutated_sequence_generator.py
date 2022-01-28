import os
import sys
import argparse
from common_libs.filters import vcf_to_bed
from common_libs.filters import filter_fasta
from common_libs.mutated_sequence import get_mutated_sequence


def check_args(args=None):
    """ Parse arguments """

    help_text = \
        """
        VCF Filter

        Description:
        This tool will take a set of SNP data (in format of VCF), then by the location
        information in the VCF files it will extract an upstream and downstream region
        from the genome in a given length.

        The tool will generate two sequence files in FASTA format. One output file will
        contain all the wild type sequences, while the other file will contain the mutated
        sequences. In the FASTA file, each sequence will be tagged with the SNP id and also
        with the sequence type (wild type vs. mutated).

        If the genomic location can not be found in the reference genome, then the sequences
        will not be added to the output FASTA files (meaning there will be no "empty" sequence
        in the output).

        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("--input_vcf",
                        help="<path to the input VCF file)> [mandatory]",
                        dest="input",
                        action="store",
                        required=True)

    parser.add_argument("--genome",
                        help="<path to the genome (in single FASTA format))> [mandatory]",
                        dest="genome",
                        action="store",
                        required=True)

    parser.add_argument("--output_wild_type",
                        help="<path to the new output FASTA file with wild type sequences)> [mandatory]",
                        dest="output_wild_type",
                        action="store",
                        required=True)

    parser.add_argument("--output_mutated",
                        help="<path to the new output FASTA file with mutated sequences> [mandatory]",
                        dest="output_mutated",
                        action="store",
                        required=True)

    parser.add_argument("--region_length",
                        help="<0..1000: using e.g. 100 will fetch region: [X-100, x, x+100] with total length of 201> [mandatory]",
                        type=int,
                        dest="region_length",
                        action="store",
                        required=True)

    # Get arguments
    results = parser.parse_args(args)
    return results.input, results.genome, results.output_wild_type, results.output_mutated, results.region_length


def check_pars(input_vcf, genome, output_wild_type, output_mutated, region_length):
    """
    Check types for every argument for mutate.
    """

    if not isinstance(input_vcf, str):
        sys.stderr.write("input_vcf must be a str!")
        sys.exit(101)

    if not isinstance(genome, str):
        sys.stderr.write("genome must be a str!")
        sys.exit(102)

    if not isinstance(output_wild_type, str):
        sys.stderr.write("output_wild_type must be a str!")
        sys.exit(103)

    if not isinstance(output_mutated, str):
        sys.stderr.write("output_mutated must be a str!")
        sys.exit(104)

    if not isinstance(region_length, int):
        sys.stderr.write("region_length must be a str!")
        sys.exit(105)

    if region_length < 1:
        sys.stderr.write("region_length must be >= 1")
        sys.exit(106)


def mutate(input_vcf, genome, output_wild_type, output_mutated, region_length):
    """
    Mutate a sequence
    """
    # Check the arguments
    check_pars(input_vcf, genome, output_wild_type, output_mutated, region_length)

    # Convert the promoter .vcf to .bed
    promoter_VCF = vcf_to_bed.ProcessVcf(input_vcf, verbose=False)
    promoter_VCF.vcf2bed(output_wild_type + ".tmp.bed", region_length, chr_decorator="")

    # Extract the FASTA genome sequences using the promoter .bed
    filter_fasta.extract_genome(genome, output_wild_type + ".tmp.bed", output_wild_type + ".old.fasta")
    filter_fasta.transform_the_wild_type_fasta(output_wild_type + ".old.fasta", input_vcf, output_wild_type)
    os.remove(output_wild_type + ".tmp.bed")
    os.remove(output_wild_type + ".old.fasta")

    # Mutate the promoter sequence
    get_mutated_sequence.generate(output_wild_type, input_vcf, output_mutated, region_length)


if __name__ == "__main__":
    input_vcf, genome, output_wild_type, output_mutated, region_length = check_args()
    mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)
