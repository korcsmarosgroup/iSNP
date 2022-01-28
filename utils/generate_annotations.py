import os
import sys
import argparse
from time import strftime
from collections import defaultdict


def parse_args(args):
    help_text = \
        """
        === Annotation files maker === 
        Generates promoter/enhancer and protein coding regions annotation files for the miRNA and Transcription factors 
        predictions. These annotation files generated from a gtf files which have the Genocde annotation format. 
        """
    parser = argparse.ArgumentParser(description=help_text)
    parser.add_argument("-i", "--input-gtf-file",
                        help="<paths to the input vcf files> [mandatory]",
                        type=str,
                        dest="input_gtf_file",
                        action="store",
                        required=True)
    parser.add_argument("-o", "--output-folder",
                        help="<path to the output folder> [mandatory]",
                        dest="output_folder",
                        action="store",
                        required=True)
    parser.add_argument("-e", "--enhancers-file",
                        help="<paths to the enhancer annotation files> [optional]",
                        type=str,
                        dest="enhancers_file",
                        action="store",
                        required=False)
    parser.add_argument("-v", "--verbose",
                        help="<verbosity level for debugging> [optional]",
                        dest="verbose",
                        action="store_true",
                        required=False)
    results = parser.parse_args(args)
    return results.input_gtf_file, results.output_folder, results.enhancers_file, results.verbose


def check_params(input_gtf_file):
    if not os.path.isfile(input_gtf_file):
        sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The specified input file does not exists: {input_gtf_file}')
        sys.exit(1)


def generate_annotations():
    input_gtf_file, output_folder, enhancers_file, verbose = parse_args(sys.argv[1:])
    check_params(input_gtf_file)
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Parameters are fine, starting...')
    protein_coding_regions = os.path.join(output_folder, "protein_coding_regions.bed")
    promoter_regions = os.path.join(output_folder, "promoter_regions.bed")
    genomic_strand = defaultdict()

    # Create the protein coding regions file
    with open(input_gtf_file, 'r') as gtf, \
            open(protein_coding_regions, 'w') as protein:
        if verbose:
            print(f'DEBUG [{strftime("%H:%M:%S")}]: Writing protein coding regions.')
        for line in gtf:
            if line.startswith("#"):
                continue
            line = line.strip().split('\t')
            informations = line[8].split("; ")
            # if 'gene_type "protein_coding"' in informations and line[2] == "gene":
            if 'gene_type "protein_coding"' in informations and line[2] == "CDS":
                chromosome = line[0]
                start = line[3]
                end = line[4]
                gene_id = informations[0].split('"')[1]
                gene_id = f"gene_id:gene;ensembl;{gene_id}"
                strand = line[6]
                protein.write(chromosome + '\t' + start + '\t' + end + '\t' + gene_id + '\n')
                genomic_strand[gene_id] = strand

    # Create the promoter regions file
    with open(protein_coding_regions, 'r') as protein_coding, \
            open(promoter_regions, 'w') as promoter:
        if verbose:
            print(f'DEBUG [{strftime("%H:%M:%S")}]: Writing promoter regions.')
        for line in protein_coding:
            line = line.strip().split('\t')
            chromosome = line[0]
            gene_id = line[3]
            if genomic_strand[gene_id] == "+" or "-":
                new_end = line[1]
                if verbose:
                    print(f'DEBUG [{strftime("%H:%M:%S")}]: {gene_id}:{genomic_strand[gene_id]} -> Subtracting 2000bp.')
                new_start = int(new_end) - 2000
            else:
                new_end = line[2]
                if verbose:
                    print(f'DEBUG [{strftime("%H:%M:%S")}]: {gene_id}:{genomic_strand[gene_id]} -> Adding 2000bp.')
                new_start = int(new_end) + 2000
            promoter.write(chromosome + '\t' + format(new_start) + '\t' + new_end + '\t' + gene_id + '\n')
        if enhancers_file:
            if verbose:
                print(f'DEBUG [{strftime("%H:%M:%S")}]: Enhancer files found. Writing enhancers.')
            with open(enhancers_file, 'r') as enhancers_annotations:
                for enhancer in enhancers_annotations:
                    promoter.write(enhancer)

    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Annotation finished successfully!')


if __name__ == "__main__":
    generate_annotations()
