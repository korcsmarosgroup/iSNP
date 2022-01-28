import sys
import os
import argparse
import subprocess


def parse_args(argv):
    help_text = \
        """
        === Generate Background Input Files ===
        This script creates the file inputs required for the background file creation
        for the transcription factor binding sites predictions.
        """
    parser = argparse.ArgumentParser(description=help_text)
    parser.add_argument("--regulatory-regions-file",
                        help="<bed file of the regulatory regions> [mandatory]",
                        type=str,
                        dest="regulatory_regions_file",
                        action="store",
                        required=True)
    parser.add_argument("--genome-folder",
                        help="<directory of fasta files holding each chromosome of the genome> [mandatory]",
                        type=str,
                        dest="genome_folder",
                        action="store",
                        required=True)
    parser.add_argument("--background-fasta",
                        help="<name of the background fasta file to be created.> [mandatory]",
                        type=str,
                        dest="background_input_fasta",
                        action="store",
                        required=True)
    parser.add_argument("--output",
                        help="<path to the output folder> [mandatory]",
                        dest="output_folder",
                        action="store",
                        required=True)
    parser.add_argument("--seq-length",
                        help="<Filter selection parameter for sequence lengths. Default: 60.> [optional]",
                        dest="seq_length",
                        type=int,
                        action="store",
                        required=False)
    parser.add_argument("--strand-type",
                        help="<Create background for single or double strands. Default: double.> [optional]",
                        dest="strand_type",
                        action="store",
                        default="double",
                        required=False)
    parser.add_argument("--order",
                        help="<Order for the background models. Default: 2.> [optional]",
                        dest="order",
                        action="store",
                        type=int,
                        default=1,
                        required=False)
    parser.add_argument("--model-name",
                        help="<Name of the prefix of the fimo and rsat model file. Default: None.> [optional]",
                        dest="model_name",
                        action="store",
                        default=None,
                        required=False)
    return parser.parse_args(argv[1:])


def get_n_lines(file):
    """ Get number of lines in the file """
    return sum(1 for _ in open(file))


def _create_genome_dict(genome_folder):
    genome_dictionary = {}
    for filepath in os.listdir(genome_folder):
        chr_file_input = os.path.join(genome_folder, filepath)
        chr_in_dict = filepath.split(".")[0]
        with open(chr_file_input, 'r') as c_file:
            for line in c_file:
                if chr_in_dict not in genome_dictionary:
                    genome_dictionary[chr_in_dict] = line
    return genome_dictionary


def generate_input_fasta(regulatory_regions_file, background_input_fasta, genome_folder):
    """ Create input fasta file  """
    genome_dictionary = _create_genome_dict(genome_folder=genome_folder)
    all_lines = get_n_lines(regulatory_regions_file)
    with open(regulatory_regions_file, 'r') as regulatory_regions, \
            open(background_input_fasta, 'w') as out_fasta:
        index = 0
        for line in regulatory_regions:
            percent = int(index / all_lines * 100)
            sys.stdout.write(f'\r{index}/{all_lines} - {percent}% done!')
            line = line.strip().split('\t')
            chromosome = line[0]
            start_p = int(line[1])
            end_p = int(line[2])
            gene = line[3].split(";")[3]
            if start_p > end_p:
                start_p = int(line[2])
                end_p = int(line[1])
            if chromosome in genome_dictionary:
                start_pos = int(start_p - 1)
                regulatory_region_sequence = genome_dictionary[chromosome][start_pos:end_p]
                out_fasta.write(f'>{gene}' + '\t' + f'size: {len(regulatory_region_sequence)}' + '\n')
                out_fasta.write(regulatory_region_sequence + '\n')
            index += 1
        sys.stdout.write(f'\r{all_lines}/{all_lines} - 100% done!')


def filter_input_fasta(input_file, output_file, n=60):
    """ Filter the input fasta file to split the sequences by the n lines """
    all_lines = get_n_lines(input_file)
    with open(input_file, 'r') as i, \
            open(output_file, 'w') as out:
        index = 0
        for sequence in i:
            percent = int(index / all_lines * 100)
            sys.stdout.write(f'\r{index}/{all_lines} - {percent}% done!')
            sequence = sequence.strip()
            if sequence.startswith(">"):
                out.write(sequence + '\n')
            if not sequence.startswith(">"):
                split_sequence = [sequence[i:i + n] for i in range(0, len(sequence), n)]
                out.write('\n'.join(split_sequence) + '\n')
            index += 1
        sys.stdout.write(f'\r{all_lines}/{all_lines} - 100% done!')


def generate_fimo_bg_file(input_file_path, ouptut_file_path, order=2):
    """ Run FIMO background file creat """
    cmd = ["sh", "generate_fimo_bg_file.sh",
           "-i", input_file_path,
           "-o", ouptut_file_path,
           "-m", order]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()


def generate_rsat_bg_file(input_file_path, output_file_path, strand_type):
    """ Run RSAT background file creation """
    cmd = ["sh", "generate_rsat_bg_file.sh",
           "-i", input_file_path,
           "-o", output_file_path,
           "-s", strand_type]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()


def generate_background(argv):
    """ Main logic to generate background files """
    args = parse_args(argv)

    # Create fasta files
    sys.stdout.write(f"Generating input fasta file...")
    generate_input_fasta(regulatory_regions_file=args.regulatory_regions_file,
                         background_input_fasta=args.background_input_fasta,
                         genome_folder=args.output_folder)

    # Reformat the created fasta files
    sys.stdout.write(f"Reformatting fasta file...")
    filter_input_fasta(input_file=args.input_file,
                       output_file=args.output_folder,
                       n=args.seq_length)

    # Generate the fimo background model
    fimo_bg_file_name = f"{args.model_name}_background_fimo.model"
    fimo_bg_file_path = os.path.join(args.output_folder, fimo_bg_file_name)
    sys.stdout.write(f"Generate FIMO background file: \t{fimo_bg_file_path}")
    generate_fimo_bg_file(input_file_path=args.input_file,
                          ouptut_file_path=fimo_bg_file_path,
                          order=2)

    # Generate the rsat background model
    rsat_bg_file_name = f"{args.model_name}_background_rsat.model"
    rsat_bg_file_path = os.path.join(args.output_folder, rsat_bg_file_name)
    sys.stdout.write(f"Generate RSAT background file: \t {rsat_bg_file_path}")
    generate_rsat_bg_file(input_file_path=args.input_file,
                          output_file_path=rsat_bg_file_path,
                          strand_type=args.strands)


if __name__ == "__main__":
    sys.exit(generate_background(sys.argv))
