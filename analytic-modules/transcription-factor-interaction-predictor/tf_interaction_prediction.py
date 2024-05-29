import subprocess
import sys
import argparse
import os
import shutil
sys.path.append("/home/liuy447/iSNP/analytic-modules/")
from rsat_prediction import scan_matrix
from rsat_prediction import write_rsat_results
from fimo_prediction import run_fimo


def parse_args(argv=None):
    help_text = \
        """
        === Transcription Factor Interaction Predictor ===
        Name of the tool: TF Interaction Predictor.
        \nDescription:

        Using the matrix-scan function from the rsat software prediction software.

        \nParameters:
        --path_to_fasta: Path to the FASTA file with the nucleotide sequences.
        --out_path: Path where the output mitab file is to be written.
        --path_to_matrix: Path to the file with transcription matrix binding profiles. By default, the JASPAR database for vertebrates (http://jaspar.genereg.net/downloads/).
        --format_matrix: Format of the file in path_to_matrix. transfac by default (recommended, as other types may cause an overheard).
        --pval_threshold: Only those interactions with a p-value lower than this value will be output.
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--fasta",
                        help="<path to an existing fasta file> [mandatory]",
                        dest="path_to_fasta",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output",
                        help="<output path> [mandatory]",
                        dest="out_path",
                        action="store",
                        required=True)

    parser.add_argument("-b", "--tf_background_rsat",
                        help="<background file> [mandatory]",
                        dest="tf_background_rsat",
                        action="store",
                        required=True)

    parser.add_argument("-bf", "--tf_background_fimo",
                        help="<background file> [mandatory]",
                        dest="tf_background_fimo",
                        action="store",
                        required=True)

    parser.add_argument("-m", "--matrix",
                        help="<path to matrix file> [Optional]",
                        dest="path_to_matrix",
                        action="store",
                        default=False,
                        required=False)

    parser.add_argument("-f", "--format",
                        help="<format of the matrix> [Optional]",
                        dest="format_matrix",
                        action="store",
                        default=False,
                        required=False)

    parser.add_argument("-t", "--threshold",
                        help="<threshold for probability of prediction (p-value)> [Optional]",
                        type=float,
                        dest="pval_threshold",
                        action="store",
                        default=False,
                        required=False)

    parser.add_argument("-pn", "--patient_folder",
                        help="<name of the patient to separate the rsat helper file> [mandatory]",
                        type=str,
                        dest="patient_folder",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)
    return results


def find_tf_sites(path_to_fasta,
                  out_path,
                  actual_patient_folder,
                  path_to_matrix=None,
                  format_matrix=None,
                  pval_threshold=None,
                  background=None):
    
    if "wt" in path_to_fasta:
        rsat_helper_file = f"{actual_patient_folder}_rsat_matrixscan_wt.txt"

    if "mut" in path_to_fasta:
        rsat_helper_file = f"{actual_patient_folder}_rsat_matrixscan_mut.txt"

    rsat_path = os.path.join(actual_patient_folder, rsat_helper_file)
    scan_matrix(path_to_fasta, rsat_path, path_to_matrix, format_matrix, background, actual_patient_folder)
    write_rsat_results(rsat_path, out_path, pval_threshold, actual_patient_folder)

    saving_command = ["arv", "keep", "put", "--project-uuid", "arkau-j7d0g-ch51898kwlrotjn", "--name", "Laurel_outputs", f"{rsat_helper_file}"]
    subprocess.run(saving_command, stderr = None, stdout = None)
    os.remove(rsat_path)


def main(argv):
    """
    Main function. Parses args and catches exit codes.
    :param argv:
    :return:
    """
    args = parse_args(argv)
    parse_fasta = args.path_to_fasta
    output_file = args.out_path
    actual_patient_folder = args.patient_folder

    if os.stat(parse_fasta).st_size == 0:

        print(f'====== The input fasta file is empty! ======')
        open(output_file, "a").close()

    else:
        new_fasta_file_path = f'{parse_fasta.split(".")[0]}_modified.fasta'
        with open(parse_fasta, 'r') as f, open(new_fasta_file_path, 'w') as new_f:

            for line in f:
                line = line.strip()

                if line[:1] == ">":
                    new_line = line.replace(" ", "")
                    new_f.write(new_line + '\n')

                else:
                    new_f.write(line + '\n')

        f.close()
        new_f.close()

        # Run FIMO
        fimo_output_folder = os.path.dirname(args.out_path)
        fimo_output_filename = os.path.join(fimo_output_folder, f"fimo_{os.path.basename(args.out_path)}")
        run_fimo(motif_file=args.path_to_matrix,
                 sequence_file=args.path_to_fasta,
                 output_folder=fimo_output_folder,
                 background_file=args.tf_background_fimo,
                 mitab_output_folder=fimo_output_filename)

        # Run RSAT
        find_tf_sites(new_fasta_file_path, args.out_path, actual_patient_folder, args.path_to_matrix, args.format_matrix, args.pval_threshold, args.tf_background_rsat)
        os.remove(new_fasta_file_path)

        # Merge the rsat and the fimo results
        new_rsat_file_name = f"rsat_{os.path.basename(args.out_path)}"
        new_rsat_file_path = os.path.join(os.path.dirname(args.out_path), new_rsat_file_name)
        shutil.copyfile(args.out_path, new_rsat_file_path)
        filenames = [fimo_output_filename, new_rsat_file_path]
        with open(args.out_path, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
