""" miRNA interaction predictor """
import argparse
import gzip
import os
import random
import string
import subprocess
import sys
import tempfile
import time
from collections import namedtuple

from Bio import SeqIO

sys.path.append("/rds/general/user/jno25/home/iSNP/analytic-modules")
from common_libs.mitab_handler import mitab_handler

# Keys for named-tuple to hold the scan results
results_keys = "Seq1, Seq2, Max_Score, Max_Energy"


class InvalidMirandaParameter(Exception):
    pass


def parse_args(argv=None):
    """ Commandline argument parser for the mirna interaction predictor module """

    help_text = \
        """
        === miRNA Interaction Predictor ===
        Name of the tool: miRNA Interaction Predictor.
        \nDescription:
        Using the miranda prediction software, where the 'file1' is a FASTA file with a microRNA query
        'file2' is a FASTA file containing the sequence(s).
        \nParameters:
        --mirna <path to an existing file> [mandatory]
        --genomic <path to the new output file> [mandatory]
        --output [mandatory]
        --score <score threshold, integer> [Optional]
        --energy <energy threshold [kcal/mol], negative integer> [Optional]
        --strict (Demand strict 5' seed pairing (default: less strict, seed region 2-8)) [Optional]
        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    # Input file path
    parser.add_argument("-m", "--mirna",
                        help="<path to an existing mirna fasta file> [mandatory]",
                        dest="mirna",
                        action="store",
                        required=True)

    # Genome data
    parser.add_argument("-g", "--genomic",
                        help="<path to an existing file> [mandatory]",
                        dest="genomic",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-o", "--output",
                        help="<path to an new file> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    # Threshold for prediction
    parser.add_argument("-sc", "--score",
                        help="<scoring threshold, integer> [Optional]",
                        type=int,
                        dest="score",
                        action="store",
                        default=False,
                        required=False)

    # Threshold for prediction
    parser.add_argument("-e", "--energy",
                        help="<energy threshold [kcal/mol], negative integer> [Optional]",
                        type=int,
                        dest="energy",
                        action="store",
                        default=False,
                        required=False)

    parser.add_argument("-s", "--strict",
                        help="<Demand strict 5' seed pairing (default: less strict, seed region 2-8)> [Optional]",
                        dest="strict",
                        action="store_const",
                        const=True,
                        default=False,
                        required=False)

    results = parser.parse_args(argv)

    return results


def _check_threshold(score, energy):
    """
    Checks that the threshold is a non-negative value if energy is being used, and
    if alignment score is being used then the value is a positive integer.

    Parameters
    ----------
    threshold: integer, the defined minimum total score/energy for filtering
    score_type: string, either 'energy' or 'alignment' for miranda tool scoring threshold

    """

    if energy:
        if not isinstance(energy, int) or energy > 0:
            raise InvalidMirandaParameter('Energy must be less than 0.')

    if score:
        if not isinstance(score, int) or score < 0:
            raise InvalidMirandaParameter('Alignment must be greater than 0.')


def check_fasta(file):
    """
    Check that the fasta file exists and that it's the correct file format. If the any of
    these statements fail, exit the program. Exit Code: 1.

    Parameters
    ----------
    file: string, file path to the fasta file

    """

    if not file.endswith('.fasta'):
        sys.exit(1)
    elif not file.endswith('.fa'):
        sys.exit(1)

    if not os.path.isfile(file):
        sys.exit(1)


def check_output(path):
    """ Ensure the output path is an actual file path, if not create the directory """
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except IsADirectoryError:
            raise IsADirectoryError('Could not make directory: {path}.')


def parse_database(database, verbose=False):
    """
    Convert the mirna database to a fasta file such that it is ready for the miranda
    tool can begin to predict interactions.

    Parameters
    ----------
    verbose: boolean, whether there is additional data (EMBL) format given
    database: string, file paths to the mirna database file

    Returns
    -------
    new_file: str, path to clean fasta file removing any non-human mirna from the
        database fasta file.
    details: dict, a dictionary holding the sequence info extracted from the
        fasta file

    """

    details = {}

    if database.endswith('fasta') or database.endswith('fa'):
        fasta = SeqIO.parse(database, 'fasta')

        cleaned = []

        for record in fasta:
            if record.id.__contains__('hsa'):
                cleaned.append(record)
                info = record.description.split()
                details[info[0]] = [info[1], info[4]]

        new_file = os.path.join(
            tempfile.gettempdir(),
            "mirna_db_%s.fasta" % ''.join(random.choices(string.ascii_lowercase, k=16))
        )

        SeqIO.write(cleaned, new_file, "fasta")
        os.sync()
        time.sleep(1)

        return new_file, details
    else:
        print(f'Converting database to fasta file...')

        new_file = os.path.splitext(database)[0].rstrip('.dat')
        new_file = f"{new_file}.fasta"

        with (gzip.open(database, 'rt') if database.endswith('.gz') else open(database)) as input_handle:
            with open(new_file, "w") as output_handle:
                embl = SeqIO.parse(input_handle, "embl")

                if verbose:
                    detailed_fasta = _parse_database(database)
                else:
                    for record in embl:
                        if record.id.__contains__('>hsa'):
                            details.append(record)

                SeqIO.write(embl, output_handle, "fasta")
                return new_file, details


def parse_sequences(sequences):
    """
    A function to parse, validate and extract information from the patient fasta files
    before mirna interactions are predicted.

    Parameters
    ----------
    sequences: str, a fasta file holding the sequences found in the patient

    Returns
    -------
    new_file: str, path to clean fasta file (remove any seqeunces that do not
        conform to the alphabet
    seq_info: dict, a dictionary holding the sequence info extracted from the
        fasta file

    """
    alphabet = 'ACTGU'
    cleaned = []
    seq_info = {}

    for seq in SeqIO.parse(sequences, 'fasta'):
        if all(i in alphabet for i in seq.seq.upper()):
            cleaned.append(seq)
            desc = seq.description.replace("|", "")
            desc = desc.split()
            seq_info[desc[0]] = desc[1:]

    new_file = os.path.join(
        tempfile.gettempdir(),
        "seq_%s.fasta" % ''.join(random.choices(string.ascii_lowercase, k=16))
    )

    SeqIO.write(cleaned, new_file, "fasta")
    os.sync()
    time.sleep(1)

    return new_file, seq_info


def _parse_database(file, spices='hsa'):
    """
    Parse the miRNA database from mirBase EMBL format. This is when the entire mirBase is given to the
    to the mirnda predictor. Here it will extract the data from the .dat file all the data possible.

    EMBL format: https://www.genomatix.de/online_help/help/sequence_formats.html

    Parameters
    ----------
    file: string, path to the full mirna database file
    spices: string, make sure the right tax id is being filtered. Default 'hsa'.

    Results
    -------
    results: list, a list of dicts holding information extracted from the full mirna database file
        which can then be used to create the new mitab file after mirna prediction.

    """

    mirna = []

    mature_mirna_id = 'FT   miRNA    '
    database_id = 'DR'

    for line in file:
        if line.startswith('ID') and line.__contains__(spices):
            setup_id = []
            mirna.append(setup_id)
        mirna[-1].append(line)

    results = []
    for group in mirna:
        name = group[0][5:23].strip()
        identifier = group[2][3:-2].strip()
        description = group[4][3:-1].strip()

        extracted_data = {
            'name': name,
            'description': description,
            'identifier': identifier
        }

        mature_mirna_lines = []
        database_lines = []

        for i, element in enumerate(group):
            if mature_mirna_id in element:
                mature_mirna_lines.append(i)
            if database_id in element:
                database_lines.append(i)

        extracted_data['products'] = [
            {
                'location': group[index][10:-1].strip(),
                'accession': group[index + 1][33:-2],
                'product': group[index + 2][31:-2],
            }
            for index in mature_mirna_lines
        ]

        for index in database_lines:
            dr_full = group[index].split()

            extracted_data['xrefs'] = [
                {
                    'database': dr_full[1].strip(';'),
                    'identifier': dr_full[2].strip(';'),
                }
            ]

        results.append(extracted_data)

    return results


def _predictor(sequences, database, score, energy, strict):
    """
    A private function to call the miranda mirna prediction tool.

    Parameters
    ----------
    sequences: str, file path to the patient sequences
    database: str, file path to the mirna database
    score: int, threshold for the scoring metric
    energy: int, threshold for the engery metric
    strict: str, strict parameter definition

    Returns
    -------
    my_stdout: str, string returned by the mirnda tool containing the mirna
        prediction results and alignments

    """

    miranda_args = ["miranda", database, sequences]

    if score:
        miranda_args.extend(['-sc', str(score)])
    if energy:
        miranda_args.extend(['-en', str(energy)])
    if strict:
        miranda_args.extend(['-strict'])

    temp_file_prefix = "%s.log" % ''.join(random.choices(string.ascii_lowercase, k=16))
    stdout_temp_file = os.path.join(tempfile.gettempdir(), 'stdout_' + temp_file_prefix)
    stderr_temp_file = os.path.join(tempfile.gettempdir(), 'stderr_' + temp_file_prefix)

    with open(stdout_temp_file, 'w+') as fout:
        with open(stderr_temp_file, 'w+') as ferr:
            child_process = subprocess.Popen(miranda_args, stderr=ferr, stdout=fout, universal_newlines=True)
            return_code = child_process.wait()

    if return_code != 0 or os.path.getsize(stderr_temp_file) != 0:
        print("ERROR - miranda args: " + str(miranda_args))
        print("ERROR - return code: " + str(return_code))
        print("ERROR - STDERR: ")
        with open(stderr_temp_file, 'r') as ferr:
            for line in ferr:
                print(line)
        # print("ERROR - STDOUT: " )
        # with open(stdout_temp_file, 'r') as fout:
        #     for line in fout:
        #         print(line)
        raise InvalidMirandaParameter()

    os.sync()
    time.sleep(1)

    mirna_connections = []
    with open(stdout_temp_file, 'r') as fout:
        for line in fout:
            if line.startswith(">>"):
                mirna_connections.append(line)

    os.remove(stdout_temp_file)
    os.remove(stderr_temp_file)
    os.sync()
    time.sleep(1)

    return mirna_connections


def extract_results(predictions):
    """
    A function to read the stdout of the miranda into a data-structure that can be used to
    create the NavigOmix output file.

    Parameters
    ----------
    predictions: str, the predictions returned by the miranda tool which need to extracted

    Returns
    -------
    mirna_pred: namedtuple, returns the top mirna prediction based on the parameters given to
        miranda tool. This are held in a namedtuple called scan where the sequence 1, sequence 2
        Max Score and Max Energy are all saved.

    """

    print(f'{len(predictions)} mirna sites found.')

    mirna_preds = []
    scan = namedtuple('scan', results_keys)

    for prediction_string in predictions:
        result = prediction_string.replace("\t", ",").strip('>>').split(",")

        mirna_preds.append(scan(Seq1=result[0],
                                Seq2=result[1],
                                Max_Score=result[4],
                                Max_Energy=result[5]))

    return mirna_preds


def create_network_file(mirna_preds, sequence_info, output):
    """
    A method to create the network file with the newly predicted interactions using the
    mitab handler from the common libs

    Parameters
    ----------
    mirna_preds: namedtuple, namedtuple (scan) holding the prediction results
    sequence_info: dict, dictionary holding any meta data about the interaction predicted
    output: str, file path location for the output mitab file

    """

    mitab = mitab_handler.MiTabHandler()
    inner_structure = {}

    for idx, mirna in enumerate(mirna_preds):

        interaction = mitab.new_interaction()

        # Clean and extract data
        mirna_interaction_score = f"score: {mirna.Max_Score}; energy: {mirna.Max_Energy}"
        mirna_interaction_score = mirna_interaction_score.rstrip(';')
        mirna_target = f'uniprotac:{mirna.Seq2.split(";")[2]}'

        # Add Interactor A and B
        interaction[mitab.uidA] = f'mirbase:{mirna.Seq1}'
        interaction[mitab.uidB] = f'{mirna_target}'
        interaction[mitab.taxA] = "taxid:9606('homo sapiens')"
        interaction[mitab.taxB] = "taxid:9906('homo sapiens')"

        # Add meta-data
        interaction[mitab.annotA] = f'start:micro rna;mirbase;{mirna.Seq1}'
        interaction[mitab.annotB] = f'end:{mirna.Seq2.split(":")[1]}'
        interaction[mitab.annotInter] = f'origin:snp;dbsnp;{sequence_info[mirna.Seq2][0].split(":")[1]}' \
                                        f' | {mirna_interaction_score} | {sequence_info[mirna.Seq2][1]}'

        inner_structure[idx] = interaction

    mitab.build_network_frame(inner_structure=inner_structure)

    # Write network file
    mitab.serialise_mitab(output, add_header=False)


def run(mirna, genomic, output, score, energy, strict):
    """
    Basic logic:
        (1) Use the miRNA sequences from mirBase (this will be an input parameter for the module)
        (2) Use tool miRanda to calculate the probability (entropy/strength) of the connection
        (3) Keep the miRNA gene connection if the probability is above a given threshold

    Parameters
    ----------
    mirna: str, file path to the fasta file holding the patient sequence
    genomic: str, file path to the fasta file holding mirna database file
    output: str, file path to the new mitab network file after mirna prediction
    score: int, the minimum alignment score allowed for an interaction
    energy: int, the negative minimum entropy of an interaction
    strict: str, defaults to 2'-8'

    """
    print(f"Starting Prediction")
    database_file_tmp, database_info = parse_database(genomic)
    genomic_file_tmp, sequences_info = parse_sequences(mirna)
    predictions = _predictor(genomic_file_tmp, database_file_tmp, score, energy, strict)
    mirna_preds = extract_results(predictions)
    create_network_file(mirna_preds, sequences_info, output)
    print(f"Finished!")

    os.remove(database_file_tmp)
    os.remove(genomic_file_tmp)
    os.sync()
    time.sleep(1)
    print(f"Temporary files deleted!")


def main(argv):
    """ Main function. Parses args and catches exit codes. """
    try:
        args = parse_args(argv)

        if os.stat(args.genomic).st_size == 0:

            print(f'====== The input fasta file is empty! ======')
            open(args.output, "a").close()

        else:

            _check_threshold(args.score, args.energy)
            run(args.mirna, args.genomic, args.output, args.score, args.energy, args.strict)

    except RuntimeError:
        sys.exit(2)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
