import os
import re
import sys
import csv
import json
import errno
import logging
import argparse
from pathlib import Path
from collections import defaultdict
from logging.handlers import QueueHandler, QueueListener
from multiprocessing import Queue
from concurrent.futures import ProcessPoolExecutor, as_completed

ID_VERSION_SEP = "."
HEADER = ["source", "target", "score", "region", "snp", "file", "mutated", "tool"]
CONVERT_FILES_TF = ["tf_gene_connections_wt.tsv", "tf_gene_connections_mut.tsv"]
CONVERT_FILES_MIRNA = ["mirna_gene_connections_wt.tsv", "mirna_gene_connections_mut.tsv"]


def setup_logging(output_dir):
    """ Setup logging to both console and a file using a queue for process safety """
    # Create the queue for logging
    log_queue = Queue()
    
    log_file = "isnp_wrangler.log"
    
    # Formatter for all logs
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    
    # Listener to handle the queue in the main process
    listener = QueueListener(log_queue, console_handler, file_handler)
    listener.start()
    
    # Configure the root logger in the main process to use the queue
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    # Remove existing handlers to avoid double logging
    for hdlr in root.handlers[:]:
        root.removeHandler(hdlr)
    
    qh = QueueHandler(log_queue)
    root.addHandler(qh)
    
    return listener, log_queue


class PatientAdapter(logging.LoggerAdapter):
    """ Custom adapter to prefix logs with patient ID """
    def process(self, msg, kwargs):
        return f"[{self.extra['patient']}] {msg}", kwargs


def parse_args(argv):
    """ VCF file lift over commands default cli interface """
    help_text = \
        """
        === ISNP Output Files Wrangler ===
        This script generates a more user-friendly output from the results of the iSNP pipeline.
        """
    args_parser = argparse.ArgumentParser(description=help_text)
    args_parser.add_argument('-i', '--input_dir', action='store', dest="input_dir", required=True)
    args_parser.add_argument('-o', '--output_dir', action='store', dest="output_dir", required=True)
    args_parser.add_argument('-m', '--mapping_file_path', action='store', dest="mapping_file_path", required=True)
    args_parser.add_argument('-t', '--target_id', action='store', dest="target_id", default="uniprotac")
    args_parser.add_argument('-v', '--verbose', action='store', dest="verbose", type=int)
    args_parser.add_argument('-c', '--compare_on', action='store', dest="compare_on", nargs='+',
                             default=["source", "target", "snp", "tool"])
    args_parser.add_argument('-p', '--processes', action='store', dest="processes", type=int, default=os.cpu_count(),
                             help="Number of parallel processes to use")
    return args_parser.parse_args(argv[1:])


def _strip_dict(d):
    """ Strip white space from interaction line dictionary """
    return {k: v.strip() for k, v in d.items()}


def _load_mapping_dict(mapping_file_path, target_id_type):
    """ Load the sherlock mapping dictionary """
    mapping_dict = {}
    with open(mapping_file_path, "r") as mapping_file:
        for line in mapping_file:
            try:
                map_line = json.loads(line.strip())
                if map_line["to_id_type"].lower() == target_id_type:
                    mapping_dict[map_line["from_id"].lower()] = map_line["to_id"].lower()
            except json.JSONDecodeError:
                continue
    return mapping_dict


def remap_ids(identifier, mapping_dict):
    """ A method to remap the ids """
    parts = re.split("[;:]", identifier)
    current_id = parts[-1]
    if ID_VERSION_SEP in current_id:
        current_id = current_id.split(ID_VERSION_SEP)[0]
    
    converted_id = mapping_dict.get(current_id.lower())
    if converted_id:
        return converted_id.upper()
    return None


def _get_sequence_type(file_path):
    """ Get mut/wt from the file name """
    return re.split("[_.]", os.path.basename(file_path))[-2]


def _generate_output_file_name(file_path):
    """ Generate the output file path """
    return f"converted_{os.path.basename(file_path)}"


def _get_mutated(dir_path, patient_type, input_dir):
    """ Get if the SNP is mutated or not from the generated fasta file """
    mutated_dict = defaultdict(set)
    if patient_type == "wt":
        fasta_file_name = "snp_in_promoter-regions_wt.fasta"
    else:
        fasta_file_name = "snp_in_promoter-regions_mut.fasta"
    fasta_file_path = os.path.join(input_dir, fasta_file_name)
    with open(fasta_file_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                line = line.split("|")
                snp = line[1].split(":")[-1].strip()
                mutated_dict[snp] = line[2].split(":")[-1].strip()
    return mutated_dict


def _check_if_complex(protein_id):
    """ Check if the id is a protein complex or not """
    protein_id = protein_id.split(";")
    if len(protein_id) > 3:
        return True


def handle_tf(tf_file_path, mapping_dict, output_path, input_dir):
    """ A method to correct formatting issues with the fimo output """
    file_mut_type = _get_sequence_type(tf_file_path)
    tf_output_file_path = os.path.join(output_path, _generate_output_file_name(tf_file_path))
    mutated_dict = _get_mutated(input_dir, file_mut_type, input_dir)
    results = []
    with open(tf_file_path, "r") as tf_file, \
            open(tf_output_file_path, "w") as tf_output_file:
        reader = csv.reader(tf_file, delimiter="\t")
        writer = csv.DictWriter(tf_output_file, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if not row: continue
            row = [x for x in row if x != "-"]
            source = remap_ids(row[0], mapping_dict)
            target = remap_ids(row[1], mapping_dict)
            if source is not None and target is not None:
                if _check_if_complex(row[0]):
                    complex_source = row[0].split(";")[2]
                    complex_source = remap_ids(complex_source, mapping_dict)
                    if complex_source:
                        source = f"{complex_source}/{source}"
                score = row[4]
                tool_name = "fimo" if "fimo" in score else "rsat"
                region = row[6].split(";")[1]
                snp = row[7].split(";")[-1]
                mutated = mutated_dict.get(snp, "Unknown")
                new_interaction = {HEADER[0]: source,
                                   HEADER[1]: target,
                                   HEADER[2]: score,
                                   HEADER[3]: region,
                                   HEADER[4]: snp,
                                   HEADER[5]: file_mut_type,
                                   HEADER[6]: mutated,
                                   HEADER[7]: tool_name}
                new_interaction = _strip_dict(new_interaction)
                writer.writerow(new_interaction)
                results.append(new_interaction)
    return results


def handle_mirna(mirna_file_path, mapping_dict, output_path):
    """ A method to correct formatting issues with the miranda output """
    file_mut_type = _get_sequence_type(mirna_file_path)
    mirna_output_file_path = os.path.join(output_path, _generate_output_file_name(mirna_file_path))
    results = []
    with open(mirna_file_path, "r") as tf_file, \
            open(mirna_output_file_path, "w") as mirna_output_file:
        reader = csv.reader(tf_file, delimiter="\t")
        writer = csv.DictWriter(mirna_output_file, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if not row: continue
            row = [x for x in row if x != "-"]
            source = row[0].split(':')[1]
            target = remap_ids(row[1], mapping_dict)
            if source is not None and target is not None:
                meta_data = row[6].split("|")
                score = meta_data[1].replace(" ", "")
                region = "protein-coding"
                snp = meta_data[0].split(";")[-1]
                mutated = meta_data[2].split(":")[1].capitalize()
                new_interaction = {HEADER[0]: source,
                                   HEADER[1]: target,
                                   HEADER[2]: score,
                                   HEADER[3]: region,
                                   HEADER[4]: snp,
                                   HEADER[5]: file_mut_type,
                                   HEADER[6]: mutated,
                                   HEADER[7]: "miranda"}
                new_interaction = _strip_dict(new_interaction)
                writer.writerow(new_interaction)
                results.append(new_interaction)
    return results


def _write_network(network, output_file_path):
    """ A method to write out the resulting network """
    with open(output_file_path, "w") as interaction_network:
        writer = csv.DictWriter(interaction_network, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(network)


def get_multi_key(my_dict, keys):
    """ A helper function to return a tuple of limited keys for hashing """
    return tuple(my_dict.get(key) for key in keys)


def _network_difference(network_a, network_b, compare_on, logger=logging):
    """ A helper method to merge the networks efficiently using sets """
    if not network_a or not network_b:
        return network_a + network_b
    
    network_a_keys = {get_multi_key(inter, compare_on) for inter in network_a}
    network_b_keys = {get_multi_key(inter, compare_on) for inter in network_b}
    
    combined = network_a + network_b
    diff = [i for i in combined
            if get_multi_key(i, compare_on) not in network_a_keys or
            get_multi_key(i, compare_on) not in network_b_keys]
    
    if diff:
        logger.info(f'There are {len(diff)} differences.')
    return diff


def calculate_differences(tf_all, mirna_all, output_dir, compare_on, logger=logging):
    """ Calculate and write differences for TF and miRNA networks """
    mirna_diff = _network_difference(
        [r for r in mirna_all if r['file'] == 'wt'],
        [r for r in mirna_all if r['file'] == 'mut'],
        compare_on,
        logger=logger
    )
    _write_network(mirna_diff, os.path.join(output_dir, "mirna_differences.tsv"))

    tf_diff = _network_difference(
        [r for r in tf_all if r['file'] == 'wt'],
        [r for r in tf_all if r['file'] == 'mut'],
        compare_on,
        logger=logger
    )
    _write_network(tf_diff, os.path.join(output_dir, "tf_differences.tsv"))

    all_diff = mirna_diff + tf_diff
    _write_network(all_diff, os.path.join(output_dir, "differences.tsv"))
    return all_diff


def process_patient(patient, args, mapping_dict):
    """ Process a single patient's data """
    # Use the adapter to prefix all logs in this process
    logger = PatientAdapter(logging.getLogger(), {'patient': patient})
    logger.info("Starting conversion")
    
    patient_output_dir_path = os.path.join(args.output_dir, f"mapped_{patient}")
    patient_input_dir_path = os.path.join(args.input_dir, patient)
    
    try:
        os.makedirs(patient_output_dir_path, exist_ok=True)
    except OSError as e:
        logger.error(f"Failed to create directory {patient_output_dir_path}: {e}")
        return False

    tf_all = []
    for file in CONVERT_FILES_TF:
        file_path = os.path.join(patient_input_dir_path, file)
        if os.path.exists(file_path):
            logger.info(f"Processing TF file: {file}")
            tf_all.extend(handle_tf(file_path, mapping_dict, patient_output_dir_path, patient_input_dir_path))
        else:
            logger.warning(f"TF file not found: {file_path}")

    mirna_all = []
    for file in CONVERT_FILES_MIRNA:
        file_path = os.path.join(patient_input_dir_path, file)
        if os.path.exists(file_path):
            logger.info(f"Processing miRNA file: {file}")
            mirna_all.extend(handle_mirna(file_path, mapping_dict, patient_output_dir_path))
        else:
            logger.warning(f"miRNA file not found: {file_path}")
    
    logger.info("Calculating network differences")
    calculate_differences(tf_all, mirna_all, patient_output_dir_path, args.compare_on, logger=logger)
    logger.info("Finished successfully")
    return True


def reformat(argv):
    """ Main logic for isnp results reformatter with parallel processing """
    args = parse_args(argv)
    listener, log_queue = setup_logging(args.output_dir)
    logging.info(f"Starting iSNP Wrangler with {args.processes} processes")
    
    try:
        mapping_dict = _load_mapping_dict(args.mapping_file_path, args.target_id)
        patients = [p for p in os.listdir(args.input_dir) if os.path.isdir(os.path.join(args.input_dir, p))]
        
        if args.processes > 1 and len(patients) > 1:
            with ProcessPoolExecutor(max_workers=args.processes) as executor:
                futures = {executor.submit(process_patient, patient, args, mapping_dict): patient for patient in patients}
                for future in as_completed(futures):
                    patient = futures[future]
                    try:
                        future.result()
                    except Exception as e:
                        logging.error(f"Error processing patient {patient}: {e}")
        else:
            for patient in patients:
                process_patient(patient, args, mapping_dict)
    
    finally:
        logging.info("iSNP Wrangler finished")
        listener.stop()


if __name__ == "__main__":
    sys.exit(reformat(sys.argv))
