import os
import re
import sys
import csv
import json
import errno
import argparse
from pathlib import Path
from collections import defaultdict

ID_VERSION_SEP = "."
HEADER = ["source", "target", "score", "region", "snp", "file", "mutated", "tool"]
CONVERT_FILES_TF = ["tf_gene_connections_wt.tsv", "tf_gene_connections_mut.tsv"]
CONVERT_FILES_MIRNA = ["mirna_gene_connections_wt.tsv", "mirna_gene_connections_mut.tsv"]


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
    return args_parser.parse_args(argv[1:])


def _strip_dict(d):
    """ Strip white space from interaction line dictionary """
    return {k: v.strip() for k, v in d.items()}


def _load_mapping_dict(mapping_file_path, target_id_type):
    """ Load the sherlock mapping dictionary """
    mapping_dict = defaultdict(set)
    with open(mapping_file_path, "r") as mapping_file:
        for line in mapping_file:
            map_line = json.loads(line.strip())
            if map_line["to_id_type"].lower() == target_id_type:
                mapping_dict[map_line["from_id"].lower()] = map_line["to_id"].lower()
    return mapping_dict


def remap_ids(identifier, mapping_dict):
    """ A method to remap the ids """
    identifier = re.split("[;:]", identifier)
    current_id = identifier[-1]
    if ID_VERSION_SEP in current_id:
        current_id = current_id.split(".")[0]
    converted_id = mapping_dict.get(current_id)
    if converted_id is not None:
        return "".join(converted_id).upper()
    else:
        return converted_id


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
    with open(tf_file_path, "r") as tf_file, \
            open(tf_output_file_path, "w") as tf_output_file:
        reader = csv.reader(tf_file, delimiter="\t")
        writer = csv.DictWriter(tf_output_file, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for row in reader:
            row = [x for x in row if x != "-"]
            source = remap_ids(row[0], mapping_dict)
            target = remap_ids(row[1], mapping_dict)
            if source is not None and target is not None:
                if _check_if_complex(row[0]):
                    complex_source = row[0].split(";")[2]
                    complex_source = remap_ids(complex_source, mapping_dict)
                    source = f"{complex_source}/{source}"
                score = row[4]
                if "fimo" in score:
                    tool_name = "fimo"
                else:
                    tool_name = "rsat"
                region = row[6].split(";")[1]
                snp = row[7].split(";")[-1]
                mutated = mutated_dict[snp]
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


def handle_mirna(mirna_file_path, mapping_dict, output_path):
    """ A method to correct formatting issues with the miranda output """
    file_mut_type = _get_sequence_type(mirna_file_path)
    mirna_output_file_path = os.path.join(output_path, _generate_output_file_name(mirna_file_path))
    with open(mirna_file_path, "r") as tf_file, \
            open(mirna_output_file_path, "w") as mirna_output_file:
        reader = csv.reader(tf_file, delimiter="\t")
        writer = csv.DictWriter(mirna_output_file, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for row in reader:
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


def _write_network(network, output_file_path):
    """ A method to write out the resulting network """
    with open(output_file_path, "w") as interaction_network:
        writer = csv.DictWriter(interaction_network, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(network)


def get_multi_key(my_dict, keys):
    """ A helper function to return a dict of limited keys """
    return {my_dict.get(key) for key in keys}


def _network_difference(interaction_network, key, compare_on):
    """ A helper method to merge the networks """
    compare_on = ["source", "target", "snp", "tool"]
    keys = [k for k in list(interaction_network.keys()) if key in k]
    network_a = interaction_network[keys[0]]
    network_b = interaction_network[keys[1]]
    combined = network_a + network_b
    network_a = [get_multi_key(inter, compare_on) for inter in network_a]
    network_b = [get_multi_key(inter, compare_on) for inter in network_b]
    diff = [i for i in combined
            if get_multi_key(i, compare_on) not in network_a or
            get_multi_key(i, compare_on) not in network_b]
    result = len(diff) == 0
    if not result:
        print(f'There are {len(diff)} differences.')
    return diff


def network_differences(output_dir, compare_on):
    """ A method to combine the corrected datasets """
    patient_networks = defaultdict(dict)
    entries = Path(output_dir)
    for interaction_file in entries.iterdir():
        patient_name = os.path.basename(interaction_file.name)
        patient_name = patient_name.split(".")[0]
        if "gene_connections" in patient_name:
            with open(interaction_file) as interaction:
                patient_name = os.path.basename(interaction_file.name)
                patient_name = patient_name.split(".")[0]
                patient_name = patient_name.split("_")
                patient_name = f"{patient_name[1]}_{patient_name[4]}"
                reader = csv.DictReader(interaction, delimiter="\t")
                patient_networks[patient_name] = list(reader)
    mirna_differences = _network_difference(patient_networks, "mirna", compare_on)
    mirna_differences_file_path = os.path.join(output_dir, "mirna_differences.tsv")
    _write_network(mirna_differences, mirna_differences_file_path)
    tf_differences = _network_difference(patient_networks, "tf", compare_on)
    tf_differences_file_path = os.path.join(output_dir, "tf_differences.tsv")
    _write_network(tf_differences, tf_differences_file_path)
    differences = mirna_differences + tf_differences
    differences_file_path = os.path.join(output_dir, "differences.tsv")
    _write_network(differences, differences_file_path)
    return differences


def reformat(argv):
    """ Main logic for isnp results reformatter """
    args = parse_args(argv)
    mapping_dict = _load_mapping_dict(args.mapping_file_path, args.target_id)
    # mapping_dict = ""
    # entries = Path(args.input_dir)
    for patient in os.listdir(args.input_dir): # entries.iterdir():
        # if not patient.name.startswith("."):
        #     if patient.name.endswith(".vcf"):
        print(f"Converting: {patient}")
        patient_output_dir_path = os.path.join(args.output_dir, f"mapped_{patient}")
        patient_input_dir_path = os.path.join(args.input_dir, patient)
        try:
            os.mkdir(patient_output_dir_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        for file in CONVERT_FILES_TF:
            handle_tf(os.path.join(patient_input_dir_path, file), mapping_dict, patient_output_dir_path, patient_input_dir_path)
        for file in CONVERT_FILES_MIRNA:
            handle_mirna(os.path.join(patient, file), mapping_dict, patient_output_dir_path)
        network_differences(patient_output_dir_path, args.compare_on)


if __name__ == "__main__":
    sys.exit(reformat(sys.argv))
