import argparse
import json
import os
import sys
sys.path.append("/rds/general/user/bbohar/home/projects/iSNP/analytic-modules")
from collections import defaultdict

import pandas as pd

from common_libs.mitab_handler import mitab_handler


def parse_args(args):
    help_text = \
        """
        === Network ID Mapper ===

        **Description:** 

        This tool will take a single network file as input, and mapping all the node identifiers to a
        specific target identifier type. Example: having an input protein interaction network with 
        nodes named with Ensembl Gene ID and Ensembl Protein ID, I can create a similar network where 
        all the node names are converted to Uniprot id.
        
        The ID mapping is not necessarily unambiguous. E.g. if there is no valid mapping, or there
        are more than one valid mappings for the given molecule, then the tool can skip or return 
        with multiple uniprot interactions for a single original interaction.
        
        The user can specify how the tools should behave if there is no mapping for a given input 
        node. The tool can either remove this node (with all its connections) from the output file
        or it can keep it without mapping.
        
        The user can also define a list of molecule type to filter the nodes for mapping. If the
        input file contains multiple kinds of molecules (MICRO RNA, PROTEIN, GENE, ...) and the user
        specify a filter like 'PROTEIN,GENE', then all the MICRO RNA nodes will be copied to the output
        network without doing any mapping on them. The valid molecule type names are given in the
        MI ontology. We use the standard names, but in a case insensitive way.
        
        Th output of the tool is a standard MITAB file (using standard NavigOmiX identifiers). The
        tool is not sensitive to self-loops, it will copy them to the output file as well, if they
        were present in the input. But there will be no duplicate links in the output file. 
        
        All the metadata from the input connections will be copied to the output file as well without 
        any change, only the fully qualified NavigOmiX identifiers will be updated with the correct 
        molecule type and id type.
        
        From the point of this tool, the "protein" and "gene" molecule types are the same. If you want
        to map the genes and the proteins at the same time, it is enough to give just one of them 
        (protein or gene). (First example!)
        
        
        **Parameters:** 
        
        -i, --input <path>                : input MITAB file path [Mandatory]
        
        -r, --remove                      : should the tool remove the nodes if no mapping can 
                                            be found. [Optional, default: keep]
        
        -f, --molecule-type-selector <str>  : comma separated list of molecule types given in standard 
                                            MI term names (case insensitive). The tool will apply
                                            id mapping only on the nodes with this type. If a node
                                            has different molecule type, then it will be copied to the
                                            output even if `--remove true` was used
                                            [Optional, default: no filtering, try to map all nodes]
        
        -t, --target-id-type <str>        : target id type given in UniProtKB DBref (case insensitive)
                                            [Mandatory]
        
        -m, --mapping-data <paths>        : comma separated list of paths, contains the ID mapping 
                                            information (it is possible to use multiple mapping data
                                            files for optimization) [Mandatory]
        
        -o, --output <path>                : output MITAB file
        
        
        **Exit codes**
        
        Exit code 1: The specified mitab input file doesn't exists!
        Exit code 2: One of the specified mapping file doesn't exists!
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
                        help="<path to a mitab (network) file> [mandatory]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=True)

    parser.add_argument("-r", "--remove",
                        help="should the tool remove the nodes if no mapping can be found [optional]",
                        dest="remove",
                        action="store_true",
                        default=None)

    parser.add_argument("-f", "--molecule-type-selector",
                        help="<comma separated list of molecule types given in standard MI term names \
                        (case insensitive)>",
                        type=str,
                        dest="molecule_type_selector",
                        action="store",
                        default='')

    parser.add_argument("-t", "--target-id-type",
                        help="<target id type given in UniProtKB DBref (case insensitive)> [mandatory]",
                        type=str,
                        dest="target_id_type",
                        action="store",
                        required=True)

    parser.add_argument("-m", "--mapping-data",
                        help="<comma separated list of paths, contains the ID mapping information (it is possible to \
                             use multiple mapping data files for optimization)>  [mandatory]",
                        type=str,
                        dest="mapping_data",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output-file",
                        help="<path to an output mitab (network) file> [mandatory]",
                        type=str,
                        dest="output_file",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_file, results.remove, results.molecule_type_selector.lower(), \
           results.target_id_type.lower(), results.mapping_data, results.output_file


def check_params(input_file, mapping_file_paths):
    if not os.path.isfile(input_file):
        sys.stderr.write(f"ERROR! the specified mitab input file doesn't exists: {input_file}")
        sys.exit(1)

    for mapping_file in mapping_file_paths:
        if not os.path.isfile(mapping_file):
            sys.stderr.write(f"ERROR! one of the specified mapping file doesn't exists: {mapping_file}")
            sys.exit(2)


def write_to_output(output_network, source_node, source_node_id, target_node, target_node_id, tax_id, source_node_molecule_type, target_node_molecule_type,
                    old_source_node, old_source_node_id, old_target_node, old_target_node_id, metadata):
    mitab = mitab_handler.MiTabHandler()

    interaction = mitab.new_interaction()

    interaction[mitab.uidA] = f'{source_node}:{source_node_id}'
    interaction[mitab.uidB] = f'{target_node}:{target_node_id}'
    interaction[mitab.altA] = f'{old_source_node}:{old_source_node_id}'
    interaction[mitab.altB] = f'{old_target_node}:{old_target_node_id}'
    interaction[mitab.taxA] = f'{tax_id}'
    interaction[mitab.taxB] = f'{tax_id}'
    interaction[mitab.annotA] = f'start:{source_node_molecule_type};{source_node};{source_node_id}'
    interaction[mitab.annotB] = f'end:{target_node_molecule_type};{target_node};{target_node_id}'
    interaction[mitab.annotInter] = f'{metadata}'

    mitab.add_interaction(interaction)

    mitab.serialise_mitab(output_network, add_header=False)


def import_mapping_data(mapping_file_paths, target_id_type):

    mapping_dictionary1 = defaultdict(set)
    mapping_dictionary2 = defaultdict(set)

    for file in mapping_file_paths:
        with open(file, "r") as mapping_file:
            for line in mapping_file:
                map_line = json.loads(line.strip())
                if map_line["to_id_type"].lower() == target_id_type:
                    mapping_dictionary1[map_line["from_id"].lower()].add(map_line["to_id"].lower())
                if map_line["from_id_type"].lower() == "uniquename":
                    mapping_dictionary2[map_line["from_id"].lower()].add(map_line["to_id"].lower())
    return mapping_dictionary1, mapping_dictionary2


def map_single_id(id, id_type, molecule_type, requested_mapped_id_type, molecule_types_list, remove, mapping_dictionary,
                  mapping_dictionary_uniquename):

    keep_original = {id}, id_type

    if len(molecule_types_list) > 0 and molecule_type not in molecule_types_list:
        return keep_original

    mapped_ids = mapping_dictionary[id]
    mapped_unique_ids = mapping_dictionary_uniquename[id]
    mapped_id_type = requested_mapped_id_type

    # first we check if it is a protein / gene name, then we will map it to the reviewed ids (if there is any)
    if len(mapped_unique_ids) > 0:
        return mapped_unique_ids, mapped_id_type

    if len(mapped_ids) > 0 or remove:
        return mapped_ids, mapped_id_type

    return keep_original


def id_mapping(input_file, remove, molecule_types_list, requested_mapped_id_type, mapping_dictionary,
               mapping_dictionary_uniquename, output_file):
    with open(output_file, "w") as output_network:

        mitab = mitab_handler.MiTabHandler()
        input_network = mitab.parse_mitab(input_file)
        rows = input_network.values

        for row in rows:
            tax_id = 'taxid:9606(Homo sapiens)'
            source_node_id_type = row[0].split(":")[0].lower()
            target_node_id_type = row[1].split(":")[0].lower()
            source_node_id = row[0].split(":")[1].lower()
            target_node_id = row[1].split(":")[1].lower()
            source_node_molecule_type = row[25].split(":")[1].split(";")[0].lower()
            target_node_molecule_type = row[26].split(":")[1].split(";")[0].lower()
            metadata = row[27]

            mapped_source_ids, mapped_source_id_type = map_single_id(source_node_id, source_node_id_type,
                                                                     source_node_molecule_type,
                                                                     requested_mapped_id_type, molecule_types_list,
                                                                     remove, mapping_dictionary,
                                                                     mapping_dictionary_uniquename)

            mapped_target_ids, mapped_target_id_type = map_single_id(target_node_id, target_node_id_type,
                                                                     target_node_molecule_type,
                                                                     requested_mapped_id_type, molecule_types_list,
                                                                     remove, mapping_dictionary,
                                                                     mapping_dictionary_uniquename)

            for new_source_node_id in mapped_source_ids:
                for new_target_node_id in mapped_target_ids:
                    write_to_output(output_network, mapped_source_id_type, new_source_node_id,
                                    mapped_target_id_type, new_target_node_id, tax_id, source_node_molecule_type,
                                    target_node_molecule_type, source_node_id_type, source_node_id, target_node_id_type,
                                    target_node_id, metadata)

    output_network.close()


def main():

    input_file, remove, molecule_type_selector, target_id_type, mapping_data, output_file = parse_args(sys.argv[1:])
    mapping_files_paths = mapping_data.split(",")

    print(f'====== Starting Network ID Mapper ======')

    check_params(input_file, mapping_files_paths)

    molecule_types = set(filter(lambda x: x != '', map(lambda x: x.strip(), molecule_type_selector.split(","))))
    if "gene" in molecule_types or "protein" in molecule_types:
        molecule_types.add("gene")
        molecule_types.add("protein")

    for mapping_file in mapping_files_paths:
        print(f'====== Loading mapping data from: {mapping_file} ======')
    mapping_dictionary, mapping_dictionary_uniquename = import_mapping_data(mapping_files_paths, target_id_type)

    print(f'====== Writing mapped connections to: {output_file} ======')
    id_mapping(input_file, remove, molecule_types, target_id_type, mapping_dictionary, mapping_dictionary_uniquename,
               output_file)
    print(f'====== Network ID mapper finished successfully! ======')


if __name__ == "__main__":
    main()
