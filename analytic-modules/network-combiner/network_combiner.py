import argparse
import sys
import os
from common_libs.mitab_handler import mitab_handler


def parse_args(args):
    help_text = \
        """
        === Network Combiner ===

        **Description:** 

        The network combiner takes a set of networks as input, and delivers a single network by 
        applying an operation on the input networks. It can either calculate the union or the 
        intersection of the networks, or the set difference of the arguments.

        At least one input network file must be provided. If all the input files are representing 
        empty networks, then an output file will be also empty. The output network will contain no 
        duplicate links. We are treating all the links in the input files as undirected.

        When the 'union' method is chosen, then all those links will appear in the output networks,
        which were present in any of the input networks. For the 'intersection' method, only the links
        present in all input networks should be added to the output. If the 'difference' method was selected, 
        then only those links will appear in the output, which were not present in all the input files. 
        (So in this definition: difference = union - intersection)

        If there is only a single network as an input, then the union and the intersection will be equal 
        to the input file, while the difference is an empty network.

        The metadata stored for the input files for the nodes or links will not be preserved in the 
        output file.


        **Parameters:**

        --input-files <comma separated list of network files from each iteration> [mandatory]

        --output-file <path to an output network set file> [mandatory]   

        --method <method name: union|intersection|difference> [mandatory]
        
        
        **Exit codes:**

        Exit code 1: One of the specified input file doesn't exists!
        Exit code 2: The method doesn't exists!
        """

    parser = argparse.ArgumentParser(description=help_text)

    # Comma separated paths of input network files
    parser.add_argument("-i", "--input-files",
                        help="<paths to the input vcf files> [mandatory]",
                        type=str,
                        dest="input_files",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-o", "--output-file",
                        help="<path to an output network set file> [mandatory]",
                        dest="output_file",
                        action="store",
                        required=True)

    # Method
    parser.add_argument("-m", "--method",
                        help="<union or intersection or difference> [mandatory]",
                        type=str,
                        dest="method",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_files, results.output_file, results.method


def check_params(input_file_list, method):
    for input_file in input_file_list:
        if not os.path.isfile(input_file):
            sys.stderr.write(f"ERROR! one of the specified input file doesn't exists: {input_file}")
            sys.exit(1)

    methods = ["union", "intersection", "difference"]
    if method not in methods:
        sys.stderr.write(f"ERROR! the method doesn't exists: {method}")
        sys.exit(2)


def one_input_file(input_file, output_file, method):

    with open(input_file, "r") as input_network, open(output_file, "w") as output_network:
        for line in input_network:
            info = line.strip().split("\t")

            mitab = mitab_handler.MiTabHandler()

            interaction = mitab.new_interaction()

            interaction[mitab.uidA] = f'{info[0]}'
            interaction[mitab.uidB] = f'{info[1]}'
            interaction[mitab.taxA] = f'taxid:9606(Homo Sapiens)'
            interaction[mitab.taxB] = f'taxid:9606(Homo Sapiens)'
            interaction[mitab.annotA] = f'{info[25]}'
            interaction[mitab.annotB] = f'{info[26]}'
            interaction[mitab.annotInter] = f'{info[27]}'

            mitab.add_interaction(interaction)

            if method == "union":
                mitab.serialise_mitab(output_network, add_header=False)

            elif method == "intersection":
                mitab.serialise_mitab(output_network, add_header=False)

            elif method == "difference":
                output_network.write("")

        output_network.close()


def parse_network_files(input_files_path, meta_data, snp_meta_data):

    results = []

    for file in input_files_path:

        output = set()

        with open(file, "r") as input_network:

            for line in input_network:

                info = line.strip().split("\t")
                interactor_a = f'{info[0]}'
                interactor_b = f'{info[1]}'
                molecule_type_a = f'{info[25].split(":")[1].split(";")[0]}'
                molecule_type_b = f'{info[26].split(":")[1].split(";")[0]}'
                metadata = f'{info[27]}'

                if (interactor_a, interactor_b) not in snp_meta_data:
                    snp_meta_data[interactor_a, interactor_b] = []
                snp_meta_data[interactor_a, interactor_b].append(metadata)

                if interactor_a not in meta_data:
                    meta_data[interactor_a] = []
                meta_data[interactor_a].append(molecule_type_a)

                if interactor_b not in meta_data:
                    meta_data[interactor_b] = []
                meta_data[interactor_b].append(molecule_type_b)

                tuples = (f'{interactor_a}', f'{interactor_b}')
                tuples_reverse = (f'{interactor_b}', f'{interactor_a}')

                if [item for item in results if tuples_reverse in item]:

                    if tuples_reverse in output:
                        continue

                    else:
                        output.add(tuples_reverse)

                else:
                    if tuples_reverse in output:
                        continue

                    else:
                        output.add(tuples)

            results.append(output)

    return results


def write_to_file(output_file, input_method_array, meta_data, snp_meta_data):
    with open(output_file, "w") as output_network:
        for interactions in input_method_array:

            mitab = mitab_handler.MiTabHandler()

            interaction = mitab.new_interaction()

            interaction[mitab.uidA] = f'{interactions[0]}'
            interaction[mitab.uidB] = f'{interactions[1]}'
            interaction[mitab.taxA] = f'taxid:9606(Homo Sapiens)'
            interaction[mitab.taxB] = f'taxid:9606(Homo Sapiens)'
            metadata_a = interactions[0].split(":")
            interaction[mitab.annotA] = f'start:{meta_data[interactions[0]][0]};{metadata_a[0]};{metadata_a[1]}'
            metadata_b = interactions[1].split(":")
            interaction[mitab.annotB] = f'end:{meta_data[interactions[1]][0]};{metadata_b[0]};{metadata_b[1]}'

            metadata_array = []
            if interactions in snp_meta_data:
                for snp_metadata in snp_meta_data[interactions]:
                    if snp_metadata not in metadata_array:
                        metadata_array.append(snp_metadata)

            delimiter = "|"
            interaction[mitab.annotInter] = delimiter.join(metadata_array)

            mitab.add_interaction(interaction)

            mitab.serialise_mitab(output_network, add_header=False)

    output_network.close()


def comparing_networks(input_files, output_file, method, meta_data, snp_meta_data):
    results = parse_network_files(input_files, meta_data, snp_meta_data)

    union = set.union(*results)
    intersection = set.intersection(*results)
    difference = union - intersection

    union_array = sorted(union)
    intersection_array = sorted(intersection)
    difference_array = sorted(difference)

    if method == "union":
        write_to_file(output_file, union_array, meta_data, snp_meta_data)

    elif method == "intersection":
        write_to_file(output_file, intersection_array, meta_data, snp_meta_data)

    elif method == "difference":
        write_to_file(output_file, difference_array, meta_data, snp_meta_data)


def main():

    input_files, output_file, method = parse_args(sys.argv[1:])
    input_file_list = input_files.split(',')

    method = method.lower()

    meta_data = {}
    snp_meta_data = {}

    check_params(input_file_list, method)

    if len(input_file_list) == 1:
        if os.stat(input_file_list[0]).st_size == 0:
            open(output_file, "a").close()

        else:
            one_input_file(input_file_list[0], output_file, method)

    else:
        empty_files = []
        for file in input_file_list:
            if os.stat(file).st_size == 0:
                empty_files.append(file)

        if len(input_file_list) == len(empty_files):
            open(output_file, "a").close()

        else:
            comparing_networks(input_file_list, output_file, method, meta_data, snp_meta_data)

    print(f"====== Network combiner finished successfully! ======")


if __name__ == "__main__":
    main()

