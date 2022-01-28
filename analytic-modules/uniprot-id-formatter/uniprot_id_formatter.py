import argparse
import sys
import os
import re


def parse_args(args):
    help_text = \
        """
        === Uniprot ID Formatter modul===
        **Description:**

        The tool will take a mitab file as an input, and writes a new output mitab file. The output file is identical to
        the input file, except that all the uniprot identifiers are transformed according to the optional parameters.
        
        You don't have to give an output file or folder, because the script automatically make an output file, according to the
        input network file. It add a "_formatted" tag to the end of the input file name.
        (input file: test.tsv, then the output will be: test_formatted.tsv)
        But if you want to specify a special output file, you can do it with the --output-network-file parameter, as well.
        (Example 4)
        
        The tool should search for uniprot ids in the following columns of the mitab file:
        - columns 1-4 (IDs using colon separator uniprotac:a6pvc2-1)
        - column 26-28 (NOX fully qualified IDs using semicolon separator: protein;uniprotac;a6pvc2-1)
        
        The tool will not change the information what in the 3rd or 4th column (alternative id's). The tool will change this
        information only if the --no-isoform parameter was given and the actual uniprot id has an izoform. In this case
        the tool save the uniprot ID with the izoform to the 3rd or 4th column.
        
        
        **Parameters**
        
        -i, --input-network-file <path>                   : path to an input network file [Mandatory]
        
        -lc, --lower-case                                 : if this parameter is given, all of the ID's will be lowercase, 
                                                            default: Null [Optional]
        
        -uc, --upper-case                                 : if this parameter is given, all of the ID's will be uppercase, 
                                                            default: Null [Optional]
        
        -ni, --no-isoform                                 : if this parameter is given, all the isoforms of the ID's will be
                                                            deleted, default: Null [optional]
        
        -o, --output-network-file <path>                  : path to an output network file [optional]
        
        
        **Exit codes**
        
        Exit code 1: The input network file does not exists!
        Exit code 2: Just one of these parameters (lower_case, upper_case) must be given!
        
        
        **Notes**
        
        1) if both --lover-case and --upper-case parameter is given, then return with non-zero exit code
        2) if neither the --lover-case nor the --upper-case parameter is given, then the case will not be changed
        3) if none of the three parameters are specified, then the output will be exactly same as the input
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-network-file",
                        help="<path to an input network file> [mandatory]",
                        type=str,
                        dest="input_network_file",
                        action="store",
                        required=True)

    parser.add_argument("-lc", "--lower-case",
                        help="<if this parameter is given, all of the ID's will be lowercase default: Null> [optional]",
                        dest="lower_case",
                        action="store_true",
                        default=None)

    parser.add_argument("-uc", "--upper-case",
                        help="<if this parameter is given, all of the ID's will be uppercase default: Null> [optional]",
                        dest="upper_case",
                        action="store_true",
                        default=None)

    parser.add_argument("-ni", "--no-isoform",
                        help="if this parameter is given, all the isoforms of the ID's will be deleted "
                             "default: Null> [optional]",
                        dest="no_isoform",
                        action="store_true",
                        default=None)

    parser.add_argument("-o", "--output-network-file",
                        help="<path to an output network file> [optional]",
                        type=str,
                        dest="output_network_file",
                        action="store",
                        required=False)

    results = parser.parse_args(args)
    return results.input_network_file, results.lower_case, results.upper_case, results.no_isoform, \
           results.output_network_file


def check_params(input_network_file, lower_case, upper_case):

    if not os.path.isfile(input_network_file):
        sys.stderr.write(f'ERROR! the input network file does not exists: {input_network_file}')
        sys.exit(1)

    if lower_case and upper_case:
        sys.stderr.write(f'ERROR! just one of these parameters (lower_case, upper_case) must be given')
        sys.exit(2)


def id_formatting(dictionary, upper_case, no_isoform):
    """
    e.g. having the uniprot ID: uniprotac:P2E7D5-54
    the matching will be:

    Full match	0-19	uniprotac:P2E7D5-54
    Group 1.	10-19	P2E7D5-54
    Group 2.	10-16	P2E7D5
    Group 3.	14-16	D5
    Group 4.	16-19	-54
    """
    pattern = "uniprotac\S(([A-Za-z][0-9]([A-Za-z0-9][A-Za-z0-9]){2,4})(-[0-9]{1,2}){0,1})"

    for index in range(1, len(dictionary)):

        uniprot_id = re.findall(pattern, dictionary[index])

        if uniprot_id != []:

            izoform = uniprot_id[0][3]
            original_id = uniprot_id[0][0]
            new_id = uniprot_id[0][1]
            upper_new_id = new_id.upper()
            upper_original_id = original_id.upper()

            if no_isoform:
                if upper_case:

                    if index == 1 or index == 2:
                        dictionary[index] = f'uniprotac:{upper_new_id}'

                        if izoform != '':
                            if dictionary[index + 2] != '-':
                                dictionary[index + 2] = f'{dictionary[index + 2]}|uniprotac:{upper_original_id}'
                            else:
                                dictionary[index + 2] = f'uniprotac:{upper_original_id}'

                    elif index == 26 or index == 27:
                        dictionary[index] = f'{dictionary[index].split(";")[0]};{dictionary[index].split(";")[1]}' \
                            f';{upper_new_id}'

                else:
                    if index == 1 or index == 2:
                        dictionary[index] = f'uniprotac:{new_id}'

                        if izoform != '':
                            if dictionary[index + 2] != '-':
                                dictionary[index + 2] = f'{dictionary[index + 2]}|uniprotac:{original_id}'
                            else:
                                dictionary[index + 2] = f'uniprotac:{original_id}'

                    elif index == 26 or index == 27:
                        dictionary[index] = f'{dictionary[index].split(";")[0]};{dictionary[index].split(";")[1]}' \
                            f';{new_id}'

            else:
                if upper_case:
                    if index == 1 or index == 2:
                        dictionary[index] = f'uniprotac:{upper_original_id}'

                    elif index == 26 or index == 27:
                        dictionary[index] = f'{dictionary[index].split(";")[0]};{dictionary[index].split(";")[1]}' \
                            f';{upper_original_id}'

    return dictionary


def main():

    input_network_file, lower_case, upper_case, no_isoform, output_network_file = parse_args(sys.argv[1:])

    check_params(input_network_file, lower_case, upper_case)
    print(f'====== Parameters are fine, starting... ======')

    output_file = f'{input_network_file.split(".")[0]}_formatted.tsv'

    if output_network_file:
        output_file = output_network_file

    with open(input_network_file, 'r') as input_network, open(output_file, "w") as output_network:

        print(f'====== Write results to file: {output_file} ======')
        for line in input_network:
            line = line.strip().split("\t")
            dictionary = {}

            for index in range(0, len(line)):
                dictionary[index + 1] = line[index]

            formatted_dictionary = id_formatting(dictionary, upper_case, no_isoform)

            for keys, value in formatted_dictionary.items():
                if keys == 42:
                    output_network.write(value)
                else:
                    output_network.write(value + '\t')
            output_network.write('\n')

    print(f'====== Uniprot ID Formatter succesfully finished ======')


if __name__ == "__main__":
    main()
