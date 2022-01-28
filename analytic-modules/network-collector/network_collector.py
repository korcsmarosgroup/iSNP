import argparse
import os
import sys
import tarfile


def collect_networks(input_file_list, output_file):
    with tarfile.open(output_file, "w:gz") as tar:
        for name in input_file_list:
            tar.add(name)


def check_params(input_file_list):
    for input_file in input_file_list:
        if not os.path.isfile(input_file):
            print("ERROR! the specified input file doesn't exists: " + input_file)
            sys.exit(1)


def parse_args(args=None):
    description = \
        """
    Description:

    The collector modules are the output points of iterations in the workflow. These modules can be used to propagate 
    the results from the iterations to the rest of the workflow. This collector takes a single network file as an input 
    from each iteration, and producing a single network set file as an output.

    The input format is a standard NavigOmiX MITAB file. The output is a NavigOmiX network set file, a tar.gz archive 
    with all the network files collected in each iteration.

    Parameters:

    --input-files <comma separated list of network files from each iteration> [mandatory]
    --output-file <path to the output tar.gz network set file> [mandatory]   
    """

    parser = argparse.ArgumentParser(description=description)

    # Input folder
    parser.add_argument("-i", "--input-files",
                        dest='input_files',
                        action='store',
                        required=True)

    # Output file path
    parser.add_argument("-o",
                        "--output-file",
                        dest='output_file',
                        action='store',
                        required=True)

    results = parser.parse_args(args)

    return results.input_files, results.output_file


def main():
    """ Run module """
    input_files, output_file = parse_args(sys.argv[1:])
    input_file_list = input_files.split(',')
    check_params(input_file_list)
    collect_networks(input_file_list, output_file)


if __name__ == "__main__":
    main()
