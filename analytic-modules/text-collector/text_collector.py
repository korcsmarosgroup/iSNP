import argparse
import os
import sys


class TextCollector:

    def __init__(self, input_files, output_file):
        self.input_files = input_files
        self.output_file = output_file

    def check_params(self):
        for input_file in self.input_files:
            if not os.path.isfile(input_file):
                print("ERROR! the specified input file doesn't exists: " + input_file)
                sys.exit(1)

    def run_tool(self):
        with open(self.output_file, 'w') as outp:
            for input_path in self.input_files:
                print("processing input file: " + input_path)
                with (open(input_path, 'r')) as inp:
                    line = inp.readline()
                    while line:
                        outp.write(line)
                        if not line.endswith("\n"):
                            outp.write("\n")
                        line = inp.readline()


def parse_args(args=None):
    description = \
        """
    Description:
    
    The collector modules are the output points of iterations in the workflow. These modules can be 
    used to propagate the results from the iterations to the rest of the workflow.
    
    This collector takes NavigOmiX primitive files (text files) as input from 
    each iteration, and combining them into a single primitive output variable.
    
    The tool will make sure that the input primitives from each iteration will be started in a 
    new line (adding extra new line characters if needed).
    
    Parameters:
    --input-files <comma separated list of primitive files from each iteration> [mandatory]
    --output-file <path to the output primitive file> [mandatory]   
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input-files",
                        dest='input_files',
                        action='store',
                        required=True)
    parser.add_argument("-o",
                        "--output-file",
                        dest='output_file',
                        action='store',
                        required=True)

    results = parser.parse_args(args)
    return results.input_files, results.output_file


def main():
    input_files, output_file = parse_args(sys.argv[1:])

    print("starting TextCollector")

    input_file_list = list(filter(lambda x: len(x) > 0, input_files.split(',')))
    module = TextCollector(input_file_list, output_file)
    module.check_params()
    module.run_tool()

    print("TextCollector finished successfully")


if __name__ == "__main__":
    main()
