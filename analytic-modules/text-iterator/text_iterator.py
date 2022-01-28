import argparse
import os.path
import sys


class TextIterator:

    def __init__(self, input_file, output_folder):
        self.input_file = input_file
        self.output_folder = output_folder

    def check_params(self):
        if not os.path.isfile(self.input_file):
            print("ERROR! the specified input file doesn't exists: " + self.input_file)
            sys.exit(1)
        if not os.path.isdir(self.output_folder):
            print("ERROR! the specified output folder doesn't exists: " + self.output_folder)
            sys.exit(2)

    def run_tool(self):
        with open(self.input_file, 'r') as inp:
            for index, line in enumerate(inp):
                output_file_path = self.get_output_file_path(index)
                with (open(output_file_path, 'w')) as outp:
                    outp.write(line)
                print("output file created: " + output_file_path)

    def get_output_file_path(self, index):
        output_file_name = "iteration-{}".format(index + 1)
        return os.path.join(self.output_folder, output_file_name)


def parse_args(args=None):
    help_text = \
        """
        === NaviTest ===
        Name of the tool: TextIterator.
        \nDescription:
        The iterator modules are the entry points of iterations in the workflow. This 
        iterator takes a NavigOmiX primitive (text file) as an input, and then 
        iterate over each line of the input file, providing single NavigOmiX primitive 
        variable for each iteration, represented as a separate text file.
        \nThe output is an ordered set of NavigOmiX primitive files in a given output folder. 
        The file names will be `iteration-1`, `iteration-2`, ... The ordering of the 
        output files is deterministic, following the order of lines in the input file.
        \nParameters:
        --input-file <path to the input file> [mandatory]   
        --output-folder [mandatory]"""

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
                        help="<path to the input file> [mandatory]",
                        dest="input",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output-folder",
                        help="<path to the output folder> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input, results.output


def main():
    input_file, output_folder = parse_args(sys.argv[1:])

    print("starting TextIterator")

    module = TextIterator(input_file, output_folder)
    module.check_params()
    module.run_tool()

    print("TextIterator finished successfully")


if __name__ == "__main__":
    main()
