""" NaviTest, Author: Madgwick """

import argparse
import os.path
import sys
import time
from shutil import copyfile


class NaviTest:
    """
    Takes the parsed file paths and copies the input file to the outpur path.
    If a wait time is given then after reading the file, NaviTest will wait
    the allotted time before writing the dictionary to a new file.
    """

    def __init__(self, input_file, exit_code):
        """ Initialise parameters """
        self.data_structure = []
        self.input_file = input_file
        self.exit_code = exit_code
        self.check_input_file

    def check_input_file(self):
        """ Die, if input file not exists """
        if not os.path.isfile(self.input_file):
            print("ERROR! the specified input file doesn't exists: " + self.input_file)
            sys.exit(1)

    def copy_to_output(self, output_file):
        """ Copies the input file to the output """
        copyfile(self.input_file, output_file)

    @staticmethod
    def wait(wait_sec):
        """ Wait for wait_sec number of seconds """
        time.sleep(wait_sec)

    def navi_exit(self):
        """ Exit with the given exit code """
        # Return exit codex
        sys.exit(self.exit_code)


def check_args(args=None):
    """ Parse arguments """

    help_text = \
        """
        === NaviTest ===
        Name of the tool: NaviTest.
        \nDescription:
        The tool will be used in integration tests, simulating the different behaviours of the
        analytical tools. It is basically repeating in the output file what it reads from the
        input file. Before creating the output it is optionally waiting for the specified time.
        After producing the output file, it will exit with the exit code given as a parameter.
        It will create the output file even if the exit code is non-zero.
        \nParameters:
        --input <path to an existing file> [mandatory]
        --output <path to the new output file> [mandatory]
        --waitsec <0..255> [optional, default is 0]
        --exitcode <0..255> [optional, default is 0]
        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    # Input file path
    parser.add_argument("-i", "--input",
                        help="<path to an existing file> [mandatory]",
                        dest="input",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-o", "--output",
                        help="<path to an new file> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    # Wait time
    parser.add_argument("-w", "--waitsec",
                        help="<0..255> [optional, default is 0]",
                        type=int,
                        dest="waitsec",
                        action="store",
                        default=0)

    # Exit code
    parser.add_argument("-e", "--exitcode",
                        help="<0..255> [optional, default is 0]",
                        type=int,
                        dest="exitcode",
                        action="store",
                        default=0)

    # Get arguments
    results = parser.parse_args(args)

    return results.input, results.output, results.waitsec, results.exitcode


def main():
    """ Main file for the commandline parser and builds new NaviTest """
    # Commandline parser
    input_file, output_file, wait_sec, exit_code = check_args(sys.argv[1:])

    # both the stdout and stderr will be captured by the executor service and logged to the docker logs
    print("executing NaviTest")

    # Build Navi-Test
    navi = NaviTest(input_file, exit_code)
    navi.wait(wait_sec)
    navi.copy_to_output(output_file)
    navi.navi_exit()


if __name__ == "__main__":
    main()
