import argparse
import sys
import os


def parse_args(args):
    help_text = \
        """
        === VCF Data Cleaner script===
        
        **Description:** 

        This script will take a single input .vcf file or an input folder, that contains one or more .vcf files.
        Then it will replace the locations of the SNP's in the .vcf file(s) to a new one according to the reference genom
        file(s). In that case you have to use just the --input-file OR just the --input-folder parameter as input and the
        --reference-genome-file parameter as reference.
        If the script does not find a new location based on the reference genom, then the SNP will be deleted and will not be
        in the output. (Example 1)
        
        The script also can take a single SNP list. Then it will make a really new .vcf file from the input SNP list according
        to the reference vcf file(s). In that case you have to use just the --input-snp-list parameter as input and the
        --reference-vcf-file parameter as reference. In this case the script always make mutated_sequence for every SNP (in the output
        vcf files of the second example, the last column will be always "1/1").
        If the script does not find a SNP in the reference vcf file(s), then that SNP will not be in the output vcf file.
        (Example 2)
        
        The --output-folder parameter is mandatory in any cases!
        
        
        **Parameters:**
        
        -i, --input-file <path>                   : path to the input .vcf file [Optional]
        
        -f, --input-folder <path>                 : path to a folder that contains all of the input .vcf files [Optional]
        
        -s, --input-snp-list <path>               : path to the an input SNP list file [Optional]
        
        -rg, --reference-genome-files <paths>     : comma separated list of paths, contains the reference genome information [Optional]
        
        -rv, --reference-vcf-files <paths>        : comma separated list of paths, contains the reference genome information [Optional]
        
        -o, --output-folder <path>                : path to an output folder that contain the changed VCF files [Mandatory]
        
        
        **Exit codes**
        
        Exit code 1: The input .vcf file does not exists!
        Exit code 2: The input .vcf file is not a .vcf file!
        Exit code 3: The specified input folder does not exists!
        Exit code 4: The specified input SNP list does not exists!
        Exit code 5: One of the reference genome file does not exists!
        Exit code 6: One of the reference vcf file does not exists!
        Exit code 7: The specified output folder does not exists!
        Exit code 8: Have to give a single vcf input file or an input folder, which contains vcf files!
        """

    parser = argparse.ArgumentParser(description=help_text)

    # Input vcf file path
    parser.add_argument("-i", "--input-file",
                        help="<path to the input .vcf file> [optional]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=False)

    # Input folder
    parser.add_argument("-f", "--input-folder",
                        help="<path to a folder that contains all of the input .vcf files> [optional]",
                        type=str,
                        dest="input_folder",
                        action="store",
                        required=False)

    # Input SNP list file path
    parser.add_argument("-s", "--input-snp-list",
                        help="<path to the an input SNP list file> [optional]",
                        type=str,
                        dest="input_snp_list",
                        action="store",
                        required=False)

    # Reference genome file(s)
    parser.add_argument("-rg", "--reference-genome-files",
                        help="<path(s) to a set of genome reference file(s)> [mandatory]",
                        type=str,
                        dest="reference_genome_files",
                        action="store",
                        required=False)

    # Reference vcf file(s)
    parser.add_argument("-rv", "--reference-vcf-files",
                        help="<path(s) to a set of reference vcf file(s)> [mandatory]",
                        type=str,
                        dest="reference_vcf_files",
                        action="store",
                        required=False)

    # Output folder path
    parser.add_argument("-o", "--output-folder",
                        help="<path to an output folder that contain the changed VCF files> [mandatory]",
                        dest="output_folder",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_file, results.input_folder, results.input_snp_list, results.reference_genome_files,\
           results.reference_vcf_files, results.output_folder


def check_params(input_file, input_folder, input_snp_list, reference_genome_paths, reference_vcf_paths, output_folder):
    if input_file:
        if not os.path.isfile(input_file):
            sys.stderr.write(f'ERROR! the input .vcf file does not exists: {input_file}')
            sys.exit(1)

        if not input_file.endswith(".vcf"):
            sys.stderr.write(f'ERROR! the input .vcf file is not a .vcf file: {input_file}')
            sys.exit(2)

    if input_folder:
        if not os.path.isdir(input_folder):
            sys.stderr.write(f'ERROR! the specified input folder does not exists: {input_folder}')
            sys.exit(3)

    if input_snp_list:
        if not os.path.isfile(input_snp_list):
            sys.stderr.write(f'ERROR! the input SNP list file does not exists: {input_snp_list}')
            sys.exit(4)

    if reference_genome_paths:
        for file in reference_genome_paths:
            if file != "":
                if not os.path.isfile(file):
                    sys.stderr.write(f'ERROR! one of the reference genome file does not exists: {file}')
                    sys.exit(5)

    if reference_vcf_paths:
        for file in reference_vcf_paths:
            if file != "":
                if not os.path.isfile(file):
                    sys.stderr.write(f'ERROR! one of the reference vcf file does not exists: {file}')
                    sys.exit(6)

    if not os.path.isdir(output_folder):
        sys.stderr.write(f'ERROR! the specified output folder does not exists: {output_folder}')
        sys.exit(7)


def parse_vcf(input_file):

    snp_dictionary = {}

    with open(input_file, "r") as input_vcf:
        for line in input_vcf:
            line = line.strip()
            if line[:1] != "#":
                snp_id = line.split("\t")[2]
                snp_dictionary[snp_id] = "None"

    return snp_dictionary


def parse_snp_list(input_snp_list):

    snp_array = []

    with open(input_snp_list, "r") as input_snp:
        for line in input_snp:
            line = line.strip()
            snp_array.append(line)

    return snp_array


def collect_snp_locations(reference_files_paths, snp_dictionary):

    for reference_file in reference_files_paths:
        if reference_file[-4:-3] == ".":
            print(f'====== Collect new SNP locations from reference file: {reference_file} ======')
            with open(reference_file, "r") as reference:
                for row in reference:
                    if row[:1] != "#":
                        row = row.strip().split("\t")
                        snp_id = row[3]
                        location = row[2]
                        if snp_id in snp_dictionary:
                            snp_dictionary[snp_id] = location

    return snp_dictionary


def collect_snp_information_and_write_to_output(input_snp_list, reference_files_paths, snp_array, output_folder):

    new_filename = f'{input_snp_list.split("/")[-1:][0].split(".")[0]}.vcf'
    output_file = os.path.join(output_folder, new_filename)

    for reference_file in reference_files_paths:
        if reference_file[-4:-3] == ".":
            print(f'====== Collect new SNP informations from reference file: {reference_file} ======')
            with open(reference_file, "r") as reference, open(output_file, "a") as output:
                output.write("##fileformat=VCF" + '\n')
                for row in reference:
                    if row[:1] != "#":
                        row = row.strip().split("\t")
                        snp_id = row[2]
                        if snp_id in snp_array:
                            chromosome = f'chr{row[0]}'
                            location = row[1]
                            ref = row[3]
                            alt = row[4]
                            mutant = f'1/1'
                            print(f'======= Write {snp_id} SNP information to output file ======')
                            output.write(chromosome + '\t' + location + '\t' + snp_id + '\t' + ref + '\t' + alt + '\t'
                                         + "." + '\t' + "." + '\t' + "." + '\t' + "." + '\t' + mutant + '\n')

            output.close()


def snp_dictionary_cleaning(snp_dictionary_with_locations):

    snp_dictionary_cleaned = {}

    for snp_id in snp_dictionary_with_locations:
        if snp_dictionary_with_locations[snp_id] != "None":
            snp_dictionary_cleaned[snp_id] = snp_dictionary_with_locations[snp_id]

    return snp_dictionary_cleaned


def data_cleaning(input_file, output_folder, snp_dictionary_cleaned):

    output_file = os.path.join(output_folder, input_file.split("/")[-1:][0])

    with open(input_file, "r") as input_vcf, open(output_file, "w") as output:
        print(f"====== Replacing old SNP identifiers with new ones in {input_file} ======")
        for line in input_vcf:
            line = line.strip()
            if line[:1] == "#":
                output.write(line + '\n')
            else:
                snp_id = line.split("\t")[2]
                snp_location = line.split("\t")[1]
                if snp_id in snp_dictionary_cleaned:
                    new_line = line.replace(snp_location, snp_dictionary_cleaned[snp_id])
                    output.write(new_line + '\n')


def input_file_given(input_file, reference_files_paths, output_folder):

    print(f'====== Collect SNP identifiers from the input vcf ======')
    snp_dictionary = parse_vcf(input_file)

    snp_dictionary_with_locations = collect_snp_locations(reference_files_paths, snp_dictionary)

    print(f"====== Remove SNP's that not have locations")
    snp_dictionary_cleaned = snp_dictionary_cleaning(snp_dictionary_with_locations)

    data_cleaning(input_file, output_folder, snp_dictionary_cleaned)


def input_folder_given(input_folder, reference_files_paths, output_folder):

    reference_vcf_file = next(os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".vcf"))

    print(f'====== Collect SNP identifiers from the input vcf ======')
    snp_dictionary = parse_vcf(reference_vcf_file)

    snp_dictionary_with_locations = collect_snp_locations(reference_files_paths, snp_dictionary)

    print(f"====== Remove SNP's that not have locations")
    snp_dictionary_cleaned = snp_dictionary_cleaning(snp_dictionary_with_locations)

    for file in os.listdir(input_folder):
        filename = os.path.join(input_folder, file)
        if filename.endswith(".vcf"):
            data_cleaning(filename, output_folder, snp_dictionary_cleaned)


def input_snp_list_given(input_snp_list, reference_files_paths, output_folder):

    print(f'====== Collect SNP identifiers from the input SNP list file ======')
    snp_array = parse_snp_list(input_snp_list)

    collect_snp_information_and_write_to_output(input_snp_list, reference_files_paths, snp_array, output_folder)


def main():

    input_file, input_folder, input_snp_list, reference_genome_files, reference_vcf_files, output_folder = parse_args(sys.argv[1:])

    if input_file and input_folder:
        sys.stderr.write(f'ERROR MESSAGE: you have to give a single vcf input file or an input folder, '
                         f'which contains vcf files!\n')
        sys.exit(8)

    if reference_genome_files:
        reference_genome_paths = reference_genome_files.split(",")
        reference_vcf_paths = None
        check_params(input_file, input_folder, input_snp_list, reference_genome_paths, reference_vcf_paths,
                     output_folder)
        print(f'====== Parameters are fine, starting... ======')

        if input_file:
            input_file_given(input_file, reference_genome_paths, output_folder)

        if input_folder:
            input_folder_given(input_folder, reference_genome_paths, output_folder)

    if reference_vcf_files:
        reference_vcf_paths = reference_vcf_files.split(",")
        reference_genome_paths = None
        check_params(input_file, input_folder, input_snp_list, reference_genome_paths, reference_vcf_paths,
                     output_folder)
        print(f'====== Parameters are fine, starting... ======')

        if input_snp_list:
            input_snp_list_given(input_snp_list, reference_vcf_paths, output_folder)

    print(f'====== VCF Cleaner finished succesfully ======')


if __name__ == "__main__":
    main()
