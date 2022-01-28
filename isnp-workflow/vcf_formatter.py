import sys
import os
import csv
import argparse


def parse_args(args):
    """ Argument parser """
    parser = argparse.ArgumentParser(description="Hg19 vcf reformatting for iSNP workflow")

    # Input folder
    parser.add_argument("-i", "--input-folder",
                        dest='input_folder',
                        action='store',
                        required=True)

    # Output file path
    parser.add_argument("-o",
                        "--output-folder",
                        dest='output_folder',
                        action='store',
                        required=True)

    results = parser.parse_args(args)
    return results.input_folder, results.output_folder


def clean_h19_vcf_files(directory, output_folder):
    """ Clean all vcf files and write them to a new directory if they are hg19 genome """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    for vcf_file_path in os.listdir(directory):
        if vcf_file_path.endswith('.vcf'):

            vcf_file_path = os.path.join(directory, vcf_file_path)
            
            # Create new file name
            print(f"VCF path: {vcf_file_path}")
            path, filename = os.path.split(vcf_file_path)
            output_vcf_file = os.path.join(output_folder, filename)
            print(f"Ouptut: {output_vcf_file}")
            
            # Write new file
            with open(vcf_file_path, 'r') as vcf_file, open(output_vcf_file, 'w') as amended_vcf:
                reader = csv.reader(vcf_file, delimiter='\t')
                writer = csv.writer(amended_vcf, delimiter='\t')
                for row in reader:
                    if row[0].__contains__('#'):
                        writer.writerow(row)
                    else:
                        row[0] = f"chr{row[0]}"
                        writer.writerow(row)
                        
                        
def main(argv):
    """ Main method to clean the h19 vcf files """
    input_folder, output_folder = parse_args(argv)
    clean_h19_vcf_files(directory=input_folder, output_folder=output_folder)
                        
                        
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
