import sys
import argparse
from collections import defaultdict

import pandas as pd


def check_args(args=None):
    """
    A CLI module for parsing arguments for iSNP stacked classifier, including file location for
    the produced results for further analysis.
    """

    help_text = \
        """
        === CLI - maf2vcf ===
        Name of the tool: maf2vcf

        Description: Convert Mutation Annotation Files (MAF) to Variant C Files (VCF) such that the 
        iSNP workflow can process the SNP information. The recommended variant caller for the best 
        results is the Mutect2 workflow. However, this tool is compatible with other variant and mutation
        calling pipelines. Null and Novel SNPs which do not have a dbSNP ID are removed from the 
        analysis. The VCF files will be created with the tumour sample as the filename.

        Parameters:
        --maf -m <path to the original maf file (.maf)> [mandatory]
        --patientvcf -vcf <path to where the patient vcf files should be saved> [mandatory]
        --clinical -c <str, a path to where the patient clinical data is being held> [optional] 
        --verbose -v <boolean, whether to display verbose information> [optional, default is False]
        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    # Input file path
    parser.add_argument("-s", "--sample-sheet",
                        help="<path to an existing file> [mandatory]",
                        dest="sample_sheet_path",
                        action="store",
                        required=True)

    # Input file path
    parser.add_argument("-m", "--maf-path",
                        help="<path to an existing file> [mandatory]",
                        dest="maf_path",
                        action="store",
                        required=True)

    # Get arguments
    results = parser.parse_args(args)

    return results.sample_sheet_path, results.maf_path


def map_maf_to_sample_sheet(sheet_path, maf_path):
    sample_sheet_cols = ['File Name', 'Case ID']
    maf_path_cols = ['Tumor_Sample_Barcode', 'Matched_Norm_Sample_UUID']
    sample_sheet = pd.read_csv(sheet_path, sep='\t')
    maf_path = pd.read_csv(maf_path, sep='\t')

    matched_tumour_id = defaultdict()
    for id in maf_path['Tumor_Sample_Barcode'].unique():
        stripped_id = id[:12]
        matched_tumour_id[stripped_id] = id

    print(matched_tumour_id)

    sample_sheet = sample_sheet[sample_sheet_cols]
    # sample_sheet = sample_sheet.set_index(['Case ID'])
    sample_sheet['maf_id'] = sample_sheet['Case ID'].map(matched_tumour_id)
    print(sample_sheet.head())


def main(argv):
    """ Main method for converting the maf file to vcf files"""
    sheet_path, maf_path = check_args(argv)
    map_maf_to_sample_sheet(sheet_path, maf_path)


if __name__ == "__main__":
    """ Run main method """
    sys.exit(main(sys.argv[1:]))
