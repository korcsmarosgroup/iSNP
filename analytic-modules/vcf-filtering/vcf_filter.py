import argparse
import shutil
import sys
import os
sys.path.append("./analytic-modules/")
from common_libs.filters import filter_vcf


def check_args(args=None):
    """ Parse arguments """

    help_text = \
        """
        VCF Filter

        Description:
        This tool will take a set of SNP data (in format of VCF), then by the location information in the VCF files it 
        will extract an upstream and downstream region from the genome in a given length.

        Parameters:
        -i, --input : obtained as a result of a sequencing experiment - (VCF file format)
        -a, --annotation : this is a .bed file, that contains the appropriate annotations: miRNA coding regions and 
            gene promoter regions - (BED file format)
        -s, --snp : this is a NavigOmix entity set file, that contains the set of SNP IDs for filtering
        -o, --output : containing just filtered mutations - (VCF/BED file format)
        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    # Input file path
    parser.add_argument("-i", "--input_file",
                        help="<obtained as a result of a sequencing experiment - (VCF file format)> [mandatory]",
                        dest="input",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-a", "--annotation",
                        help="<this is a .bed file, that contains the appropriate annotations: miRNA coding "
                             "regions and gene promoter regions - (BED file format)> [mandatory]",
                        dest="annotation",
                        action="store",
                        default=None)

    # Wait time
    parser.add_argument("-s", "--snp",
                        help="<this is a NavigOmix entity set file, that contains the set of "
                             "SNP IDs for filtering> [optional, default is 0]",
                        dest="snp",
                        action="store",
                        default=None)

    # Exit code
    parser.add_argument("-o", "--output",
                        help="<containing just filtered mutations - (VCF/BED file format)>",
                        dest="output",
                        action="store",
                        required=True)

    # Get arguments
    results = parser.parse_args(args)
    return results.input, results.annotation, results.snp, results.output


def vcf_filter(input_vcf_file, output, annotation=None, snp=None):
    """
    Filter of SNPs in VCF files based on bed files and/or snp lists.
    """
    if snp:
        temp_out = output + ".temp"
    else:
        temp_out = output
    if annotation:
        print(input_vcf_file, " : ", annotation, " : ", temp_out)
        filter_vcf.filter_vcf(input_vcf_file, annotation, temp_out, bed_info=True)
    else:
        shutil.copyfile(input_vcf_file, temp_out)
    if snp:
        filter_vcf.filter_vcf_by_id(temp_out, snp, output)
    if snp:
        os.remove(temp_out)


def main(argv):
    input_file, annotation, snp, output = check_args(argv)
    if os.stat(input_file).st_size == 0:
        print(f'====== The input VCF file is empty! ======')
        open(output, "a").close()
    else:
        vcf_filter(input_vcf_file=input_file,
                   output=output,
                   annotation=annotation,
                   snp=snp)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
