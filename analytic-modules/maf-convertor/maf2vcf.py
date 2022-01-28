""" Convert MAF files to VCF files """
import argparse
import os
import sys
import pandas as pd
import numpy as np
from collections import Counter


# VCF Header information
vcf_header = \
    """##fileformat=VCFv4.3
##fileDate=20181211
##source=PLINKv2.00
##contig=<ID=1,length=247177331>
##contig=<ID=2,length=242626207>
##contig=<ID=3,length=199071699>
##contig=<ID=4,length=190980781>
##contig=<ID=5,length=180512666>
##contig=<ID=6,length=170723976>
##contig=<ID=7,length=158710966>
##contig=<ID=8,length=146264219>
##contig=<ID=9,length=140111181>
##contig=<ID=10,length=135040420>
##contig=<ID=11,length=134268352>
##contig=<ID=12,length=132228484>
##contig=<ID=13,length=114106016>
##contig=<ID=14,length=106296207>
##contig=<ID=15,length=100194179>
##contig=<ID=16,length=88655878>
##contig=<ID=17,length=78634367>
##contig=<ID=18,length=76115294>
##contig=<ID=19,length=63556292>
##contig=<ID=20,length=62322843>
##contig=<ID=21,length=46840115>
##contig=<ID=22,length=49503533>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n"""

# MAF file column headings
MAF_COLUMNS = ["case_id", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "dbSNP_RS", "Reference_Allele",
               "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]

# VCF file column headings
VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


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
    parser.add_argument("-m", "--maf",
                        help="<path to an existing file> [mandatory]",
                        dest="maf",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-o", "--output",
                        help="<path to an new file> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    # Clinical file path
    parser.add_argument("-c", "--clinical",
                        help="<path to clinical file> [optional]",
                        dest="clinical_file",
                        action="store",
                        required=False)

    # Threshold
    parser.add_argument("-t", "--threshold",
                        help="<integer value for threshold to filter somatic mutations> [optional]",
                        dest="threshold",
                        type=int,
                        action="store",
                        required=False)

    # Reference
    parser.add_argument("-r", "--reference",
                        help="<path to reference file> [optional]",
                        dest="reference",
                        action="store",
                        required=False)

    # Verbosity level
    parser.add_argument("-v", "--verbose",
                        help="Verbosity level (Default: False)",
                        dest="verbose",
                        action="store",
                        default=False)

    # Get arguments
    results = parser.parse_args(args)

    return results.maf, results.output, results.clinical_file, results.threshold, results.reference, results.verbose


def map_maf_files_to_clinical(patient_id, clinical_file, output_path):
    """
    A function to clean and process the clinical data to the patients annonated by the MAF file.
    If a column has more than half of the attributes as a NaN then the entire column is dropped
    and will not be considered in any downstream analysis. The dataframe is also reindexed for the
    case_id.

    Parameters
    ----------
    patient_id : list, list of a string containing the tumour-id showing the case (patient)
    clinical_file : str, path to the currently unprocessed clinical file
    output_path : str, path to the created directory for the cleaned clinical data

    Returns
    -------
    cleaned_clinical_file_path : str, file path to the newly wrangled and mapped clinical file path
    """
    cleaned_file = f"cleaned-{os.path.basename(clinical_file)}"
    output_path = os.path.join(output_path, cleaned_file)
    raw_clinical_file = pd.read_csv(clinical_file, sep='\t', engine='python')
    col_threshold = len(raw_clinical_file) / 2
    row_threshold = len(list(raw_clinical_file.columns.values)) / 2
    raw_clinical_file.set_index('case_id', inplace=True)
    raw_clinical_file = raw_clinical_file.replace('--', np.NaN)
    raw_clinical_file = raw_clinical_file.dropna(thresh=col_threshold, axis=1, how="all")
    raw_clinical_file = raw_clinical_file.dropna(thresh=row_threshold, axis=0)
    cleaned_patient_id = [case_id for case_id in list(raw_clinical_file.index.values) if case_id in patient_id]
    raw_clinical_file = raw_clinical_file.loc[cleaned_patient_id]
    raw_clinical_file.to_csv(output_path, sep='\t')


def _create_patient_vcf_file(patient_vcf_data, patient_folder, patient_id):
    """
    From the pre-indexed patient maf file create a vcf file for the patient. The SNPs are
    filtered and only point mutations are taken. Frame-shift and large deletions are
    ignored from the somatic mutations listed in the MAF files. The case-id is used as the
    patient_id and then QUAL and FILTER are placed a blank as default of a VCF file. To
    determine the genotype information the reference allele is compared against the tumour
    allele 1 and tumour allele 2. If they match a 0 is placed down and if they are different
    a 1 is put placed in the genotype string.

    Parameters
    ----------
    patient_vcf_data : Pandas DataFrame, a dataframe holding the
    patient_folder : str, path to the folder where the patient vcf files should be written too
    patient_id : list, list of a string containing the tumour-id showing the case (patient)

    Returns
    -------
    snp_identifiers : list, a  list of SNPs (dbSNP) found in the patient
    """
    patient_file_name = os.path.join(patient_folder, f"{patient_id}.vcf")
    tumour_id = f"{patient_id}_{patient_id}"
    vcf_cols = VCF_COLUMNS + [tumour_id]
    snp_identifiers = []

    with open(patient_file_name, "w") as vcf_file:
        # Write header and contig information
        vcf_file.writelines(vcf_header)
        vcf_file.writelines([col + "\t" for col in vcf_cols])
        vcf_file.write("\n")

        # Create new vcf columns
        new_row = dict.fromkeys(vcf_cols)
        for row in patient_vcf_data.itertuples():
            if str(row.dbSNP_RS).startswith("rs") and row.Chromosome != "chrX":
                if len(row.Reference_Allele) == 1 and len(row.Tumor_Seq_Allele2) == 1 and row.Reference_Allele != "-":
                    # Process the genotype
                    genotype_value = ["0", "/", "0"]
                    if row.Reference_Allele != row.Tumor_Seq_Allele2:
                        genotype_value[0] = "1"
                    if row.Tumor_Seq_Allele1 != row.Tumor_Seq_Allele2:
                        genotype_value[2] = "1"
                    genotype_value = "".join(genotype_value)

                    # Create new entry
                    new_row["CHROM"] = row.Chromosome
                    new_row["POS"] = row.Start_Position
                    new_row["ID"] = row.dbSNP_RS
                    new_row["REF"] = row.Reference_Allele
                    if row.Tumor_Seq_Allele2 == "-":
                        new_row["ALT"] = "*"
                    else:
                        new_row["ALT"] = row.Tumor_Seq_Allele2
                    new_row["QUAL"] = "."
                    new_row["FILTER"] = "."
                    new_row["INFO"] = "PR"
                    new_row["FORMAT"] = "GT"
                    new_row[tumour_id] = genotype_value

                    for item in new_row.items():
                        vcf_file.write(f"{item[1]}\t")

                    vcf_file.write("\n")
                    snp_identifiers.append(row.dbSNP_RS)

    return snp_identifiers


def convert_maf_to_vcf(maf, patient_folder, verbose):
    """
    For each of the tumour samples and therefore cases found in the MAF file, iterate over these
    and create a vcf from the subset indexed on the tumour id. Then call the internal function
    _create_patient_vcf_file() to create the individual vcf files. This will return a list of
    SNPs found in the patient.

    Parameters
    ----------
    maf : str, path to the maf file
    patient_folder : str, to the folder of the patients
    verbose : boolean, to display verbose details or not

    Returns
    -------
    maf_frame : Pandas dataframe, containing a tabulated MAF file
    """
    snp_identifiers = []
    maf_frame = pd.read_csv(maf, skiprows=5, delimiter="\t", usecols=MAF_COLUMNS, index_col=False)
    maf_frame = maf_frame.dropna()
    maf_frame = maf_frame[~maf_frame.dbSNP_RS.str.contains("novel")]
    patient_ids = pd.unique(maf_frame[MAF_COLUMNS[0]])

    for case in patient_ids:
        patient_frame = maf_frame[MAF_COLUMNS[0]].isin([case])
        patient_frame = maf_frame[patient_frame]
        snp_list = _create_patient_vcf_file(patient_frame, patient_folder, case)
        snp_identifiers.extend(snp_list)

    return snp_identifiers, patient_ids


def filter_snp_identifiers(snp_identifiers, reference, threshold):
    """
    Filter unique patient files based on the reference (census) data set. There is also a further option to
    instead use a threshold for variance control of these SNPs. If a SNP doesn't appear within the cutoff
    value it will not be added to the final SNP-identifier file.
    Parameters
    ----------
    snp_identifiers : list, a list of all the snp identifiers found in the patient cohort
    reference : str, file path to the reference genome
    threshold : int, a value for the number of occurrence of a SNP throughout the patient cohort

    Returns
    -------
    snp_identifiers : list, a list of SNP names from all the patient VCF files after being filtered
    """
    snp_filtered = []
    snp_occur = Counter(snp_identifiers)
    for snp in snp_occur.items():
        if int(snp[1]) >= threshold:
            snp_filtered.append(snp[0])
    return snp_filtered


def create_snp_identifiers_list(snp_identifiers, patient_folder, reference, threshold):
    """
    Creates a snp-identifiers file for the iSNP workflow from the all the processed MAF files. This means
    that no SNPs that are invalid to the workflow end up in the pipeline.

    Parameters
    ----------
    snp_identifiers : list, a list of SNP names from all the patient VCF files
    patient_folder : str, path to the folder where the patient vcf files should be written too
    threshold : int, an integer threshold for finding the snp list from the snps in the patient cohort
    reference : str, a file path to a reference disease-related snps to filter the patient cohort
    """
    if reference or threshold:
        snp_identifiers = filter_snp_identifiers(snp_identifiers, reference, threshold)

    file_name = "SNP-identifiers.txt"
    patient_path = os.path.join(patient_folder, file_name)
    with open(patient_path, "w") as snp_file:
        for snp in snp_identifiers:
            snp_file.write(f"{snp}\n")


def check_file_format(maf, patient_folder):
    """
    A method to check the file type of the maf input, if it is not a maf file then a TypeError will be
    raised.

    Parameters
    ----------
    maf : str, path to the maf file
    patient_folder : str, path to output folder, here termed patient_folder
    """
    # Check file is a maf file
    if not maf.lower().endswith(('.txt', '.maf')):
        sys.exit(2)

    # Check that base folder exists
    if not os.path.exists(os.path.dirname(patient_folder)):
        try:
            os.mkdir(os.path.dirname(patient_folder))
        except OSError as exc:
            raise FileNotFoundError(f"Directory {exc} failed to generate.")


def _create_output_dir(output_path):
    """
    A method to create the create input directory structure.
    root-path
    | - cleaned-cohort-vcf
    | - cleaned-clinical
    | - cleaned-snp-list

    Returns
    -------
    dir_list : list, a list of directories and created folders
    """
    dir_list = ["cleaned-cohort-vcf/", "cleaned-clinical/", "cleaned-snp-list/"]
    for idx, folder in enumerate(dir_list):
        folder_path = os.path.join(output_path, folder)
        print(folder_path)
        dir_list[idx] = folder_path
        os.mkdir(os.path.dirname(folder_path))

    return dir_list


def main(argv):
    """ Main method for converting the maf file to vcf files"""
    try:
        maf, output_folder, clinical_file, threshold, reference, verbose = check_args(argv)
        check_file_format(maf, output_folder)
        cohort_path, clinical_path, snp_path = _create_output_dir(output_folder)
        snp_identifiers, patient_ids = convert_maf_to_vcf(maf, cohort_path, verbose)
        create_snp_identifiers_list(snp_identifiers, snp_path, reference, threshold)
        if clinical_file:
            map_maf_files_to_clinical(patient_ids, clinical_file, clinical_path)
    except TypeError as e:
        print(f"maf2vcf conversion failed: {e}")


if __name__ == "__main__":
    """ Run main method """
    sys.exit(main(sys.argv[1:]))
