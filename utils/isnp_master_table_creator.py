import pandas as pd
import argparse
import sys
import os
import re


def parse_args(argv):
    """ Command line interface for the the module """
    help_text = \
        """
        === ISNP Master Table Creator ===
        This script generates a table, which contains all of the TF-target gene and miRNA-target gene interactions per patients.
        """
    parser = argparse.ArgumentParser()

    # Input folder path to learn from
    parser.add_argument("-pf", "--patient_folder",
                        help="<path to the folder where the patient specific folders are> [mandatory]",
                        dest="patient_folder",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)

    return results


def read_in_folder(patient_folder):
    """
    - param folder: patient folder which we use
    - return:

    There should be two different type of interaction.
    # 1: The TF with the uniprot interactions and their target genes
    # 2. The miRNA and their target genes
    """

    # Defining final dictionary with key: patient value: set of SNP affected proteins / interactions
    dictionary_line_out = {}
    dictionary_affected_protein_out = {}
    dictionary_snps_out = {}

    # Defining sets which contain all possible data
    all_affected_proteins = set()
    all_outlines = set()
    all_snps = set()

    subfolders = [f.path for f in os.scandir(patient_folder) if f.is_dir()]

    for folder in subfolders:
        patient = re.findall(r"(\w*).vcf", folder)
        patient = patient[0]
        df = pd.read_csv(os.path.join(folder, "differences.tsv"), sep = "\t", header = 0)
        affected_proteins = set()
        outset_line = set()
        rs_set = set()

        for row in df:

            if row["file"] == "mut" and row["mutated"] == False:
                continue

            else:
                source = row["source"]
                target = row["target"]
                snp = row["snp"]
                rs_set.add(str(snp).strip())
                all_snps.add(str(snp).strip())
                outset_line.add(str(str(source).strip() + "\t" + str(target).strip() + "\t" + str(snp).strip() + "\t" +
                                    str(row["file"]).strip() + "\t" + str(row["tool"]).strip()))
                all_outlines.add(str(str(source).strip() + "\t" + str(target).strip() + "\t" + str(snp).strip() + "\t" +
                                     str(row["file"]).strip() + "\t" + str(row["tool"]).strip()))
                affected_proteins.add(target)
                all_affected_proteins.add(target)

        dictionary_affected_protein_out[patient] = affected_proteins
        dictionary_line_out[patient] = outset_line
        dictionary_snps_out[patient] = rs_set

    final_affected_proteins = {}
    final_outlines = {}
    final_rs = {}

    for patient in dictionary_line_out:
        final_rs[patient] = {}

        for rs in all_snps:

            if rs in dictionary_snps_out[patient]:
                final_rs[patient][rs] = 1

            else:
                final_rs[patient][rs] = 0

        final_affected_proteins[patient] = {}
        for protein in all_affected_proteins:

            if protein in dictionary_affected_protein_out[patient]:
                final_affected_proteins[patient][protein] = 1

            else:
                final_affected_proteins[patient][protein] = 0

        final_outlines[patient] = {}
        for outline in all_outlines:

            if outline in dictionary_line_out[patient]:
                final_outlines[patient][outline] = 1

            else:
                final_outlines[patient][outline] = 0

    # Creating dataframes from the dictionaries
    proteins_df = pd.DataFrame.from_dict(final_affected_proteins)
    outline_df = pd.DataFrame.from_dict(final_outlines)
    rs_df = pd.DataFrame.from_dict(final_rs)

    # Indexing the dataframes
    proteins_df.index = proteins_df.index.str.upper()
    outline_df.index = outline_df.index.str.upper()
    rs_df.index = rs_df.index.str.upper()

    # Set the names of the axises
    proteins_df = proteins_df.rename_axis("id")
    outline_df = outline_df.rename_axis("Source\tTarget\tSNP\tMutated\tInteraction_source")
    rs_df = rs_df.rename_axis("SNP")

    # Write out the results to the specific output files
    proteins_df.to_csv(os.path.join(patient_folder, "affected_proteins.tsv"), sep = "\t")
    outline_df.to_csv(os.path.join(patient_folder, "affected_proteins_TFs_mirs.tsv"), sep = "\t")
    rs_df.to_csv(os.path.join(patient_folder, "SNPs.tsv"), sep = "\t")


def main(argv):
    """Main function for the script"""
    args = parse_args(argv)
    read_in_folder(args.patient_folder)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
