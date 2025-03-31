
import pandas as pd
import argparse
import sys
import os
import re


def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()

    # Input file path to learn from
    parser.add_argument("-pf", "--patient_folder",
                        help="<path to the folder where the patient specific folders are> [mandatory]",
                        dest="patient_folder",
                        action="store",
                        required=True)

    parser.add_argument("-of", "--output_file_suffix",
                        help="<output_file_suffix>, [optional]",
                        default="",
                        dest="output_file_suffix",
                        required=False)

    parser.add_argument("-asnps","--accessible_snp_file_list",
                        help="<path to the SNP list files> [optional]",
                        dest="accessible_snp_file_list",
                        action="store",
                        required=False,
                        nargs="*")

    parser.add_argument("-tfsnps", "--transcription_factor_snps",
                        help="<path to the transcription factor specifc SNPs> [optional]",
                        dest="transcription_factor_snps",
                        action= "store",
                        required=False)

    parser.add_argument("-mirnasnps", "--mirna_snp_file",
                        help="<path to the mirNa sNP files> [optional]",
                        dest = "mirna_snp_file",
                        action="store",
                        required=False)

    results = parser.parse_args(argv)
    return results


def read_in_folder(
        patient_folder,
        accessible_snp_file_list = False,
        transcription_factor_snps = False,
        mirna_snp_file = False,
        suffix_ = ""
    ):
    """
    :param folder: patient folder which we use
    :return:
    There should be two different type of interaction.
    # 1: The TF with the uniprot interactions and their target genes.
    # 2. The miRNA and the target genes
    # Theese should be joint and will be good.
    """

    # Defining final dictionary with key: patient value: set of SNP affected proteins / interactions
    dictionary_line_out = {}
    dictionary_affected_protein_out = {}
    dictionary_snps_out = {}

    # Defining sets which contain all possible data
    all_affected_proteins = set()
    all_outlines = set()
    all_snps = set()

    # If we have filtering SNP list files we use them to work with
    # Defining accessible SNPs from the SNP file list
    accessible_snps = set()

    if not accessible_snp_file_list == False:

        for accessible_snp_file in accessible_snp_file_list:

            with open(accessible_snp_file, 'r') as inp:

                for line in inp:
                    line = line.strip()
                    accessible_snps.add(line)

    # Secondary check for miRNA SNPs
    tf_snps =set()

    if transcription_factor_snps:

        with open (transcription_factor_snps, "r") as inp:

            for line in inp:
                line= line.strip()
                tf_snps.add(line)

    # Secondary check for TF SNPs
    mirna_snps = set()

    if mirna_snp_file:

        with open(mirna_snp_file, "r") as inp:

            for line in inp:
                line = line.strip()
                mirna_snps.add(line)

    subfolders = [f.path for f in os.scandir(patient_folder) if f.is_dir()]

    for folder in subfolders:

        patient = re.findall(r"(\w*).vcf", folder)
        print(patient)

        if len(patient) > 0:
            patient = patient[0]
            df = pd.read_csv(os.path.join(folder, "differences.tsv"), sep="\t", header=0)
            outset_line = set()
            rs_set = set()
            affected_proteins = set()

            for index, row in df.iterrows():

                if row["file"] == "mut" and row["mutated"] == False:
                    print ("This is what I am talking about!")
                    print (row)

                else:
                    source = row["source"].strip()
                    target = row["target"].strip()
                    snp = row["snp"].strip()
                    tool = str(row["tool"]).strip()

                    if accessible_snp_file_list == False: # If we do not have filtering parameter then we keep every SNP.
                        accessible_snps.add(snp)

                    if transcription_factor_snps == False:
                        tf_snps.add(snp)

                    if mirna_snp_file == False:
                        mirna_snps.add(snp)

                    if snp in accessible_snps and ((tool == "miranda" and snp in mirna_snps) or
                                                   ((tool == "rsat" or tool == "fimo") and snp in tf_snps)):

                        rs_set.add(str(snp).strip())
                        all_snps.add(str(snp).strip())
                        outset_line.add(str(str(source).strip() + "\t" + str(target).strip() + "\t" + str(snp).strip() + "\t" +
                                            str(row["file"]).strip() + "\t" + str(row["tool"]).strip()))
                        all_outlines.add(str(str(source).strip() + "\t" + str(target).strip() + "\t" + str(snp).strip() + "\t" +
                                             str(row["file"]).strip() + "\t" + str(row["tool"]).strip()))
                        affected_proteins.add(target)
                        all_affected_proteins.add(target)
            dictionary_line_out[patient] = outset_line
            dictionary_affected_protein_out[patient] = affected_proteins
            dictionary_snps_out[patient] = rs_set

    final_rs = {}
    final_affected_proteins = {}
    final_outlines = {}

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

    outline_df = pd.DataFrame.from_dict(final_outlines)
    proteins_df = pd.DataFrame.from_dict(final_affected_proteins)
    rs_df = pd.DataFrame.from_dict(final_rs)

    outline_df.index = outline_df.index.str.upper()
    proteins_df.index = proteins_df.index.str.upper()
    rs_df.index = rs_df.index.str.upper()

    outline_df = outline_df.rename_axis("Source\tTarget\tSNP\tMutated\tInteraction_source")
    proteins_df = proteins_df.rename_axis("id")
    rs_df = rs_df.rename_axis("SNP")

    # pandas recognises from the csv module the quoting. It requires the escape character
    outline_df.to_csv(os.path.join(patient_folder, "affected_proteins_TFs_mirs" + suffix_ + ".txt"), sep = '\t')
    inp = open(os.path.join(patient_folder,"affected_proteins_TFs_mirs" + suffix_ + ".txt"))
    text = inp.read()
    text = text.replace('"', '')
    out = open(os.path.join(patient_folder, "affected_proteins_TFs_mirs" + suffix_ + ".txt"), "w")
    out.write(text)
    out.close()
    proteins_df.to_csv(os.path.join(patient_folder, "affected_proteins" + suffix_ + ".txt"), sep = '\t')
    rs_df.to_csv(os.path.join(patient_folder, "SNPs" + suffix_ + ".txt"), sep = '\t')


def main(argv):
    args = parse_args(argv)

    if args.accessible_snp_file_list is not None:

        try:
            accessible_snp_file_list = list(args.accessible_snp_file_list)
            print(accessible_snp_file_list)

            if len(accessible_snp_file_list[0]) == 1:
                raise FileNotFoundError("The SNP file is not correctly given. Please give it as a python list.")

        except:
            raise FileNotFoundError("The SNP file is not correctly given. Please give it as a python list.")

    else:
        accessible_snp_file_list = False

    if args.transcription_factor_snps is not None:
        transcription_factor_snps = args.transcription_factor_snps

    else:
        transcription_factor_snps = False

    if args.mirna_snp_file is not None:
        mirna_snp_file = args.mirna_snp_file

    else:
        mirna_snp_file = False

    read_in_folder(
        args.patient_folder,
        accessible_snp_file_list,
        transcription_factor_snps,
        mirna_snp_file,
        suffix_ = args.output_file_suffix
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
