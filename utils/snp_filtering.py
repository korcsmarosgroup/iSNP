from typing import List
import pandas as pd
import argparse
import sys


def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()

    # Input file path to learn from
    parser.add_argument("-i", "--inputfile",
                        help="<Inputfile of the SNPs and interactions> [madatory]",
                        dest="inputfile",
                        action="store",
                        required=True)

    parser.add_argument("-a", "--SNP_annotation_file",
                        help="<path to the SNP annotation file or SNP list where we think there are miRNAs> [optional]",
                        dest="SNP_annotation_file",
                        action="store",
                        required=False)

    parser.add_argument("-o", "--outfile",
                        help="<path to the outfile name. The script makes three of them.> [mandatory]",
                        dest="outfile",
                        action="store",
                        required=True)

    parser.add_argument("-t", "--type",
                        help= "<Analysis type by default TF only.\n"
                              "ga : gene annotation file. Searching only in exons.\n"
                              "tf: TF only\n"
                              "snp: SNP list as input for miRNAs> [optional]",
                        dest="analysis_type",
                        action="store",
                        required=False)

    parser.add_argument("-log", "--double_snp_out_file",
                        help="<Log file for the double store RNAs> [optional]",
                        dest="double_snp_out_file",
                        action="store",
                        required=False)

    results = parser.parse_args(argv)
    return results


def selecting_SNP_affected_proteins(SNP_annotation_file, double_snp_out_file, affected_protein_interactions, outfile):
    genome_coordinates = pd.read_csv(SNP_annotation_file, header = 0, sep = "\t")
    genome_coordinates.head()

    """ 
    From this table we need only two columns: The SNP ID and whether certain SNP has a cDNA. 
    We can have problem if the SNP has a gene at the positive and at the negative strand. First better to check that.
    We first map every SNP to uniqe ensmeble gene ids to avoid the transcript variants.
    Then we group and sum the transcript chromosome strands. If they are 0 then it means that there is a positive (1) 
    and a negative (-1) strand. 
    We write them out with the gene information for check most of them are LRGs (Locus Reference Genomic). They are 
    references of the genome. 
    You need to change the final file and delete the interaction if the SNP only in one of the stranded genes.
    """
    only_one_splice_per_gene = genome_coordinates[["refsnp_id", "ensembl_transcript_chrom_strand", "ensembl_gene_stable_id"]].groupby(by = "ensembl_gene_stable_id").first()

    snps_id_splice = only_one_splice_per_gene.groupby(by = "refsnp_id").sum()
    double_gene_snps = snps_id_splice[snps_id_splice["ensembl_transcript_chrom_strand"] == 0]
    genome_coordinates.merge(double_gene_snps, on="refsnp_id").to_csv(double_snp_out_file, sep = "\t")

    important_information = genome_coordinates[["refsnp_id", "cdna_start"]]
    potetially_affected_miRNA_SNPs = important_information[pd.notna(important_information["cdna_start"])]["refsnp_id"].unique()
    snps = list(potetially_affected_miRNA_SNPs)
    snps_upper = [x.upper() for x in snps]

    TF_miRNA_interactions = pd.read_csv(affected_protein_interactions, sep = "\t")
    TF_miRNA_interactions.head()

    filtered_interactions = TF_miRNA_interactions[TF_miRNA_interactions.Interaction_source.isin(["RSAT", "FIMO"]) |
                                                  (TF_miRNA_interactions.Interaction_source.isin(["MIRANDA"]) &
                                                   TF_miRNA_interactions.SNP.isin(snps_upper))]

    filtered_interactions.to_csv(outfile, sep = "\t")
    patients = list(filtered_interactions.columns)[5:]

    patients_plus_target = ["Target"] + patients
    affected_proteins = filtered_interactions[patients_plus_target].groupby(by=["Target"]).max()
    affected_proteins.to_csv(outfile.replace(".txt", "affected_proteins.txt"), sep = "\t")

    patients_plus_snps = ["SNP"] + patients
    affecting_snps = filtered_interactions[patients_plus_snps].groupby(by=["SNP"]).max()
    affecting_snps.to_csv(outfile.replace(".txt", "SNPs.txt"), sep = "\t")


def selecting_SNP_affected_proteins_tf_only(affected_protein_interactions, outfile):
    TF_miRNA_interactions = pd.read_csv(affected_protein_interactions, sep = "\t")
    TF_miRNA_interactions.head()

    filtered_interactions = TF_miRNA_interactions[TF_miRNA_interactions.Interaction_source.isin(["RSAT", "FIMO"])]
    filtered_interactions.to_csv(outfile, sep = "\t")
    patients = list(filtered_interactions.columns)[5:]

    patients_plus_target = ["Target"] + patients
    affected_proteins = filtered_interactions[patients_plus_target].groupby(by=["Target"]).max()
    affected_proteins = affected_proteins.rename_axis("id")
    affected_proteins.to_csv(outfile.replace(".txt", "affected_proteins.txt"), sep = "\t")

    patients_plus_snps = ["SNP"] + patients
    affecting_snps = filtered_interactions[patients_plus_snps].groupby(by=["SNP"]).max()
    affecting_snps = affecting_snps.rename_axis("SNP")
    affecting_snps.to_csv(outfile.replace(".txt", "SNPs.txt"), sep = "\t")


def selecting_SNP_affected_proteins_snp_list(affected_protein_interactions, snp_list_file, outfile):
    TF_miRNA_interactions = pd.read_csv(affected_protein_interactions, sep = "\t")
    snp_list = []
    file_ = open(snp_list_file, "r")

    for line in file_:
        line=line.strip()
        line=line.split("\t")
        snp_list.append(line[0])

    file_.close()

    snp_list = [x.upper() for x in snp_list]

    filtered_interactions = TF_miRNA_interactions[TF_miRNA_interactions.Interaction_source.isin(["RSAT", "FIMO"]) |
                                                  (TF_miRNA_interactions.Interaction_source.isin(["MIRANDA"]) &
                                                   TF_miRNA_interactions.SNP.isin(snp_list))]

    filtered_interactions.to_csv(outfile, sep = "\t")
    patients = list(filtered_interactions.columns)[5:]

    patients_plus_target = ["Target"] + patients
    affected_proteins = filtered_interactions[patients_plus_target].groupby(by=["Target"]).max()
    affected_proteins = affected_proteins.rename_axis("id")
    affected_proteins.to_csv(outfile.replace(".txt", "affected_proteins.txt"), sep = "\t")

    patients_plus_snps = ["SNP"] + patients
    affecting_snps = filtered_interactions[patients_plus_snps].groupby(by=["SNP"]).max()
    affecting_snps = affecting_snps.rename_axis("SNP")
    affecting_snps.to_csv(outfile.replace(".txt", "SNPs.txt"), sep = "\t")


def main(argv):
    args = parse_args(argv)
    affected_protein_interactions = args.inputfile
    outfile = args.outfile

    if args.double_snp_out_file:
        double_snp_out_file = args.double_snp_out_file

    else:
        double_snp_out_file = "log.txt"

    a = outfile.split(".")

    if len(a) > 2:
        print("The path has some problems.")
        raise

    try:
        out_file_extension = a[1]
        if out_file_extension != "tsv":
            print("The output format needs to be tsv")

    except:
        print("The file has no extension")
        raise

    if args.analysis_type != "snp" and args.analysis_type != "ga":
        print ("Transcription Factor based filtering. miRNA-TS connections are ignored. ")
        selecting_SNP_affected_proteins_tf_only(affected_protein_interactions, outfile)

    elif args.analysis_type == "ga":
        try:
            annotation = args.SNP_annotation_file
        except:
            print ("No annotation file given")
            raise
        print("mirNA annotations are considered in the gene annotation file. Please check the annotation file if there are errors.")
        selecting_SNP_affected_proteins(annotation, double_snp_out_file, affected_protein_interactions, outfile)

    elif args.analysis_type == "snp":
        try:
            annotation = args.SNP_annotation_file
        except:
            print ("No annotation file given")
            raise
        print ("The miRNAs are in a prefiltered SNP list. The SNP list must be tab delimited file. \n"
               "The first column is considered as the SNP list.")
        selecting_SNP_affected_proteins_snp_list(affected_protein_interactions, annotation, outfile)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
