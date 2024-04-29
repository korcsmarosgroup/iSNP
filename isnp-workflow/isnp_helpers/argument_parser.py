import argparse
import os
from argparse import RawTextHelpFormatter


class ArgumentParser:

    def __init__(self, navigomix_path, args=None):
        self.navigomix_path = navigomix_path

        help_text = "This tool will execute the iSNP workflow, using the following parameters."

        parser = argparse.ArgumentParser(description=help_text, formatter_class=RawTextHelpFormatter)

        parser.add_argument("-d", "--debug",
                            help="print out debug message",
                            action="store_true")

        parser.add_argument("--no-docker-build",
                            dest="no_docker_build",
                            help="don't build the analytical docker images, but use the current navi-analytics:isnp image",
                            action="store_true")

        parser.add_argument("--only-build-docker",
                            dest="only_build_docker",
                            help="only building the docker images (navi-base and navi-analytics:isnp), but don't run anything else",
                            action="store_true")

        parser.add_argument("--separate",
                            dest="separate",
                            help="execute each analytical step in a separated docker container (default: start a single container and exec each step)",
                            action="store_true")

        default_input_folder_path = os.path.join(navigomix_path, 'doc', 'iSNP-dummy-data')
        parser.add_argument("-i", "--input-folder",
                            help="path to the folder where all the input files are located, it must be under the navigomix git repo\ndefault: " + default_input_folder_path,
                            dest="input_folder",
                            action="store",
                            default=default_input_folder_path)

        default_output_folder_path = os.path.join(navigomix_path, 'output-data')
        parser.add_argument("-o", "--output-folder",
                            help="path to the output folder where all the generated files will be placed, it must be under the navigomix git repo\ndefault: " + default_output_folder_path,
                            dest="output_folder",
                            action="store",
                            default=default_output_folder_path)

        parser.add_argument("--patient_vcf",
                            help="name of the VCF file stores all the input SNP for the patient (default: 0001-patient.vcf)",
                            dest="patient_vcf",
                            action="store",
                            default="0001-patient.vcf")

        parser.add_argument("--disease_snp_list",
                            help="name of the text file stores all the SNP IDs we filter for the disease (default: SNP-identifier.txt)",
                            dest="snp_id_list",
                            action="store",
                            default="uc-snp-identifier.txt")

        parser.add_argument("--promoter_regions_bed",
                            help="name of the annotation file for promoter regions (default: annotation-promoter-regions.bed)",
                            dest="promoter_regions",
                            action="store",
                            default="annotation-promoter-regions-5k-extended-hgnc.bed")

        parser.add_argument("--protein_coding_regions_bed",
                            help="name of the annotation file for protein coding regions (default: annotation-protein-coding-regions.bed)",
                            dest="protein_coding_regions",
                            action="store",
                            default="annotation-protein-coding-regions.bed")

        parser.add_argument("--genome",
                            help="genome in a single fasta file (default: human-genome.fasta)",
                            dest="genome",
                            action="store",
                            default="human-genome.fasta")

        parser.add_argument("--snp_genome_region_radius_protein_coding",
                            help="length of region to fetch from the genome around SNPs in protein coding regions (default: 21, resulting 43 long sequences)",
                            dest="snp_genome_region_radius_protein_coding",
                            action="store",
                            type=int,
                            default=21)

        parser.add_argument("--snp_genome_region_radius_promoter",
                            help="length of region to fetch from the genome around SNPs in promoter regions (default: 50, resulting 101 long sequences)",
                            dest="snp_genome_region_radius_promoter",
                            action="store",
                            type=int,
                            default=100)

        parser.add_argument("--mirna_fasta",
                            help="miRNA sequences in fasta format (default: mirna.fasta)",
                            dest="mirna_fasta",
                            action="store",
                            default="mirna.fasta")

        parser.add_argument("--miranda_score_threshold",
                            help="integer score threshold to use for filtering the strong mirna-gene interactions (default: 90)",
                            dest="miranda_score_threshold",
                            action="store",
                            type=int,
                            default=90)

        parser.add_argument("--miranda_energy_threshold",
                            help="energy threshold to use for filtering the strong mirna-gene interactions (negative integer, default: -20 kcal/mol threshold)",
                            dest="miranda_energy_threshold",
                            action="store",
                            type=int,
                            default=-20)

        parser.add_argument("--tf_binding_matrices",
                            help="matrix file for TF binding simulation (default: jaspar_matrices.txt)",
                            dest="tf_binding_matrices",
                            action="store",
                            default="jaspar_matrices.txt")

        parser.add_argument("--tf_background_rsat",
                            help="background file for the TF binding simulation (default: background.freq)",
                            dest="tf_background_rsat",
                            action="store",
                            default="background_rsat_own5k.freq")

        parser.add_argument("--tf_background_fimo",
                            help="background file for the TF binding simulation (default: fimo.model)",
                            dest="tf_background_fimo",
                            action="store",
                            default="background_fimo_own5k.model")

        parser.add_argument("--tf_score_threshold",
                            help="threshold to use for filtering the strong tf-gene interactions (default: 0.0001)",
                            dest="tf_score_threshold",
                            action="store",
                            type=float,
                            default=0.001)

        parser.add_argument("--id_mapping_json",
                            help="json file name contains ID mapping data (default: uniprot_id_mapping.json)",
                            dest="id_mapping_json",
                            action="store",
                            default="uniprot_id_mapping.json")

        parser.add_argument("--reference_interactions_for_enrichment_tsv",
                            help="mitab file name for interaction data we use for enrichment (default: omnipath.tsv)",
                            dest="reference_interactions_for_enrichment_tsv",
                            action="store",
                            default="omnipath.tsv")

        parser.add_argument("-p", "--patient-files",
                            help="<list of the paths of the specific patient files> [mandatory]",
                            type=str,
                            dest="patient_files",
                            action="store",
                            default="testlist")

        parser.add_argument("-pf", "--patients-folder",
                            help="<folder of the input patient specific VCF files> [mandatory]",
                            type=str,
                            dest="patients_folder",
                            action="store",
                            default="testlist")

        parser.add_argument("-ci", "--counting-index",
                            help="<counting index for the paralell runs> [mandatory]",
                            type=str,
                            dest="counting_index",
                            action="store",
                            default="testlist")

        results = parser.parse_args(args)

        self.debug = results.debug
        self.no_docker_build = results.no_docker_build
        self.only_build_docker = results.only_build_docker
        self.separate = results.separate
        self.input_folder = results.input_folder
        self.output_folder = results.output_folder
        self.patient_vcf = results.patient_vcf
        self.snp_id_list = results.snp_id_list
        self.promoter_regions = results.promoter_regions
        self.protein_coding_regions = results.protein_coding_regions
        self.genome = results.genome
        self.snp_genome_region_radius_protein_coding = results.snp_genome_region_radius_protein_coding
        self.snp_genome_region_radius_promoter = results.snp_genome_region_radius_promoter
        self.mirna_fasta = results.mirna_fasta
        self.miranda_score_threshold = results.miranda_score_threshold
        self.miranda_energy_threshold = results.miranda_energy_threshold
        self.tf_binding_matrices = results.tf_binding_matrices
        self.tf_background_rsat = results.tf_background_rsat
        self.tf_background_fimo = results.tf_background_fimo
        self.tf_score_threshold = results.tf_score_threshold
        self.id_mapping_json_files = list(map(lambda x: x.strip(), results.id_mapping_json.split(',')))
        self.reference_interactions_for_enrichment_tsv = results.reference_interactions_for_enrichment_tsv
        self.patient_files = results.patient_files
        self.patients_folder = results.patients_folder
        self.counting_index = results.counting_index
