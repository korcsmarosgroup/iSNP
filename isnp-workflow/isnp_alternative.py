from isnp_helpers.argument_parser import ArgumentParser
from time import strftime
import subprocess
import logging
import shutil
import sys
import os


def run_pipeline(params, input_folder, output_folder, patient_folder, patient_file):
    actual_patient = patient_file.split(".")[0]
    actual_patient_folder = os.path.join(output_folder, actual_patient)

    if os.path.isdir(actual_patient_folder):
        shutil.rmtree(actual_patient_folder)

    if not os.path.isdir(actual_patient_folder):
        os.mkdir(actual_patient_folder)

    actual_patient_log_file_name = f"{actual_patient}.log"
    actual_patient_log_file = os.path.join(actual_patient_folder, actual_patient_log_file_name)

    if os.path.isfile(actual_patient_log_file):
        os.remove(actual_patient_log_file)

    logging.basicConfig(filename = actual_patient_log_file, level = logging.INFO)
    logging.info(f"### [{strftime('%H:%M:%S')}] Starting the pipeline on the patient: {actual_patient}")

    module_0_command = ["python3", "../analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"{patient_folder}" + patient_file,
                        "--snp", f"{input_folder}" + params.snp_id_list,
                        "--output", f"{output_folder}/{actual_patient}/disease_filtered.vcf"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 1/16 ======= running analytical task with command: {module_0_command}")
    subprocess.run(module_0_command, check = True)

    module_1_command = ["python3", "../analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"{output_folder}/{actual_patient}/disease_filtered.vcf",
                        "--annotation", f"{input_folder}/" + params.promoter_regions,
                        "--output", f"{output_folder}/{actual_patient}/promoter-regions.vcf"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 2/16 ======= running analytical task with command: {module_1_command}")
    subprocess.run(module_1_command, check = True)

    module_2_command = ["python3", "../analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"{output_folder}/{actual_patient}/disease_filtered.vcf",
                        "--annotation", f"{input_folder}/" + params.protein_coding_regions,
                        "--output", f"{output_folder}/{actual_patient}/protein-coding-regions.vcf"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 3/16 ======= running analytical task with command: {module_2_command}")
    subprocess.run(module_2_command, check = True)

    module_3_command = ["python3", "../analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"{output_folder}/{actual_patient}/protein-coding-regions.vcf",
                        "--genome", f"{input_folder}/" + params.genome,
                        "--output_wild_type", f"{output_folder}/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--output_mutated", f"{output_folder}/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_protein_coding)]
    logging.info(f"### [{strftime('%H:%M:%S')}] 4/16 ======= running analytical task with command: {module_3_command}")
    subprocess.run(module_3_command, check = True)

    module_4_command = ["python3", "../analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"{output_folder}/{actual_patient}/promoter-regions.vcf",
                        "--genome", f"{input_folder}/" + params.genome,
                        "--output_wild_type", f"{output_folder}/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--output_mutated", f"{output_folder}/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_promoter)]
    logging.info(f"### [{strftime('%H:%M:%S')}] 5/16 ======= running analytical task with command: {module_4_command}")
    subprocess.run(module_4_command, check = True)

    module_5_command = ["python3", "../analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--mirna", f"{output_folder}/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--genomic", f"{input_folder}/" + params.mirna_fasta,
                        "--output", f"{output_folder}/{actual_patient}/mirna_gene_connections_mut.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    # logging.info(f"### [{strftime('%H:%M:%S')}] 6/16 ======= running analytical task with command: {module_5_command}")
    # subprocess.run(module_5_command, check = True)
    logging.info(f"### [{strftime('%H:%M:%S')}] 6/16 ======= Skipping the miranda module with the mutant region")

    module_6_command = ["python3", "../analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--mirna", f"{output_folder}/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--genomic", f"{input_folder}/" + params.mirna_fasta,
                        "--output", f"{output_folder}/{actual_patient}/mirna_gene_connections_wt.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    # logging.info(f"### [{strftime('%H:%M:%S')}] 7/16 ======= running analytical task with command: {module_6_command}")
    # subprocess.run(module_6_command, check = True)
    logging.info(f"### [{strftime('%H:%M:%S')}] 7/16 ======= Skipping the miranda module with the wild region")

    module_7_command = ["python3", "../analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"{output_folder}/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--matrix", f"{input_folder}/" + params.tf_binding_matrices,
                        "--tf_background_rsat", f"{input_folder}/" + params.tf_background_rsat,
                        "--tf_background_fimo", f"{input_folder}/" + params.tf_background_fimo,
                        "--output", f"{output_folder}/{actual_patient}/tf_gene_connections_mut.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    logging.info(f"### [{strftime('%H:%M:%S')}] 8/16 ======= running analytical task with command: {module_7_command}")
    subprocess.run(module_7_command, check = True)

    module_8_command = ["python3", "../analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"{output_folder}/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--matrix", f"{input_folder}/" + params.tf_binding_matrices,
                        "--tf_background_rsat", f"{input_folder}/" + params.tf_background_rsat,
                        "--tf_background_fimo", f"{input_folder}/" + params.tf_background_fimo,
                        "--output", f"{output_folder}/{actual_patient}/tf_gene_connections_wt.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    logging.info(f"### [{strftime('%H:%M:%S')}] 9/16 ======= running analytical task with command: {module_8_command}")
    subprocess.run(module_8_command, check = True)

    module_9_command = ["python3", "../analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"{output_folder}/{actual_patient}/mirna_gene_connections_mut.tsv,{output_folder}/{actual_patient}/tf_gene_connections_mut.tsv",
                        "--output-file", f"{output_folder}/{actual_patient}/combined_connections_mut.tsv",
                        "--method", "union"]
    # logging.info(f"### [{strftime('%H:%M:%S')}] 10/16 ======= running analytical task with command: {module_1_command}")
    # subprocess.run(module_9_command, check = True)
    logging.info(f"### [{strftime('%H:%M:%S')}] 10/16 ======= Skipping this network combiner module, because miranda is missing")

    module_10_command = ["python3", "../analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"{output_folder}/{actual_patient}/mirna_gene_connections_wt.tsv,{output_folder}/{actual_patient}/tf_gene_connections_wt.tsv",
                        "--output-file", f"{output_folder}/{actual_patient}/combined_connections_wt.tsv",
                        "--method", "union"]
    # logging.info(f"### [{strftime('%H:%M:%S')}] 11/16 ======= running analytical task with command: {module_1_command}")
    # subprocess.run(module_10_command, check = True)
    logging.info(f"### [{strftime('%H:%M:%S')}] 11/16 ======= Skipping this network combiner module, because miranda is missing")

    # module_11_command = ["python3", "../analytic-modules/network-combiner/network_combiner.py",
    #                     "--input-files", f"{output_folder}/{actual_patient}/combined_connections_mut.tsv,{output_folder}/{actual_patient}/combined_connections_wt.tsv",
    #                     "--output-file", f"{output_folder}/{actual_patient}/differences_between_mut_wt_networks.tsv",
    #                     "--method", "difference"]
    # logging.info(f"### [{strftime('%H:%M:%S')}] 12/16 ======= running analytical task with command: {module_1_command}")
    # subprocess.run(module_11_command, check = True)

    module_11_command = ["python3", "../analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"{output_folder}/{actual_patient}/tf_gene_connections_mut.tsv,{output_folder}/{actual_patient}/tf_gene_connections_wt.tsv",
                        "--output-file", f"{output_folder}/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--method", "difference"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 12/16 ======= running analytical task with command: {module_11_command}")
    subprocess.run(module_11_command, check = True)

    module_12_command = ["python3", "../analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"{output_folder}/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--lower-case",
                        "--no-isoform"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 13/16 ======= running analytical task with command: {module_12_command}")
    subprocess.run(module_12_command, check = True)

    module_13_command = ["python3", "../analytic-modules/network-id-mapper/network_id_mapper.py",
                        "--input", f"{output_folder}/{actual_patient}/differences_between_mut_wt_networks_formatted.tsv",
                        "--target-id-type", "uniprotac",
                        # "--molecule-type-filter", "gene",  # do we need this ???
                        # "--remove",
                        "--mapping-data", ",".join(map(lambda x: f"{input_folder}/" + x, params.id_mapping_json_files)),
                        "--output", f"{output_folder}/{actual_patient}/uniprot_differences.tsv"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 14/16 ======= running analytical task with command: {module_13_command}")
    subprocess.run(module_13_command, check = True)

    module_14_command = ["python3", "../analytic-modules/network-enrichment/network_enrichment.py",
                        "--input", f"{output_folder}/{actual_patient}/uniprot_differences.tsv",
                        "--output", f"{output_folder}/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--reference-net", f"{input_folder}/" + params.reference_interactions_for_enrichment_tsv,
                        "--distance", "1"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 15/16 ======= running analytical task with command: {module_14_command}")
    subprocess.run(module_14_command, check = True)

    module_15_command = ["python3", "../analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"{output_folder}/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--upper-case",
                        "--no-isoform"]
    logging.info(f"### [{strftime('%H:%M:%S')}] 16/16 ======= running analytical task with command: {module_15_command}")
    subprocess.run(module_15_command, check = True)

    logging.info(f"### [{strftime('%H:%M:%S')}] Finished on the patient: {actual_patient}")


def main():
    navigomix_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
    params = ArgumentParser(navigomix_path, sys.argv[1:])

    for patient_file in params.patient_files.split(","):
        run_pipeline(params, params.input_folder, params.output_folder, params.patients_folder, patient_file)


if __name__ == "__main__":
    main()
