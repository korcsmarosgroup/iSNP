from isnp_helpers.argument_parser import ArgumentParser
import subprocess
import logging
import sys
import os


def run_pipeline(params, patient_file, container_name):
    actual_patient = patient_file.split(".")[0].split("/")[-1]

    module_0_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", "/patient_specific_VCF_files/" + patient_file,
                        "--snp", "/input/" + params.snp_id_list,
                        "--output", f"/output/{actual_patient}/disease_filtered.vcf"]
    subprocess.run(module_0_command, check = True)
    
    module_1_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"/output/{actual_patient}/disease_filtered.vcf",
                        "--annotation", "/input/" + params.promoter_regions,
                        "--output", f"/output/{actual_patient}/promoter-regions.vcf"]
    subprocess.run(module_1_command, check = True)

    
    module_2_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"/output/{actual_patient}/disease_filtered.vcf",
                        "--annotation", "/input/" + params.protein_coding_regions,
                        "--output", f"/output/{actual_patient}/protein-coding-regions.vcf"]
    subprocess.run(module_2_command, check = True)
    
    module_3_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"/output/{actual_patient}/protein-coding-regions.vcf",
                        "--genome", "/input/" + params.genome,
                        "--output_wild_type", f"/output/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--output_mutated", f"/output/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_protein_coding)]
    subprocess.run(module_3_command, check = True)
    
    module_4_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"/output/{actual_patient}/promoter-regions.vcf",
                        "--genome", "/input/" + params.genome,
                        "--output_wild_type", f"/output/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--output_mutated", f"/output/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_promoter)]
    subprocess.run(module_4_command, check = True)
    
    module_5_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--genomic", f"/output/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--mirna", "/input/" + params.mirna_fasta,
                        "--output", f"/output/{actual_patient}/mirna_gene_connections_mut.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    subprocess.run(module_5_command, check = True)
    
    module_6_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--genomic", f"/output/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--mirna", "/input/" + params.mirna_fasta,
                        "--output", f"/output/{actual_patient}/mirna_gene_connections_wt.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    subprocess.run(module_6_command, check = True)
    
    module_7_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"/output/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--matrix", "/input/" + params.tf_binding_matrices,
                        "--tf_background_rsat", "/input/" + params.tf_background_rsat,
                        "--tf_background_fimo", "/input/" + params.tf_background_fimo,
                        "--output", f"/output/{actual_patient}/tf_gene_connections_mut.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    subprocess.run(module_7_command, check = True)
    
    module_8_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"/output/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--matrix", "/input/" + params.tf_binding_matrices,
                        "--tf_background_rsat", "/input/" + params.tf_background_rsat,
                        "--tf_background_fimo", "/input/" + params.tf_background_fimo,
                        "--output", f"/output/{actual_patient}/tf_gene_connections_wt.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    subprocess.run(module_8_command, check = True)
    
    module_9_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/mirna_gene_connections_mut.tsv,/output/tf_gene_connections_mut.tsv",
                        "--output-file", f"/output/{actual_patient}/combined_connections_mut.tsv",
                        "--method", "union"]
    subprocess.run(module_9_command, check = True)
    
    module_10_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/mirna_gene_connections_wt.tsv,/output/tf_gene_connections_wt.tsv",
                        "--output-file", f"/output/{actual_patient}/combined_connections_wt.tsv",
                        "--method", "union"]
    subprocess.run(module_10_command, check = True)
    
    module_11_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/combined_connections_mut.tsv,/output/combined_connections_wt.tsv",
                        "--output-file", f"/output/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--method", "difference"]
    subprocess.run(module_11_command, check = True)
    
    module_12_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"/output/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--lower-case",
                        "--no-isoform"]
    subprocess.run(module_12_command, check = True)
    
    module_13_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/network-id-mapper/network_id_mapper.py",
                        "--input", f"/output/{actual_patient}/differences_between_mut_wt_networks_formatted.tsv",
                        "--target-id-type", "uniprotac",
                        # "--molecule-type-filter", "gene",  # do we need this ???
                        # "--remove",
                        "--mapping-data", ",".join(map(lambda x: "/input/" + x, params.id_mapping_json_files)),
                        "--output", f"/output/{actual_patient}/uniprot_differences.tsv"]
    subprocess.run(module_13_command, check = True)
    
    module_14_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/network-enrichment/network_enrichment.py",
                        "--input", f"/output/{actual_patient}/uniprot_differences.tsv",
                        "--output", f"/output/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--reference-net", "/input/" + params.reference_interactions_for_enrichment_tsv,
                        "--distance", "1"]
    subprocess.run(module_14_command, check = True)
    
    module_15_command = ["docker", "exec", f"{container_name}",
                        "python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"/output/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--upper-case",
                        "--no-isoform"]
    subprocess.run(module_15_command, check = True)


def main():
    navigomix_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
    params = ArgumentParser(navigomix_path, sys.argv[1:])

    for patient_file in params.patient_files.split(","):
        run_pipeline(params, patient_file, params.container_name)


if __name__ == "__main__":
    main()
