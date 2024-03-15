from isnp_helpers.argument_parser import ArgumentParser
import subprocess
import sys
import os


def run_pipeline(params, input_folder, output_folder, patient_folder, patient_file):
    image_name = "navi-analytics:isnp"
    actual_patient = patient_file.split(".")[0].split("/")[-1]
    container_name = actual_patient

    module_0_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_0", f"{image_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", "/patient_specific_VCF_files/" + patient_file,
                        "--snp", "/input/" + params.snp_id_list,
                        "--output", f"/output/{actual_patient}/disease_filtered.vcf"]
    print(f"\n\n1/16 ======= running analytical task in a single docker container with command: {module_0_command}")
    subprocess.run(module_0_command, check = True)
    
    module_1_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_1", f"{image_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"/output/{actual_patient}/disease_filtered.vcf",
                        "--annotation", "/input/" + params.promoter_regions,
                        "--output", f"/output/{actual_patient}/promoter-regions.vcf"]
    print(f"\n\n2/16 ======= running analytical task in a single docker container with command: {module_1_command}")
    subprocess.run(module_1_command, check = True)

    module_2_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_2", f"{image_name}",
                        "python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                        "--input", f"/output/{actual_patient}/disease_filtered.vcf",
                        "--annotation", "/input/" + params.protein_coding_regions,
                        "--output", f"/output/{actual_patient}/protein-coding-regions.vcf"]
    print(f"\n\n3/16 ======= running analytical task in a single docker container with command: {module_2_command}")
    subprocess.run(module_2_command, check = True)
    
    module_3_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_3", f"{image_name}",
                        "python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"/output/{actual_patient}/protein-coding-regions.vcf",
                        "--genome", "/input/" + params.genome,
                        "--output_wild_type", f"/output/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--output_mutated", f"/output/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_protein_coding)]
    print(f"\n\n4/16 ======= running analytical task in a single docker container with command: {module_3_command}")
    subprocess.run(module_3_command, check = True)
    
    module_4_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_4", f"{image_name}",
                        "python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                        "--input_vcf", f"/output/{actual_patient}/promoter-regions.vcf",
                        "--genome", "/input/" + params.genome,
                        "--output_wild_type", f"/output/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--output_mutated", f"/output/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--region_length", str(params.snp_genome_region_radius_promoter)]
    print(f"\n\n5/16 ======= running analytical task in a single docker container with command: {module_4_command}")
    subprocess.run(module_4_command, check = True)
    
    module_5_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_5", f"{image_name}",
                        "python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--mirna", f"/output/{actual_patient}/snp_in_protein-coding-regions_mut.fasta",
                        "--genomic", "/input/" + params.mirna_fasta,
                        "--output", f"/output/{actual_patient}/mirna_gene_connections_mut.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    print(f"\n\n6/16 ======= running analytical task in a single docker container with command: {module_5_command}")
    subprocess.run(module_5_command, check = True)
    
    module_6_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_6", f"{image_name}",
                        "python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                        "--mirna", f"/output/{actual_patient}/snp_in_protein-coding-regions_wt.fasta",
                        "--genomic", "/input/" + params.mirna_fasta,
                        "--output", f"/output/{actual_patient}/mirna_gene_connections_wt.tsv",
                        "--score", str(params.miranda_score_threshold),
                        "--energy", str(params.miranda_energy_threshold)]
    print(f"\n\n7/16 ======= running analytical task in a single docker container with command: {module_6_command}")
    subprocess.run(module_6_command, check = True)
    
    module_7_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_7", f"{image_name}",
                        "python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"/output/{actual_patient}/snp_in_promoter-regions_mut.fasta",
                        "--matrix", "/input/" + params.tf_binding_matrices,
                        "--tf_background_rsat", "/input/" + params.tf_background_rsat,
                        "--tf_background_fimo", "/input/" + params.tf_background_fimo,
                        "--output", f"/output/{actual_patient}/tf_gene_connections_mut.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    print(f"\n\n8/16 ======= running analytical task in a single docker container with command: {module_7_command}")
    subprocess.run(module_7_command, check = True)
    
    module_8_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_8", f"{image_name}",
                        "python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                        "--fasta", f"/output/{actual_patient}/snp_in_promoter-regions_wt.fasta",
                        "--matrix", "/input/" + params.tf_binding_matrices,
                        "--tf_background_rsat", "/input/" + params.tf_background_rsat,
                        "--tf_background_fimo", "/input/" + params.tf_background_fimo,
                        "--output", f"/output/{actual_patient}/tf_gene_connections_wt.tsv",
                        "--format", "transfac",
                        "--threshold", str(params.tf_score_threshold)]
    print(f"\n\n9/16 ======= running analytical task in a single docker container with command: {module_8_command}")
    subprocess.run(module_8_command, check = True)
    
    module_9_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_9", f"{image_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/mirna_gene_connections_mut.tsv,/output/tf_gene_connections_mut.tsv",
                        "--output-file", f"/output/{actual_patient}/combined_connections_mut.tsv",
                        "--method", "union"]
    print(f"\n\n10/16 ======= running analytical task in a single docker container with command: {module_9_command}")
    subprocess.run(module_9_command, check = True)
    
    module_10_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_10", f"{image_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/mirna_gene_connections_wt.tsv,/output/tf_gene_connections_wt.tsv",
                        "--output-file", f"/output/{actual_patient}/combined_connections_wt.tsv",
                        "--method", "union"]
    print(f"\n\n11/16 ======= running analytical task in a single docker container with command: {module_10_command}")
    subprocess.run(module_10_command, check = True)
    
    module_11_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_11", f"{image_name}",
                        "python3", "/analytic-modules/network-combiner/network_combiner.py",
                        "--input-files", f"/output/{actual_patient}/combined_connections_mut.tsv,/output/combined_connections_wt.tsv",
                        "--output-file", f"/output/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--method", "difference"]
    print(f"\n\n12/16 ======= running analytical task in a single docker container with command: {module_11_command}")
    subprocess.run(module_11_command, check = True)
    
    module_12_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_12", f"{image_name}",
                        "python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"/output/{actual_patient}/differences_between_mut_wt_networks.tsv",
                        "--lower-case",
                        "--no-isoform"]
    print(f"\n\n13/16 ======= running analytical task in a single docker container with command: {module_12_command}")
    subprocess.run(module_12_command, check = True)
    
    module_13_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_13", f"{image_name}",
                        "python3", "/analytic-modules/network-id-mapper/network_id_mapper.py",
                        "--input", f"/output/{actual_patient}/differences_between_mut_wt_networks_formatted.tsv",
                        "--target-id-type", "uniprotac",
                        # "--molecule-type-filter", "gene",  # do we need this ???
                        # "--remove",
                        "--mapping-data", ",".join(map(lambda x: "/input/" + x, params.id_mapping_json_files)),
                        "--output", f"/output/{actual_patient}/uniprot_differences.tsv"]
    print(f"\n\n14/16 ======= running analytical task in a single docker container with command: {module_13_command}")
    subprocess.run(module_13_command, check = True)
    
    module_14_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_14", f"{image_name}",
                        "python3", "/analytic-modules/network-enrichment/network_enrichment.py",
                        "--input", f"/output/{actual_patient}/uniprot_differences.tsv",
                        "--output", f"/output/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--reference-net", "/input/" + params.reference_interactions_for_enrichment_tsv,
                        "--distance", "1"]
    print(f"\n\n15/16 ======= running analytical task in a single docker container with command: {module_14_command}")
    subprocess.run(module_14_command, check = True)
    
    module_15_command = ["podman", "run", "--rm", "-v", f"{input_folder}:/input", "-v", f"{output_folder}:/output", "-v" , f"{patient_folder}:/patient_specific_VCF_files",
                        "--name", f"{container_name}_15", f"{image_name}",
                        "python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                        "--input-network-file", f"/output/{actual_patient}/enriched_uniprot_differences.tsv",
                        "--upper-case",
                        "--no-isoform"]
    print(f"\n\n16/16 ======= running analytical task in a single docker container with command: {module_15_command}")
    subprocess.run(module_15_command, check = True)


def main():
    navigomix_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
    params = ArgumentParser(navigomix_path, sys.argv[1:])

    for patient_file in params.patient_files.split(","):
        run_pipeline(params, params.input_folder, params.output_folder, params.patients_folder, patient_file)


if __name__ == "__main__":
    main()
