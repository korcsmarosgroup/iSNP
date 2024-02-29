#!/usr/bin/env python3

import os
import shutil
import signal
import sys
from curses import wrapper

from isnp_helpers.argument_parser import ArgumentParser
from isnp_helpers.display import IsnpDisplay
from isnp_helpers.docker_helper import DockerHelper

image_name_with_tag = 'navi-analytics:isnp'
docker_helper = None


def build_analytics_docker_image(docker_helper, navigomix_path):
    # first copy all analytic module source code to the docker build dir:
    abs_path_src = os.path.join(navigomix_path, 'analytic-modules')
    abs_path_trg = os.path.join(navigomix_path, 'deploy/containers/navi-analytics/analytic-modules')
    shutil.rmtree(abs_path_trg, ignore_errors=True)
    shutil.copytree(abs_path_src, abs_path_trg)

    # now building the base docker image
    docker_helper.build_docker(os.path.join(navigomix_path, 'deploy/containers/navi-base'), 'navi-base:1.0.0')

    # then build the analytic container with tag 'isnp'
    docker_helper.build_docker(os.path.join(navigomix_path, 'deploy/containers/navi-analytics'), 'navi-analytics:isnp')


def execute_command(docker_helper, command_id, display, args):
    display.print("running modules {:05.2f}% ({})".format(float(command_id) / 17 * 100, args[1]))
    display.set_running_status(command_id)
    exit_code = docker_helper.execute_analytical_task(args)
    if exit_code == 0:
        display.set_success_status(command_id)
    else:
        display.set_error_status(command_id)
        display.print("analytical module execution failed! use --debug option to see the details")
        docker_helper.kill_running_container()
        display.exit(exit_code)


def main(stdscr, params, navigomix_path):
    global docker_helper
    display = IsnpDisplay(stdscr, params.debug)
    docker_helper = DockerHelper(navigomix_path, display, params.input_folder, params.output_folder, image_name_with_tag, not params.debug, params.separate)

    if not params.no_docker_build:
        display.print("building docker images...")
        build_analytics_docker_image(docker_helper, navigomix_path)
        display.print("docker images built successfully :)")

    if params.only_build_docker:
        display.print("'--only-build-docker' parameter was used, exiting now")
        display.exit(0)

    shutil.rmtree(params.output_folder, ignore_errors=True)
    os.makedirs(params.output_folder)

    if not params.separate:
        docker_helper.start_long_term_container()
        # time.sleep(1000)

    execute_command(docker_helper, 0, display,
                    ["python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                     "--input", params.input_folder + params.patient_vcf,
                     "--snp", params.input_folder + params.snp_id_list,
                     "--output", params.output_folder + "disease_filtered.vcf"])

    execute_command(docker_helper, 1, display,
                    ["python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                     "--input", params.output_folder + "disease_filtered.vcf",
                     "--annotation", params.input_folder + params.promoter_regions,
                     "--output", params.output_folder + "promoter-regions.vcf"])

    execute_command(docker_helper, 2, display,
                    ["python3", "/analytic-modules/vcf-filtering/vcf_filter.py",
                     "--input", params.output_folder + "disease_filtered.vcf",
                     "--annotation", params.input_folder + params.protein_coding_regions,
                     "--output", params.output_folder + "protein-coding-regions.vcf"])

    execute_command(docker_helper, 3, display,
                    ["python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                     "--input_vcf", params.output_folder + "protein-coding-regions.vcf",
                     "--genome", params.input_folder + params.genome,
                     "--output_wild_type", params.output_folder + "snp_in_protein-coding-regions_wt.fasta",
                     "--output_mutated", params.output_folder + "snp_in_protein-coding-regions_mut.fasta",
                     "--region_length", str(params.snp_genome_region_radius_protein_coding)])

    execute_command(docker_helper, 4, display,
                    ["python3", "/analytic-modules/mutated-sequence-generator/mutated_sequence_generator.py",
                     "--input_vcf", params.output_folder + "promoter-regions.vcf",
                     "--genome", params.input_folder + params.genome,
                     "--output_wild_type", params.output_folder + "snp_in_promoter-regions_wt.fasta",
                     "--output_mutated", params.output_folder + "snp_in_promoter-regions_mut.fasta",
                     "--region_length", str(params.snp_genome_region_radius_promoter)])

    execute_command(docker_helper, 5, display,
                    ["python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                     "--genomic", params.output_folder + "snp_in_protein-coding-regions_mut.fasta",
                     "--mirna", params.input_folder + params.mirna_fasta,
                     "--output", params.output_folder + "mirna_gene_connections_mut.tsv",
                     "--score", str(params.miranda_score_threshold),
                     "--energy", str(params.miranda_energy_threshold)])

    execute_command(docker_helper, 6, display,
                    ["python3", "/analytic-modules/mirna-interaction-predictor/mirna_interaction_predictor.py",
                     "--genomic", params.output_folder + "snp_in_protein-coding-regions_wt.fasta",
                     "--mirna", params.input_folder + params.mirna_fasta,
                     "--output", params.output_folder + "mirna_gene_connections_wt.tsv",
                     "--score", str(params.miranda_score_threshold),
                     "--energy", str(params.miranda_energy_threshold)])

    execute_command(docker_helper, 7, display,
                    ["python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                     "--fasta", params.output_folder + "snp_in_promoter-regions_mut.fasta",
                     "--matrix", params.input_folder + params.tf_binding_matrices,
                     "--tf_background_rsat", params.input_folder + params.tf_background_rsat,
                     "--tf_background_fimo", params.input_folder + params.tf_background_fimo,
                     "--output", params.output_folder + "tf_gene_connections_mut.tsv",
                     "--format", "transfac",
                     "--threshold", str(params.tf_score_threshold)])

    execute_command(docker_helper, 8, display,
                    ["python3", "/analytic-modules/transcription-factor-interaction-predictor/tf_interaction_prediction.py",
                     "--fasta", params.output_folder + "snp_in_promoter-regions_wt.fasta",
                     "--matrix", params.input_folder + params.tf_binding_matrices,
                     "--tf_background_rsat", params.input_folder + params.tf_background_rsat,
                     "--tf_background_fimo", params.input_folder + params.tf_background_fimo,
                     "--output", params.output_folder + "tf_gene_connections_wt.tsv",
                     "--format", "transfac",
                     "--threshold", str(params.tf_score_threshold)])

    execute_command(docker_helper, 9, display,
                    ["python3", "/analytic-modules/network-combiner/network_combiner.py",
                     "--input-files", params.output_folder + "mirna_gene_connections_mut.tsv,/output/tf_gene_connections_mut.tsv",
                     "--output-file", params.output_folder + "combined_connections_mut.tsv",
                     "--method", "union"])

    execute_command(docker_helper, 10, display,
                    ["python3", "/analytic-modules/network-combiner/network_combiner.py",
                     "--input-files", params.output_folder + "mirna_gene_connections_wt.tsv,/output/tf_gene_connections_wt.tsv",
                     "--output-file", params.output_folder + "combined_connections_wt.tsv",
                     "--method", "union"])

    execute_command(docker_helper, 11, display,
                    ["python3", "/analytic-modules/network-combiner/network_combiner.py",
                     "--input-files", params.output_folder + "combined_connections_mut.tsv,/output/combined_connections_wt.tsv",
                     "--output-file", params.output_folder + "differences_between_mut_wt_networks.tsv",
                     "--method", "difference"])

    execute_command(docker_helper, 12, display,
                    ["python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                     "--input-network-file", params.output_folder + "differences_between_mut_wt_networks.tsv",
                     "--lower-case",
                     "--no-isoform"])

    execute_command(docker_helper, 13, display,
                    ["python3", "/analytic-modules/network-id-mapper/network_id_mapper.py",
                     "--input", params.output_folder + "differences_between_mut_wt_networks_formatted.tsv",
                     "--target-id-type", "uniprotac",
                     # "--molecule-type-filter", "gene",  # do we need this ???
                     # "--remove",
                     "--mapping-data", ",".join(map(lambda x: params.input_folder + x, params.id_mapping_json_files)),
                     "--output", params.output_folder + "uniprot_differences.tsv"])

    execute_command(docker_helper, 14, display,
                    ["python3", "/analytic-modules/network-enrichment/network_enrichment.py",
                     "--input", params.output_folder + "uniprot_differences.tsv",
                     "--output", params.output_folder + "enriched_uniprot_differences.tsv",
                     "--reference-net", params.input_folder + params.reference_interactions_for_enrichment_tsv,
                     "--distance", "1"])

    execute_command(docker_helper, 15, display,
                    ["python3", "/analytic-modules/uniprot-id-formatter/uniprot_id_formatter.py",
                     "--input-network-file", params.output_folder + "enriched_uniprot_differences.tsv",
                     "--upper-case",
                     "--no-isoform"])

    if not params.separate:
        docker_helper.kill_running_container()

    display.print("workflow finished successfully! :)")
    display.wait_for_key()


def signal_handler(sig, frame):
    print('\n\nYou pressed Ctrl+C! try to kill and remove the running docker container if possible')
    if docker_helper:
        docker_helper.kill_running_container()
    sys.exit(1)

if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    navigomix_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
    params = ArgumentParser(navigomix_path, sys.argv[1:])
    if params.debug:
        main(None, params, navigomix_path)
    else:
        wrapper(main, params, navigomix_path)
