""" Fimo TF interaction predictor """
import csv
import os
import sys
import subprocess
import tempfile
sys.path.append("/rds/general/user/jno25/home/iSNP/analytic-modules")
from common_libs.mitab_handler import mitab_handler


def _create_interaction(pred, mitab, snp, uniprot_id, motif_id):
    """ Add a new interaction line to the mitab dict builder """
    fimo_confidence_score = f"fimo-p-value:{pred[6]};fimo-q-value:{pred[6]}"
    interaction = mitab.new_interaction()
    interaction[mitab.uidA] = f'entity:tf;{uniprot_id.replace(" ", "")}'
    interaction[mitab.uidB] = f'{pred[2]}'
    interaction[mitab.taxA] = "taxid:9606('homo sapiens')"
    interaction[mitab.taxB] = "taxid:9606('homo sapiens')"
    interaction[mitab.confidence] = fimo_confidence_score
    interaction[mitab.annotA] = f'motif_id:{motif_id}'
    interaction[mitab.annotB] = f'sequence_name:{pred[2]}'
    interaction[mitab.annotInter] = f'origin:snp;dbsnp;{snp.split(":")[1]}'
    return interaction


def create_network_file(fimo_output_predictions, uniprot_motif_mapping_dict, dbsnp_gene_dict, output):
    """ A method to create the network file with the newly predicted interactions using the mitab handler """
    inner_structure = {}
    mitab = mitab_handler.MiTabHandler()

    with open(fimo_output_predictions, "r") as fimo_preds:
        reader = csv.reader(fimo_preds, delimiter='\t')
        next(reader, None)
        uniprot_id_complex = None
        idx = 0
        for fimo_prediction in reader:
            if len(fimo_prediction) == 0:
                break

            snp = dbsnp_gene_dict[fimo_prediction[2]]
            uniprot_id = uniprot_motif_mapping_dict[fimo_prediction[0]]
            motif_id = fimo_prediction[0]
            
            if len(motif_id.split("::")) > 2:
                temp_uniprot_id = uniprot_id.split(";")
                temp_motif_id = fimo_prediction[0].split("::")
                uniprot_id = f"{temp_uniprot_id[0]};{temp_uniprot_id[1]}"
                motif_id = temp_motif_id[0]
                uniprot_id_complex = f"{temp_uniprot_id[0]};{temp_uniprot_id[2]}"
                motif_id_complex = temp_motif_id[1]

            interaction = _create_interaction(fimo_prediction, mitab, snp, uniprot_id, motif_id)
            inner_structure[idx] = interaction
            if uniprot_id_complex:
                idx += 1
                interaction = _create_interaction(fimo_prediction, mitab, snp, uniprot_id_complex, motif_id_complex)
                inner_structure[idx] = interaction
                uniprot_id_complex = None
            idx += 1

    mitab.build_network_frame(inner_structure=inner_structure)
    mitab.serialise_mitab(output, add_header=False)


def extract_and_map_motif_ids(motif_file):
    """ Extract the uniprot accessions from the transfec file """
    accession_list = []
    uniprot_list = []
    with open(motif_file, "r") as transfec_matrix:
        for line in transfec_matrix:
            if line.startswith("ID"):
                accession_list.append(line.split(" ", 1)[1].strip("\n"))
            if "uniprot_ids" in line:
                line = line.replace("uniprot_ids", "uniprotac")
                line = line.replace(":", ";")
                uniprot_list.append(line.split(" ", 1)[1].lower().strip("\n"))
    return dict(zip(accession_list, uniprot_list))


def extract_and_map_dbsnp(fasta_file):
    """ Extract the dbsnp id from the patient fasta file """
    uniprot_gene_id_list = []
    snp_id_list = []
    with open(fasta_file, "r") as patient_fasta:
        for line in patient_fasta:
            if line.startswith(">"):
                line = line.strip(">")
                uniprot_gene_id_list.append(line.split("|")[0].strip(" "))
                snp_id_list.append(line.split("|")[1].strip(" "))
    return dict(zip(uniprot_gene_id_list, snp_id_list))


def run_fimo(motif_file, sequence_file, output_folder, background_file, mitab_output_folder):
    """ Fimo prediction module using the FIMO tool to predict the TF interactions """
    with tempfile.TemporaryDirectory() as tmpdirname:
        print(f"Creating temp directory for fimo results: {tmpdirname}")
        fimo_output = os.path.join(output_folder, 'fimo.log')
        meme_motif_file_path = os.path.join(tmpdirname, 'temp_pfms_meme.txt')
        with open(fimo_output, "w") as fimo, open(meme_motif_file_path, "w") as meme_motif_file:
            convertion_args = ["transfac2meme", motif_file]
            fimo_prediction_args = ["fimo",
                                    "--bgfile", background_file,
                                    "--oc", output_folder,
                                    "--thresh", "1e-3",
                                    meme_motif_file_path,
                                    sequence_file]
            convert_motif_file = subprocess.Popen(convertion_args, stdout=meme_motif_file, universal_newlines=True)
            _ = convert_motif_file.wait()
            fimo_preds = subprocess.Popen(fimo_prediction_args, stdout=fimo, universal_newlines=True)
            return_code = fimo_preds.wait()
            if return_code != 0:
                _display_return_code(return_code)

        uniprot_motif_mapping_dict = extract_and_map_motif_ids(motif_file)
        dbsnp_gene_mapping_dict = extract_and_map_dbsnp(sequence_file)
        fimo_output_preds = os.path.join(output_folder, "fimo.tsv")
        create_network_file(fimo_output_preds, uniprot_motif_mapping_dict, dbsnp_gene_mapping_dict, mitab_output_folder)


def _display_return_code(return_code):
    """ Displays the error code from the Fimos modules for debugging purposes """
    print("There was an error running FIMO to predict the TFBS. The following error number was raised...\n")
    print(f"Return Code: {return_code}")
    print("Exit status of the child process. Typically, an exit status of 0 indicates that it ran successfully.")
    print("A negative value -N indicates that the child was terminated by signal N (POSIX only).")
