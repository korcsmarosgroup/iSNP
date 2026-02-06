from time import strftime
import multiprocessing
import numpy as np
import argparse
import logging
import sys
import os


def ExecuteProcess(command):
    os.system(command)


def split_list_np(lst, x):
    result = np.array_split(lst, x)
    return result


def parse_args(args):
    help_text = \
        """
        === Create Containers ===

        **Description:** 

        This tool is creating all the necessary docker containers for the paralell iSNP pipeline run.
        
        
        **Parameters:** 
        
        -i, --input-folder <path>                 : input folder [Mandatory]
        
        -o, --output-folder <path>                : output folder [Mandatory]

        -p, --patient-folder <path>               : folder, where all the patient specific VCF files are [Mandatory]
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-folder",
                        help="<path to the input folder> [mandatory]",
                        type=str,
                        dest="input_folder",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output-folder",
                        help="<path to the output folder> [mandatory]",
                        type=str,
                        dest="output_folder",
                        action="store",
                        required=True)

    parser.add_argument("-p", "--patient-folder",
                        help="<path to patient specific input VCF files> [mandatory]",
                        type=str,
                        dest="patient_folder",
                        action="store",
                        required=True)

    parser.add_argument("-n", "--number-of-runs",
                        help="<number of the paralell runs> [mandatory]",
                        type=int,
                        dest="number_of_runs",
                        action="store",
                        required=True)

    parser.add_argument("-id", "--identifier",
                        help="<identifier of the run> [mandatory]",
                        type=str,
                        dest="identifier",
                        action="store",
                        required=True)

    parser.add_argument("-dp", "--done-patients",
                        help="<path to the done patients file> [mandatory]",
                        type=str,
                        dest="done_patients",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_folder, results.output_folder, results.patient_folder, results.number_of_runs, results.identifier, results.done_patients


def main():

    input_folder, output_folder, patient_folder, number_of_runs, identifier, done_patients = parse_args(sys.argv[1:])

    if os.path.isfile(f"SubprocessRun_{identifier}.log"):
        os.remove(f"SubprocessRun_{identifier}.log")

    logging.basicConfig(filename = f"SubprocessRun_{identifier}.log", level = logging.INFO)
    logging.info(f"### [{strftime('%H:%M:%S')}] Starting the subprocesses!")
    multiprocessing_tuple = tuple()
    all_patient_files = []
    done_patients_list = []

    logging.info(f"### [{strftime('%H:%M:%S')}] Reading the done patients file")
    with open(done_patients, "r") as dp:
        for line in dp:
            done_patients_list.append(line.strip())

    logging.info(f"### [{strftime('%H:%M:%S')}] The number of done patients: {len(done_patients)}\nThey will be excluded from the subprocesses")

    logging.info(f"### [{strftime('%H:%M:%S')}] Creating the multiprocessing tuple")
    for actual_patient in os.listdir(patient_folder):

        if ".vcf" not in actual_patient:
            continue

        if actual_patient in done_patients:
            logging.info(f"### [{strftime('%H:%M:%S')}] The patient {actual_patient} is already done, skipping it")
            continue

        if actual_patient.split(".")[0] not in all_patient_files:
            all_patient_files.append(actual_patient)
    
    all_lists = split_list_np(all_patient_files, number_of_runs)

    counting_index = 1
    for list in all_lists:
        actual_list = ",".join(list)
        multiprocessing_tuple = multiprocessing_tuple + (f"python3 isnp_alternative.py -i {input_folder} -o {output_folder} -p {actual_list} -pf {patient_folder} -ci {counting_index} -ri {identifier}",)
        counting_index += 1

    logging.info(f"### [{strftime('%H:%M:%S')}] This is the multiprocessing tuple: {multiprocessing_tuple}")

    logging.info(f"### [{strftime('%H:%M:%S')}] Run the processes paralell; the number of the paralell processes: {number_of_runs}")
    process_pool = multiprocessing.Pool(processes = number_of_runs)
    process_pool.map(ExecuteProcess, multiprocessing_tuple)

    logging.info(f"### [{strftime('%H:%M:%S')}] The subprocesses finished successfully!")


if __name__ == "__main__":
    main()
