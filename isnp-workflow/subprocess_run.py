import multiprocessing
import numpy as np
import argparse
import shutil
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

    results = parser.parse_args(args)
    return results.input_folder, results.output_folder, results.patient_folder, results.number_of_runs


def main():

    print(f'====== Starting! ======')
    input_folder, output_folder, patient_folder, number_of_runs = parse_args(sys.argv[1:])
    multiprocessing_tuple = tuple()
    all_patient_files = []

    for actual_patient in os.listdir(patient_folder):

        if ".vcf" not in actual_patient:
            continue

        actual_patient_name = actual_patient.split(".")[0]
        actual_patient_folder = os.path.join(output_folder, actual_patient_name)

        if actual_patient not in all_patient_files:
            all_patient_files.append(actual_patient)

        if os.path.isdir(actual_patient_folder):
            shutil.rmtree(actual_patient_folder)

        if not os.path.isdir(actual_patient_folder):
            os.mkdir(actual_patient_folder)
    
    all_lists = split_list_np(all_patient_files, number_of_runs)

    for list in all_lists:
        actual_list = ",".join(list)
        multiprocessing_tuple = multiprocessing_tuple + (f"python3 isnp_alternative.py -i {input_folder} -o {output_folder} -p {actual_list} -pf {patient_folder}",)

    print(multiprocessing_tuple)

    process_pool = multiprocessing.Pool(processes = number_of_runs)
    process_pool.map(ExecuteProcess, multiprocessing_tuple)

    print(f'====== Create containers finished successfully! ======')


if __name__ == "__main__":
    main()