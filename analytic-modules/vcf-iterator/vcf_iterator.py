import argparse
import os
import shutil
import subprocess
import sys
import tarfile


def parse_args(args):
    help_text = \
        """
        === VCF-Iterator ===
        
        **Description:**

        The iterator modules are the entry points of iterations in the workflow.
        This iterator takes a set of VCF files as an input, and then iterate over them,
        providing single VCF variable for each iteration.
        
        The tool supports three input formats. The input can be specified as
        - a paths to VCF files
        - a path to a tar.gz archive file, contains a set of VCF files
        - a path to a standard plink files
        
        The output is an ordered set of VCF files in a given output folder. The file names
        will be `iteration-1`, `iteration-2`, ... The ordering of the output files is
        deterministic, following order how the user defined the vcf file paths, and the order
        of the vcf files in the archive, and the order of columns in the plink files.
        
        At least one of the input parameters must be specified. If multiple input parameters
        are used, then the tool will first iterate over the vcf file paths given by the user,
        then processes the archive file and finally takes the plink files.
        
        
        **Parameters:**
        
        --input-vcf-files <comma separated list of VCF file paths> [optional]
        
        --input-vcf-archive <path to the input tar.gz file, contains VCF files> [optional]

        --input-plink-files <comma separeted list to plink file paths> [optional]
        
        --output-folder [mandatory]
        
        
        **Exit codes**
        
        Exit code 1: One of the specified .vcf file doesn't exists!
        Exit code 2: The specified tar.gz file doesn't exists!
        Exit code 3: One of the specified plink file doesn't exists!
        Exit code 4: The specified output folder doesn't exists!
        Exit code 5: It's not a .vcf file!
        Exit code 6: It's not a tar.gz file!
        Exit code 7: The plink found an error during individual extraction!
        Exit code 8: The plink found an error during making vcf files!
        """

    parser = argparse.ArgumentParser(description=help_text)

    # Comma separated paths of input vcf files
    parser.add_argument("-i", "--input-vcf-files",
                        help="<paths to the input vcf files> [optional]",
                        type=str,
                        dest="input_files",
                        action="store",
                        default=False)

    # Input tar.gz file
    parser.add_argument("-t", "--input-vcf-archive",
                        help="<path to the tar.gz input file> [optional]",
                        dest="targz_file",
                        action="store",
                        default=False)

    # Comma separated paths of plink files (.bed, .fam, .bim)
    parser.add_argument("-p", "--input-plink-files",
                        help="<paths to the plink files> [optional]",
                        type=str,
                        dest="plink_files",
                        action="store",
                        default=False)

    # Output file path
    parser.add_argument("-o", "--output-folder",
                        help="<path to the output folder> [mandatory]",
                        dest="output_folder",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_files, results.targz_file, results.plink_files, results.output_folder


def check_input_paths(input_paths):
    input_vcf_files = input_paths.split(",")

    for file in input_vcf_files:
        if not os.path.isfile(file):
            sys.stderr.write(f"ERROR! One of the specified .vcf file doesn't exists: {file}")
            sys.exit(1)  # Exit code 1: One of the specified .vcf file doesn't exists!


def check_targz_file(targz_file):
    if not os.path.isfile(targz_file):
        sys.stderr.write(f"ERROR! The specified tar.gz file doesn't exists: {targz_file}")
        sys.exit(2)  # Exit code 2: The specified tar.gz file doesn't exists!


def check_plink_files(plink_file_paths):
    plink_files = plink_file_paths.split(",")

    for file in plink_files:
        if not os.path.isfile(file):
            sys.stderr.write(f"ERROR! One of the specified plink file doesn't exists: {file}")
            sys.exit(3)  # Exit code 3: One of the specified plink file doesn't exists!


def check_output_folder(output_folder):
    if not os.path.isdir(output_folder):
        sys.stderr.write(f"ERROR! The specified output folder doesn't exists: {output_folder}")
        sys.exit(4)  # Exit code 4: The specified output folder doesn't exists!


def input_paths_iteration(input_paths, output_folder):
    input_vcf_files = input_paths.split(",")
    index_input = 1

    for files in input_vcf_files:
        shutil.copy(files, output_folder)

        for file in os.listdir(output_folder):
            word = "iteration-"
            if word in file:
                continue

            else:
                old_file = os.path.join(output_folder, file)
                new_file = os.path.join(output_folder, f"iteration-{index_input}")
                os.rename(old_file, new_file)
                index_input = index_input + 1


def extract_targz_file(targz_file, output_folder):
    tar = tarfile.open(targz_file, "r:gz")

    for members in tar:
        tar.extract(members, output_folder)
    tar.close()


def tar_file_iteration(output_folder, index_targz):
    for dirpath, dirnames, filenames in os.walk(output_folder):
        for filename in filenames:
            word = "iteration-"
            if word in filename:
                continue
            else:
                old_file = os.path.join(dirpath, filename)
                new_file = os.path.join(output_folder, f"iteration-{index_targz}")
                os.rename(old_file, new_file)
                index_targz = index_targz + 1


def plink_extraction(plink_file_paths, output_folder, index_plink):
    in_dataset = plink_file_paths.split(",")[0].split("/")[-1][:-4]
    plink_files = plink_file_paths.split(",")

    shutil.copy(plink_files[0], output_folder)
    shutil.copy(plink_files[1], output_folder)
    shutil.copy(plink_files[2], output_folder)

    fam_format = f"{in_dataset}.fam"
    file = os.path.join(output_folder, fam_format)

    with open(file, 'r') as fam_file:

        for line in fam_file:
            line = line.strip()
            cells = line.split()
            id = cells[0]

            tmp_file = os.path.join(output_folder, "tmp.keep")
            with open(tmp_file, "w") as tmp_keep:
                tmp_keep.write("%s %s\n" % (id, id))

            plinkCMD = ["plink2", "--bfile", in_dataset, "--keep", "tmp.keep", "--make-bed", "--out", id]

            proc = subprocess.Popen(plinkCMD, cwd=output_folder)
            return_code = proc.wait()

            if proc.returncode == 0:
                print("====== analytical tool executed successfully! ======")

            else:
                sys.stderr.write(f"ERROR! Plink exited with error code: {return_code}")
                sys.exit(7)  # Exit code 7: The plink found an error during individual extraction!

            iteration = os.path.join(f"iteration-{index_plink}")

            plinkProc = ["plink2", "--bfile", id, "--recode", "vcf", "--out", iteration]

            proc2 = subprocess.Popen(plinkProc, cwd=output_folder)
            return_code = proc2.wait()

            if proc.returncode == 0:
                print("====== analytical tool executed successfully! ======")

            else:
                sys.stderr.write(f"ERROR! Plink exited with error code: {return_code}")
                sys.exit(8)  # Exit code 8: The plink found an error during making vcf files!

            index_plink = index_plink + 1

    for files in os.listdir(output_folder):
        filenames = os.path.join(output_folder, files)

        if filenames.endswith(".bed"):
            os.remove(filenames)
        elif filenames.endswith(".fam"):
            os.remove(filenames)
        elif filenames.endswith(".bim"):
            os.remove(filenames)
        elif filenames.endswith(".log"):
            os.remove(filenames)
        elif filenames.endswith(".keep"):
            os.remove(filenames)
        else:
            continue


def main():
    input_files, targz_file, plink_files, output_folder = parse_args(sys.argv[1:])

    check_output_folder(output_folder)

    if input_files and targz_file and plink_files:

        check_input_paths(input_files)
        check_targz_file(targz_file)
        check_plink_files(plink_files)
        input_paths_iteration(input_files, output_folder)
        index_targz = len(input_files.split(",")) + 1
        extract_targz_file(targz_file, output_folder)
        tar_file_iteration(output_folder, index_targz)
        index_plink = len([name for name in os.listdir(output_folder) if name.endswith(".vcf")]) + 1
        plink_extraction(plink_files, output_folder, index_plink)

    elif input_files and targz_file:

        check_input_paths(input_files)
        check_targz_file(targz_file)
        input_paths_iteration(input_files, output_folder)
        index_targz = len(input_files.split(",")) + 1
        extract_targz_file(targz_file, output_folder)
        tar_file_iteration(output_folder, index_targz)

    elif input_files and plink_files:

        check_input_paths(input_files)
        check_plink_files(plink_files)
        input_paths_iteration(input_files, output_folder)
        index_plink = len(input_files.split(",")) + 1
        plink_extraction(plink_files, output_folder, index_plink)

    elif targz_file and plink_files:

        check_targz_file(targz_file)
        check_plink_files(plink_files)
        extract_targz_file(targz_file, output_folder)
        index_targz = 1
        tar_file_iteration(output_folder, index_targz)
        index_plink = len([name for name in os.listdir(output_folder) if name.endswith(".vcf")]) + 1
        plink_extraction(plink_files, output_folder, index_plink)

    elif input_files:

        check_input_paths(input_files)
        input_paths_iteration(input_files, output_folder)

    elif targz_file:

        check_targz_file(targz_file)
        extract_targz_file(targz_file, output_folder)
        index_targz = 1
        tar_file_iteration(output_folder, index_targz)

    elif plink_files:

        check_plink_files(plink_files)
        index_plink = 1
        plink_extraction(plink_files, output_folder, index_plink)

    print(f"====== VCF-Iterator finished successfully! ======")


if __name__ == "__main__":
    main()
