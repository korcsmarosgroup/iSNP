import os
import re
import subprocess
import sys


def _create_bedtools_command(path_to_vcf,
                             path_to_bed,
                             header,
                             bed_info):
    """
    https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

    Parameters
    ----------
    - path_to_vcf: character path to the .vcf file to filter.
    - path_to_bed: character path to the .bed file to define the filtering.
    - out_vcf: character path where to save the filtered .vcf.
    - header: boolean. Should the header of the .vcf file be included in out_vcf? default = True
    - bed_info: boolean. Should the information in the bed file be included (e.g. names)? Note that
    the output file may not be .vcf compliant. default = False.

    Returns
    -------
    bed_tools_command: call for bedtools
    """
    bed_tools_command = ['bedtools', 'intersect', '-a', path_to_vcf, '-b', path_to_bed, '-wa']
    if header:
        bed_tools_command.append("-header")
    if bed_info:
        bed_tools_command.append("-wb")
    else:
        bed_tools_command.append("-u")
    return bed_tools_command


def filter_vcf(path_to_vcf, path_to_bed, out_vcf, header=True, bed_info=False):
    """
    Filters the contents of a .vcf file based on the contents of a .bed file using the intersect
    function from bedtools. The result is written as a .vcf file in out_vcf_file.

    Note that, by default, the information in the bed file (chrom, start, end, name) is included as
    additional fields in the .vcf. Therefore, the output file is not compliant with the .vcf format.
    I am doing this just because it is convenient to propagate the NOXid through the workflow (using 
    the name of the .bdf). By setting bed_info = FALSE, these fields are not output, so the file is
    .vcf compliant (although the information from the bed is lost).
    
    Parameters
    ----------
    - path_to_vcf: character path to the .vcf file to filter.
    - path_to_bed: character path to the .bed file to define the filtering.
    - out_vcf: character path where to save the filtered .vcf.
    - header: boolean. Should the header of the .vcf file be included in out_vcf? default = True
    - bed_info: boolean. Should the information in the bed file be included (e.g. names)? Note that
    the output file may not be .vcf compliant. default = False.
    """

    if not os.path.exists(path_to_vcf):
        sys.stderr.write("vcf file not found: " + path_to_vcf)
        sys.exit(201)

    if not os.path.exists(path_to_bed):
        sys.stderr.write("Could not find bed file: " + path_to_bed)
        sys.exit(203)

    bed_tools_command = _create_bedtools_command(path_to_vcf,
                                                 path_to_bed,
                                                 header,
                                                 bed_info)
    p = subprocess.Popen(bed_tools_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    my_stdout, my_stderr = p.communicate()
    print(my_stderr.decode("utf-8"))
    with open(out_vcf, mode="w") as out_vcf_file:
        out_vcf_file.write(my_stdout.decode("utf-8"))


def filter_vcf_by_id(path_to_vcf, path_to_snp_list, out_vcf, header=True):
    """
    Filters the contents of a .vcf file based on the contents of a .bed file using the intersect
    function from bedtools. The result is written as a .vcf file in out_vcf_file.
    Parameters
    ----------
    - path_to_vcf: character path to the .vcf file to filter.
    - path_to_bed: character path to the .bed file to define the filtering.
    - out_vcf: character path where to save the filtered .vcf.
    - header: boolean. Should the header of the .vcf file be included in out_vcf? default = True
    """
    snp_list = []
    with open(path_to_snp_list) as fin:
        for each_line in fin:
            snp_list.append(each_line.strip())
    nox_id_regex = re.compile('[^;]+;[^;]+;([^;]+)')
    snp_list = list(map(lambda x: nox_id_regex.match(x).group(1) if nox_id_regex.match(x) else x, snp_list))
    with open(path_to_vcf) as in_vcf:
        with open(out_vcf, mode = "w") as out_vcf_file:
            for each_line in in_vcf:
                if each_line.startswith("#"):
                    if header:
                        out_vcf_file.write(each_line)
                else:
                    my_data = each_line.split("\t")
                    if len(my_data) < 3:
                        print("SNP had length " + str(len(my_data)) + ". A zero SNP will be generated.")
                    else:
                        snp_id = my_data[2].strip()
                        if snp_id in snp_list:
                            out_vcf_file.write(each_line)
