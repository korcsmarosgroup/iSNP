import os
import subprocess
from common_libs.mitab_handler import mitab_handler


def process_rsat_results(in_path, pval_threshold=None):
    """
    Converts the output of matrix-scan to a dictionary

    Parameters
    ----------
    in_path: Path where the output of RSAT is located.
    pval_threshold: Only connections with pval < pval_threshold are output. None by default (no filter).

    Output
    ------
    rsat_results: dictionary of resuls from rsat
    """
    rsat_results = {}
    with open(in_path) as fin:
        for each_line in fin:
            if each_line.startswith(";"):
                continue
            if not each_line:
                continue
            if each_line.startswith("#"):
                each_line = each_line.replace("#", "")
                each_line = each_line.strip()
                header_data = each_line.split("\t")
                try:
                    seq_index = header_data.index("seq_id")
                except ValueError:
                    print("Could not find index for seq_id")
                try:
                    prot_index = header_data.index("ft_name")
                except ValueError:
                    print("Could not find index for ft_name")
                try:
                    pval_index = header_data.index("Pval")
                except ValueError:
                    print("Could not find index for Pval")
            else:
                each_line = each_line.rstrip()
                line_data = each_line.split("\t")
                seq = line_data[seq_index]
                prot = line_data[prot_index]
                pval = float(line_data[pval_index])
                if pval_threshold:
                    if pval > pval_threshold:
                        continue
                if (seq, prot) in rsat_results:
                    if pval < rsat_results[(seq, prot)]:
                        rsat_results[(seq, prot)] = pval
                else:
                    rsat_results[(seq, prot)] = pval
    return rsat_results


def write_rsat_results(in_path, out_path, pval_threshold=None):
    """
    Write output dictionary to a MITAB file.

    Parameters
    ----------
    in_path: Path where the output of RSAT is located.
    out_path: Path where the output MITAB file is to be written.
    pval_threshold: Only connections with pval < pval_threshold are output. None by default (no filter).

    """
    rsat_results = process_rsat_results(in_path, pval_threshold)
    idx = 0
    mitab = mitab_handler.MiTabHandler()
    inner_structure = {}
    for tseq, tfprot in rsat_results:
        pval = rsat_results[(tseq, tfprot)]
        for tf in tfprot.split("::"):
            tf = tf.split("(")[0]
            interaction = mitab.new_interaction()
            interaction[mitab.uidA] = "name:%s" % tf
            interaction[mitab.uidB] = "uniprotac:%s" % tseq.split(";")[3].split("|")[0]
            interaction[mitab.taxA] = "taxid:9606('homo sapiens')"
            interaction[mitab.taxB] = "taxid:9906('homo sapiens')"
            interaction[mitab.confidence] = "rsat_pvalue:%.16f" % pval
            interaction[mitab.annotA] = "start:protein;name;%s" % tf
            interaction[mitab.annotB] = f'end:{tseq.split(":")[1].split("|")[0]}'
            interaction[mitab.annotInter] = f'origin:snp;dbsnp;{tseq.split("|")[1].split(":")[1]}'
            inner_structure[idx] = interaction
            idx += 1
    mitab.build_network_frame(inner_structure=inner_structure)
    mitab.serialise_mitab(out_path, add_header=False)


def scan_matrix(path_to_fasta, out_path, path_to_matrix=None, format_matrix=None, path_to_background=None):
    """
    Finds TF binding sites in a (set of) secuence using the matrix-scan function of RSAT.

    Parameters
    ----------
    path_to_fasta: Path to the FASTA file with the sequences.
    out_path: Path to write the output.
    path_to_matrix: Path to the file with the matrix data. None by default (the default matrix).
    format_matrix: Format of the matrix.
    path_to_background: Path to the background file.

    Output
    ------
    A text file with a table with 11 columns (eq_id	ft_type, ft_name, strand, start, end, sequence,
        weight, Pval, ln_Pval, sig, normw).

    """
    default_matrix_path = "test.transfac"
    if path_to_matrix is None: 
        path_to_matrix = default_matrix_path
    if format_matrix is None:
        format_matrix = "transfac"
    if not os.path.exists(path_to_fasta):
        raise FileNotFoundError("Could not find fasta file: " + path_to_fasta)
    my_call = ['matrix-scan',
               '-matrix_format', format_matrix,
               '-m', path_to_matrix,
               '-i', path_to_fasta,
               '-bgfile', path_to_background,
               '-quick',
               '-return', 'pval',
               '-return', 'normw', '-2str',
               '-v', '1',
               '-seq_format', 'fasta',
               '-o', out_path]
    p = subprocess.Popen(my_call, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    my_stdout, my_stderr = p.communicate()
    print(my_stderr.decode("utf-8"))


if __name__ == "__main__":
    scan_matrix("test.fasta", "rsat_matrixscan.txt")
