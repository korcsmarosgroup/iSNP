import os


def generate(path_to_fasta, path_to_vcf, out_path, read_length):
    """
    Generate a mutated sequence based on a reference one

    Details
    -------
    The functions does not perform any testing and  relies quite heavily that the input files
    are in the right format. That is:

    - There is one entry in the FASTA file per SNP.
    - The sequences in the FASTA file are in the same order as in the VCF.
    - The read_length is the same one used in vcf2bed.
    - There is no mismatch between the FASTA and the VCF.

    This should not be a problem if the function is called within the workflow, because the input
    files for this function are output of our functions, so we have quite a lot of control. But this
    function should be used with LOTS OF CARE WHEN IMPLEMENTED IN OTHER WORKFLOWS.

    Parameters
    ----------
    path_to_fasta : str
        path to a FASTA file defining the reference sequence. There must be one entry
        per SNP. It is recommended to use extract_genome().

    path_to_vcf : str
        path to a VCF file defining the SNPs.

    out_path : str
        path where to write the output sequence.

    read_length : int
        read length used to extract the genome with extract_genome().

    Returns
    -------
    A FASTA file identical to path_to_fasta, but with the mutations defined in path_to_vcf.

    """

    if not os.path.exists(path_to_fasta):
        raise FileNotFoundError("FASTA file not found: " + path_to_fasta)

    if not os.path.exists(path_to_vcf):
        raise FileNotFoundError("VCF file not found: " + path_to_vcf)

    with open(path_to_fasta) as my_fasta:
        with open(path_to_vcf) as my_vcf:
            with open(out_path, mode="w") as fout:
                my_iterator = fasta_iterator(my_fasta)
                for vcf_line in my_vcf:
                    vcf_line = vcf_line.strip().split('\t')
                    if vcf_line[0].startswith("#") or "N" in vcf_line[4] or "n" in vcf_line[4]:
                        continue
                    else:
                        head, seq = my_iterator.__next__()
                        ref = vcf_line[3]
                        bases = ['A', 'C', 'G', 'T']
                        if "1" in vcf_line[9]:
                            alt_bases = vcf_line[4].split(",")
                            for alt_base in alt_bases:
                                if alt_base == '-' or alt_base == '*':
                                    mutated_seq = seq[0:read_length] + seq[read_length + len(ref):]
                                    mutated_seq = mutated_seq.strip()
                                elif alt_base in bases:
                                    alt = alt_base.replace(".", "")
                                    mutated_seq = seq[0:read_length] + alt + seq[read_length + len(ref):]
                                    mutated_seq = mutated_seq.strip()
                                else:
                                    continue
                                out_head = head.replace("False", "True")
                                fout.write(out_head)
                                fout.write(mutated_seq)
                                fout.write("\n")
                        else:
                            fout.write(head)
                            fout.write(seq)
                            fout.write("\n")


def fasta_iterator(file_object):
    """
    An iterator for a FASTA file separating the header and nucleotide sequence of each entry.

    Parameters
    ----------
    - file_object : a file object (e.g. _io.TextIOWrapper)

    Returns
    -------
    An iterator. On each call, it yields a tuple with two elements:
    """
    head = ""
    seq = ""
    for line in file_object:
        if line.startswith(">"):
            out_head = head
            out_seq = seq
            head = line
            seq = ""
            if out_head:
                yield out_head, out_seq.replace("\n", "")
        else:
            seq = seq + line
    yield head, seq.replace("\n", "")
