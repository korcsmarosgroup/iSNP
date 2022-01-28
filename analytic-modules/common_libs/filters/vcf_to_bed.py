import os
import sys


class ProcessVcf:
    """
    A class for getting info from vcf files and output it to other formats.
    """

    def __init__(self, path, verbose=False):
        """
        Object for processing VCF files.
        The __init__ only initializes some attributes, does not read the vcf file. Methods have
        to be called for that.

        Parameters
        ----------
        path : str
            path to the vcf file to read
        verbose : boolean
            whether to print extra info to stdout.
        """
        if not os.path.exists(path):
            sys.stderr.write("vcf file not found: " + path)
            sys.exit(201)
        self.path = path
        self.verbose = verbose

    def vcf2bed(self, out_path, region_length, chr_decorator="chr"):
        with open(out_path, mode="w") as fout:
            with open(self.path) as fin:
                for this_line in fin:
                    if not this_line.startswith("#"):
                        my_snp = self._read_snp_line(this_line)
                        out_line = my_snp.get_bed_line(region_length, chr_decorator)
                        fout.write(out_line)

    @staticmethod
    def _read_snp_line(my_line, generate_nox_id=True):
        """
        Reads the information in a SNP line of a .vcf file and returns an instance of SNP.
        """
        this_data = my_line.split("\t")  # Handle possible errors here
        if len(this_data) < 5:
            print("SNP had length " + str(len(this_data)) + ". A zero SNP will be generated.")
            this_snp = SNP("0", "0", "ShortSNPError", "N", "N")
        else:
            chromosome = this_data[0].strip()
            pos = this_data[1].strip()
            this_id = this_data[2].strip()
            ref = this_data[3].strip()
            alt = this_data[4].strip()
            meta_from_bed = this_data[-1].strip().replace("gene_id:", "entity:")
            if generate_nox_id:
                this_id = f"{meta_from_bed} | origin:{this_id} | mutated:False"
            this_snp = SNP(chromosome, pos, this_id, ref, alt)
        return this_snp


class SNP:
    """
    Class for handling SNPs.
    """

    def __init__(self, chrom, pos, this_id, ref, alt):
        """
        Parameters
        ----------
        chrom : something that can be converted to str describing the chromosome position.
        pos : something that can be converted to integer describing the position of the SNP.
              if not, the position is set to int(0) and '_ValueError' appended to self.ID.
        this_id : something that can be converted to str describing the SNP ID.
        ref : something that can be converted to str describing the reference nucleotide sequence.
        alt : something that can be converted to str describing the alternative nucleotide sequence.

        Notes
        -----
        alt fields with ',' are not supported. If input, a warning message is printed and only the first
        SNP (the 1st before the coma) is saved.

        """
        self.chrom = str(chrom)
        self.this_id = str(this_id)
        try:
            self.pos = int(pos)
        except ValueError:
            print(pos + " was input as position of SNP with ID " + this_id +
                  ". Position '0' will be set and _ValueError appended to the ID.")
            self.pos = 0
            self.this_id = self.this_id + "_ValueError"
        self.ref = str(ref)
        alt = alt.replace(".", "")
        if "," in alt:
            print("WARNING: SNP with ',' are not yet supported. Only the 1st SNP is saved.")
            alt, _, _ = alt.partition(",")
        self.alt = alt.strip()

    def get_bed_line(self, region_length, chr_decorator=""):
        """
        Entry for a BED file describing the SNP position, so that it can be later extracted from
        a fasta file using bedtools.

        Parameters
        ----------
        region_length : int
            The length of the region.
        chr_decorator : chr
            Character to be added before the chr information saved in the SNP object. Empty string by default.
        """
        my_id = self.this_id
        if not type(region_length) is int:
            print("WARNING: The type of region_length is not integer for SNP " + self.this_id)
            try:
                region_length = int(region_length)
            except ValueError:
                print("WARNING: The region_length could not be converted to integer. It will be set to 1 \t"
                      "and the _BadLengthError appended to the SNP id in the bed file")
                region_length = 1
                my_id = my_id + "_BadLengthError"
        start = max(self.pos - region_length - 1, 0)
        end = self.pos + region_length
        try:
            my_chr = str(chr_decorator) + self.chrom
        except TypeError:
            print("WARNING: chr_decorator could not be combined with self.chrom. It will be omitted.")
            my_chr = self.chrom
        out_line = "\t".join([my_chr, str(start), str(end), my_id, "\n"])
        return out_line
