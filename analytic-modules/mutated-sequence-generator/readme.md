# Mutated Sequence Generator

**Description:** 

This tool will take a set of SNP data (in format of VCF), then by the location 
information in the VCF files it will extract an upstream and downstream region
from the genome in a given length. 

The tool will generate two sequence files in FASTA format. One output file will 
contain all the wild type sequences, while the other file will contain the mutated 
sequences. In the FASTA file, each sequence will be tagged with the SNP id, with the 
sequence type (mutated: true or false) and with all the tags in the input VCF file.
(the tag 'gene_id' will be renamed to 'entity')

If the genomic location can not be found in the reference genome, then the sequences 
will not be added to the output FASTA files (meaning there will be no "empty" sequence
in the output).

If the given SNP is 'inactive' in the input VCF file (having 0/0 in the 10th column of
th VCF file), then the same wild type sequence will be generated to both output files.

If the alternative allel in the input vcf file is a '-' or '*' symbol, then this is a deletion,
which will be shown in the mutated output fasta file as an empty character! So for this SNP,
there will be a one base shorter sequence in the mutant fasta file. For example: the region-length is 21, then there
will be 43 base long sequence in the wild type fasta file and 42 base long sequence in the mutant fasta file.

If the alternative allel in the input vcf file is an 'n' r 'N' character, then there will be no sequence for that SNP
in the output files.

**Parameters:**

--input_vcf <path to the input VCF file> [mandatory]

--genome <path to the genome (in single FASTA format)> [mandatory]   

--output_wild_type <path to the new output FASTA file with wild type sequences> [mandatory]

--output_mutated <path to the new output FASTA file with mutated sequences> [mandatory]

--region_length <0..1000: using e.g. 100 will fetch region: [X-100, x, x+100] with total length of 201> [mandatory]

## Error codes

In case of an exception, the following error codes are raised:

- No error: 0

### Wrong types 1XX

- Wrong type for input_vcf (not str): 101
- Wrong type for genome (not str): 102
- Wrong type for output_wild_type (not str): 103
- Wrong type for output_mutate (not str): 104
- Wrong type for region_length (not int): 105
- Region length smaller than 1: 106

### Files not found 2XX

- Not found inpup_vcf: 201
- Not found genome: 202
- Not found bed file (should not happen when mutate is called): 203

