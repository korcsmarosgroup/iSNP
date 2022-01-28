# VCF Filtering

** Description **

VCF filtering takes SNPs (e.g. from the patient's genome) and intersect these with any other 
regions defined in input annotation (.bed) files. You can filter SNPs inside miRNA coding 
region or in gene promoter region, etc. depending on what input annotation file you use. 

Alternatively you can also specify a list of SNP identifiers in a file, and the tool will also
use that for filtering.

You can specify both annotation file and both SNP identifier list file in the same time, the tool
will filter on both of these in this case (keeping the SNPs which are matching with both approaches).

The inputs are SNPs in VCF format. The output will be represented as VCF file. 

We expect for the input annotation file to contain a fully qualified NavigOmiX entity id. This 
id will be preserved in the output file. 


Example for the metadata (when using annotation files):
- in the input VCF: you have the snp id for each snp
- input annotation we expect the name column to be set with fully qualified name 
(e.g. for gene promoter regions, the name of the gene): gene;uniprotkb;P31946
- in the output bed file the tool will set the following name: gene;uniprotkb;P31946|snp;dbsnp;rs987423


Example for the SNP identifier list:
snp;dbsnp;rs987423
snp;dbsnp;rs18363
snp;dbsnp;rs338276

(the tool also accepts a simple text file, without fully qualified IDs, just a SNP ID in each line)

** Parameters: **

-i, --input <file path> : obtained as a result of a sequencing experiment - (VCF file format)

-a, --annotation <file path> : this is a .bed file, that contains the appropriate annotations: miRNA coding regions and gene promoter regions - (BED file format)

-s, --snp <file path> : this is a NavigOmix entity set file, that contains the set of SNP IDs for filtering

-o, --output <file path> : containing just filtered mutations - (VCF file format)

