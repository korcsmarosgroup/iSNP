# VCF iterator

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