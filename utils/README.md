# iSNP Pipeline Pre/post-processing Scripts
A repo to hold the various scripts created for the processing and wrangling data for input to the iSNP pipeline and output of the iSNP pipeline.

## Preprocessing Data Scripts
Input data wrangling and preprocessing scripts for the iSNP pipeline. 

### Background File Generation
The background file generation for the two transcription factor binding tools. The background model normalizes for biased distribution of individual letters in the sequences. The background model generation scripts have been wrapped in `generate_background_input.py` which also includes additional processing and data extract required for generating these files.

The arguments for this wrapping script are as follows:
```
optional arguments:
  -h, --help            show this help message and exit
  --regulatory-regions-file REGULATORY_REGIONS_FILE
                        <bed file of the regulatory regions> [mandatory]
  --genome-folder GENOME_FOLDER
                        <directory of fasta files holding each chromosome of the genome> [mandatory]
  --background-fasta BACKGROUND_INPUT_FASTA
                        <name of the background fasta file to be created.> [mandatory]
  --output OUTPUT_FOLDER
                        <path to the output folder> [mandatory]
  --seq-length SEQ_LENGTH
                        <Filter selection parameter for sequence lengths. Default: 60.> [optional]
  --strand-type STRAND_TYPE
                        <Create background for single or double strands. Default: double.> [optional]
  --order ORDER         <Order for the background models. Default: 2.> [optional]
  --model-name MODEL_NAME
                        <Name of the prefix of the fimo and rsat model file. Default: None.> [optional]
```

The following scripts generated a model for both tools (Fimo and RSAT) using the following scripts. These scripts can be run independently if you wish not to use the wrapped script.
* Fimo model: `generate_rsat_bg_file.sh`
* RSAT model: `generate_fimo_bg_file.sh`

### Generate Annotations
Generates promoter/enhancer and protein coding regions annotation files for the miRNA and Transcription factors predictions. These annotation files generated from a gtf files which have the Genocde annotation format.

The arguments for this script are as follows:
```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_GTF_FILE, --input-gtf-file INPUT_GTF_FILE
                        <paths to the input vcf files> [mandatory]
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        <path to the output folder> [mandatory]
  -e ENHANCERS_FILE, --enhancers-file ENHANCERS_FILE
                        <paths to the enhancer annotation files> [optional]
  -v, --verbose         <verbosity level for debugging> [optional]
```

### VCF Wrangling scripts
There are several scripts to pre-processing the input VCF files for the iSNP pipeline. Depending on the current format of the VCF files you have.

#### VCF Splitter
This script takes a main VCF files and splits into sample specific VCF files. This is required for the parallelization of the iSNP pipeline. This way by creating a single VCF file for each sample/individual which will then be used in each run of the iSNP pipeline.

The arguments for this script are as follows:
```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_VCF, --input_vcf INPUT_VCF
                        <the input vcf file to be split> [mandatory]
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        <output directory for each individual vcf file> [mandatory]
  -p PREFIX, --prefix PREFIX
                        <the prefix to add for chromosome field [optional]>
```

#### VCF Header Wrangler
Reformat and fixes the header of the VCF files. The individual VCFs you may need to add back in the original VCF file header information. This requires a reference header file which only contains the wanted VCF file header.

The arguments for this script are as follows:
```
optional arguments:
  -i INPUT_FOLDER, <Input folder where the VCF files are hled> [mandatory]
  -o OUTPUT_DIR, <Ouput folder to write the formatted VCF files> [mandatory]
  -r REFERENCE, <Reference VCF file header to add> [mandatory]>
```

#### VCF Formatter 
This script will reform the VCF file to include to required/matching chromosome annotation. I.e. will add `chr` as the prefix of all values in the chr column.

The arguments for this script are as follows:
```
optional arguments:
  -i INPUT_FOLDER, <Input folder where the VCF files are hled> [mandatory]
  -o OUTPUT_DIR, <Ouput folder to write the formatted VCF files> [mandatory]
```

## Postprocessing data scripts
Output data wrangling scripts for the iSNP pipeline results.

### iSNP Wrangler
This script generates a more user-friendly output from the results of the iSNP pipeline.

The arguments for this script are as follows:
```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  -m MAPPING_FILE_PATH, --mapping_file_path MAPPING_FILE_PATH
  -t TARGET_ID, --target_id TARGET_ID
  -v VERBOSE, --verbose VERBOSE
  -c COMPARE_ON [COMPARE_ON ...], --compare_on COMPARE_ON [COMPARE_ON ...] 
```

The mapped outputs which are in directories with the patient vcf files as the identifier (e.g. mapped_<patient_name>.vcf) holds the intermediates files and the combined differences which include:
* `converted_mirna_gene_connections_mut.tsv`: uniport converted files for the miRNAs of the mutant sequences for this patient with header format seen below.
* `converted_mirna_gene_connections_wt.tsv`: uniport converted files for the miRNAs of the mutant sequences for this patient with header format seen below.
* `converted_tf_gene_connections_mut.tsv`: uniport converted files for the TF of the mutant sequences for this patient with header format seen below.
* `converted_tf_gene_connections_wt.tsv`: uniport converted files  for the TF of the mutant sequences for this patient with header format seen below.
* `differences.tsv`: all differences between the combined uniport converted TF and uniport converted miRNAs form the tf_differences.tsv and the  mirna_differences.tsv  files respectively. Calculated using the given key identifiers (e.g. source, target, snp and tool)
* `mirna_differences.tsv`: all differences between the combined miRNAs from the converted_mirna_gene_connections_wt.tsv  and the converted_mirna_gene_connections_mut.tsv files. Calculated using the given key identifiers (e.g. source, target, snp and tool)
* `tf_differences.tsv`: all differences between the combined TF from the converted_tf_gene_connections_wt.tsv  and the converted_tf_gene_connections_mut.tsv files. Calculated using the given key identifiers (e.g. source, target, snp and tool)

Each file holds the following headers now which are much easier to understand and the uniprot ids have been capitalised and converted. Complex are identified by / in the name and the
```
source	target	score	region	snp	file	mutated    tool
```
* `source`: is always either the mirna or tf.
* `target`: gene it interacts with
* `score`: prediction score from tool
* `region`: where annotation came from
* `snp`: rs snp id extracted from the vcf file of the patient 
* `file`: which file it came from i.e. if it was from the wt or mut patient (which genotype the interaction originated from)
* `mutated`: if the snp was mutated as per the vcf file saying is the snp is mutated or not (if there SNP is present or not)
* `tool`: tool which the prediction came from

Example of what the output would look like for each tool:
```
source	target	score	region	snp	file	mutated	tool
hsa-mir-425-5p	Q8WUI4	score:97.00;energy:-20.61	protein-coding	rs11168249 	mut	True	miranda
Q14938	Q92665	fimo-p-value:9.06122;fimo-q-value:9.06122	enhancer	rs941823	wt	False	fimo
O00570	Q86TN4	rsat_pvalue:0.0000220000000000	enhancer	rs559928	wt	False	rsat
Q13887	P33241	rsat_pvalue:0.0001000000000000	promoter	rs907611	wt	False	rsat
```

### iSNP Master Table Creator
This script generates a table, which contains all of the TF-target gene and miRNA-target gene interactions per patients.

The argument(s) for this script are as follows:
```
optional arguments:
  -h --help            show this help message and exit
  -pf PATIENT_FOLDER, --patient_folder PATIENT_FOLDER
```

The output files of this script are the follows:
* `affected_proteins_TFs_mirs.tsv`: The master table, which contains all of the TF-target gene and miRNA-target gene interactions for the patients.
* `affected_proteins.tsv`: .
* `SNPs.tsv`: This file contains the information, that which SNP can be found in the patients, separately.

Example of what the output would look like for each file:
* *`affected_proteins.tsv`*
-- `Source`: is always either the mirna
-- `Target`: gene it interacts with or tf
-- `SNP`: rs snp identifier extracted from the vcf file of the patient
-- `Mutated`: if the snp was mutated as per the vcf file saying is the snp is mutated or not (if there SNP is present or not)
-- `Interaction_source`: the source of the tool, where the interaction was found with
-- `the rest of the columns`: patient identifiers; 1 represents that the SNP is in the patient's genome (if `Mutated` is MUT, then the interaction exists; if `Mutated` is WT, then the interaction does not exist), 0 represents that the SNP is not there (if `Mutated` is WT, then the interaction exists; if `Mutated` is MUT, then the interaction does not exist)

```
Source	Target	SNP	Mutated	Interaction_source	patientID	patientID	patientID...
O60548	Q9HAV4	RS943072	WT	RSAT 	0	0	0
HSA-MIR-7843-5P	Q3MIT2	RS4560096	WT	MIRANDA	1	0	1
HSA-MIR-3972	Q5TA45	RS12103	WT	MIRANDA	0	0	0
Q9UBR4	P15036	RS4817986	MUT	RSAT	1	0	0
```
