

# Executing the iSNP pipeline using command-line interface


##### Table of Contents  
[Overview](#overview)  
[Setting up execution environment](#setup)  
[Executing the iSNP workflow (single patient)](#singlepatient)  
[Executing the iSNP workflow (multiple patients)](#multipatient)
[Citation](#citation)


<a name="overview"/>

## Overview


The iSNP pipeline contains a set of analytical steps, each implemeted as a separate python command line tool.
You have three ways to executing the iSNP pipeline currently:
1) Executing the pipeline on your local machine for a single patient, using the `isnp.py` command line script. This script will run the analytical steps sequentially, using a single CPU core.
1) Executing the pipeline on your machine for multiple patients, using the `isnp_bunch_run.sh` script. This script  will call the `isnp.py` script on a set of patient specific input files. The 
script will use multiple CPU cores, executing many single-patient pipelines parallel.
1) Executing the pipeline in NavigOmiX. This way you can distribute the load among multiple machines, parallelizing on analytical step level. The down-side of this way is that it requires the 
maintenance and deployment of a full NavigOmiX solution. 

<a name="setup"/>

## Setting up the environment for local execution

In order to run the pipeline outside NavigOmiX, you have several prerequisite to set up.

### Set up docker environment (MacOS)
Create a python env which cannot talk to any external libraries:
Example (virtualenv):

`python3.7 -m virtualenv  --no-site-packages isnpdocker`

`source isnpdocker/bin/activate`

Example(conda):

`conda create --name isnpdocker`

`source activate isnpdocker`

Install dependencies:

`pip3 install docker==3.4.0` or `conda install -c conda-forge docker-py==3.4.0`

Set up docker host (make sense to put this to your `.bashrc` file, so you don't need to export the environment variable each time you open a new shell to execute iSNP):

```
export DOCKER_HOST=unix:///var/run/docker.sock
```

### preparing the patient data and all the reference databases

```
# the folder isnp-workflow/input-data is added to .gitignore, you can add all your files there
# it won't be committed to git

mkdir input-data
cd input-data

# ======== patient data =========
# copy the patient vcf file and the SNP-identifier.txt file


# ======== miRNA sequences (http://www.mirbase.org/ftp.shtml) ========
wget ftp://mirbase.org/pub/mirbase/22/mature.fa.gz
gunzip ./mature.fa.gz
mv ./mature.fa ./miRNA-sequences-mature-mirbase-v22.fasta
ln -sf ./miRNA-sequences-mature-mirbase-v22.fasta mirna.fasta

# ======== Jaspar TF-TFBS matrices (http://jaspar.genereg.net/downloads/) ========
JASPAR=JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt
curl http://jaspar.genereg.net/download/CORE/$JASPAR --output $JASPAR
ln -sf $JASPAR jaspar_matrices.txt


# ======== Human genoms (http://jaspar.genereg.net/downloads/) ========
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ./hg38.fa.gz
ln -sf ./hg38.fa human-genome.fasta


# ======== uniprot id mapping data ========
s3cmd get s3://sherlock-korcsmaros-group/project_zone/uniprot_id_mapping_2018_11_target_uniprotac/20181128_132321_00009_tic8w_ce91e363-12e7-4fb0-9cc1-cc635e1821b6.gz
gunzip 20181128_132321_00009_tic8w_ce91e363-12e7-4fb0-9cc1-cc635e1821b6.gz
mv 20181128_132321_00009_tic8w_ce91e363-12e7-4fb0-9cc1-cc635e1821b6 uniprot_id_mapping_2018_11_target_uniprotac.json
ln -sf ./uniprot_id_mapping_2018_11_target_uniprotac.json uniprot_id_mapping.json


# ======== omnipath database for PPI enrichment ========
s3cmd get s3://sherlock-korcsmaros-group/landing_zone/omnipath_0.7.111/interactor_a_tax_id=9606/omnipath.json
# TODO: convert it to mitab, named omnipath.tsv


# TODO: ======== Annotation files for protein coding and promoter regions ========


```

<a name="singlepatient"/>

## Executing the iSNP workflow for a single patient

You can execute the iSNP workflow using the `isnp.py` script. 
The script provides you help if you start it like: `isnp.py --help`.

You have to specify a folder where all the input files are uploaded and an other folder where all the outputs will be generated. You can also override all the default settings and input file names 
using optional parameters.

An example:

```
cd isnp-workflow
python3 ./isnp.py -i `pwd`/input-data -o `pwd`/output-data --tf_binding_matrices JASPAR_small_demo.txt --patient_vcf demo_short_2000_lines.vcf
```


<a name="multipatient"/>

## Executing the command line iSNP workflow for multiple patients 

The `isnp_bunch_run.sh` script is executing the iSNP workflow parallel on multiple patients. 
You can get the up-to-date usage information for the script: `/bin/bash ./isnp_bunch_run.sh --help`

```
Usage:

 isnp_bunch_run.sh [OPTIONS] [<relative path to directory contains vcf files, without input-data dir>]

options:
-t|--thread <number of threads>      : specify the number of isnp.py scripts will be executed parallel 

--no-docker-build                    : will not build the docker images (expects you already built them)  
                                       If you don't provide this option, the docker images will be built
                                       once before the starting of parallel execution.

--vcf-files                          : you can manually give a set of paths to VCF files separated by ','
                                       which will be executed by the script (the default behavior is to 
                                       provide a folder where all the VCF files will be taken)

-s|--separate-tasks-in-containers    : using this option, all the analytical tasks will be executed in 
                                       separate docker containers (default behavior: execute single 
                                       container for each patient and run all the steps for a patient in 
                                       it's container)


```

<a name="citation"/>

## Citation

To cite this work please use the following reference (the paper can be found here: https://doi.org/10.1101/692269):

```
A systems genomics approach to uncover patient-specific pathogenic pathways and proteins in a complex disease

Johanne Brooks, Dezso Modos, Padhmanand Sudhakar, David Fazekas, Azedine Zoufir, Orsolya Kapuy, Mate Szalay-Beko, 
Matthew Madgwick, Bram Verstockt, Lindsay Hall, Alastair Watson, Mark Tremelling, Miles Parkes, Severine Vermeire, 
Andreas Bender, Simon R. Carding, Tamas Korcsmaros

bioRxiv 692269; doi: https://doi.org/10.1101/692269
```
