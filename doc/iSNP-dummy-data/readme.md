# Dummy files format and some description


## SNP-identifier.txt ##

input for the "keep only UC (disease) specific SNPs from the VCF filtering"


## 0001-patient.vcf ##

input for the "keep only UC (disease) specific SNPs from the VCF filtering"


## 0001-patient-ucfilter.vcf ##

output of the "keep only UC (disease) specific SNPs from the VCF filtering"
input for the "filter mutations in protein coding regions" (VCF-filtering modul)
input for the "filter mutations in promoter regions" (VCF-filtering modul)


## anotation-protein-coding-regions.bed ##

input for the "filter mutations in protein coding regions" (VCF-filtering modul)

* Note: metadata is in the 4th "name" column
* Note: in the first column, the chromosome is an integer because of the bedtools-requirements


## annotation-promoter-regions.bed ##

input for the "filter mutations in promoter regions" (VCF-filtering modul)

* Note: metadata is in the 4th "name" column
* Note: in the first column, the chromosome is an integer because of the bedtools-requirements


## human-genome.fasta ##

input for the "mutated sequence generator"


## 0001-patient-protein-coding-regions.vcf ##

output of the "filter mutations in protein coding regions" (VCF-filtering modul)
input for the "mutated sequence generator"

* Note: metadata is in a comment at the end of the lines


## 0001-patient-promoter-regions.vcf ##

output of the "filter mutations in protein coding regions" (VCF-filtering modul)
input for the "mutated sequence generator"

* Note: metadata is in a comment at the end of the lines


## 0001-patient-msg-protein-coding-mutant.fasta ##

output of the "mutated sequence generator"
input for the "generate miRNA - gene connections"

* Note: metadata is in the header of the fasta file
* Note: region_lenght = 5 bases + snp + 5 bases


## 0001-patient-msg-protein-coding-wild.fasta ##

output of the "mutated sequence generator"
input for the "generate miRNA - gene connections"

* Note: metadata is in the header of the fasta file
* Note: region_lenght = 5 bases + snp + 5 bases


## 0001-patient-msg-promoter-mutant.fasta ##

output of the "mutated sequence generator"
input for the "generate TF - target gene connections"

* Note: metadata is in the header of the fasta file
* Note: region_lenght = 5 bases + snp + 5 bases


## 0001-patient-msg-promoter-wild.fasta ##

output of the "mutated sequence generator"
input for the "generate TF - target gene connections"

* Note: metadata is in the header of the fasta file
* Note: region_lenght = 5 bases + snp + 5 bases


## mirna.fasta ##

input for the "generate miRNA - gene connections"


## matrices.jaspar ##

input for the "generate TF - target gene connections"


## 0001-patient-miRNA-gene-connections-mutant.tsv ##

input for the "merge and enrich networks"
output of the "generate miRNA - gene connections"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html


## 0001-patient-miRNA-gene-connections-wild.tsv ##

input for the "merge and enrich networks"
output of the "generate miRNA - gene connections"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html


## 0001-patient-TF-targetgene-connections-mutant.tsv ##

input for the "merge and enrich networks"
output of the "generate TF - target gene connections"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html


## 0001-patient-TF-targetgene-connections-wild.tsv ##

input for the "merge and enrich networks"
output of the "generate TF - target gene connections"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html


## 0001-patient-enriched-network-mutant.tsv ##

output of the "merge and enrich networks"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html
* Note: the origin database name in the 28th column is "DataBaseName" for the time being, because I don't know what databases will we use.


## 0001-patient-enriched-network-wild.tsv ##

output of the "merge and enrich networks"

* Note: metadata is in the 26th, 27th, 28th columns | https://psicquic.github.io/MITAB27Format.html
* Note: the origin database name in the 28th column is "DataBaseName" for the time being, because I don't know what databases will we use.

