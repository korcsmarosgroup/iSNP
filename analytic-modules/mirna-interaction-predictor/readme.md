# miRNA interaction predictor for any potential target sequence of a mRNA

**Description:** 

This tool will take a set of sequences in FASTA format. Each sequence contains a gene
sequence which can be a potential binding site for miRNAs. For each input sequence we 
expect the following meta information to be stored: 
- mandatory: gene id (we expect the sequence to be part of a gene coding region)
- optional: SNP id (if the sequence was generated 'around' a SNP location)
 
The tool will go through all the known miRNA sequences and calculate the probability
of the potential connection with the input sequence. For this it will 
(1) use the miRNA sequences from mirBase (this will be an input parameter for the module), and then
(2) use tool miRanda to calculate the probability (entropy / strength) of the connection
between the miRNA and the input sequence
(3) keep the miRNA - gene connection if the probability is above a given threshold

The output will be a NavigOmiX network file (MITAB format). For each connections in 
the network we will store the following attributes:
- start node: miRNA id (in standard NavigOmiX format, MIRNA;mirbase;mir-121)
- end node: gene id (in standard NavigOmiX format, GENE;EnsembleGeneId;ENSG00000139618)
- link attribute: entropy score (higher score means higher probability for the connection)
- optional origin link attribute: SNP id (in standard NavigOmiX format, like origin:SNP;rsID;rs338276)

Note: the output network can contain duplicated links in the unlikely case of having 
multiple sequences belonging to the same gene in the input files (this is handled by the MiTabHandler). 


**Useful links:**
- mirTarBase: http://mirtarbase.mbc.nctu.edu.tw/php/download.php
- miranda manual: http://cbio.mskcc.org/microrna_data/manual.html
- miranda download: http://www.microrna.org/microrna/getDownloads.do
- http://www.mirbase.org/index.shtml (in particular mature.fa file)


**Parameters:**

-m --mirna <miRNA sequences in fasta format> [mandatory]
-g --genomic <genomic sequences in fasta format> [mandatory]
-o --output <path to an output network set file> [mandatory]
-sc --score <score threshold, positive integer> [Optional]>
-e --energy <energy threshold [kcal/mol], negative integer> [Optional]
-s --strict <Demand strict 5' seed pairing (default: less strict, seed region 2-8) [Optional]>

