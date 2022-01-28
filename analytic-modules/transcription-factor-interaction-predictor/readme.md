# TF interaction predictor for any potential target sequence of a gene

**Description:** 

This tool will take a set of DNA sequences in FASTA format. Each sequence is a region inside 
the coding region of a gene. For each input sequence we expect the following meta information 
to be stored, in standard NavigOmiX entity id format: 
- mandatory: gene id (we expect the sequence to be part of a promoter region)
- optional: SNP id (if the sequence was generated 'around' a SNP location)
 
The tool will go through all the known transcription factor binding matrices and calculate 
the probability of the connections. For this the module will use the RSAT tool (matrix-scan), which also
return the p-value (probability, strength) for the TF-gene interaction. The user can also 
specify a cut-off value parameter to skip the TF-gene interactions having p-value score 
below the cut-off parameter value.

The output will be a NavigOmiX network file (MITAB format). For each connections in 
the network we will store the following attributes:
- start node: TF id (in standard NavigOmix id format)
- end node: gene id (in standard NavigOmix id format)
- link attribute: p-value (higher p-value means higher probability for the connection)
- optional origin link attribute: SNP id (in standard NavigOmiX format, like origin:SNP;rsID;rs338276)

Note: the output network can contain duplicated links in the unlikely case of having 
multiple sequences belonging to the same mRNA in the input files

**Parameters**

--path_to_fasta: Path to the FASTA file with the nucleotide sequences.
--out_path: Path where the output mitab file is to be written.
--path_to_matrix: Path to the file with transcription matrix binding profiles. By default, the JASPAR database for vertebrates (http://jaspar.genereg.net/downloads/).
--format_matrix: Format of the file in path_to_matrix. transfac by default (recommended, as other types may cause an overheard).
--pval_threshold: Only those interactions with a p-value lower than this value will be output.


**Useful links:**
- http://rsat.sb-roscoff.fr/
- http://jaspar.genereg.net
- http://rsat.sb-roscoff.fr/matrix-scan_form.cgi
