# Mapping molecule IDs in molecular interaction networks

**Description:** 

This tool will take a single network file as input, and mapping all the node identifiers to a
specific target identifier type. Example: having an input protein interaction network with 
nodes named with Ensembl Gene ID and Ensembl Protein ID, I can create a similar network where 
all the node names are converted to Uniprot id.

The ID mapping is not necessarily unambiguous. E.g. if there is no valid mapping, or there
are more than one valid mappings for the given molecule, then the tool can skip or return 
with multiple uniprot interactions for a single original interaction.

The user can specify how the tools should behave if there is no mapping for a given input 
node. The tool can either remove this node (with all its connections) from the output file
or it can keep it without mapping.

The user can also define a list of molecule type to filter the nodes for mapping. If the
input file contains multiple kinds of molecules (MICRO RNA, PROTEIN, GENE, ...) and the user
specify a filter like 'PROTEIN,GENE', then all the MICRO RNA nodes will be copied to the output
network without doing any mapping on them. The valid molecule type names are given in the
MI ontology. We use the standard names, but in a case insensitive way.

Th output of the tool is a standard MITAB file (using standard NavigOmiX identifiers). The
tool is not sensitive to self-loops, it will copy them to the output file as well, if they
were present in the input. But there will be no duplicate links in the output file. 

All the metadata from the input connections will be copied to the output file as well without 
any change, only the fully qualified NavigOmiX identifiers will be updated with the correct 
molecule type and id type.

From the point of this tool, the "protein" and "gene" molecule types are the same. If you want
to map the genes and the proteins at the same time, it is enough to give just one of them (protein or gene).
(First example!)


**Parameters:** 

-i, --input <path>                : input MITAB file path [Mandatory]

-r, --remove <true|false>         : should the tool remove the nodes if no mapping can 
                                    be found. [Optional, default: true]

-f, --molecule-type-selector <str>  : comma separated list of molecule types given in standard 
                                    MI term names (case insensitive). The tool will apply
                                    id mapping only on the nodes with this type. If a node
                                    has different molecule type, then it will be copied to the
                                    output even if `--remove true` was used
                                    [Optional, default: no filtering, try to map all nodes]

-t, --target-id-type <str>        : target id type given in UniProtKB DBref (case insensitive)
                                    [Mandatory]

-m, --mapping-data <paths>        : comma separated list of paths, contains the ID mapping 
                                    information (it is possible to use multiple mapping data
                                    files for optimization) [Mandatory]

-n, --uniquename <true|false>      : should the tool map just the unique ID's or all of the ID's for the protein names 
                                     [Optional, default: false]

-o, --output <path>                : output MITAB file


**Exit codes**

Exit code 1: The specified mitab input file doesn't exists!
Exit code 2: Remove must be true or false!
Exit code 3: One of the specified mapping file doesn't exists!


**Additional info:**

1) The list of valid MI term names for molecule types: 
https://www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0313

2) The list of valid ID type names (we use UniProtKB DBref, basically the name of the ID type is 
the standard name of the database): https://www.uniprot.org/database/

3) The mapping data file format: a data lake compatible JSON file, contains a single JSON mapping 
data in each line. Example:

4) The tool also retains old identifiers stored in columns 3rd (Alternatives A) and 4th (Alternatives B)!

```json
{ “from_id_type”: “ensembl”, “to_id_type”: “uniprotac“, “from_id”: “ensg11867234247”, “to_id”: “ph39f7”}
{ “from_id_type”: “ensembl”, “to_id_type”: “uniprotac“, “from_id”: “ensg98236453842”, “to_id”: “ph39f7”}
{ “from_id_type”: “ensembl”, “to_id_type”: “uniprotac“, “from_id”: “ensg00034783463”, “to_id”: “phd734”}
{ “from_id_type”: “ensembl”, “to_id_type”: “uniprotac“, “from_id”: “ensg00034783463”, “to_id”: “ph863s”}
```

**Examples**

Example mapping file:
```
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010401", "to_id": "p00001"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010406", "to_id": "p00006"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010407", "to_id": "p00007"}
{"from_id_type": "uniprotac", "to_id_type": "uniprotac", "from_id": "g00001", "to_id": "p11111"}
{"from_id_type": "uniprotac", "to_id_type": "uniprotac", "from_id": "g00001", "to_id": "p22222"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010406", "to_id": "p10006"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010407", "to_id": "p10007"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010407", "to_id": "p11007"}
{"from_id_type": "genename", "to_id_type": "uniprotac", "from_id": "GAIR", "to_id": "p11007"}
{"from_id_type": "protein_name_recommended_short", "to_id_type": "uniprotac", "from_id": "GAIR", "to_id": "p11007"}
{"from_id_type": "protein_name", "to_id_type": "uniprotac", "from_id": "GAIR", "to_id": "p11007"}
{"from_id_type": "uniprotac", "to_id_type": "uniprotac", "from_id": "q29999", "to_id": "q29999"}
{"from_id_type": "uniprotac", "to_id_type": "uniprotac", "from_id": "gair29999", "to_id": "g29999"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010408", "to_id": "p00008"}
{"from_id_type": "ensembl", "to_id_type": "uniprotac", "from_id": "ensg00000010408", "to_id": "p10008"}
```

Example input network file:
```
ensembl:ensg00000010401	mirbase:mi0000002	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;ensembl;ensg00000010401	end:micro rna;mirna;mi0000002	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:q23004	ensembl:ensg00000010404	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;q23004	end:gene;ensembl;ensg00000010404	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
ensembl:ensg00000010406	uniprotac:g00001	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;ensembl;ensg00000010406	end:protein;uniprotac;g00001	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
ensembl:ensg00000010407	ensembl:ensg00000010408	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;ensembl;ensg00000010407	end:gene;ensembl;ensg00000010408	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
genename:GAIR	uniprotac:q29999	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;genename;GAIR	end:protein;uniprotac;q29999	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:gair29999	genename:GAIR	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;gair29999	end:gene;genename;GAIR	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

1) First example parameters:
-i --input: example input network file
-r --remove: true
-f --molecule-type-selector: protein
-t --target-id-type: uniprotac
-m --mapping data: example mapping file
-o --output: output network

The output network file will be:
```
uniprotac:p00001	mirbase:mi0000002	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00001	end:micro rna;mirbase;mi0000002	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00006	uniprotac:p11111	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00006	end:protein;uniprotac;p11111	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00006	uniprotac:p22222	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00006	end:protein;uniprotac;p22222	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10006	uniprotac:p11111	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10006	end:protein;uniprotac;p11111	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10006	uniprotac:p22222	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10006	end:protein;uniprotac;p22222	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:q29999	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:protein;uniprotac;q29999	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:g29999	uniprotac:p11007	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;g29999	end:gene;uniprotac;p11007	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

2) Second example parameters:
-i --input: example input network file
-r --remove: true
-f --molecule-type-selector: "micro rna"
-t --target-id-type: uniprotac
-m --mapping data: example mapping file
-o --output: output network

The output network file will be:
```
uniprotac:q23004	ensembl:ensg00000010404	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;q23004	end:gene;ensembl;ensg00000010404	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
ensembl:ensg00000010406	uniprotac:g00001	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;ensembl;ensg00000010406	end:protein;uniprotac;g00001	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
ensembl:ensg00000010407	ensembl:ensg00000010408	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;ensembl;ensg00000010407	end:gene;ensembl;ensg00000010408	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
genename:gair	uniprotac:q29999	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;genename;gair	end:protein;uniprotac;q29999	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:gair29999	genename:gair	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;gair29999	end:gene;genename;gair	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```
*Note: Of course there is no mapping data for `micro rna` molecules with `uniprotac` type, so we expect the mapping
module to remove all the micro rna interactions, since the `remove` parameter was set to ‘true’.


3) Third example parameters:
-i --input: example input network file
-r --remove: false
-t --target-id-type: uniprotac
-m --mapping data: example mapping file
-o --output: output network

The output network file will be:
```
uniprotac:p00001	mirbase:mi0000002	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00001	end:micro rna;mirbase;mi0000002	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:q23004	ensembl:ensg00000010404	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;q23004	end:gene;ensembl;ensg00000010404	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10006	uniprotac:p22222	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10006	end:protein;uniprotac;p22222	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10006	uniprotac:p11111	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10006	end:protein;uniprotac;p11111	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00006	uniprotac:p22222	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00006	end:protein;uniprotac;p22222	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00006	uniprotac:p11111	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00006	end:protein;uniprotac;p11111	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p00007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p00007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p10007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p10007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:p00008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:gene;uniprotac;p00008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:p10008	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:gene;uniprotac;p10008	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p11007	uniprotac:q29999	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:gene;uniprotac;p11007	end:protein;uniprotac;q29999	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:g29999	uniprotac:p11007	-	-	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;g29999	end:gene;uniprotac;p11007	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```
