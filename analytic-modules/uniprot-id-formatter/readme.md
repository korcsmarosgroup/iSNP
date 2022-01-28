# Uniprot ID Formatter modul


**Description:**

The tool will take a mitab file as an input, and writes a new output mitab file. The output file is identical to
the input file, except that all the uniprot identifiers are transformed according to the optional parameters.

You don't have to give an output file or folder, because the script automatically make an output file, according to the
input network file. It add a "_formatted" tag to the end of the input file name.
(input file: test.tsv, then the output will be: test_formatted.tsv)
But if you want to specify a special output file, you can do it with the --output-network-file parameter, as well.
(Example 4)

The tool should search for uniprot ids in the following columns of the mitab file:
- columns 1-4 (IDs using colon separator uniprotac:a6pvc2-1)
- column 26-28 (NOX fully qualified IDs using semicolon separator: protein;uniprotac;a6pvc2-1)

The tool will not change the information what in the 3rd or 4th column (alternative ID's). The tool will change this
information only if the --no-isoform parameter was given and the actual uniprot id has an izoform. In this case
the tool save the uniprot ID with the izoform to the 3rd or 4th column as an alternative ID.


**Parameters**

-i, --input-network-file <path>                   : path to an input network file [Mandatory]

-lc, --lower-case                                 : if this parameter is given, all of the ID's will be lowercase, 
                                                    default: Null [Optional]

-uc, --upper-case                                 : if this parameter is given, all of the ID's will be uppercase, 
                                                    default: Null [Optional]

-ni, --no-isoform                                 : if this parameter is given, all the isoforms of the ID's will be
                                                    deleted, default: Null [Optional]

-o, --output-network-file <path>                  : path to an output network file [Optional]


**Exit codes**

Exit code 1: The input network file does not exists!
Exit code 2: Just one of these parameters (lower_case, upper_case) must be given!


**Notes**

1) if both --lover-case and --upper-case parameter is given, then return with non-zero exit code
2) if neither the --lover-case nor the --upper-case parameter is given, then the case will not be changed
3) if none of the three parameters are specified, then the output will be exactly same as the input


**Example**

Example input network file:
```
mirbase:hsa-mir-99b-5p	uniprotac:o95407	mirbase:hsa-mir-99b-5p	uniprotac:o95407	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:micro rna;mirbase;hsa-mir-99b-5p	end:gene;uniprotac;o95407	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:o95076	uniprotac:o95150-1	name:alx3	uniprotac:o95150-1	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;o95076	end:gene;uniprotac;o95150-1	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p15172-1	uniprotac:q6ste5	-	-	-	-	-	-	pudmed:15870273;pudmed:22068056	-	taxid:9606 ('homo sapiens')	psi-mi:'mi:407'(unknown)	psi-mi:'mi:463'(unknown)|psi-mi:'mi:2214'(unknown)	-	-	-	-	-	-	-	psi-mi:'mi:326'(unknown)	psi-mi:'mi:326'(unknown)	-	-	-	start:protein;uniprotac;p15172-1	end:protein;uniprotac;q6ste5	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

Example 1 parameters:
-i --input-network-file: example network file
-lc --lower-case: not given
-uc --upper-case: given
-ni --no-isoform: given
-o --output-network-file: not given

The output will be:
```
mirbase:hsa-mir-99b-5p	uniprotac:O95407	mirbase:hsa-mir-99b-5p	uniprotac:o95407	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:micro rna;mirbase;hsa-mir-99b-5p	end:gene;uniprotac;O95407	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:O95076	uniprotac:O95150	name:alx3	uniprotac:o95150-1|uniprotac:O95150-1	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;O95076	end:gene;uniprotac;O95150	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:P15172	uniprotac:Q6S8E5	uniprotac:P15172-1	-	-	-	-	-	pudmed:15870273;pudmed:22068056	-	taxid:9606 ('homo sapiens')	psi-mi:'mi:407'(unknown)	psi-mi:'mi:463'(unknown)|psi-mi:'mi:2214'(unknown)	-	-	-	-	-	-	-	psi-mi:'mi:326'(unknown)	psi-mi:'mi:326'(unknown)	-	-	-	start:protein;uniprotac;P15172	end:protein;uniprotac;Q6S8E5	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

Example 2 parameters:
-i --input-network-file: example network file
-lc --lower-case: not given
-uc --upper-case: given
-ni --no-isoform: not given
-o --output-network-file: not given

The output will be:
```
mirbase:hsa-mir-99b-5p	uniprotac:O95407	mirbase:hsa-mir-99b-5p	uniprotac:o95407	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:micro rna;mirbase;hsa-mir-99b-5p	end:gene;uniprotac;O95407	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:O95076	uniprotac:O95150-1	name:alx3	uniprotac:o95150-1	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;O95076	end:gene;uniprotac;O95150-1	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:P15172-1	uniprotac:Q6S8E5	-	-	-	-	-	-	pudmed:15870273;pudmed:22068056	-	taxid:9606 ('homo sapiens')	psi-mi:'mi:407'(unknown)	psi-mi:'mi:463'(unknown)|psi-mi:'mi:2214'(unknown)	-	-	-	-	-	-	-	psi-mi:'mi:326'(unknown)	psi-mi:'mi:326'(unknown)	-	-	-	start:protein;uniprotac;P15172-1	end:protein;uniprotac;Q6S8E5	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

```

Example 3 parameters:
-i --input-network-file: example network file
-lc --lower-case: not given
-uc --upper-case: not given
-ni --no-isoform: given
-o --output-network-file: not given

The output will be:
```
mirbase:hsa-mir-99b-5p	uniprotac:o95407	mirbase:hsa-mir-99b-5p	uniprotac:o95407	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:micro rna;mirbase;hsa-mir-99b-5p	end:gene;uniprotac;o95407	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:o95076	uniprotac:o95150	name:alx3	uniprotac:o95150-1|uniprotac:o95150-1	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;o95076	end:gene;uniprotac;o95150	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:p15172	uniprotac:q6s8e5	uniprotac:p15172-1	-	-	-	-	-	pudmed:15870273;pudmed:22068056	-	taxid:9606 ('homo sapiens')	psi-mi:'mi:407'(unknown)	psi-mi:'mi:463'(unknown)|psi-mi:'mi:2214'(unknown)	-	-	-	-	-	-	-	psi-mi:'mi:326'(unknown)	psi-mi:'mi:326'(unknown)	-	-	-	start:protein;uniprotac;p15172	end:protein;uniprotac;q6s8e5	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

Example 4 parameters:
-i --input-network-file: example network file
-lc --lower-case: not given
-uc --upper-case: given
-ni --no-isoform: given
-o --output-network-file: "results.tsv"

The output will be (results.tsv):
```
mirbase:hsa-mir-99b-5p	uniprotac:O95407	mirbase:hsa-mir-99b-5p	uniprotac:o95407	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:micro rna;mirbase;hsa-mir-99b-5p	end:gene;uniprotac;O95407	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:O95076	uniprotac:O95150	name:alx3	uniprotac:o95150-1|uniprotac:O95150-1	-	-	-	-	-	taxid:9606(homo sapiens)	taxid:9606(homo sapiens)	-	-	-	-	-	-	-	-	-	-	-	-	-	-	start:protein;uniprotac;O95076	end:gene;uniprotac;O95150	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
uniprotac:P15172	uniprotac:Q6S8E5	uniprotac:P15172-1	-	-	-	-	-	pudmed:15870273;pudmed:22068056	-	taxid:9606 ('homo sapiens')	psi-mi:'mi:407'(unknown)	psi-mi:'mi:463'(unknown)|psi-mi:'mi:2214'(unknown)	-	-	-	-	-	-	-	psi-mi:'mi:326'(unknown)	psi-mi:'mi:326'(unknown)	-	-	-	start:protein;uniprotac;P15172	end:protein;uniprotac;Q6S8E5	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```
