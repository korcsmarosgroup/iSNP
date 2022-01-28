
# Naming requirements of bedtools for filtering

Bedtools implements the unctions *intersect* for filtering between common files format (VCF/bed/GFF) (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html). Also, it has the function *getfasta* that extracts the nucleotide sequence from a FASTA file according to the positions defined in a .bed file (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html). Both seem quite efficient and easy to integrate in our workflow using *subprocess*.

## *getfasta*

The *getfasta* function requires as input the path to a .bed file and the one to a .fasta file:

```
> $ bedtools getfasta [OPTIONS] -fi *path_to_FASTA* -bed *path_to_BED*
```
Then, the function extracts the nucleotide sequence according to the positions defined in the .bed file, creating a new FASTA file (that can be sent to stdout or written).  It has several options, maybe the most relevant for us is *-name*, that allows to name each nucleotide sequence in the output FASTA according to the *name* field in the BED file (4th column). This will allow us to propagate IDs.

However, in order for the function to work, it must know where to look in the FASTA file. For that reason, the *chrm* field of the .BED file (first column) must match the description line of the FASTA file (the text after the ">").

For instance, I have this bed file (*short_bed.bed*):
```
chr19	80	121	.
chr19	14	53	rs6054257
```

And this FASTA file (the one uploaded by Balazs):

```
>human | reference genom | 1000bp
AAAGACTGAAATACAACCGCTAAGCCTTGCCCTATTTAGTAAGTTCTAGAGCGAACGGCGGTTTTTACAACATGAAAAAC
CGTTACATTCCTCGATAAGTAGAGTACGAAGGTGCTCTCCCTATTACGTGCCGGGTCACCTTTAGAGCGGACCGTGTATA
...
```

If I try to call *getfasta*, I get the following warning:

```
my_test_files agarre$ bedtools getfasta -fi human-annotation.fasta -bed short_bed.bed -name -fo stdout

WARNING. chromosome (chr19) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr19) was not found in the FASTA file. Skipping.
```

The reason for this is that the *chrm* field in the .bed file ("chr19") does not match the the one in the .fasta. I have modified the header of the FASTA so that it matches the bed (*short.fasta*):

```
>chr19
AAAGACTGAAATACAACCGCTAAGCCTTGCCCTATTTAGTAAGTTCTAGAGCGAACGGCGGTTTTTACAACATGAAAAAC
CGTTACATTCCTCGATAAGTAGAGTACGAAGGTGCTCTCCCTATTACGTGCCGGGTCACCTTTAGAGCGGACCGTGTATA
CTCAAGCGTTCAGGGGGTTCTCCCTTCGTAAAGATCATCGTTCCGAGTGTAGC
...
```

Now it works:

```
bedtools getfasta -fi short.fasta -bed short_bed.bed -name -fo stdout
>.
CGTTACATTCCTCGATAAGTAGAGTACGAAGGTGCTCTCCC
>rs6054257
AACCGCTAAGCCTTGCCCTATTTAGTAAGTTCTAGAGCG
```

Note that each sequence in the output FASTA is named according to the *name* in the .bed file.

So, **the name field in the .bed file must be identical to the descriptor in the FASTA**. 

## intersect

The *intersect* function requires the path to one reference file and one (or more) file(s) defining the filter:

```
bedtools intersect [OPTIONS] -a <FILE> -b <FILE1, FILE2, ..., FILEN>
```

The function extracts the information of the original file (-a) according to the filter defined (-b) and outputs a file in the same file as -a. 

Again, in order to do this, the namings between all the files must match. For instance, we have this .vcf file (0001-patient.vcf):

```
...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	7	snp1	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
1	12	snp2	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
1	35	snp3	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
1	53	snp4	G	A	31	.	AC=1;AN=2;DP=196;NS=63
1	182	snp5	C	T	25	.	AC=1;AN=2;DP=275;NS=66
...
```

that we want to filter according to a .bed file (annotation-promoter-regions.bed):

```
chrom   chromStart  chromEnd    gene-ID
chr1	31	41	ENSG00000010401
chr1	148	158	ENSG00000010402
chr1	284	294	ENSG00000010403
chr1	399	409	ENSG00000010404
chr1	514	524	ENSG00000010405
chr1	648	658	ENSG00000010406
chr1	751	761	ENSG00000010407
chr1	880	890	ENSG00000010408
```

If we call *bedtools intersect*...

```
bedtools intersect -a 0001-patient.vcf -b annotation-promoter-regions.bed 

***** WARNING: File 0001-patient.vcf has inconsistent naming convention for record:
1	7	snp1	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65

***** WARNING: File 0001-patient.vcf has inconsistent naming convention for record:
1	7	snp1	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
```

*bedtools* is not able to do the filter because the names of the chromosomes do not match in both files. This can be corrected by just removing the "chr" in the .bed file. Note that this is the only valid modification, because the #CHROM field in the .vcf must be an integer. I have made a file (annotation-short.bed) with this modification:

```
chrom   chromStart  chromEnd    gene-ID
1	31	41	ENSG00000010401
1	148	158	ENSG00000010402
1	284	294	ENSG00000010403
1	399	409	ENSG00000010404
1	514	524	ENSG00000010405
1	648	658	ENSG00000010406
1	751	761	ENSG00000010407
1	880	890	ENSG00000010408
```

And now...

```
bedtools intersect -a 0001-patient.vcf -b annotation-short.bed 

1	35	snp3	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
1	292	snp7	T	A	67	.	AC=1;AN=2;DP=225;NS=67
1	521	snp11	T	C	30	.	AC=1;AN=2;DP=205;NS=64
1	753	snp16	C	T	51	.	AC=1;AN=2;DP=206;NS=68
1	885	snp17	T	C	35	.	AC=1;AN=2;DP=244;NS=61
```

It works :D (If we don't want the header of the .vcf file, we can add -header)

However, the gene-ID that was defined in the .bed file is lost. If we want to keep it, we need to add -wa -wb to the call to *intersect*:

```
bedtools intersect -a 0001-patient.vcf -b annotation-short.bed -wa -wb

1	35	snp3	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65	1	31	41	ENSG00000010401
1	292	snp7	T	A	67	.	AC=1;AN=2;DP=225;NS=67	1	284	294	ENSG00000010403
1	521	snp11	T	C	30	.	AC=1;AN=2;DP=205;NS=64	1	514	524	ENSG00000010405
1	753	snp16	C	T	51	.	AC=1;AN=2;DP=206;NS=68	1	751	761	ENSG00000010407
1	885	snp17	T	C	35	.	AC=1;AN=2;DP=244;NS=61	1	880	890	ENSG00000010408
```

On the one hand, this outputs additional columns with the information from the .bed file, retaining the ID. On the other hand, I don't think this file is now compliant with the .vcf format, so we should be careful with how we use it in later steps. When we have a better defined way of propagating IDs through the workflow, I will figure this out.


