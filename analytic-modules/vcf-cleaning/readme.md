# VCF Data Cleaner script

**Description:** 

This script will take a single input .vcf file or an input folder, that contains one or more .vcf files.
Then it will replace the locations of the SNP's in the .vcf file(s) to a new one according to the reference genom
file(s). In that case you have to use just the --input-file OR just the --input-folder parameter as input and the
--reference-genome-file parameter as reference.
If the script does not find a new location based on the reference genom, then the SNP will be deleted and will not be
in the output. (Example 1)

The script also can take a single SNP list. Then it will make a really new .vcf file from the input SNP list according
to the reference vcf file(s). In that case you have to use just the --input-snp-list parameter as input and the
--reference-vcf-file parameter as reference. In this case the script always make mutants for every SNP (in the output
vcf files of the second example, the last column will be always "1/1").
If the script does not find a SNP in the reference vcf file(s), then that SNP will not be in the output vcf file.
(Example 2)

The --output-folder parameter is mandatory in any cases!


**Parameters:**

-i, --input-file <path>                   : path to the input .vcf file [Optional]

-f, --input-folder <path>                 : path to a folder that contains all of the input .vcf files [Optional]

-s, --input-snp-list <path>               : path to the an input SNP list file [Optional]

-rg, --reference-genome-files <paths>     : comma separated list of paths, contains the reference genome information [Optional]

-rv, --reference-vcf-files <paths>        : comma separated list of paths, contains the reference genome information [Optional]

-o, --output-folder <path>                : path to an output folder that contain the changed VCF files [Mandatory]


**Exit codes**

Exit code 1: The input .vcf file does not exists!
Exit code 2: The input .vcf file is not a .vcf file!
Exit code 3: The specified input folder does not exists!
Exit code 4: The specified input SNP list does not exists!
Exit code 5: One of the reference genome file does not exists!
Exit code 6: One of the reference vcf file does not exists!
Exit code 7: The specified output folder does not exists!
Exit code 8: Have to give a single vcf input file or an input folder, which contains vcf files!


**Example**

Example 1: reference genome file (.bed)
```
chr21	40972191	40972192	rs2837891	0	+
chr21	41019094	41019095	rs8127758	0	+
chr21	41132428	41132429	rs1571731	0	-
```

Example 1: input .vcf file
```
##fileformat=VCFv4.3
##fileDate=20181211
##source=PLINKv2.00
##contig=<ID=1,length=247177331>
##contig=<ID=2,length=242626207>
##contig=<ID=3,length=199071699>
##contig=<ID=4,length=190980781>
##contig=<ID=5,length=180512666>
##contig=<ID=6,length=170723976>
##contig=<ID=7,length=158710966>
##contig=<ID=8,length=146264219>
##contig=<ID=9,length=140111181>
##contig=<ID=10,length=135040420>
##contig=<ID=11,length=134268352>
##contig=<ID=12,length=132228484>
##contig=<ID=13,length=114106016>
##contig=<ID=14,length=106296207>
##contig=<ID=15,length=100194179>
##contig=<ID=16,length=88655878>
##contig=<ID=17,length=78634367>
##contig=<ID=18,length=76115294>
##contig=<ID=19,length=63556292>
##contig=<ID=20,length=62322843>
##contig=<ID=21,length=46840115>
##contig=<ID=22,length=49503533>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	IW3134268_IW3134268
chr21	41265988	rs2837891	G	A	.	.	PR	GT	1/1
chr21	41312891	rs8127758	A	C	.	.	PR	GT	1/1
chr21	41426226	rs1571731	G	A	.	.	PR	GT	1/1
chr21	39373482	ccc-21-39373482-G-A	G	A	.	.	PR	GT	0/1
chr21	39373657	ccc-21-39373657-C-T	G	A	.	.	PR	GT	0/0
```

Example 1 parameters:
-i --input-file: example input .vcf file
-f --input-folder: -
-s --input-snp-list: -
-rg --reference-genome-files: example reference genome file
-rv --reference-vcf-files: -
-o --output: an output folder

The output .vcf file will be:
```
##fileformat=VCFv4.3
##fileDate=20181211
##source=PLINKv2.00
##contig=<ID=1,length=247177331>
##contig=<ID=2,length=242626207>
##contig=<ID=3,length=199071699>
##contig=<ID=4,length=190980781>
##contig=<ID=5,length=180512666>
##contig=<ID=6,length=170723976>
##contig=<ID=7,length=158710966>
##contig=<ID=8,length=146264219>
##contig=<ID=9,length=140111181>
##contig=<ID=10,length=135040420>
##contig=<ID=11,length=134268352>
##contig=<ID=12,length=132228484>
##contig=<ID=13,length=114106016>
##contig=<ID=14,length=106296207>
##contig=<ID=15,length=100194179>
##contig=<ID=16,length=88655878>
##contig=<ID=17,length=78634367>
##contig=<ID=18,length=76115294>
##contig=<ID=19,length=63556292>
##contig=<ID=20,length=62322843>
##contig=<ID=21,length=46840115>
##contig=<ID=22,length=49503533>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	IW3134268_IW3134268
chr21	40972192	rs2837891	G	A	.	.	PR	GT	1/1
chr21	41019095	rs8127758	A	C	.	.	PR	GT	1/1
chr21	41132429	rs1571731	G	A	.	.	PR	GT	1/1
```


Example 2: reference vcf file (.vcf)
```
1	10177	rs367896724	A	AC	.	.	RS=367896724;RSPOS=10177;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5747,0.4253;COMMON=1;TOPMED=0.76728147298674821,0.23271852701325178
1	10352	rs555500075	T	TA	.	.	RS=555500075;RSPOS=10352;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5625,0.4375;COMMON=1;TOPMED=0.86356396534148827,0.13643603465851172
1	10616	rs376342519	CCGCCGTTGCAAAGGCGCGCCG	C	.	.	RS=376342519;RSPOS=10617;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005040026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;KGPhase3;CAF=0.006989,0.993;COMMON=1
```

Example 2: input SNP list
```
rs367896724
rs32456
rs555500075
```

Example 2 parameters:
-i --input-file: -
-f --input-folder: -
-s --input-snp-list: example input SNP list
-rg --reference-genome-files: -
-rv --reference-vcf-files: example reference vcf file
-o --output: an output folder

The output .vcf file will be:
```
##fileformat=VCF
chr1	10177	rs367896724	A	AC	.	.	.	.	1/1
chr1	10352	rs555500075	T	TA	.	.	.	.	1/1
```
