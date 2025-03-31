library('rtracklayer')
library('dplyr')
library('biomaRt')
library('stringr')

log <- file("SNP_check_local.log", open = "w")
sink(log, type = "message", append = TRUE)
message("\nStarting SNP_check_local.R script\n")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Wrong number of command line input parameters. Please check!")
}

snps <- read.table(args[1], sep = '\t', header = TRUE)
genome_version <- args[2]

if (genome_version == 'grch38') {
  genome_version <- ''
  gtf_file_URL <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v47.annotation.gtf.gz"
}
if (genome_version == 'grch37') {
  genome_version <- paste0(genome_version, '.')
  gtf_file_URL <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
}

destination_file <- args[3]
download.file(gtf_file_URL, destination_file)

message("GTF file downloaded\n")

variation = useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = paste0("https://", genome_version, "ensembl.org"))

message("Done with Ensembl\n")

rses <- snps$SNP

snps_and_place_of_the_SNPs <- getBM(
  attributes = c(
    'ensembl_gene_stable_id',
    "associated_gene",
    "chr_name",
    "chrom_start",
    "chrom_end",
    "chrom_strand",
    "minor_allele",
    "allele_1",
    "refsnp_id"
  ),
  filters = "snp_filter",
  values = rses,
  mart = variation
)

message("Done with getBM\n")

my_obj <- import(destination_file)
gene_annotation_original <- as.data.frame(my_obj)

gene_annotation_replace <- str_replace(gene_annotation_original$gene_id, "\\.\\d", "")
gene_annotation <- transform(gene_annotation_original, exon_number = as.numeric(exon_number)) 
gene_annotation$gene_id <- gene_annotation_replace
ensg_genes_in_SNPs <- unique(snps_and_place_of_the_SNPs$ensembl_gene_stable_id)
need_to_check_genes <- filter(gene_annotation, gene_id %in% ensg_genes_in_SNPs)

message("Starting with introns\n")

introns <-data.frame(
  "start" = 0 ,
  "end" = 0,
  "chromosome" = "",
  "gene_name" = "",
  "strand"  =  "+",
  "gene_id"  =  "",
  "type" = "first_intron"
)

for (gene in ensg_genes_in_SNPs) {
  
  if (length(gene_annotation$gene_id == gene) > 0) {
    strand_here = unique(as.character(gene_annotation[gene_annotation$gene_id == gene,"strand"]))
    
    if (length(strand_here) > 0) {
      #Checking positive strand
      if (strand_here == "+") {

        for (start_nucleotide in filter(gene_annotation, gene_id == gene & exon_number == 1)$end) {

          for (end_nucleotide in filter(gene_annotation, gene_id == gene & exon_number == 2)$start) {
            gene_name = unique(as.character(gene_annotation[gene_annotation$gene_id == gene, "gene_name"]))
            chromosome = unique(as.character(gene_annotation[gene_annotation$gene_id == gene, "seqnames"]))
            new_intron_data <- data.frame("start" = start_nucleotide, 
                                          "end" = end_nucleotide, 
                                          "chromosome" = chromosome,
                                          "strand"= strand_here,
                                          "gene_id" = gene,
                                          "gene_name" = gene_name,
                                          "type" = "first_intron")
            introns <- rbind(introns, new_intron_data)
          }

        }
        
      } else if (strand_here == "-") {
        #checking negative strand
        for (start_nucleotide in filter(gene_annotation, gene_id == gene & exon_number == 1)$start) {

          for (end_nucleotide in filter(gene_annotation, gene_id == gene & exon_number == 2)$end) {

            gene_name = unique(as.character(gene_annotation[gene_annotation$gene_id == gene, "gene_name"]))
            chromosome = unique(as.character(gene_annotation[gene_annotation$gene_id == gene, "seqnames"]))
            new_intron_data <- data.frame("start" = start_nucleotide, 
                                          "end" = end_nucleotide, 
                                          "chromosome" = chromosome,
                                          "strand"= strand_here,
                                          "gene_id" = gene,
                                          "gene_name" = gene_name,
                                          "type" = "first_intron")
            introns <- rbind(introns, new_intron_data)
          }

        }

      }

    }
    
  }

}

message("Starting with exons\n")

exons <- filter(need_to_check_genes, type == "exon")
needed_data <- exons[, c("start", "end", "seqnames", "strand", "gene_id", "gene_name", "type")]
colnames(needed_data) <- c("start", "end", "chromosome", "strand", "gene_id", "gene_name", "type")

potential_miRNA_location <- rbind(introns, needed_data)
our_snps <- data.frame("snp" = "" , "type" = "")

for (snp in unique(snps_and_place_of_the_SNPs$refsnp_id)) {
  pos = unique(snps_and_place_of_the_SNPs[snps_and_place_of_the_SNPs$refsnp_id == snp, "chrom_start"])
  chromosome = unique(snps_and_place_of_the_SNPs[snps_and_place_of_the_SNPs$refsnp_id == snp, "chr_name"])
  miRNA_locations <- filter(potential_miRNA_location, chromosome == chromosome & start <= pos & end >= pos)
  
  if (length(miRNA_locations) > 0) {
    types = unique(miRNA_locations$type)
    
    for (type in types) {
      row = data.frame("snp" = snp, "type" = type) 
      our_snps <- rbind(our_snps, row)
    }

  } 

}

message("Writing out the outputs\n")

write.table(
  our_snps,
  "output_SNPs_exons_introns.tsv",
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
)

write.table(
  our_snps$snp,
  "miRNA_filter.tsv",
  sep = '\t',
  row.names = FALSE,
  quote = FALSE,
  col.names = FALSE
)

message("\nThe SNP_check_local.R script finished\n")
