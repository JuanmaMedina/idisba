library("ExomeDepth")

# Exonic positions for the hg19 assembly
data("exons.hg19")
# print(head(exons.hg19))

# Prepare aggregate reference set from "controls":
# 1. Assumption: CNVs of interest are absent from the reference (EXCLUDE CASE SAMPLES)
# 2. Aggregate reference made from 5-10 exome samples (7 in HS_4000 cluster case)
# 3. Miss common CNVs if the calls are present in the aggregate reference (GOOD!)
# 4. 150-280 CNV calls per exome sample (2/3 are deletions)

# GenomicRanges object with the counts data per sample & per exon for 7 HS control exomes
HS_7CON_counts <- getBamCounts(
  
  # Define the exonic regions
  bed.frame = exons.hg19,
  
  # List of .bams (and .bais) as vector of characters
  bam.files = c(as.character(read.table("Almacen/seq/familias_cardio_SEQ/HS_4000_ED_7CON_SORTED.txt", header = FALSE)[, 1])),
  index.files = c(as.character(read.table("Almacen/seq/familias_cardio_SEQ/HS_4000_ED_7CON_SORTED.BAI.txt", header = FALSE)[, 1])),
  
  min.mapq = 20,
  
  # Remove the string ’chr’ from the  chromosome names of the target .bed file.
  include.chr = TRUE,
  
  # To compute GC content
  referenceFasta = "/home/genomica/Almacen/ref/hg19/genome.fa")
  
## Examine this warning ##
## WARNING: Each of the 2 combined objects has sequence levels not in the other ##

# This function reads a case sample exome and creates a GenomicRanges object with the counts
# data per sample & per exon with the getBamCounts function form ExomeDepth pacjage
case_HS_counts <- function(sorted_bam, sorted_bai){
  
  # Path to the input sorted.bam
  case_bam <- paste("Almacen/seq/familias_cardio_SEQ/", 
                     sorted_bam, 
                     "_sorted.bam", 
                     sep = "")
  
  # Path to the input sorted.bam.bai
  case_bai <- paste("Almacen/seq/familias_cardio_SEQ/", 
                    sorted_bai, 
                    "_sorted.bam.bai", 
                    sep = "")
  
  # Read count data calculation
  counts <- getBamCounts(
  
  # Define the exonic regions
  bed.frame = exons.hg19, 
    
  # .bam / .bai as character vectors
  bam.files = as.character(case_bam), 
  index.files = as.character(case_bai), 
    
  min.mapq = 20, 
    
  # Remove 'chr' string from the chromosome names of the target .bed file    
  include.chr = TRUE, 
    
  # To compute GC content
  referenceFasta = "/home/genomica/Almacen/ref/hg19/genome.fa")
  
  return(counts)
  
}

# List of affected cases
fam2_51609 <- case_HS_counts(sorted_bam = "familia_2/51609", 
                             sorted_bai = "familia_2/51609")
