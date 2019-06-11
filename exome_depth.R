library("ExomeDepth")

# Exonic positions for the hg19 assembly (no 'chr' string)
data("exons.hg19")

# Notions about the preparation of the aggregate reference set from the "control" exomes:
# a) Assumption: CNVs of interest are absent from the reference (EXCLUDE CASE SAMPLES)
# b) Aggregate reference made from 5-10 exome samples (7 in HS_4000 cluster case)
# c) Miss common CNVs if the calls are present in the aggregate reference (GOOD!)
# d) 150-280 CNV calls per exome sample (2/3 are deletions)

## 1. Depth of coverage calculations (based on samtools, takes a bit)
# GenomicRanges object with the counts data per sample & per exon for 7 HS control exomes
HS_7CON_counts <- getBamCounts(
  
  # Define the exonic regions
  bed.frame = exons.hg19,
  
  # List of .bams and .bais as vectors of characters
  bam.files = c(as.character(read.table("Almacen/seq/familias_cardio_SEQ/HS_4000_ED_7CON_SORTED.txt", header = FALSE)[, 1])),
  index.files = c(as.character(read.table("Almacen/seq/familias_cardio_SEQ/HS_4000_ED_7CON_SORTED.BAI.txt", header = FALSE)[, 1])),
  
  min.mapq = 20,
  
  # Setting TRUE here removes the string ’chr’ from the chromosome names of the target .bed
  # file (which they do not have, what is this for exactly?)
  include.chr = T,
  
  # Inpute the reference genome to compute GC content
  referenceFasta = "/home/genomica/Almacen/ref/hg19/genome.fa")

## WARNING: Each of the 2 combined objects has sequence levels not in the other ##
## This warning seems to be affecting at some nomenclature problem between chroms X,Y,M
## and as I am not comparing sexual chromosomes, I am going to let it as is for now

# Convert into DF, removing again the 'chr' string
HS_7CON_counts <- as(HS_7CON_counts[, colnames(HS_7CON_counts)], 'data.frame')
HS_7CON_counts$chromosome <- gsub(as.character(HS_7CON_counts$space), 
                                     pattern = 'chr', 
                                     replacement = '')
  
# This function reads a case sample exome and returns a GenomicRanges object with the counts
# data per sample & per exon with the getBamCounts function from ExomeDepth package
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
    
  # Remove 'chr' string from the chromosome names of the target .bed file (again, ?)    
  include.chr = TRUE, 
    
  # To compute GC content
  referenceFasta = "/home/genomica/Almacen/ref/hg19/genome.fa")
  
  # Convert into DF
  counts_df <- as(counts[, colnames(counts)], 'data.frame')
  
  counts_df$chromosome <- gsub(as.character(counts_df$space), 
                               pattern = 'chr', 
                               replacement = '')
  
  return(counts_df)
  
}

## List of affected cases ## 
# Family 2, case 51609
fam2_51609 <- case_HS_counts(sorted_bam = "familia_2/51609", 
                             sorted_bai = "familia_2/51609")
# Family 2, case 51610
fam2_51610 <- case_HS_counts(sorted_bam = "familia_2/51610", 
                             sorted_bai = "familia_2/51610")
# Family 3, case 20592
fam3_20592 <- case_HS_counts(sorted_bam = "familia_3/20592", 
                             sorted_bai = "familia_3/20592")
# Family 3, case 46632
fam3_46632 <- case_HS_counts(sorted_bam = "familia_3/46632", 
                             sorted_bai = "familia_3/46632")
# Family 4, case 52718
fam4_52718 <- case_HS_counts(sorted_bam = "familia_4/52718", 
                             sorted_bai = "familia_4/52718")
# Family 5, case 23281
fam5_23281 <- case_HS_counts(sorted_bam = "familia_5/23281", 
                             sorted_bai = "familia_5/23281")
# Family 5, case 22010
fam5_22010 <- case_HS_counts(sorted_bam = "familia_5/22010", 
                             sorted_bai = "familia_5/22010")

## 2. Construction of the best reference dataset and CNV calling

# List of reference exomes IDs contained in the DF
reference_exomes <- c("X1754_sorted.bam", "X1755_sorted.bam", "X51611_sorted.bam", 
                      "X46631_sorted.bam", "X52716_sorted.bam", "X23279_sorted.bam", 
                      "X23280_sorted.bam")

# This function takes the counts of a case and reference DFs, combines the reference 
# exomes to build the best reference dataset and performs CNV calling based on depths LHRs
# and an underlying HMM
CNV_caller <- function(case_exome, reference){
  
  # Select the most appropriate reference sample
  case_test <- select.reference.set(test.counts = case_exome, 
                                    reference.counts = as.matrix(reference[, reference_exomes]), 
                                    bin.length = (reference$end - reference$start)/1000,
                                    n.bins.reduced = 10000)
  
  # Build the reference dataset
  case_matrix <- as.matrix(reference[, case_test$reference.choice])
  
  # ?
  selected <- apply(X = case_matrix, MARGIN = 1, FUN = sum)

  # Apply beta-binomial model to the full set of exons to create an ExomeDepth object
  # including the computed LHs
  comparison <- new('ExomeDepth',
                    test = case_exome,
                    reference = selected,
                    formula = 'cbind(test, reference) ~ 1')
  
  # CNV calling
  # This function fits a HMM to the read depth data with three hidden states (Nor, Del, Dup)
  comparison <- CallCNVs(x = comparison,
                         chromosome = reference$space,
                         start = reference$start,
                         end = reference$end,
                         name = reference$names,
                         
                         # Transition p. from the normal CN state to a deletion/duplication
                         transition.probability = 10^-4, 
                         
                         # Used in Viterbi to compute the transition between states
                         expected.CNV.length = 50000)

  return(comparison)
      
}

# For meaningful results, the correlation between reference and test should be > 0.97
# Sort by Bayes factor = log10 LHR (CNV call (ALT. HIP.) / null call (NULL HIP.))
# e.g. BF = 38 means the presence of a CNV is 10**38 times more likely than no CNV
# This can be unconvincing, specially for short CNVs
# BF should be used to gain specificity over the aggressive sensitivity of this method

# fam2_51609 case: 137
calls_51609 <- CNV_caller(case_exome = fam2_51609$X51609_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_51609@CNV.calls)
head(calls_51609@CNV.calls[order(calls_51609@CNV.calls$BF, decreasing = T),])

# fam2_51610 case: 168
calls_51610 <- CNV_caller(case_exome = fam2_51610$X51610_sorted.bam,  
                          reference = HS_7CON_counts)
nrow(calls_51610@CNV.calls)
head(calls_51610@CNV.calls[order(calls_51610@CNV.calls$BF, decreasing = T),])

# fam3_20592 case: 180
calls_20592 <- CNV_caller(case_exome = fam3_20592$X20592_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_20592@CNV.calls)
head(calls_20592@CNV.calls[order(calls_20592@CNV.calls$BF, decreasing = T),])

# fam3_46632 case: 168
calls_46632 <- CNV_caller(case_exome = fam3_46632$X46632_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_46632@CNV.calls)
head(calls_46632@CNV.calls[order(calls_46632@CNV.calls$BF, decreasing = T),])

# fam4_52718 case: 151
calls_52718 <- CNV_caller(case_exome = fam4_52718$X52718_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_52718@CNV.calls)
head(calls_52718@CNV.calls[order(calls_52718@CNV.calls$BF, decreasing = T),])

# fam5_23281 case: 148
calls_23281 <- CNV_caller(case_exome = fam5_23281$X23281_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_23281@CNV.calls)
head(calls_23281@CNV.calls[order(calls_23281@CNV.calls$BF, decreasing = T),])

# fam5_51609 case: 131
calls_22010 <- CNV_caller(case_exome = fam5_22010$X22010_sorted.bam, 
                          reference = HS_7CON_counts)
nrow(calls_22010@CNV.calls)
head(calls_22010@CNV.calls[order(calls_22010@CNV.calls$BF, decreasing = T),])


# calls <- c(calls_51609@CNV.calls, calls_51610@CNV.calls, calls_20592@CNV.calls, 
#            calls_46632@CNV.calls, calls_52718@CNV.calls, calls_23281@CNV.calls, 
#            calls_22010@CNV.calls)
# 
# 
# for (i in calls) {
#   
#   #   i$length <- (i$end - i$start)
# }





                         

      ####


# # Further annotation of CNVs with Conrad et al 2010 Nature identified CNVs
# # (not suitable to identify rare CNVs)
# data(Conrad.hg19)
# 
# calls_51609_conrad <- AnnotateExtra(x = calls_51609,
#                                     reference.annotation = Conrad.hg19.common.CNVs,
#                                     min.overlap = 0.5,
#                                     column.name = 'Conrad.hg19')

# # Further annotation with the exonic information
# exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
#                                              IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
#                                              names = exons.hg19$name)
# 
# calls_51609_exons <- AnnotateExtra(x = calls_51609,
#                                    reference.annotation = exons.hg19.GRanges,
#                                    min.overlap = 0.0001,
#                                    column.name = 'exons.hg19')

# # This plot represents the ratio between observed and expected read depth, with the 95% CI
# # marked by a grey shaded area
# plot(calls_51609,
#       sequence = '1',
#       xlim = c(110230497 - 100000, 110235917 + 100000),
#       count.threshold = 10,
#       main = 'Gene' ,
#       cex.lab = 0.8,
#       with.gene = TRUE)

