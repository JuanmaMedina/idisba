#### CNMOPS PIPELINE (CN estimation by a Mixture Of PoissonS) #### 
# This pipeline uses R package cn.mops to predict CNVs from WGS data

## Theoretical premises

# 1. Models depth across samples at each genomic segment
# Hence, cnmops predictions remain unaffected by read count variations across chrs

# 2. Bayesian framework
# Guided decomposition of depth variation across samples into integer CN and noise
# Maximization of posterior (CN + noise) through EM
# LH assumes a linear dependency between average read counts in segment and its CN
# Dirichlet prior of a CN = 2 for all samples

# Bayesian paradigm:
# "The more the read count data drives the posterior away from the prior, the more likely
#                           the data is caused by a CNV"

library("cn.mops")

setwd("Almacen/seq/familias_cardio_SEQ/WGS")

# Input data from BAM files
BAMFiles <- list.files(pattern = ".bam$")

# Chrs to analyze (exclusion of chrM reference sequence)
chrs <- readLines("~/Juanma/CNV_analysis/WGS/chr_List.txt")

# Time 1
matrix_build_start_time <- Sys.time()

# Generate read count matrix (~ 7 hours)
# Default window length = 25.000
bamDataRanges <- getReadCountsFromBAM(BAMFiles = BAMFiles,
                                      sampleNames = c("06690","06692","06693","1752","1756","1765"),
                                      refSeqNames = chrs,
                                      WL = 25000,
                                      parallel = 4)

# Time 2
matrix_build_end_time <- Sys.time()

# Compute time
print(matrix_build_end_time - matrix_build_start_time)

# CNV calling
res <- cn.mops(input = bamDataRanges,
               # Expected FC of the CN classes
               I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
               # CN classes
               classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8"),
               # Strength of the prior on the results (LH?)
               # The higher, the more samples will be assumed to have CN=2
               priorImpact = 1, 
               cyc = 20, 
               parallel = 0, # not necessary    
               # Normalization of read counts through a Poisson model: GOOD
               norm = 1, normType = "poisson",
               
               # CNVs above and below these thresholds will be called as 
               # amplifications and deletions, respectively (in log2 scale). I.e.,
               # CNVs expressing >= 1.41 times, or <= 0.53 will be called as
               # amplifications and deletions, respectively (in linear scale)
               upperThreshold = 0.5,  # 1.41 in linear scale 
               lowerThreshold = -0.9, # 0.53 in linear scale
               
               # Default faster segmentation algorithm. NB:
               # We can perform CBS here with "DNAcopy"
               segAlgorithm = "fast",
               # Number of segments the CNV should span to be called 
               # Each segment = 25.000
               # CNV call > 50.000 bps
               minWidth = 2, 
               # If all samples are below this value in a segment, the algorithm
               # will not be applied over it
               minReadCount = 5,
               returnPosterior = T)

# Calculate integer CN (CNVdetectionResult object, split in regions and individual CNVs (361))
CNVs <- calcIntegerCopyNumbers(res)

# GenomicRanges objects
regions <- cnvr(CNVs)
individual <- cnvs(CNVs)


# Theoretical note about median/mean columns in the GenomicRanges metadata
# These columns correspond to the median/mean individual I/NI (informative/non-informative) 
# call for this CN segment. The higher, the more informative ("significant") is the call

# Build DF
x <- as.data.frame(individual)[, c(1,2,3,6,4,9)]

# Modify sample names
x$sampleName <- paste(x$sampleName, ".bam", sep = "")

# Further specify CNV type in new column
x$type <- ifelse(x$CN == "CN0" | x$CN == "CN1", "DEL", "AMP")


# Writing output
write.table(x, "~/Juanma/CNV_analysis/WGS/cnmops_CNVs.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")

# BASH bedtools: overlap with SD/SA

# bedtools intersect -wa -f 0.8 -v -a cnmops_CNVs.bed -b 
#   ~/Juanma/UCSC_tracks_SVs/segmental_duplications_hg19.bed | 
#   bedtools intersect -wa -f 0.8 -v -a - -b ~/Juanma/UCSC_tracks_SVs/self_alignments_hg19.bed 
#   > calls_cnmops_HS_noSDSA_50k.bed


### Plots ###

# Error in .External.graphics --> https://support.bioconductor.org/p/101191/
dev.off()

# Plot the log normalized read counts and the detected segments as a segmentation plot
# x-axis represents genomic location, and y-axis represents the LR of the read counts
# (black) and the CN of each segment (red)
segplot(CNVs, 
        # Index of sample to be plotted
        sampleIdx = 1, 
        # Chromosomes to be plotted
        seqnames = "chr1")

# Plot segments, specified by index 
# Blue lines mark samples with a CN loss, and red lines mark samples with a CN gain
res@cnvr@elementMetadata@nrows
plot(res, which = 140)
