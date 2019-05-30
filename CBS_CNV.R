# Circular binary segmentation (CBS) of DNA into regions of estimated 
# equal CN to detect regions with abnormal CN

# Start pipeline by running VarScan between the control and the case, i.e.

# java -jar ~/bin/VarScan.v2.3.9.jar copynumber 
# ~/Almacen/seq/familias_cardio_SEQ/HS_4000.mpileup 
# ~/Almacen/seq/familias_cardio_SEQ/familia_2/51609.mpileup 
# controlVS51609 --min-coverage 20 --min-base-qual 20 --min-map-qual 20

# Libraries
library("DNAcopy")
library("tibble")   # For column reordering 
library("ggplot2")

# Determine chromosomes (will be used repeatedly across the script)
chroms <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
            'chr9','chr10','chr11','chr12','chr13','chr14','chr15',
            'chr16','chr17','chr18','chr19','chr20','chr21','chr22',
            'chrX','chrY')


# This function reads a .copynumber file, builds a CNA object and computes the CBS
# after removing outliers to compute the p-vals for the change points found, 
# writing the output to a .txt file

CBS_detection <- function(copynumber_in, CBS_out){
  
  # Path to the input .copynumber file
  case_file <- paste("~/Juanma/CNV_analysis/", copynumber_in, 
                     ".copynumber", sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
  
  # Copy Number Array object (matrix with log ratio changes)
  CNA.object <- CNA(genomdat = input[,7], chrom = input[,1],
                    maploc = input[,2], data.type = 'logratio')
  
  # Detect outliers and smooth the data prior to CBS
  CNA.object.smoothed <- smooth.CNA(CNA.object)
  
  # CBS
  CNA.object.smoothed.seg <- segment(CNA.object.smoothed, verbose = 0)
  
  # Plot with data ordered by chromosome and map positions.
  
  # Red dots correspond to the segment (CNV) mean LR change values of the case when 
  # compared to the control. I.e., if the CNV corresponds to a deletion in the affected
  # compared to the normal case, it will be marked with a negative mean LR change value 
  # and it will be shown below the grey line (which will be plotted as a CYAN dot/line
  # later on). On the contrary, if the CNV corresponds to an amplification, it will be
  # marked with a positive mean LR change value and it will be shown above the grey line 
  # (SALMON dot)
  
  # print(plot(CNA.object.smoothed.seg, plot.type="w"))
  
  # Compute p-vals and CIs for the change points found by the CBS. The data
  # are permuted to the left and right of the identified change point and the
  # location of the maximal binary segmentation is computed. The CI is given
  # by the quantiles of the permutation distribution of the locations
  
  # IMPROVE: Insert correction for multiple testing here?
  seg.pvalue <- segments.p(CNA.object.smoothed.seg, 
                           ngrid = 100, 
                           tol = 1e-6, 
                           alpha = 0.05, 
                           search.range = 100, 
                           nperm = 1000)
  
  # 2 additional fixes for posterior mergeSegment.pl script
  
  # a) Add a second column with "chrom"
  seg.pvalue <- add_column(seg.pvalue, fixed = "chrom", .after = 1)
  
  # b)  Keep only true chromosomes, not seq-generated artifact denominations 
  seg.pvalue.clean <- seg.pvalue[seg.pvalue$chrom %in% chroms, ]
  
  # Path to the output .txt file
  filename_out <- paste("~/Juanma/CNV_analysis/", CBS_out, ".txt", sep = "")
  
  # Write out file
  write.table(seg.pvalue.clean, filename_out,
              row.names = F, col.names = F, quote = F, sep = "\t")

}


# This function points a CBS file path to mergeSegments.pl script, calls it and
# writes the output in a pointed path

seg_merger <- function(CBS_in, merged_out){
  
  # Path to the input CBS file
  CBS_file <- paste("~/Juanma/CNV_analysis/", 
                    CBS_in, 
                    ".txt", 
                    sep = "")
  
  # Path to the output merged file
  merged_file <- paste("~/Juanma/CNV_analysis/", 
                       merged_out, 
                       sep = "")
  
  # bash:
  # perl mergeSegments.pl CBS_fam2_11vs9_pvals.txt 
  # --ref-arm-sizes refArmSize.txt 
  # --output-basename merged_11vs9_pvalue
  
  # MergeSegments.pl call
  system(paste("perl mergeSegments.pl", 
               CBS_file,
               "--ref-arm-sizes",
               "~/Juanma/CNV_analysis/refArmSize.txt",
               "--output-basename",
               merged_file,
               sep=" "))
  
}


# This function reads a mergeSegments.pl output file (.events.tsv), removes 
# identified CNVs non classified as deletions/amplifications and plots them after 
# some filtering by size / p-value, presenting their coordinate along the x-axis 
# and their CNV LR mean along the y-axis, coloring them by their SV nature

CNV_dotplot_LRmean <- function(merged_file){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
  
  # Add column with indicative middle coordinate
  input$chr_mid <- round((input$chr_start + 
                                    input$chr_stop) / 2, 0)
    
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                 
                 # Filter by CNV size
                 & input$event_size > 1e5
                               
                 # & merged_file$event_size < 1e5
                               
                 # Filter by CNV p-value
                 & input$p_value <= 0.01, ]
    
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$chrom, 
                                      levels=chroms)
    
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input) + geom_point(aes(x = chr_mid / 10e6, y = seg_mean, 
                                                     color = event_type)) +
      scale_y_continuous(name = 'CNV LR mean' , 
                         labels = scales::comma) +  
      facet_grid(.~chrom_order, scales="free") +
      scale_x_continuous(name = 'Genomic mid-coordinate (Mb)', 
                         labels = scales::comma) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
      theme(legend.position = "none")
    
    return(p)
}


# This function reads a mergeSegments.pl output file (.events.tsv), removes 
# identified CNVs non classified as deletions/amplifications and plots them after 
# some filtering by size / p-value, presenting their coordinate along the x-axis 
# and their size along the y-axis, coloring them by their SV nature
  
CNV_dotplot_size <- function(merged_file){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
    
  # Add column with indicative middle coordinate
  input$chr_mid <- round((input$chr_start + 
                                    input$chr_stop) / 2, 0)
    
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                               
                               # Filter by CNV size
                               & input$event_size > 1e5
                               # & input$event_size < 1e5
                               
                               # Filter by CNV p-value
                               & input$p_value <= 0.01, ]
    
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$chrom, levels=chroms)
    
  # Apply log10 scale along y axis
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input) + geom_point(aes(x = chr_mid / 10e6, y = event_size, 
                                                     color = event_type)) +
    scale_y_continuous(name = 'CNV size', 
                         # trans = 'log10', 
                         labels = scales::comma) +
    facet_grid(.~chrom_order, scales="free") +
    scale_x_continuous(name = 'Genomic mid-coordinate (Mb)', 
                         labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
    theme(legend.position = "none")
    
  return(p)
  }

# The same as before, but as a barplot

CNV_barplot_size <- function(merged_file){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
  
  # Add column with indicative middle coordinate
  input$chr_mid <- round((input$chr_start + 
                            input$chr_stop) / 2, 0)
  
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                             
                             # Filter by CNV size
                             & input$event_size > 1e5
                             # & merged_file$event_size < 1e5
                             
                             # Filter by CNV p-value
                             & input$p_value <= 0.01, ]
  
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$chrom, levels=chroms)
  
  # Apply log10 scale along y axis
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input) + geom_col(aes(x = chr_mid / 10e6, y = event_size, 
                                                 color = event_type)) +
    scale_y_continuous(name = 'CNV size',
                       # trans = 'log10',
                       labels = scales::comma) +
    facet_grid(.~chrom_order, scales="free") +
    scale_x_continuous(name = 'Genomic mid-coordinate (Mb)', 
                       labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
    theme(legend.position = "none")
  
  return(p)
}

setwd("~/Juanma/CNV_analysis/")


# CBS runs
# Family 2: CONTROL VS 51609
CBS_detection(copynumber_in = "controlVS51609", 
              CBS_out = "CBS_fam2_controlVS51609_pvals")
CBS_detection(copynumber_in = "controlVS51610", 
              CBS_out = "CBS_fam2_controlVS51610_pvals")

# mergeSegments.pl runs
# Family 2: CONTROL VS 51609
seg_merger(CBS_in = "CBS_fam2_controlVS51609_pvals", 
           merged_out = "merged_controlVS51609_pvalue")

# Plots
# Family 2: CONTROL VS 51609
CNV_dotplot_LRmean(merged_file = "merged_controlVS51609_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS51609_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS51609_pvalue")
    

############ DIRT, dont look below, at your own risk ###########

# # Family 2, 51611 VS 51610
# fam2_11vs10 <- read.table ("11vs10.copynumber", header=TRUE)
# fam2_11vs10_pvals <- CBS_detection(case = fam2_11vs9)
# write.table (fam2_11vs10_pvals, file="CBS_fam2_11vs10_pvals.txt", 
#              row.names=F, col.names=F, quote=F, sep="\t")

## Family 2 ##
# N11vsA9 <- read.table("merged_11vs9_pvalue.events.tsv", header=TRUE)
# N11vsA10 <- read.table("merged_11vs10_pvalue.events.tsv", header=TRUE)

# This function takes a .copynumber R file, builds a CNA object and 
# computes and returns the CBS after removing outliers to compute the 
# p-vals for the change points found

# CBS_detection <- function(CN_file){
#   
#   # Copy Number Array object (matrix with log ratio changes)
#   CNA.object <- CNA(genomdat = CN_file[,7], chrom = CN_file[,1],
#                     maploc = CN_file[,2], data.type = 'logratio')
#   
#   # Detect outliers and smooth the data prior to CBS
#   CNA.object.smoothed <- smooth.CNA(CNA.object)
#   
#   # CBS
#   CNA.object.smoothed.seg <- segment(CNA.object.smoothed, verbose = 0)
#   
#   # Plot with data ordered by chromosome and map positions.
#   
#   # Red dots correspond to the segment (CNV) mean LR change values of the case when 
#   # compared to the control. I.e., if the CNV corresponds to a deletion in the affected
#   # compared to the normal case, it will be marked with a negative mean LR change value 
#   # and it will be shown below the grey line (which will be plotted as adot/line
#   # later on). On the contrary, if the CNV corresponds to an amplification, it will be
#   # marked with a positive mean LR change value and it will be shown above the grey line 
#   # (SALMON dot)
#   
#   # print(plot(CNA.object.smoothed.seg, plot.type="w"))
#   
#   # Compute p-vals and CIs for the change points found by the CBS. The data
#   # are permuted to the left and right of the identified change point and the
#   # location of the maximal binary segmentation is computed. The CI is given
#   # by the quantiles of the permutation distribution of the locations
#   seg.pvalue <- segments.p(CNA.object.smoothed.seg, 
#                            ngrid = 100, 
#                            tol = 1e-6, 
#                            # A bit more restrictive than default
#                            alpha = 0.01, 
#                            search.range = 100, 
#                            nperm = 1000)
#   
#   # 2 additional fixes for posterior mergeSegment.pl script
#   
#   # Add a second column with "chrom"
#   seg.pvalue <- add_column(seg.pvalue, fixed = "chrom", .after = 1)
#   
#   # Keep only true chromosomes, not seq-generated artifact denominations 
#   seg.pvalue.clean <- seg.pvalue[seg.pvalue$chrom %in% chroms, ]
#   
#   return(seg.pvalue.clean)
# }

    