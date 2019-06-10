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
library("ggrepel")  # To automatically separate the labels in scatterplots

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
                           # A bit more restrictive than default
                           # alpha = 0.05, 
                           alpha = 0.01,
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


# This function points a merged file path to Bedtools, calls it to perform the
# intersection with a .bed file containing every unique gene from hg19 assembly
# and its coordinate (sorted by left-most coordinate), and writes the output 
# in a pointed path

overlapper <- function(merged_file, overlapped_out){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Path to the output overlapping .bed file
  case_out <- paste("~/Juanma/CNV_analysis/", 
                    overlapped_out, 
                    ".bed", 
                    sep = "")
  
  # bash:
  # bedtools intersect -wo 
  
  ## Minimum overlap required as a fraction of B
  ## i.e. 10% of the gene must be overlapping the CNV to be included
  # -F 0.1
  
  # -a merged_controlVS20592_pvalue.events.tsv 
  # -b ../ucsc_genes/unique_genes_hg19_ucsc.bed > test_20592_inter_genes.bed

  # Bedtools call
  system(paste("/home/genomica/bin/bedtools2/bin/bedtools", 
               "intersect -wo -F 0.1",
               "-a", 
               case_file,
               "-b",
               "~/Juanma/ucsc_genes/unique_sorted_genes_hg19_ucsc.bed",
               ">",
               case_out,
               sep=" "))
  
}


# This function reads a mergeSegments.pl output file (.events.tsv), removes 
# identified CNVs non classified as deletions/amplifications and plots them after 
# some filtering by size / p-value / LR change, presenting their coordinate along 
# the x-axis and their CNV LR mean along the y-axis, coloring them by their SV 
# nature

CNV_dotplot_LRmean <- function(merged_file){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
  
  # Add column with indicative CNV middle coordinate
  input$chr_mid <- round((input$chr_start + 
                                    input$chr_stop) / 2, 0)
    
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                 
                 # Remove sex chromosomes
                 & !(input$chrom == "chrX")
                 & !(input$chrom == "chrY")
                 
                 # Filter by -1 <= mean LR change >= 1
                 # NB! Only keep segments with a 2-fold increase in coverage in case
                 # compared to control
                 & (input$seg_mean <= -1 | input$seg_mean >= 1)
                 
                 # Filter by CNV size
                 & input$event_size > 1e5
                               
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
# some filtering by size / p-value / LR change, presenting their coordinate along 
# the x-axis and their size along the y-axis, coloring them by their SV nature
  
CNV_dotplot_size <- function(merged_file){
  
  # Path to the input merged file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     merged_file, 
                     ".events.tsv", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=TRUE)
    
  # Add column with indicative CNV middle coordinate
  input$chr_mid <- round((input$chr_start + 
                                    input$chr_stop) / 2, 0)
    
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                               
                 # Remove sex chromosomes
                 & !(input$chrom == "chrX")
                  & !(input$chrom == "chrY")
                 
                 # Filter by -1 <= mean LR change >= 1 
                 & (input$seg_mean <= -0.5 | input$seg_mean >= 0.5)
                 
                 # Filter by CNV size
                 & input$event_size > 1e5
 
                 # Filter by CNV p-value
                 & input$p_value <= 0.01, ]
    
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$chrom, levels=chroms)
  
  # To print out the number of reported CNVs
  print(nrow(input))
    
  # Apply log10 scale along y axis
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input) + geom_point(aes(x = chr_mid / 10e6, y = event_size, 
                                                     color = event_type)) +
    scale_y_continuous(name = 'CNV size (log10)', 
                         trans = 'log10',
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
  
  # Add column with indicative CNV middle coordinate
  input$chr_mid <- round((input$chr_start + 
                            input$chr_stop) / 2, 0)
  
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$event_type == "neutral")
                 
                 # Remove sex chromosomes
                 & !(input$chrom == "chrX")
                 & !(input$chrom == "chrY")
                 
                 # Filter by -1 <= mean LR change >= 1 
                 & (input$seg_mean <= -1 | input$seg_mean >= 1)
                             
                 # Filter by CNV size
                 & input$event_size > 1e5
                            
                 # Filter by CNV p-value
                 & input$p_value <= 0.01, ]
  
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$chrom, levels=chroms)
  
  # Apply log10 scale along y axis
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input) + geom_col(aes(x = chr_mid / 10e6, y = event_size, 
                                                 color = event_type, fill = event_type)) +
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

# This function reads a .bed output file resulting from the intersection between
# the predicted CNVs and the UCSC unique sorted genes. It further removes CNVs
# non classified as deletions/amplifications and plots them after some filtering
# by size / p-value / LR change, presenting their coordinate along the x-axis 
# and their size along the y-axis, coloring them by their SV nature

CNV_I_genes_dotplot <- function(input, bed_out){
  
  # Path to the .bed overlapping file
  case_file <- paste("~/Juanma/CNV_analysis/", 
                     input, 
                     ".bed", 
                     sep = "")
  
  # Read in file
  input <- read.table(file = case_file, header=FALSE)

  # cols 1:13  --> CNV info
  # cols 14:17 --> gene info
  # col 18     --> overlapping info
  
  colnames(input) <- c("CNV_chrom", "CNV_start", "CNV_stop", "CNV_seg_mean", 
                       "CNV_num_segments", "CNV_num_markers", "CNV_pval", 
                       "CNV_event_type", "CNV_event_size", "CNV_class", "CNV_arm", 
                       "CNV_arm_fraction", "CNV_chrom_fraction", "gene_chrom", 
                       "gene_start", "gene_stop", "gene_name", "overlapping_fract")
  
  # Add column with indicative CNV middle coordinate
  input$CNV_chr_mid <- round((input$CNV_start + input$CNV_stop) / 2, 0)
  
  # Remove neutral events (non deletions/amplifications) 
  input <- input[!(input$CNV_event_type == "neutral")
                 
                 # Remove sex chromosomes
                 & !(input$CNV_chrom == "chrX")
                 & !(input$CNV_chrom == "chrY")
                 
                 # Filter by -1 <= mean LR change >= 1 
                 & (input$CNV_seg_mean <= -1 | input$CNV_seg_mean >= 1)
                 
                 # Filter by CNV size
                 & input$CNV_event_size > 1e5
                 
                 # Filter by CNV p-value
                 & input$CNV_pval <= 0.01, ]
  
  # To order the chromosomic facets in ascending order
  input$chrom_order <- factor(input$CNV_chrom, levels=chroms)
  
  print(nrow(input))
  
  # Path to the output .txt file
  filename_out <- paste("~/Juanma/CNV_analysis/", bed_out, 
                        "_CNVs_genes_filtered.bed", sep = "")
  
  # Write out file (CNV coords, type CNV, gene coords)
  write.table(input[, c(1,2,3,8,9,15,16,17)], filename_out,
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  # Apply log10 scale along y axis
  # Cyan: deletions -- salmon: amplifications
  p <- ggplot(data = input, aes(x = CNV_chr_mid / 10e6, y = CNV_event_size, 
                                color = CNV_event_type, label = gene_name)) + 
    geom_point() +
    
    geom_text_repel(aes(label = gene_name), 
                    box.padding = unit(0.25, "lines")) +
    
    scale_y_continuous(name = 'CNVs with >= 1 overlapping gene', 
                       trans = 'log10',
                       labels = scales::comma) +
    facet_grid(.~chrom_order, scales="free") +
    scale_x_continuous(name = 'Genomic mid-coordinate (Mb)', 
                       labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
    theme(legend.position = "none")
  
  return(p)
}


############ ANALYSIS ############

#### To polish-up: ####
  # mean LR change (-0.5 <= mean LR change >= 0.5)
# Segment size parameter of Varscan (longer segments better? 100-10000)
# Number of markers
# Number of windows 

setwd("~/Juanma/CNV_analysis/")
    



### NEW CBS runs, with segment size 100-10000
## In average, 0.5 kk segments found in each

# Family 2: 51609
CBS_detection(copynumber_in = "51609_NEW", 
              CBS_out = "CBS_NEW_fam2_51609_pvals")
  
# Family 2: 51610
CBS_detection(copynumber_in = "51610_NEW", 
              CBS_out = "CBS_NEW_fam2_51610_pvals")

# Family 3: 20592 (DONE)
CBS_detection(copynumber_in = "20592_NEW", 
              CBS_out = "CBS_NEW_fam3_20592_pvals")

# Family 3: 46632 (DONE)
CBS_detection(copynumber_in = "46632_NEW", 
              CBS_out = "CBS_NEW_fam3_46632_pvals")

# Family 4: 52718 (DONE)
CBS_detection(copynumber_in = "52718_NEW", 
              CBS_out = "CBS_NEW_fam4_52718_pvals")

# Family 5: 22010 (DONE) 
CBS_detection(copynumber_in = "22010_NEW", 
              CBS_out = "CBS_NEW_fam5_22010_pvals")

# Family 5: 23281 (DONE)
CBS_detection(copynumber_in = "23281_NEW", 
              CBS_out = "CBS_NEW_fam5_23281_pvals")


### NEW mergesegments runs, with segment size 100-10000

seg_merger(CBS_in = "CBS_NEW_fam2_51609_pvals", 
           merged_out = "merged_NEW_fam2_51609_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam2_51610_pvals", 
           merged_out = "merged_NEW_fam2_51610_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam3_20592_pvals", 
           merged_out = "merged_NEW_fam3_20592_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam3_46632_pvals", 
           merged_out = "merged_NEW_fam3_46632_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam4_52718_pvals", 
           merged_out = "merged_NEW_fam4_52718_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam5_22010_pvals", 
           merged_out = "merged_NEW_fam5_22010_pvalue")

seg_merger(CBS_in = "CBS_NEW_fam5_23281_pvals", 
           merged_out = "merged_NEW_fam5_23281_pvalue")






### bedtools runs ### 
# Family 2: 51609 (DONE)
overlapper(merged_file = "merged_controlVS51609_pvalue", 
           overlapped_out = "51609_inter_genes")

# test 51609 (segment size 100-10000 instead of 10-100): 
# 0.5 kk segments found instead of 3 kk
# 18 amplifications // 57 deletions
overlapper(merged_file = "merged_NEW_fam2_51609_pvalue", 
           overlapped_out = "51609_NEW_inter_genes")

# Family 2: 51610 (DONE)
overlapper(merged_file = "merged_controlVS51610_pvalue", 
           overlapped_out = "51610_inter_genes")

# test 51610 (segment size 100-10000 instead of 10-100): 
# 0.5 kk segments found instead of 3.5 kk
# 23 amplifications // 70 deletions
overlapper(merged_file = "merged_NEW_fam2_51610_pvalue", 
           overlapped_out = "51610_NEW_inter_genes")

# Family 3: 20592 (DONE)
overlapper(merged_file = "merged_controlVS20592_pvalue", 
           overlapped_out = "20592_inter_genes")

# Family 3: 46632 (DONE)
overlapper(merged_file = "merged_controlVS46632_pvalue", 
           overlapped_out = "46632_inter_genes")

# Family 4: 52718 (DONE)
overlapper(merged_file = "merged_controlVS52718_pvalue", 
           overlapped_out = "52718_inter_genes")

# Family 5: 22010 (DONE)
overlapper(merged_file = "merged_controlVS22010_pvalue", 
           overlapped_out = "22010_inter_genes")

# Family 5: 23281 (DONE)
overlapper(merged_file = "merged_controlVS23281_pvalue", 
           overlapped_out = "23281_inter_genes")

### Plots ### 
# Family 2: 51609
CNV_dotplot_LRmean(merged_file = "merged_controlVS51609_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS51609_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS51609_pvalue")
CNV_I_genes_dotplot(input = "51609_inter_genes", bed_out = "51609")

# test 51609 (segment size 100-10000 instead of 10-100): 
# 0.5 kk segments found instead of 3 kk
# 18 amplifications // 57 deletions
CNV_dotplot_LRmean(merged_file = "merged_NEW_fam2_51609_pvalue")
CNV_dotplot_size(merged_file = "merged_NEW_fam2_51609_pvalue")
CNV_I_genes_dotplot(input = "51609_NEW_inter_genes", bed_out = "51609_NEW")

# Family 2: 51610
CNV_dotplot_LRmean(merged_file = "merged_controlVS51610_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS51610_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS51610_pvalue")
CNV_I_genes_dotplot(input = "51610_inter_genes", bed_out = "51610")

# test 51610 (segment size 100-10000 instead of 10-100): 
# 0.5 kk segments found instead of 3.5 kk
# 23 amplifications // 70 deletions
CNV_dotplot_LRmean(merged_file = "merged_NEW_fam2_51610_pvalue")
CNV_dotplot_size(merged_file = "merged_NEW_fam2_51610_pvalue")
CNV_I_genes_dotplot(input = "51610_NEW_inter_genes", bed_out = "51610_NEW")

# Family 3: 20592
CNV_dotplot_LRmean(merged_file = "merged_controlVS20592_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS20592_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS20592_pvalue")
CNV_I_genes_dotplot(input = "20592_inter_genes", bed_out = "20592")

# Family 3: 46632
CNV_dotplot_LRmean(merged_file = "merged_controlVS46632_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS46632_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS46632_pvalue")
CNV_I_genes_dotplot(input = "46632_inter_genes", bed_out = "46632")

# Family 4: 52718
CNV_dotplot_LRmean(merged_file = "merged_controlVS52718_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS52718_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS52718_pvalue")
CNV_I_genes_dotplot(input = "52718_inter_genes", bed_out = "52718")

# Family 5: 22010
CNV_dotplot_LRmean(merged_file = "merged_controlVS22010_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS22010_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS22010_pvalue")
CNV_I_genes_dotplot(input = "22010_inter_genes", bed_out = "22010")

# Family 5: 23281
CNV_dotplot_LRmean(merged_file = "merged_controlVS23281_pvalue")
CNV_dotplot_size(merged_file = "merged_controlVS23281_pvalue")
CNV_barplot_size(merged_file = "merged_controlVS23281_pvalue")
CNV_I_genes_dotplot(input = "23281_inter_genes", bed_out = "23281")




### OLD CBS runs, with segment size 10-100 ###
## In average, 3kk segments found in each

# # Family 2: 51609 (DONE)
# CBS_detection(copynumber_in = "controlVS51609", 
#               CBS_out = "CBS_fam2_controlVS51609_pvals")
# 
# # Family 2: 51610 (DONE)
# CBS_detection(copynumber_in = "controlVS51610", 
#               CBS_out = "CBS_fam2_controlVS51610_pvals")

# ...

# ### OLD mergeSegments.pl runs, with segment size 10-100 ###
# # Family 2: 51609 (DONE):
# # 7725 amplifications // 2905 deletions
# seg_merger(CBS_in = "CBS_fam2_controlVS51609_pvals", 
#            merged_out = "merged_controlVS51609_pvalue")
# 
# # Family 2: 51610 (DONE)
# # 8000 amplifications // 3600 deletions
# seg_merger(CBS_in = "CBS_fam2_controlVS51610_pvals", 
#            merged_out = "merged_controlVS51610_pvalue")

# ...

