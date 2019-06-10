
#!/bin/bash

# Define paths
REF_GENOME=/media/genomica/ref/hg19/genome.fa
PICARD_DIR=/media/genomica/software/SeqMule/exe/picard


#CreateSequenceDictionary: needed to reorderSam --> OK
#java -jar picard.jar CreateSequenceDictionary R=reference.fasta O=reference.dict


for FILE in *[0-9].bam; do

	echo 'Starting GATK massage for file:' $FILE

	echo 'Starting SAMTOOLS procedure for file:' $FILE

	# Sort alignment by left-most coordinate
	echo 'Sorting file ... ' $FILE
	samtools sort -@ 4 -o "$(basename "$FILE" .bam)_sorted.bam" $FILE
	echo 'Sorting done'

	# Generate .bam.bai index of sorted .bam for fast random access
	echo 'Generating index ... ' $FILE
	samtools index "$(basename "$FILE" .bam)_sorted.bam"
	echo 'Index built'

	echo 'SAMTOOLS done for file:' $FILE

	# FASTER PICARD PROCEDURE WITH INDEX GENERATED
	echo 'Starting PICARD procedure for file:' $FILE

	# Identifies duplicated reads from e.g. library construction using PCR
	echo 'PICARD MarkDuplicates ... '
	java -Xmx10g -jar $PICARD_DIR/MarkDuplicates.jar I="$(basename "$FILE" .bam)_sorted.bam" O="$(basename "$FILE" .bam)_markedups.bam" M=$FILE.metrics.txt
	echo 'Duplicates marked'

	# Replace all RG in .bam with a single new read group and assign all reads to it in the output .bam
	echo 'PICARD AddOrReplaceGroups ... ' $FILE
        java -Xmx10g -jar $PICARD_DIR/AddOrReplaceReadGroups.jar I="$(basename "$FILE" .bam)_markedups.bam" O="$(basename "$FILE" .bam)_RG.bam" RGLB=dna RGPL=illumina RGPU=unknown RGSM="$(basename "$FILE" .bam)sample"

	# RGLB: DNA/RNA
	# RGPL: Platform
	# RGPU: Platform
	# RGSM: Library name

	echo 'Read groups added'

	# Re-orders reads to match the contig ordering in the reference file, as determined by exact name matching of contigs
	echo 'Picard ReorderSam ... ' $FILE
	java -Xmx10g -jar $PICARD_DIR/ReorderSam.jar I="$(basename "$FILE" .bam)_RG.bam" O="$(basename "$FILE" .bam)_reordered.bam" REFERENCE=$REF_GENOME
	echo 'Reordering completed'

	rm "$(basename "$FILE" .bam)_markedups.bam"
	rm "$(basename "$FILE" .bam)_RG.bam"
	rm "$(basename "$FILE" .bam)_sorted.bam.bai"

	echo 'PICARD done for file:' $FILE

	# Generate a .bam .bai index of reordered .bam
	echo 'Generating index ... ' $FILE
	samtools index "$(basename "$FILE" .bam)_reordered.bam"
	echo 'Index of .reordered built'
 

	echo 'GATK MASSAGE done for file:' $FILE
done

