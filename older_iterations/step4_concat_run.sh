#!/usr/bin/env bash
# This is a script to pass files to featureCounts for the Roach et al. dataset

# First the dirs:
OUTPUT="/data16/marcus/working/210528_NanoporeRun_0639_L3s"
DTSTR=$(date +"%y%m%d_%H:%M")
READABLE_DT=$(date +"%D %H:%M")

# First lets merge a huge bam file, index it, and make a sam from it:
samtools merge -f $OUTPUT/cat_allReads.bam $OUTPUT/bam_files/*/*.sorted.bam
samtools index $OUTPUT/cat_allReads.bam
samtools view $OUTPUT/cat_allReads.bam > $OUTPUT/cat_allReads.sam

# Merge all the polya calls:
head -n 1 $OUTPUT/polya_files/000*.tsv > $OUTPUT/cat_polya.tsv
tail -n +1 $OUTPUT/polya_files/*.tsv >> $OUTPUT/cat_polya.tsv
head -n 1 $OUTPUT/cat_polya.tsv > $OUTPUT/cat_polya_passed.tsv
grep PASS $OUTPUT/cat_polya.tsv >> $OUTPUT/cat_polya_passed.tsv

# Merge featureCounts
cat $OUTPUT/featureCounts/*/*.bam.featureCounts >> $OUTPUT/cat_featureCounts.tsv
grep Assigned $OUTPUT/cat_featureCounts.tsv >> $OUTPUT/cat_featureCounts_assigned.tsv
