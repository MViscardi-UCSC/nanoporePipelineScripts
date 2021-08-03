#!/usr/bin/env bash
# This is a script to pass files to featureCounts for the Roach et al. dataset

# First the dirs:
GENOME="/data16/marcus/genomes/elegansRelease100"
OUTPUT="/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir"
DTSTR=$(date +"%y%m%d_%H:%M")
READABLE_DT=$(date +"%D %H:%M")
LOG=$OUTPUT/"$DTSTR"_featureCounts.log
touch $LOG

mkdir $OUTPUT/featureCounts

echo ; echo ; echo "featureCounts starting at: " $READABLE_DT 2>&1 | tee -a $LOG ; echo
for bam_folder in $OUTPUT/bam_files/* ; do
    name="$(basename $bam_folder)"
    bam=$OUTPUT/bam_files/$name/$name.bam
    echo ; echo "featureCounts for "$name" at: " $READABLE_DT 2>&1 | tee -a $LOG
    mkdir $OUTPUT/featureCounts/$name
    featureCounts -L -R CORE -a $GENOME/Caenorhabditis_elegans.WBcel235.100.gtf \
    -o $OUTPUT/featureCounts/$name/$DTSTR $bam 2>&1 | tee -a $LOG
done
echo ; echo ; echo "featureCounts finished at: " $READABLE_DT 2>&1 | tee -a $LOG ; echo
