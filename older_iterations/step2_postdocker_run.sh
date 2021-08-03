#!/usr/bin/env bash
# This is a script to pass files to and run Minimap2, Samtools and Nanopolish

# First the dirs:
GENOME="/data16/marcus/genomes/elegansRelease100"
DATA="/data16/marcus/minknow/data/210525_PolyA_WT_L3/210525_PolyA_WT_L3/20210526_0115_MN36576_FAP47925_42184fdb/fast5"
OUTPUT="/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir"
DATESTR=$(date +"%y%m%d")

mkdir $OUTPUT/bam_files
mkdir $OUTPUT/polya_files

echo ; echo ; echo "mapping and polya calling starting at: $(date +"%D %H:%M ")" ; echo
for fastq in $OUTPUT/cat_fastqs/*.fastq ; do
    name="$(basename $fastq)"
    name="${name:0:3}"
    echo ; echo "Nanopolish, Minimap2, & Samtools for "$name" at: $(date +"%D %H:%M")"
    echo $DATA/fast5/$name
    nanopolish index --directory=$DATA/fast5/$name $fastq
    mkdir $OUTPUT/bam_files/$name
    minimap2 -a -x map-ont $GENOME/200430_allChrs.fa $fastq | samtools view -b - -o $OUTPUT/bam_files/$name/$name.bam
    samtools sort -T tmp -o $OUTPUT/bam_files/$name/$name.sorted.bam $OUTPUT/bam_files/$name/$name.bam \
    && samtools index $OUTPUT/bam_files/$name/$name.sorted.bam
    nanopolish polya --threads=36 --reads=$OUTPUT/cat_fastqs/$name.fastq --bam=$OUTPUT/bam_files/$name/$name.sorted.bam \
    --genome=$GENOME/200430_allChrs.fa > $OUTPUT/polya_files/$name.polya.tsv
done
echo ; echo ; echo "mapping and polya calling finished at: $(date +"%D %H:%M")" ; echo
