#!/usr/bin/env bash
# This is a script to pass files to and start up the docker to run guppy_basecaller and other scripts

# First the scripts:
GENOME="/data16/marcus/genomes/elegansRelease100"
DATA="/data16/marcus/minknow/data/210525_PolyA_WT_L3/210525_PolyA_WT_L3/20210526_0115_MN36576_FAP47925_42184fdb"
OUTPUT="/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir"
DATETIMESTR=$(date +"%y%m%d_%H:%M")

sudo docker run -v "${GENOME}":/usr/src/working/genome_dir \
                -v "${DATA}":/usr/src/working/data_dir \
                -v "${OUTPUT}":/usr/src/working/output_dir \
                -it nanopore_run 2>&1 | tee $OUTPUT/"$DATETIMESTR".log
