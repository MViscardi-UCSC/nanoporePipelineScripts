#!/usr/bin/env bash

echo ; echo "guppy_basecaller and nanopolish index starting at: $(date)" ; echo
OUTDIR="/usr/src/working/output_dir"
mkdir $OUTDIR/fastqs
mkdir $OUTDIR/cat_fastqs
mkdir $OUTDIR/seq_summaries
echo; echo ; echo "directories build and ready for basecalling" ; echo
guppy_basecaller --num_callers 3 --cpu_threads_per_caller 10 -c rna_r9.4.1_70bps_hac.cfg -i /usr/src/working/data_dir/fast5 -s $OUTDIR/fastqs
echo ; echo "basecalling completed" ; echo
cat $OUTDIR/fastqs/*.fastq > $OUTDIR/cat_fastqs/cat.fastq
echo "done"
echo ; echo "guppy_basecaller and nanopolish index done at: $(date)" ; echo
