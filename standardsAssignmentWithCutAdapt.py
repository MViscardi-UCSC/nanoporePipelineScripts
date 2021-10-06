"""
standardsAssignmentWithCutAdapt.py
Marcus Viscardi,    October 05, 2021

Goal here is to try to use cutadapt to assign
    standard reads based on their 3' UMIs.

This script will mostly just be a wrapper for the
    cutadapt call. I will eventually also add the
    steps to process the cutadapt info file and
    combine that information with the bam files.

Going to use multiple named adapters as noted here:
    https://cutadapt.readthedocs.io/en/stable/guide.html#named-adapters

Going to try and pull the assigned adapters with this:
    https://cutadapt.readthedocs.io/en/stable/guide.html#info-file
"""
import subprocess
import pandas as pd


def cutadapt_call(adaptor_file, output_file, input_file):
    subprocess.run(f"cutadapt "
                   f"-a file:{adaptor_file} "
                   f"-o {output_file} "
                   f"--info-file=deleteme.tsv "
                   f"-j 20 -e 0.3 "
                   f"--no-indels "
                   # f"--overlap 16 "
                   f"--action lowercase "
                   "--rename '{header} adapter={adapter_name}' "
                   f"{input_file}",
                   shell=True)


if __name__ == '__main__':
    path_to_fa = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
                   "output_dir/cat_files/cat.fastq"
    cutadapt_call(adaptor_file="/data16/marcus/scripts/nanoporePipelineScripts/standardUMIs.fa",
                  output_file="./deleteme.fa",
                  input_file=path_to_fa)
