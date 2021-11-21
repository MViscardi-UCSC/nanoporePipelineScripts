"""
samParsingFor5TERA.py
Marcus Viscardi,    November 19, 2021

Goal here is to test out parsing sam files and adding a new TAG that I can use
    for looking at 5TERA reads in IGV
"""
import simplesam as ssam
from step0_nanopore_pipeline import live_cmd_call
from nanoporePipelineCommon import FastqFile


def convert_read_id_tagged_sam(in_sam_path="/data16/marcus/working/"
                                           "211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/"
                                           "output_dir/cat_files/deleteme.sam"):
    # The following loop will take stuff from a SAM file that was made by running minimap2 on
    #   a fastq that was run through cutadapt with the parameter: --rename '{id}:{adapter_name} {comment}'
    # This is a fine solution but is going to make an absolute MESS when it comes to the other things that
    #   read the fastq file (nanopolish and featureCounts, I think?)
    # For this reason, I'll need to very selectively keep or drop this other information as we go!
    #   Or maybe just make a copy of the fastq with the only change being the inclusion of the cutadapt
    #   tag to the name. Then I just feed the original to everything besides minimap2? then immediately
    #   after minimap2 I do the following to remove the cutadapt tag and make it a bam/sam tag?!
    with ssam.Reader(open(in_sam_path, 'r')) as in_sam:
        with ssam.Writer(open(in_sam_path.strip(".sam") + ".customComment.sam", 'w'),
                         in_sam.header) as out_sam:
            for read in in_sam:
                print(read)
                read["ZA"] = read.qname.split(":")[1]
                read.qname = read.qname.split(":")[0]
                print(read.tags["ZA"])
                out_sam.write(read)


if __name__ == '__main__':
    threads = 20
    output_dir_path = "/data16/marcus/working/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/output_dir"
    fastq_file_path = f"{output_dir_path}/cat_files/cat.fastq"
    output_fastq_path = f'{output_dir_path}/cat_files/cat.teraAdapters.fastq'
    # tera5_adapter = "AAUGAUACGGCGACCACCGAGAUCUACACUCUUUCCCUACACGACGCUCUUCCGANNN"  # Add X to front
    # cutadapt_call = f"cutadapt -g 5TERA={'X' + tera5_adapter} --action=trim -j {threads} " \
    #                f"--overlap 31 --error-rate 0.29 " \
    #                "--rename '{id}:{adapter_name} {comment}' " \
    #                f"{fastq_file_path} > {output_fastq_path}"
    # live_cmd_call(cutadapt_call)
    adapter_modified_fastq = FastqFile(output_fastq_path)
    