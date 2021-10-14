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
import pprint
import subprocess
import pandas as pd
from Bio import SeqIO

pd.set_option("display.max_columns", None)


def cutadapt_call(adaptor_file, output_file, input_file,
                  error_fraction=0.3, cores=20, overlap=16, allow_indels=False, head=None):
    tsv_name = output_file.rstrip(".fa") + ".tsv"
    call = f"cutadapt -a file:{adaptor_file} -o {output_file} --info-file={tsv_name} " \
           f"-j {cores} -e {error_fraction} --overlap {overlap} --action lowercase " \
           "--rename '{header} adapter={adapter_name}' "
    if allow_indels:
        call += "--no-indels "
    if isinstance(head, int):
        buffer = subprocess.check_output(f"head -n{head * 4} {input_file}", shell=True)
        # print(buffer.decode('utf-8'))
        subprocess.run(call+f"-",
                       shell=True,
                       input=buffer)
    else:
        subprocess.run(call+f"{input_file}",
                       shell=True)


def fa_to_df(fa_file) -> pd.DataFrame:
    with open(fa_file) as fa_file:  # Will close handle cleanly
        list_of_read_dicts = []
        for record in SeqIO.QualityIO.FastqGeneralIterator(fa_file):
            list_of_read_dicts.append(dict(zip(["header", "seq", "qual"], record)))
    fastq_df = pd.DataFrame(list_of_read_dicts)
    header_columns = ["read_id",
                      "run_id",
                      "sample_id",
                      "read_num",
                      "channel",
                      "start_time",
                      "adapter"]
    fastq_df[header_columns] = fastq_df["header"].str.split(expand=True)
    for column, prefix in zip(header_columns[1:], ["runid", "sampleid", "read", "ch", "start_time", "adapter"]):
        fastq_df[column] = fastq_df[column].str.lstrip(f"{prefix}=")
    fastq_df.drop(columns="header", inplace=True)
    fastq_df = fastq_df[["read_id",
                         "run_id",
                         "sample_id",
                         "read_num",
                         "channel",
                         "start_time",
                         "adapter",
                         "seq",
                         "qual",
                         ]]
    return fastq_df


if __name__ == '__main__':
    path_to_fa = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
                   "output_dir/cat_files/cat.fastq"
    subset_size = 1000
    error_dict = {}
    heatmap_df = pd.DataFrame([])
    for error in range(30, 74, 2):
        error /= 100
        for overlap in range(2, 18, 2):
            cutadapt_call(adaptor_file="./standardUMIs.fa",
                          output_file="./deleteme.fa",
                          input_file=path_to_fa,
                          head=subset_size,
                          error_fraction=error,
                          overlap=overlap)
            df = fa_to_df("deleteme.fa")
            # print("\n", error)
            # print(df["adapter"].value_counts())
            test_dict = df["adapter"].value_counts().to_dict()
            test_dict["error"] = error
            test_dict["overlap"] = overlap
            heatmap_df = heatmap_df.append(test_dict, ignore_index=True)
            error_dict[f"e:{error}, o:{overlap}"] = df["adapter"].value_counts().to_dict()
    pprint.pprint(error_dict)
    print()
    import seaborn as sea
    import matplotlib.pyplot as plt
    
    # Great heatmap stackoverflow:
    #   https://stackoverflow.com/questions/28356359/one-colorbar-for-seaborn-heatmaps-in-subplot
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    cbar_ax = fig.add_axes([.91, .1, .03, .8])
    fig.subplots_adjust(top=0.9)
    fig.suptitle(f"Adjusting allowed error and overlap values for cutadapt (subset={subset_size})")
    heatmap_df.fillna(0, inplace=True)
    for i, (ax, value) in enumerate(zip(axs.flat, ["05", "10", "15", "30"])):
        column = f"{value}A_adapter"
        heatmap_df[f"{column}_norm"] = heatmap_df[column].div(subset_size)
        heatmap_df_pivot = heatmap_df.pivot(index="error", columns="overlap", values=f"{column}_norm")
        ax = sea.heatmap(heatmap_df_pivot, ax=ax,
                         cbar=i == 0,
                         cmap="rocket_r",
                         vmin=0.001,
                         vmax=0.5,
                         cbar_ax=None if i else cbar_ax)
        ax.set_title(f"Heatmap of {column} calls")
        ax.set(xlabel='minimum required overlap', ylabel='maximum fractional error')
    plt.tight_layout(rect=[0, 0, .9, 1])
    plt.show()
    print("Done!")
