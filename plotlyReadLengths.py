"""
plotlyReadLengths.py
Marcus Viscardi,    October 19, 2021

Plotting read lengths, specifically genes from MtDNA
"""
import pandas as pd
pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

import plotly.express as px
from nanoporePipelineCommon import find_newest_matching_file

def load_df():
    pass


if __name__ == '__main__':
    pathdict = {
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/"
                 "merge_files/*_mergedOnReads.tsv",
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                    "output_dir/merge_files/*_mergedOnReads.tsv",
        "totalRNA2": "/data16/marcus/working/"
                     "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir/"
                     "merge_files/*_mergedOnReads.tsv",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/"
                 "merge_files/*_mergedOnReads.tsv",
        "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/"
                  "output_dir/merge_files/*_mergedOnReads.tsv",
        "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
                 "output_dir/merge_files/*_mergedOnReads.tsv"
    }
    
    run_with = ["polyA2"]
    
    pathdict = {name: find_newest_matching_file(pathdict[name]) for name in run_with}
    
    for name, lib in pathdict.items():
        df = pd.read_csv(lib, sep="\t")
        df["read_len"] = df["sequence"].str.len()
        print(df)
        mt_genes = df[df["chr_id"] == "MtDNA"]
        mt_genes["gene_name"] = mt_genes["gene_name"].fillna("Unmapped")
        fig = px.histogram(mt_genes, x="read_len",
                           color="gene_name",
                           )
        fig.show()
        
        fig = px.box(mt_genes, y="read_len", x="gene_name")
        fig.show()
        print("Done..")
