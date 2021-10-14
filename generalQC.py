"""
generalQC.py
Marcus Viscardi,    August 30, 2021

Trying to write a simple QC script to spit out metrics from
    my nanopore runs. Nothing crazy for now

Initial goal is to spit out a real estimate mean of read lengths,
    as this is a metric that MinKnow is trying to estimate during
    Nanopore runs.
"""
import pandas as pd
from step0_nanopore_pipeline import find_newest_matching_file


def pick_file(lib_key: str, file_suffix: str) -> str:
    working_dir_dict = {
        "polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
        "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
        "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
        "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
        "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
        "xrn-1": "210905_nanoporeRun_totalRNA_5108_xrn-1-KD",  # First XRN-1 Knockdown
    }
    if lib_key in working_dir_dict.keys():
        working_dir = f"/data16/marcus/working/{working_dir_dict[lib_key]}"
        path_to_file = find_newest_matching_file(f"{working_dir}"
                                                 f"/output_dir"
                                                 f"/merge_files"
                                                 f"/*{file_suffix}")
        print(f"Found matching file @  {path_to_file}")
        return path_to_file
    else:
        raise NotImplementedError(f"Libary key \"{lib_key}\" not found!\n"
                                  f"Please use one of these: "
                                  f"{list(working_dir_dict.keys())}")


def load_mergedOnReads(lib_key: str) -> pd.DataFrame:
    suffix = "_mergedOnReads.tsv"
    path = pick_file(lib_key, suffix)
    print(f"Loading dataframe from {path}")
    df = pd.read_csv(path, sep="\t")
    print(f"Done!\n")
    return df


def get_read_lengths(lib_key: str) -> pd.Series:
    df = load_mergedOnReads(lib_key)
    read_lengths = df.sequence.apply(len)
    # read_lengths = read_lengths.add(df.polya_length)
    return read_lengths


def plotly_read_lengths(lib_keys: [str, str]):
    import plotly.graph_objects as go
    series1 = get_read_lengths(lib_keys[0])
    series2 = get_read_lengths(lib_keys[1])
    max1 = series1.max()
    max2 = series2.max()
    fig = go.Figure()
    bin_size = 10
    fig.add_trace(go.Histogram(x=series1.to_list(),
                               histnorm="probability",
                               name=lib_keys[0],
                               nbinsx=int((max1+1) / bin_size)))
    fig.add_trace(go.Histogram(x=series2.to_list(),
                               histnorm="probability",
                               name=lib_keys[1],
                               nbinsx=int((max2+1) / bin_size)))
    # TODO: get below working!!
    # fig.add_trace(px.ecdf(series1))
    # fig.add_trace(px.ecdf(series2))
    
    # Overlay both histograms
    fig.update_layout(barmode='overlay')
    # Reduce opacity to see both histograms
    fig.update_traces(opacity=0.75)
    fig.update_xaxes(title_text="Called Read Length")
    fig.update_yaxes(title_text="Probability",
                     fixedrange=True)
    fig.show()


if __name__ == '__main__':
    libraries = [
        "xrn-1",
        "totalRNA2",
    ]
    plotly_read_lengths(libraries)
