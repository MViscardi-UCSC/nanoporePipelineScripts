#!/usr/bin/python3
"""
roach_bulkTailPlotting.py
Marcus Viscardi     Feb 9, 2021

Going to plot all called tail lengths to try an compare the distribution to that of 
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sea
import numpy as np
import pandas as pd

from roach_compareReplicates_run import load_comp_on_genes

pd.set_option("display.max_columns", None)


def plot_by_genes(compressed_on_genes_filepath):
    # Thinking about it: this is probably not the best move as it
    #   skews the overall tail distribution towards the lowly
    #   expressed genes' tail lengths.
    df = load_comp_on_genes(compressed_on_genes_filepath)
    print(df.head(5))
    pass


def plot_all_tails(compressed_on_genes_filepath, max_value=250, window_size=3, subsample=None):
    def bin_and_plot_df(plotting_df, num_bins, rolling_window_size, title=None, plot=True):
        # Tuple with the range of numbers that are basically the indexes of our bins (+1)
        index_list = range(1, num_bins + 1)
        # bin our values into num_bins bins
        binned_series = pd.cut(plotting_df.plot_values, bins=num_bins, labels=index_list).value_counts()
        # reorder binned series based on indexes (which correspond to Tail length here)
        binned_ordered_series = binned_series.sort_index()
        # find moving average for each bin
        rolling_average = binned_ordered_series.rolling(rolling_window_size,
                                                        center=True,
                                                        min_periods=int(rolling_window_size / 2)
                                                        ).mean()
        # # Converting to a list makes it easier to pass into MatPlotLib
        # binned_ordered_list = binned_ordered_series.values
        # bin_factor allows us to translate the bin indexes back to the UTR_length values
        bin_factor = int(plotting_df.plot_values.max() / num_bins)

        fig, axs = plt.subplots()
        # # plot binned data
        # axs.plot([x * bin_factor for x in index_list], binned_ordered_list)

        # plot moving average
        x = [x * bin_factor for x in index_list]
        if plot:
            axs.plot(x,
                     rolling_average,
                     label=f"Simple Moving Average, window={rolling_window_size}",
                     linewidth=2.3, color=(0 / 255, 114 / 255, 178 / 255, 255 / 255))

            # plot histogram with data binned independently by matplotlib (slow)
            axs.hist(plotting_df.plot_values.values,
                     bins=num_bins,
                     label=f"Binned Histogram, bins={num_bins}",
                     color=(213 / 255, 94 / 255, 0 / 255, 255 / 255))

            axs.legend()
            if not title:
                axs.set_title("Tail Lengths from Roach et al. 2020")
            else:
                axs.set_title(title)
            plt.xlabel(f"Tail Lengths (nt)")  # TODO: Change this to 'fraction of reads'
            plt.ylabel(f"Number of UTRs Counted,\nn={plotting_df.shape[0]}")

            plt.xticks(ticks=(30, 45, 60), labels=("30nt", "45nt", "~60nt"))
            plt.xlim((25, 70))

            plt.tight_layout()
            plt.show()
        return x, rolling_average, binned_ordered_series

    df = load_comp_on_genes(compressed_on_genes_filepath)
    np_tails = np.hstack(df["polya_lengths"].to_numpy()).astype(float)
    print(f"Total called tail count ({compressed_on_genes_filepath[68:72]}): {len(np_tails)}")

    plot_df = pd.Series(np_tails).to_frame(name="plot_values")
    if subsample:
        plot_df = plot_df.sample(subsample)
    # Change all values above max_value down to max_value
    plot_df.loc[plot_df.plot_values >= max_value, 'plot_values'] = max_value

    plot_title = f"Tail Lengths from Roach et al. 2020 ({compressed_on_genes_filepath[68:72]})"

    x, y, raw_y = bin_and_plot_df(plot_df, num_bins=int(plot_df.plot_values.max()), rolling_window_size=window_size,
                                  title=plot_title,
                                  plot=False,
                                  )
    return x, y, raw_y


def roach_main():
    ABS_PATH = f"/data16/marcus/working/210106_nanopolish_wRoachData"
    OUTDIR1 = f"{ABS_PATH}/output_dir_adultrep1"
    OUTDIR2 = f"{ABS_PATH}/output_dir_adultrep2"
    cutoff_length, rolling_window, subsample = 250, 3, 2000
    xy1 = plot_all_tails(f"{OUTDIR1}/210124_compressed_on_genes.tsv",
                         max_value=cutoff_length, window_size=rolling_window,
                         subsample=subsample)
    xy2 = plot_all_tails(f"{OUTDIR2}/210128_compressed_on_genes.tsv",
                         max_value=cutoff_length, window_size=rolling_window,
                         subsample=subsample)
    for i, (x, y, raw_y) in enumerate((xy1, xy2)):
        y = y.to_list()
        y = [i / sum(y) for i in y]
        # cur_count = 0
        # cdf = []
        # for value in y:
        #     cur_count += value
        #     cdf.append(cur_count)
        # y = [i/max(cdf) for i in cdf]
        p = sea.lineplot(x=x, y=y, label=f"Adult Rep{i + 1}")
    p.set(title=f"Moving Average of All Tail Lengths from Roach et al. 2020 (samp: {subsample})",
          xlabel="Tail length (nt)", ylabel="Fraction of calls at 'x' length")
    p.text(50, 0.019, f"Moving average window size = {rolling_window}nt")
    p.text(cutoff_length + 5, 0, f"(Values greater than {cutoff_length} rounded down)", rotation=90)
    plt.show()


if __name__ == '__main__':
    outputDir = "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir"
    compressed_path = f"{outputDir}/merge_files/210601_05:05:54PM_compressedOnGenes.tsv"
    cutoff_length, rolling_window, subsample = 200, 3, 2000

    x, y, raw_y = plot_all_tails(compressed_path, # subsample=subsample,
                                 max_value=cutoff_length, window_size=rolling_window)
    y = y.to_list()
    fractional_y = [i / sum(y) for i in y]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    p = sea.lineplot(x=x, y=y, ax=ax1)
    p.set(title=f"Moving Average of All Tail Lengths (210528_polyA-dRNA)",
          xlabel="Tail length (nt)", ylabel="Number of calls at 'x' length")
    p.text(50, 0.019, f"Moving average window size = {rolling_window}nt")
    p.text(cutoff_length + 5, 0, f"(Values greater than {cutoff_length} rounded down)", rotation=90)

    # Create a Rectangle patch
    rect = patches.Rectangle((20, 5000), 50, 10000, linewidth=1, edgecolor='grey', facecolor='none')
    # Add the patch to the Axes
    ax1.add_patch(rect)

    ax2 = plt.axes([.55, .55, .3, .3])
    p2 = sea.lineplot(x=x, y=y, ax=ax2)
    for line_x in [30, 45, 60]:
        ax2.axvline(x=line_x, c='gray')
    ax2.set_xlim(20, 70)
    ax2.set_ylim(5000, 15000)
    plt.setp(ax2, xticks=[30, 45, 60], yticks=[5000, 10000, 15000])
    
    for file_format in ["svg", "png"]:
        plt.savefig(f"{outputDir}/plots/210603_tailLengths.{file_format}", format=file_format)
    plt.show()
