## From JA:
* <s>(i) figure or subpanel containing poly(A) tail information, and reproducibility of measurements at various read count cutoffs</s>
* (ii) supp figure or subpanel reproducibility of abundance (avg count, or SD/avg) at various read count cutoffs, to justify the counts we chose
* (iii) supp figure or subpanel  with “saturation analysis”, showing how many genes we contain reads for.

## Planning:
### Target (i): Tail information and reproducibility
I have done this analysis at some point. I think it predated the poly(A) paper?! I will need to find the code where I decided on the cutoffs for reads.
3 plots:
1. Standards: violin with 10, 15, and 60 bp tails
    - I produced these recently for Josh, just polish them and pop them in!
2. Subsampled tails: x-axis = number_of_test_samples, y-axis = mean_tail_length (and SEM)
    - This will be a plot for each standard.
3. Reproducibility of tail uniqueness (for CDFs we use):
    - CDF for subsets of reads per standard (eg. 5, 10, 25, 50)
    - Perform each subsetting 100 times or so, this will give us a bunch of CDFs
    - Then take the middle 95% of the CDFs and the mean of the middle 95% of the CDFs for each tail length
    - Plot the mean CDF and the 95% +/- CDFs for each tail standard (10, 15, & 60)
    - This will give us a good way to show that we can differentiate between tails and that we can do so reproducibly!
### Target (ii): Reproducibility of abundance
X-axis = RPM/read counts for per gene
Y-axis = SD/avg RPM/read counts per gene

This could be an overall scatter, or instead take the average of the SD/avg for each read_count bin and plot that as a line. This should show that at some point the data doesn't keep getting "better" with more depth

Another way to target the same idea is just a rocket plot with a vertical and horizontal cutoff showing what we identified as out cutoff!


### Target (iii): Saturation analysis
I did this for the Poly(A) paper in [subsamplingReadsVsProteinCoding.ipynb](../polyA_manuscriptPostReviewScripts/subsamplingReadsVsProteinCoding.ipynb). I just adjusted some inputs and ran it again... Worked great!

Josh thinks we could additionally have a plot with:
X-axis = reads_per_gene cutoff (1 to 200?)
Y-axis = cumulative number of genes above X cutoff