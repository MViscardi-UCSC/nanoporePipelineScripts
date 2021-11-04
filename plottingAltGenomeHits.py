"""
plottingAltGenomeHits.py
Marcus Viscardi,    November 03, 2021

General goal of this is going to be to plot what the yeast genome
hits look like between the 210706 and the 210709 libraries

Going to just use bash word count -line to get a quick answer to this question:

wc -l /data16/marcus/working/211101_nanoporeSoftLinks/21070*/output_dir/cat_files/*.sam
   3503015 210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/cat_files/cat.altGenome.sam
   49878 210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/cat_files/cat.sorted.mappedAndPrimary.sam
   107864 210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/cat_files/cat.sorted.sam

   10266 210709_NanoporeRun_totalRNA_0639_L3/output_dir/cat_files/cat.altGenome.sam
   89636 210709_NanoporeRun_totalRNA_0639_L3/output_dir/cat_files/cat.sorted.mappedAndPrimary.sam
   94095 210709_NanoporeRun_totalRNA_0639_L3/output_dir/cat_files/cat.sorted.sam
"""
import pandas as pd
import seaborn as sea
import matplotlib.pyplot as plt
from nanoporePipelineCommon import get_dt

nested_dict = [
    {"lib": "riboD",
     "sam": "elegans",
     "normed": True,
     "count": 107864 / (3503015 + 107864),
     },
    {"lib": "riboD",
     "sam": "yeast",
     "normed": True,
     "count": 3503015 / (3503015 + 107864),
     },
    {"lib": "totalRNA",
     "sam": "yeast",
     "normed": True,
     "count": 10266 / (94095 + 10266),
     },
    {"lib": "totalRNA",
     "sam": "elegans",
     "normed": True,
     "count": 94095 / (94095 + 10266),
     },
    {"lib": "riboD",
     "sam": "elegans",
     "normed": False,
     "count": 107864,
     },
    {"lib": "riboD",
     "sam": "yeast",
     "normed": False,
     "count": 3503015,
     },
    {"lib": "totalRNA",
     "sam": "yeast",
     "normed": False,
     "count": 10266,
     },
    {"lib": "totalRNA",
     "sam": "elegans",
     "normed": False,
     "count": 94095,
     },
]
df = pd.DataFrame(nested_dict)
# sea.barplot(data=df, x="lib", y="count", hue="sam")
# plt.show()

# sea.histplot(data=df,
#              x="lib",
#              weights="count",
#              hue="sam",
#              multiple="stack",
#              shrink=0.8,
#              )
# g = sea.FacetGrid(df, row="normed", sharey=False)
# g.map(sea.histplot,
#       data=df,
#       x="lib",
#       weights="count",
#       hue="sam",
#       multiple="stack",
#       shrink=0.8,
#       )
fig, axes = plt.subplots(2, 1, figsize=(5, 5),
                         sharey=False, sharex=True)
sea.histplot(data=df[~df["normed"]],
             ax=axes[0],
             x="lib",
             weights="count",
             hue="sam",
             multiple="stack",
             shrink=0.8,
             legend=True,
             )
sea.histplot(data=df[df["normed"]],
             ax=axes[1],
             x="lib",
             weights="count",
             hue="sam",
             multiple="stack",
             shrink=0.8,
             legend=False,
             )
save_path = f"/home/marcus/211025_polyA-vs-totalRNA_plotsAndThings/" \
            f"{get_dt(for_file=True)}_yeastContaminationInTotalRNA1.svg"
plt.savefig(save_path)
plt.show()
print("yay?")
