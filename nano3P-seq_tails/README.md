# nano3P-seq_tails/README.md
## Marcus Viscardi     Dec 6, 2022
***
As of guppy_basecaller v6.4.2 the *--fast5_out* param is no longer supported!!
This is a major issue b/c tailfindr is reliant on the Events/move table that
is included inside the outputted fast5s from guppy.

The crappy workaround that I found was to download the last version that supported
this parameter from the ONT website (by messing with the download URL*) **[v6.3.8]**,
and unpacking that in my scripts folder. In the future we can specify if we want
to use the old version, or whatever the newest version apt has installed. We can
do this by specifying the path of the guppy we call:

- *guppy_basecaller* - this will run whatever version apt has installed to the PATH
- */data16/marcus/scripts/ont-guppy/bin/guppy_basecaller* - will always be v6.3.8!!

Final *TODO* here will be adding a flag to the pipeline to allow choice of which\
basecaller version to use!
***
*The URL to the v6.3.8 for download:
https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.3.8_linux64.tar.gz

You can find the newest version through the ONT website or apt (b/c we added ONT
stuff to the PATH).

