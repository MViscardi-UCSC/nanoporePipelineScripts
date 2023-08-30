#!/bin/bash

# I want to have this script run the ambiguousReads.py script on a short list of libraries:
for LIBRARY in newerN2 newerS5 newerS6 oldN2 oldS6
do
    python ambiguousReads.py $LIBRARY input
done
