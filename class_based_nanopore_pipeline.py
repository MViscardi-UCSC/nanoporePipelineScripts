"""
class_based_nanopore_pipeline.py
Marcus Viscardi,    March 24, 2023

This is a semi temporary file. The goal is to move the majority
of my pipeline functionality to a class based system. The main
reason for this is that there are many flags that should be carried
between the separate steps/functions of the pipeline.

I am not sure if this is the right approach to the problem at hand,
but it is the best way I can think of!

This will also provide an opportunity to fix some shape edges
that the current pipeline has, including:

 - Lack of comprehensive logging functionality...
    
    - A large hangup in the past has been getting this to
      work with the called shell functions, but I should be
      able to do with now with a better understanding of how
      this works.
 
 - Inability to have some steps function on their own
    
    - For example, the basecaller doesn't concatenate to
      a single fastq because the "concatenate" step is its
      own thing!
 
 - Ability to diagnose common errors that occur...
   
   - The problem here is that most (all) of my pathways to
     troubleshoot errors that come up are rattling around in
     my head. These should really be integrated into the code
     so that future users don't have to rely on me. At least,
     I can make an effort to solve the easy ones!
   
   - An easy example of this is an empty fastq being passed to
     early steps of pipeline. This error will not be caught until
     WAY down the line, and should be easy to identify and notify
     the user of.
 
 - Better path handling with pathlib. Half of the current errors
   that occur are due to missing files or lack of ownership/access.

Additionally, this will be an opportunity to update with some disparate
pieces of code I have around this repo.

 - Integrate new standards' functionality, drop old standards.
 
 - Integrate Nano3P functionality!!
 
 - Add some general QC steps to produce print-outs and plots.
"""

