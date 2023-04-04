"""
loggingForPipeline.py
Marcus Viscardi,    February 08, 2023

So I think I need to set up a better form of logging for the pipeline...

Currently, I have a few of the steps store their run information into the logging folder in the output_dir,
but this is massively incomplete and does not retain overall run information!

Ideally I would have both:
- console outputs (so that I know what's going on)
AND
- overall logging to track what ran, how it ran, and when it ran

I think utilizing the python logging module is the obvious solution for this.
But it might take some messing with to get it RIGHT right.

Additionally, it would be nice to store the Defaults, SettingsFile, and CLI options into the top of the storage dict,
but IDK how I'll do that. The problem is that I don't know where the output_dir/logs directory is until I have parsed
all of these things out!!
"""

import logging
from step0_nanopore_pipeline import live_cmd_call


if __name__ == '__main__':
    # logging_example.py

    import logging

    # Create a custom logger
    logger = logging.getLogger("Marcus")

    # Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler('file.log')
    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.ERROR)

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    logger.warning('This is a warning')
    logger.error('This is an error')
    logger.info('This is some nice info')
