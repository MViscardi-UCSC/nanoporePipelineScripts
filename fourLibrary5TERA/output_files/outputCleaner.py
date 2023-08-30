"""
outputCleaner.py
Marcus Viscardi,    June 26, 2023

This is a quick script to reorganize outputFiles based on their dates
"""
from pathlib import Path
from glob import glob
from sys import argv

if len(argv) == 1:
    # Just use the current directory
    output_dir = Path.cwd()
elif len(argv) == 2:
    # Use the directory given
    output_dir = Path(argv[1])
else:
    raise RuntimeError("Too many arguments given!")

# Get all the files in the directory
all_files = glob(str(output_dir / "*"))
# Only keep files that start with 6 digits (the date)
all_files = [Path(file) for file in all_files if file.split("/")[-1][:6].isdigit()]

all_file_dates = set([file.name[:6] for file in all_files])

for date in all_file_dates:
    # Make a new directory for each date
    date_dir = output_dir / f"{date}_outputs"
    date_dir.mkdir(exist_ok=True)
    # Move all files with that date into the new directory
    for file in all_files:
        if file.name[:6] == date and file.is_file():
            file.rename(date_dir / file.name)
