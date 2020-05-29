"""proc.py
"""

# This file has been replaced by jupyter/reading_andi.ipynb

import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

from pyms.GCMS.IO.ANDI import ANDI_reader

# read the raw data
andi_file = data_directory / "gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# print info
data.info()

# write data to output file. This will create
# two ascii data tables, data.I.csv and data.mz.csv
# with intensities and m/z values
data.write("output/data")
